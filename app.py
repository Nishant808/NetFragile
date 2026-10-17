"""
NetFragile 2.0: Multi-Network Comparative Analysis Tool

Interactive Dash application for analyzing biological network robustness
with support for multiple EdgeR DEG files, comparative analysis, and
normalizing target discovery.
"""

import dash
from dash import dcc, html, Input, Output, State, callback, ctx
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import plotly.express as px
import networkx as nx
import pandas as pd
import numpy as np
import base64
import io
import json
import time
from typing import List, Dict, Any, Optional, Tuple
import logging

# Import utilities
from utils.file_parser import EdgeRParser, parse_deg_file
from utils.network_builder import NetworkBuilder, build_network_from_genes
from utils.robustness import RobustnessAnalyzer, simulate_removal
from utils.comparison import NetworkComparator, calculate_similarity_matrix
from utils.normalizing import NormalizingTargetFinder

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
HUMAN_TAXID = 9606
DEFAULT_REQUIRED_SCORE = 700
DEFAULT_ADD_NODES = 200


class NetworkData:
    """Container for network data and analysis results."""
    
    def __init__(self, factor_name: str):
        self.factor_name = factor_name
        self.gene_list = []
        self.graph = None
        self.degree_cent = None
        self.betweenness_cent = None
        self.rand_curve = None  # (x, y)
        self.targ_curve = None  # (x, y)
        self.metrics = {}
        self.metadata = {}
        self.centrality_results = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'factor_name': self.factor_name,
            'gene_list': self.gene_list,
            'metrics': self.metrics,
            'metadata': self.metadata,
            'rand_curve': {'x': self.rand_curve[0], 'y': self.rand_curve[1]} if self.rand_curve else None,
            'targ_curve': {'x': self.targ_curve[0], 'y': self.targ_curve[1]} if self.targ_curve else None,
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'NetworkData':
        """Create from dictionary."""
        nd = cls(data['factor_name'])
        nd.gene_list = data.get('gene_list', [])
        nd.metrics = data.get('metrics', {})
        nd.metadata = data.get('metadata', {})
        
        if data.get('rand_curve'):
            nd.rand_curve = (data['rand_curve']['x'], data['rand_curve']['y'])
        if data.get('targ_curve'):
            nd.targ_curve = (data['targ_curve']['x'], data['targ_curve']['y'])
        
        return nd


class AppState:
    """Global application state."""
    
    def __init__(self):
        self.networks = {}  # {factor_name: NetworkData}
        self.uploaded_files = {}  # {factor_name: file_content}
        self.similarity_matrix = None
        self.normalizing_results = None
        self.comparison_metrics = None
        self.last_analysis_time = None


# Initialize app state
app_state = AppState()

# Initialize Dash app
app = dash.Dash(__name__, suppress_callback_exceptions=True)
app.title = "NetFragile 2.0"


# ============================================================================
# LAYOUT
# ============================================================================

app.layout = html.Div(
    style={
        'fontFamily': '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif',
        'padding': '0',
        'margin': '0',
        'minHeight': '100vh',
        'background': 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
    },
    children=[
        # Header with shadow
        html.Div(
            style={
                'textAlign': 'center',
                'padding': '40px 20px 20px',
                'background': 'rgba(255, 255, 255, 0.95)',
                'boxShadow': '0 4px 6px rgba(0,0,0,0.1)',
                'marginBottom': '30px',
                'borderRadius': '0 0 20px 20px',
            },
            children=[
                html.Div(
                    style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'justifyContent': 'center',
                        'gap': '20px',
                        'marginBottom': '10px',
                    },
                    children=[
                        html.Div(
                            style={
                                'fontSize': '48px',
                                'color': '#667eea',
                            },
                            children=['🧬']
                        ),
                        html.H1("NetFragile 2.0", 
                               style={
                                   'color': '#2C3E50',
                                   'margin': '0',
                                   'fontSize': '42px',
                                   'fontWeight': '700',
                                   'textShadow': '2px 2px 4px rgba(0,0,0,0.1)',
                               }),
                    ]
                ),
                html.H3("Multi-Network Comparative Analysis for Drug Target Discovery", 
                       style={
                           'color': '#4A5568',
                           'margin': '10px 0 0 0',
                           'fontSize': '20px',
                           'fontWeight': '500',
                       }),
                html.P(
                    "Upload EdgeR DEG files → Build Networks → Find Normalizing Targets",
                    style={
                        'color': '#718096',
                        'margin': '10px 0 0 0',
                        'fontSize': '14px',
                        'fontStyle': 'italic',
                    }
                ),
            ]
        ),
        
        # Main content container
        html.Div(
            style={
                'maxWidth': '1400px',
                'margin': '0 auto',
                'padding': '0 20px 40px',
            },
            children=[
                # Tabs with modern styling
                html.Div(
                    style={
                        'background': 'rgba(255, 255, 255, 0.9)',
                        'borderRadius': '12px',
                        'boxShadow': '0 4px 12px rgba(0,0,0,0.1)',
                        'padding': '0 10px',
                        'marginBottom': '20px',
                    },
                    children=[
                        dcc.Tabs(
                            id='main-tabs',
                            value='upload',
                            style={'marginBottom': '0'},
                            colors={
                                'border': '#E2E8F0',
                                'primary': '#667eea',
                                'background': 'transparent',
                            },
                            children=[
                                dcc.Tab(
                                    label='📤 Upload',
                                    value='upload',
                                    style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                    },
                                    selected_style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                        'backgroundColor': '#667eea',
                                        'color': 'white',
                                        'borderTop': '3px solid #667eea',
                                    }
                                ),
                                dcc.Tab(
                                    label='🔍 Individual',
                                    value='individual',
                                    style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                    },
                                    selected_style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                        'backgroundColor': '#667eea',
                                        'color': 'white',
                                        'borderTop': '3px solid #667eea',
                                    }
                                ),
                                dcc.Tab(
                                    label='📊 Compare',
                                    value='comparison',
                                    style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                    },
                                    selected_style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                        'backgroundColor': '#667eea',
                                        'color': 'white',
                                        'borderTop': '3px solid #667eea',
                                    }
                                ),
                                dcc.Tab(
                                    label='🎯 Targets',
                                    value='normalizing',
                                    style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                    },
                                    selected_style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                        'backgroundColor': '#667eea',
                                        'color': 'white',
                                        'borderTop': '3px solid #667eea',
                                    }
                                ),
                                dcc.Tab(
                                    label='🌐 Graphs',
                                    value='graphs',
                                    style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                    },
                                    selected_style={
                                        'padding': '12px 24px',
                                        'fontWeight': '600',
                                        'fontSize': '14px',
                                        'borderRadius': '8px 8px 0 0',
                                        'backgroundColor': '#667eea',
                                        'color': 'white',
                                        'borderTop': '3px solid #667eea',
                                    }
                                ),
                            ]
                        ),
                    ]
                ),
                
                # Tab content
                html.Div(id='tab-content'),
                
                # Store components
                dcc.Store(id='networks-store'),
                dcc.Store(id='analysis-store'),
                dcc.Store(id='normalizing-results-store'),
                dcc.Download(id='download-csv'),
                dcc.Download(id='download-json'),
                dcc.Download(id='download-normalizing-csv'),
            ]
        ),
        
        # Footer
        html.Div(
            style={
                'textAlign': 'center',
                'padding': '20px',
                'color': 'rgba(255, 255, 255, 0.8)',
                'fontSize': '14px',
            },
            children=[
                html.P("NetFragile 2.0 | Built with ❤️ using Dash & Plotly | ", 
                       style={'margin': '0'}),
                html.A("View on GitHub", 
                       href="https://github.com/Nishant808/NetFragile",
                       style={'color': 'white', 'textDecoration': 'underline'}),
            ]
        ),
    ]
)


# ============================================================================
# TAB LAYOUTS
# ============================================================================

def upload_tab_layout():
    """Layout for the upload tab."""
    return html.Div([
        html.Div(
            style={
                'background': 'rgba(255, 255, 255, 0.95)',
                'padding': '30px',
                'borderRadius': '12px',
                'boxShadow': '0 4px 12px rgba(0,0,0,0.15)',
                'marginBottom': '20px',
                'backdropFilter': 'blur(10px)',
            },
            children=[
                html.Div(
                    style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'gap': '15px',
                        'marginBottom': '20px',
                    },
                    children=[
                        html.Div('📤', style={'fontSize': '32px'}),
                        html.H3("Upload EdgeR DEG Files", 
                               style={
                                   'color': '#2C3E50',
                                   'margin': '0',
                                   'fontSize': '24px',
                                   'fontWeight': '600',
                               }),
                    ]
                ),
                html.P(
                    "Upload your EdgeR differential expression output files (CSV/TSV). "
                    "Each file will be used to construct a network for comparative analysis.",
                    style={
                        'color': '#4A5568',
                        'marginBottom': '25px',
                        'fontSize': '15px',
                        'lineHeight': '1.6',
                    }
                ),
                
                # File upload section
                html.Div(id='file-upload-section'),
                
                html.Button(
                    '➕ Add Another File',
                    id='add-file-button',
                    n_clicks=0,
                    style={
                        'backgroundColor': 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                        'color': 'white',
                        'border': 'none',
                        'padding': '12px 24px',
                        'borderRadius': '8px',
                        'cursor': 'pointer',
                        'marginTop': '15px',
                        'fontSize': '14px',
                        'fontWeight': '600',
                        'boxShadow': '0 4px 8px rgba(102, 126, 234, 0.3)',
                        'transition': 'all 0.3s ease',
                    }
                ),
                
                # Settings
                html.Hr(style={'margin': '30px 0', 'border': 'none', 'borderTop': '2px solid #E2E8F0'}),
                html.Div(
                    style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'gap': '12px',
                        'marginBottom': '20px',
                    },
                    children=[
                        html.Div('⚙️', style={'fontSize': '24px'}),
                        html.H4("Analysis Settings", 
                               style={
                                   'color': '#2C3E50',
                                   'margin': '0',
                                   'fontSize': '18px',
                                   'fontWeight': '600',
                               }),
                    ]
                ),
                
                html.Div(
                    style={'display': 'grid', 'gridTemplateColumns': '1fr 1fr', 'gap': '20px', 'marginTop': '15px'},
                    children=[
                        html.Div([
                            html.Label("STRING Confidence Score:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                            dcc.Dropdown(
                                id='required-score',
                                options=[
                                    {'label': 'Low (400)', 'value': 400},
                                    {'label': 'Medium (700)', 'value': 700},
                                    {'label': 'High (900)', 'value': 900},
                                ],
                                value=DEFAULT_REQUIRED_SCORE,
                                style={'width': '100%'}
                            ),
                        ]),
                        html.Div([
                            html.Label("Additional Nodes:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                            dcc.Dropdown(
                                id='add-nodes',
                                options=[
                                    {'label': '0', 'value': 0},
                                    {'label': '50', 'value': 50},
                                    {'label': '100', 'value': 100},
                                    {'label': '200', 'value': 200},
                                    {'label': '500', 'value': 500},
                                ],
                                value=DEFAULT_ADD_NODES,
                                style={'width': '100%'}
                            ),
                        ]),
                    ]
                ),
                
                html.Div(
                    style={'display': 'grid', 'gridTemplateColumns': '1fr 1fr', 'gap': '20px', 'marginTop': '15px'},
                    children=[
                        html.Div([
                            html.Label("Filter by p-value:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                            dcc.Input(
                                id='pvalue-threshold',
                                type='number',
                                value=0.05,
                                step=0.01,
                                min=0,
                                max=1,
                                style={'width': '100%', 'padding': '8px', 'borderRadius': '5px', 'border': '1px solid #ddd'}
                            ),
                        ]),
                        html.Div([
                            html.Label("Filter by |logFC|:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                            dcc.Input(
                                id='logfc-threshold',
                                type='number',
                                value=1.0,
                                step=0.5,
                                min=0,
                                style={'width': '100%', 'padding': '8px', 'borderRadius': '5px', 'border': '1px solid #ddd'}
                            ),
                        ]),
                    ]
                ),
                
                # Filter type
                html.Div(
                    style={'marginTop': '15px'},
                    children=[
                        html.Label("Apply Filters:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                        dcc.Checklist(
                            id='apply-filters',
                            options=[
                                {'label': ' p-value', 'value': 'pvalue'},
                                {'label': ' |logFC|', 'value': 'logfc'},
                            ],
                            value=['pvalue', 'logfc'],
                            inline=True,
                            style={'fontSize': '14px'}
                        ),
                    ]
                ),
                
                # Run analysis button
                html.Hr(style={'margin': '30px 0', 'border': 'none', 'borderTop': '2px solid #E2E8F0'}),
                html.Div(
                    style={'textAlign': 'center'},
                    children=[
                        html.Button(
                            '🚀 RUN COMPARATIVE ANALYSIS',
                            id='run-analysis-button',
                            n_clicks=0,
                            style={
                                'backgroundColor': 'linear-gradient(135deg, #27AE60 0%, #229954 100%)',
                                'color': 'white',
                                'border': 'none',
                                'padding': '16px 36px',
                                'borderRadius': '10px',
                                'cursor': 'pointer',
                                'fontSize': '16px',
                                'fontWeight': '700',
                                'boxShadow': '0 6px 12px rgba(39, 174, 96, 0.3)',
                                'transition': 'all 0.3s ease',
                                'letterSpacing': '1px',
                            }
                        ),
                    ]
                ),
                
                # Progress and status
                html.Div(id='upload-status', style={'marginTop': '20px', 'textAlign': 'center'}),
                dcc.Loading(
                    id='analysis-loading',
                    type='circle',
                    children=html.Div(id='analysis-output')
                ),
            ]
        ),
    ])


def individual_tab_layout():
    """Layout for the individual networks tab."""
    return html.Div([
        html.Div(
            style={
                'background': 'rgba(255, 255, 255, 0.95)',
                'padding': '30px',
                'borderRadius': '12px',
                'boxShadow': '0 4px 12px rgba(0,0,0,0.15)',
                'backdropFilter': 'blur(10px)',
            },
            children=[
                html.Div(
                    style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'gap': '15px',
                        'marginBottom': '20px',
                    },
                    children=[
                        html.Div('🔍', style={'fontSize': '32px'}),
                        html.H3("Individual Network Analysis", 
                               style={
                                   'color': '#2C3E50',
                                   'margin': '0',
                                   'fontSize': '24px',
                                   'fontWeight': '600',
                               }),
                    ]
                ),
                
                # Network selector
                html.Div(
                    style={'marginBottom': '25px'},
                    children=[
                        html.Label("Select Network:", 
                                  style={
                                      'fontWeight': '600',
                                      'marginBottom': '8px',
                                      'display': 'block',
                                      'color': '#4A5568',
                                      'fontSize': '15px',
                                  }),
                        dcc.Dropdown(
                            id='individual-network-selector',
                            options=[],
                            value=None,
                            style={'width': '300px'},
                            placeholder="Choose a network..."
                        ),
                    ]
                ),
                
                # Network info
                html.Div(id='individual-network-info', style={'marginBottom': '30px'}),
                
                # Tabs for individual network
                html.Div(
                    style={
                        'background': 'rgba(255, 255, 255, 0.8)',
                        'borderRadius': '8px',
                        'padding': '0 10px',
                        'boxShadow': '0 2px 8px rgba(0,0,0,0.1)',
                    },
                    children=[
                        dcc.Tabs(
                            id='individual-subtabs',
                            value='robustness',
                            colors={
                                'border': '#E2E8F0',
                                'primary': '#667eea',
                                'background': 'transparent',
                            },
                            children=[
                                dcc.Tab(
                                    label='Robustness Curves',
                                    value='robustness',
                                    style={'padding': '10px 20px', 'fontWeight': '500'},
                                    selected_style={
                                        'padding': '10px 20px',
                                        'fontWeight': '600',
                                        'backgroundColor': '#667eea',
                                        'color': 'white',
                                    }
                                ),
                                dcc.Tab(
                                    label='Network Graph',
                                    value='graph',
                                    style={'padding': '10px 20px', 'fontWeight': '500'},
                                    selected_style={
                                        'padding': '10px 20px',
                                        'fontWeight': '600',
                                        'backgroundColor': '#667eea',
                                        'color': 'white',
                                    }
                                ),
                                dcc.Tab(
                                    label='Centrality Analysis',
                                    value='centrality',
                                    style={'padding': '10px 20px', 'fontWeight': '500'},
                                    selected_style={
                                        'padding': '10px 20px',
                                        'fontWeight': '600',
                                        'backgroundColor': '#667eea',
                                        'color': 'white',
                                    }
                                ),
                            ]
                        ),
                    ]
                ),
                
                html.Div(id='individual-subtab-content', style={'marginTop': '20px'}),
            ]
        ),
    ])


def comparison_tab_layout():
    """Layout for the comparison tab."""
    return html.Div([
        html.Div(
            style={
                'background': 'rgba(255, 255, 255, 0.95)',
                'padding': '30px',
                'borderRadius': '12px',
                'boxShadow': '0 4px 12px rgba(0,0,0,0.15)',
                'backdropFilter': 'blur(10px)',
            },
            children=[
                html.Div(
                    style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'gap': '15px',
                        'marginBottom': '20px',
                    },
                    children=[
                        html.Div('📊', style={'fontSize': '32px'}),
                        html.H3("Network Comparison", 
                               style={
                                   'color': '#2C3E50',
                                   'margin': '0',
                                   'fontSize': '24px',
                                   'fontWeight': '600',
                               }),
                    ]
                ),
                
                # Curve type selector
                html.Div(
                    style={'marginBottom': '25px'},
                    children=[
                        html.Label("Curve Type:", 
                                  style={
                                      'fontWeight': '600',
                                      'marginBottom': '10px',
                                      'display': 'block',
                                      'color': '#4A5568',
                                      'fontSize': '15px',
                                  }),
                        dcc.RadioItems(
                            id='comparison-curve-type',
                            options=[
                                {'label': ' Targeted Removal', 'value': 'targ'},
                                {'label': ' Random Removal', 'value': 'rand'},
                            ],
                            value='targ',
                            inline=True,
                            style={'fontSize': '14px'},
                            labelStyle={'marginRight': '20px'}
                        ),
                    ]
                ),
                
                # Overlay plot
                html.Div(
                    style={
                        'background': 'rgba(255, 255, 255, 0.8)',
                        'padding': '20px',
                        'borderRadius': '10px',
                        'marginBottom': '25px',
                        'boxShadow': '0 2px 8px rgba(0,0,0,0.05)',
                    },
                    children=[
                        html.H4("📈 Overlay Robustness Curves", 
                               style={'color': '#2C3E50', 'marginBottom': '15px', 'fontSize': '18px'}),
                        dcc.Graph(id='comparison-overlay-plot', config={'displayModeBar': True}),
                    ]
                ),
                
                # Similarity matrix
                html.Div(
                    style={
                        'background': 'rgba(255, 255, 255, 0.8)',
                        'padding': '20px',
                        'borderRadius': '10px',
                        'marginBottom': '25px',
                        'boxShadow': '0 2px 8px rgba(0,0,0,0.05)',
                    },
                    children=[
                        html.H4("🔗 Network Similarity Matrix", 
                               style={'color': '#2C3E50', 'marginBottom': '15px', 'fontSize': '18px'}),
                        html.Div(id='similarity-matrix'),
                    ]
                ),
                
                # Metrics table
                html.Div(
                    style={
                        'background': 'rgba(255, 255, 255, 0.8)',
                        'padding': '20px',
                        'borderRadius': '10px',
                        'boxShadow': '0 2px 8px rgba(0,0,0,0.05)',
                    },
                    children=[
                        html.H4("📊 Network Metrics", 
                               style={'color': '#2C3E50', 'marginBottom': '15px', 'fontSize': '18px'}),
                        html.Div(id='comparison-metrics'),
                    ]
                ),
            ]
        ),
    ])


def normalizing_tab_layout():
    """Layout for the normalizing targets tab."""
    return html.Div([
        html.Div(
            style={
                'background': 'rgba(255, 255, 255, 0.95)',
                'padding': '30px',
                'borderRadius': '12px',
                'boxShadow': '0 4px 12px rgba(0,0,0,0.15)',
                'backdropFilter': 'blur(10px)',
            },
            children=[
                html.Div(
                    style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'gap': '15px',
                        'marginBottom': '15px',
                    },
                    children=[
                        html.Div('🎯', style={'fontSize': '32px'}),
                        html.H3("Normalizing Target Discovery", 
                               style={
                                   'color': '#2C3E50',
                                   'margin': '0',
                                   'fontSize': '24px',
                                   'fontWeight': '600',
                               }),
                    ]
                ),
                html.P(
                    "Find genes in the target network whose removal makes its robustness "
                    "profile resemble the reference network.",
                    style={
                        'color': '#4A5568',
                        'marginBottom': '20px',
                        'fontSize': '15px',
                        'lineHeight': '1.6',
                    }
                ),
                
                # Network selectors
                html.Div(
                    style={'display': 'grid', 'gridTemplateColumns': '1fr 1fr', 'gap': '20px', 'marginBottom': '25px'},
                    children=[
                        html.Div([
                            html.Label("Target Network:", 
                                      style={
                                          'fontWeight': '600',
                                          'marginBottom': '8px',
                                          'display': 'block',
                                          'color': '#4A5568',
                                          'fontSize': '15px',
                                      }),
                            dcc.Dropdown(
                                id='normalizing-target-network',
                                options=[],
                                value=None,
                                style={'width': '100%'},
                                placeholder="Select target network..."
                            ),
                        ]),
                        html.Div([
                            html.Label("Reference Network:", 
                                      style={
                                          'fontWeight': '600',
                                          'marginBottom': '8px',
                                          'display': 'block',
                                          'color': '#4A5568',
                                          'fontSize': '15px',
                                      }),
                            dcc.Dropdown(
                                id='normalizing-reference-network',
                                options=[],
                                value=None,
                                style={'width': '100%'},
                                placeholder="Select reference network..."
                            ),
                        ]),
                    ]
                ),
                
                # Run button
                html.Div(
                    style={'textAlign': 'center', 'marginBottom': '25px'},
                    children=[
                        html.Button(
                            '🔍 FIND NORMALIZING TARGETS',
                            id='find-targets-button',
                            n_clicks=0,
                            style={
                                'backgroundColor': 'linear-gradient(135deg, #E74C3C 0%, #C0392B 100%)',
                                'color': 'white',
                                'border': 'none',
                                'padding': '14px 30px',
                                'borderRadius': '10px',
                                'cursor': 'pointer',
                                'fontSize': '15px',
                                'fontWeight': '700',
                                'boxShadow': '0 6px 12px rgba(231, 76, 60, 0.3)',
                                'transition': 'all 0.3s ease',
                                'letterSpacing': '1px',
                            }
                        ),
                    ]
                ),
                
                # Progress
                dcc.Loading(
                    id='normalizing-loading',
                    type='circle',
                    children=html.Div(id='normalizing-status')
                ),
                
                # Results
                html.Div(id='normalizing-results'),
            ]
        ),
    ])


def graphs_tab_layout():
    """Layout for the network graphs tab."""
    return html.Div([
        html.Div(
            style={
                'background': 'rgba(255, 255, 255, 0.95)',
                'padding': '30px',
                'borderRadius': '12px',
                'boxShadow': '0 4px 12px rgba(0,0,0,0.15)',
                'backdropFilter': 'blur(10px)',
            },
            children=[
                html.Div(
                    style={
                        'display': 'flex',
                        'alignItems': 'center',
                        'gap': '15px',
                        'marginBottom': '20px',
                    },
                    children=[
                        html.Div('🌐', style={'fontSize': '32px'}),
                        html.H3("Network Visualizations", 
                               style={
                                   'color': '#2C3E50',
                                   'margin': '0',
                                   'fontSize': '24px',
                                   'fontWeight': '600',
                               }),
                    ]
                ),
                
                # Network selector
                html.Div(
                    style={'marginBottom': '25px'},
                    children=[
                        html.Label("Select Network:", 
                                  style={
                                      'fontWeight': '600',
                                      'marginBottom': '8px',
                                      'display': 'block',
                                      'color': '#4A5568',
                                      'fontSize': '15px',
                                  }),
                        dcc.Dropdown(
                            id='graph-network-selector',
                            options=[],
                            value=None,
                            style={'width': '300px'},
                            placeholder="Choose a network..."
                        ),
                    ]
                ),
                
                # Graph display
                html.Div(
                    style={
                        'background': 'rgba(255, 255, 255, 0.8)',
                        'padding': '20px',
                        'borderRadius': '10px',
                        'boxShadow': '0 2px 8px rgba(0,0,0,0.05)',
                    },
                    children=[
                        dcc.Graph(id='network-graph-display', config={'displayModeBar': True}),
                    ]
                ),
                
                # Graph info
                html.Div(id='graph-info', 
                        style={
                            'marginTop': '20px',
                            'padding': '15px',
                            'background': 'rgba(255, 255, 255, 0.6)',
                            'borderRadius': '8px',
                            'textAlign': 'center',
                            'fontSize': '15px',
                            'color': '#4A5568',
                        }),
            ]
        ),
    ])


# ============================================================================
# CALLBACKS
# ============================================================================

# File upload callbacks
@callback(
    Output('file-upload-section', 'children'),
    Input('add-file-button', 'n_clicks'),
    State('file-upload-section', 'children'),
    prevent_initial_call=True
)
def add_file_upload(n_clicks, children):
    """Add a new file upload component."""
    if children is None:
        children = []
    
    # Count existing uploads
    num_uploads = len([c for c in children if isinstance(c, dict) and 'props' in c and 'id' in c['props'] and 'file-upload' in c['props']['id']])
    
    new_upload = html.Div(
        id=f'file-upload-{num_uploads}',
        style={
            'background': '#F8F9FA',
            'padding': '15px',
            'borderRadius': '8px',
            'marginBottom': '15px',
            'border': '2px dashed #BDC3C7'
        },
        children=[
            html.Div(
                style={'display': 'flex', 'gap': '15px', 'alignItems': 'center'},
                children=[
                    dcc.Upload(
                        id=f'file-upload-{num_uploads}',
                        children=html.Button(
                            f'Select File {num_uploads + 1}',
                            style={
                                'backgroundColor': '#3498DB',
                                'color': 'white',
                                'border': 'none',
                                'padding': '10px 20px',
                                'borderRadius': '5px',
                                'cursor': 'pointer'
                            }
                        ),
                        multiple=False,
                        style={'flex': '1'}
                    ),
                    dcc.Input(
                        id=f'factor-name-{num_uploads}',
                        type='text',
                        placeholder=f'Factor {num_uploads + 1} Name',
                        value=f'Condition_{num_uploads + 1}',
                        style={
                            'flex': '1',
                            'padding': '10px',
                            'borderRadius': '5px',
                            'border': '1px solid #ddd'
                        }
                    ),
                    html.Button(
                        '✕',
                        id=f'remove-file-{num_uploads}',
                        n_clicks=0,
                        style={
                            'backgroundColor': '#E74C3C',
                            'color': 'white',
                            'border': 'none',
                            'padding': '10px',
                            'borderRadius': '5px',
                            'cursor': 'pointer',
                            'width': '40px'
                        }
                    ),
                ]
            ),
            html.Div(id=f'file-info-{num_uploads}', style={'marginTop': '10px', 'fontSize': '12px', 'color': '#7F8C8D'})
        ]
    )
    
    children.append(new_upload)
    return children


@callback(
    [Output(f'file-info-{i}', 'children') for i in range(20)],
    [Input(f'file-upload-{i}', 'contents') for i in range(20)],
    [State(f'file-upload-{i}', 'filename') for i in range(20)],
    prevent_initial_call=True
)
def update_file_info(*args):
    """Update file info displays."""
    outputs = []
    
    for i in range(0, len(args), 2):
        content = args[i]
        filename = args[i + 1]
        
        if content is not None:
            # Parse file to get info
            try:
                parser = EdgeRParser()
                info = parser.get_file_info(filename)
                if 'error' in info:
                    outputs.append(f"Error: {info['error']}")
                else:
                    outputs.append(f"✓ {filename} | {info['row_count']} rows | Gene column: {info['detected_gene_column']}")
            except Exception as e:
                outputs.append(f"Error: {str(e)}")
        else:
            outputs.append("No file selected")
    
    return outputs


# Main analysis callback
@callback(
    [Output('networks-store', 'data'),
     Output('upload-status', 'children'),
     Output('analysis-output', 'children')],
    Input('run-analysis-button', 'n_clicks'),
    [State(f'file-upload-{i}', 'contents') for i in range(20)],
    [State(f'factor-name-{i}', 'value') for i in range(20)],
    [State(f'file-upload-{i}', 'filename') for i in range(20)],
    State('required-score', 'value'),
    State('add-nodes', 'value'),
    State('pvalue-threshold', 'value'),
    State('logfc-threshold', 'value'),
    State('apply-filters', 'value'),
    prevent_initial_call=True
)
def run_comparative_analysis(n_clicks, *args):
    """Run the complete comparative analysis pipeline."""
    if n_clicks == 0:
        raise PreventUpdate
    
    # Parse arguments
    file_contents = args[0:20]
    factor_names = args[20:40]
    filenames = args[40:60]
    required_score = args[60]
    add_nodes = args[61]
    pvalue_threshold = args[62]
    logfc_threshold = args[63]
    apply_filters = args[64]
    
    # Get valid uploads
    uploads = []
    for i in range(20):
        if file_contents[i] is not None and factor_names[i]:
            uploads.append({
                'content': file_contents[i],
                'filename': filenames[i],
                'factor_name': factor_names[i]
            })
    
    if not uploads:
        return {}, "Please upload at least one file.", ""
    
    # Clear previous state
    app_state.networks = {}
    app_state.uploaded_files = {}
    
    status_messages = []
    analysis_outputs = []
    
    # Process each file
    for upload in uploads:
        try:
            # Parse file
            content = upload['content']
            filename = upload['filename']
            factor_name = upload['factor_name']
            
            # Decode content
            content_type, content_string = content.split(',')
            decoded = base64.b64decode(content_string)
            
            # Save to temporary file
            temp_path = f"/tmp/{filename}"
            with open(temp_path, 'wb') as f:
                f.write(decoded)
            
            # Parse DEG file
            pvalue_filter = pvalue_threshold if 'pvalue' in apply_filters else None
            logfc_filter = logfc_threshold if 'logfc' in apply_filters else None
            
            gene_list = parse_deg_file(
                temp_path,
                pvalue_threshold=pvalue_filter,
                logfc_threshold=logfc_filter,
                min_genes=10
            )
            
            status_messages.append(f"✓ {factor_name}: {len(gene_list)} genes parsed")
            
            # Build network
            builder = NetworkBuilder(HUMAN_TAXID)
            graph, metadata = builder.build_network(
                gene_list,
                required_score=required_score,
                add_nodes=add_nodes
            )
            
            # Create NetworkData
            nd = NetworkData(factor_name)
            nd.gene_list = gene_list
            nd.graph = graph
            nd.metadata = metadata
            
            # Run robustness analysis
            analyzer = RobustnessAnalyzer()
            curves = analyzer.analyze_network(graph, ['random', 'degree'])
            nd.rand_curve = curves['random']
            nd.targ_curve = curves['degree']
            
            # Calculate centrality
            if graph.number_of_nodes() > 0:
                nd.degree_cent = nx.degree_centrality(graph)
                nd.betweenness_cent = nx.betweenness_centrality(graph)
                
                # Get top hubs
                top_hubs = sorted(nd.degree_cent.items(), key=lambda x: x[1], reverse=True)[:5]
                nd.centrality_results = {
                    'top_hubs': [n for n, v in top_hubs],
                    'top_bottlenecks': []
                }
                
                if graph.number_of_nodes() > 1:
                    top_bottlenecks = sorted(nd.betweenness_cent.items(), key=lambda x: x[1], reverse=True)[:5]
                    nd.centrality_results['top_bottlenecks'] = [n for n, v in top_bottlenecks]
            
            # Calculate metrics
            nd.metrics = {
                'nodes': graph.number_of_nodes(),
                'edges': graph.number_of_edges(),
                'density': nx.density(graph) if graph.number_of_nodes() > 1 else 0,
                'avg_degree': sum(dict(graph.degree()).values()) / graph.number_of_nodes() if graph.number_of_nodes() > 0 else 0,
            }
            
            # Store
            app_state.networks[factor_name] = nd
            app_state.uploaded_files[factor_name] = upload
            
            status_messages.append(f"✓ {factor_name}: Network built ({metadata['final_node_count']} nodes, {metadata['final_edge_count']} edges)")
            
        except Exception as e:
            status_messages.append(f"✗ {upload.get('factor_name', 'Unknown')}: Error - {str(e)}")
            logger.error(f"Error processing {upload.get('filename', 'unknown')}: {e}")
    
    # Calculate similarity matrix
    if len(app_state.networks) > 1:
        networks_data = {
            name: {
                'targ_curve': nd.targ_curve,
                'rand_curve': nd.rand_curve,
                'graph': nd.graph,
                'metrics': nd.metrics
            }
            for name, nd in app_state.networks.items()
        }
        
        comparator = NetworkComparator()
        app_state.similarity_matrix = comparator.calculate_similarity_matrix(networks_data)
        app_state.comparison_metrics = comparator.calculate_network_metrics(networks_data)
        
        status_messages.append("✓ Similarity matrix calculated")
    
    # Serialize networks for storage
    networks_serialized = {
        name: nd.to_dict()
        for name, nd in app_state.networks.items()
    }
    
    # Format status
    status_html = html.Div([
        html.H4("Analysis Complete!", style={'color': '#27AE60'}),
        html.Ul([html.Li(msg) for msg in status_messages])
    ])
    
    # Analysis summary
    analysis_summary = html.Div([
        html.H4(f"Successfully analyzed {len(app_state.networks)} networks", style={'color': '#2C3E50'}),
        html.P(f"Networks: {', '.join(app_state.networks.keys())}")
    ])
    
    return networks_serialized, status_html, analysis_summary


# Individual network tab callbacks
@callback(
    Output('individual-network-selector', 'options'),
    Input('networks-store', 'data'),
    prevent_initial_call=True
)
def update_individual_selector(networks_data):
    """Update network selector in individual tab."""
    if networks_data:
        options = [{'label': name, 'value': name} for name in networks_data.keys()]
        return options
    return []


@callback(
    Output('individual-network-info', 'children'),
    Input('individual-network-selector', 'value'),
    State('networks-store', 'data'),
    prevent_initial_call=True
)
def update_individual_info(selected_network, networks_data):
    """Update network info display."""
    if not selected_network or not networks_data:
        return html.Div("Select a network to view details.")
    
    nd = NetworkData.from_dict(networks_data[selected_network])
    
    return html.Div([
        html.H4(f"Network: {selected_network}", style={'color': '#2C3E50'}),
        html.Div(style={'display': 'grid', 'gridTemplateColumns': '1fr 1fr', 'gap': '15px'}, children=[
            html.Div([
                html.Strong("Nodes: "),
                html.Span(str(nd.metrics.get('nodes', 'N/A')))
            ]),
            html.Div([
                html.Strong("Edges: "),
                html.Span(str(nd.metrics.get('edges', 'N/A')))
            ]),
            html.Div([
                html.Strong("Density: "),
                html.Span(f"{nd.metrics.get('density', 0):.4f}")
            ]),
            html.Div([
                html.Strong("Avg Degree: "),
                html.Span(f"{nd.metrics.get('avg_degree', 0):.2f}")
            ]),
        ]),
        html.Hr(style={'margin': '15px 0'}),
        html.Div([
            html.Strong("Top Hubs: "),
            html.Span(", ".join(nd.centrality_results.get('top_hubs', [])[:5]))
        ]),
        html.Div([
            html.Strong("Top Bottlenecks: "),
            html.Span(", ".join(nd.centrality_results.get('top_bottlenecks', [])[:5]))
        ]),
    ])


@callback(
    Output('individual-subtab-content', 'children'),
    [Input('individual-subtabs', 'value'),
     Input('individual-network-selector', 'value')],
    State('networks-store', 'data'),
    prevent_initial_call=True
)
def update_individual_subtab(subtab, selected_network, networks_data):
    """Update content of individual network subtabs."""
    if not selected_network or not networks_data:
        return html.Div("Select a network to view details.")
    
    nd = NetworkData.from_dict(networks_data[selected_network])
    
    if subtab == 'robustness':
        return html.Div([
            dcc.Graph(
                id='individual-robustness-plot',
                figure=self._create_robustness_figure(nd)
            )
        ])
    elif subtab == 'graph':
        return html.Div([
            dcc.Graph(
                id='individual-network-graph',
                figure=self._create_network_figure(nd)
            )
        ])
    elif subtab == 'centrality':
        return html.Div([
            html.H4("Centrality Analysis", style={'color': '#2C3E50'}),
            html.Div([
                html.H5("Top Hubs (by Degree Centrality)", style={'color': '#2980B9'}),
                html.Ul([html.Li(gene) for gene in nd.centrality_results.get('top_hubs', [])]),
            ]),
            html.Div([
                html.H5("Top Bottlenecks (by Betweenness Centrality)", style={'color': '#2980B9'}),
                html.Ul([html.Li(gene) for gene in nd.centrality_results.get('top_bottlenecks', [])]),
            ]),
        ])
    
    return html.Div()


@callback(
    Output('comparison-overlay-plot', 'figure'),
    Input('comparison-curve-type', 'value'),
    State('networks-store', 'data'),
    prevent_initial_call=True
)
def update_comparison_overlay(curve_type, networks_data):
    """Update overlay plot in comparison tab."""
    if not networks_data:
        return go.Figure()
    
    fig = go.Figure()
    
    colors = px.colors.qualitative.Set2
    
    for i, (name, data) in enumerate(networks_data.items()):
        curve = data.get(f'{curve_type}_curve')
        if curve:
            x = curve['x']
            y = curve['y']
            color = colors[i % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode='lines+markers',
                name=name,
                line=dict(width=3, color=color),
                marker=dict(size=6),
                hovertemplate=f'{name}<br>Fraction Removed: %{{x:.2f}}<br>Integrity: %{{y:.2f}}<extra></extra>'
            ))
    
    fig.update_layout(
        title=f'Overlay Robustness Curves ({curve_type.capitalize()} Removal)',
        xaxis_title='Fraction of Network Removed',
        yaxis_title='Network Integrity (LCC Size)',
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        plot_bgcolor='white',
        height=600,
        hovermode='closest'
    )
    
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    
    return fig


@callback(
    Output('similarity-matrix', 'children'),
    Input('networks-store', 'data'),
    prevent_initial_call=True
)
def update_similarity_matrix(networks_data):
    """Update similarity matrix display."""
    if not networks_data or len(networks_data) < 2:
        return html.Div("Need at least 2 networks to calculate similarity.")
    
    # Recalculate similarity matrix
    networks_for_comparison = {
        name: {'targ_curve': (data['targ_curve']['x'], data['targ_curve']['y'])}
        for name, data in networks_data.items()
    }
    
    comparator = NetworkComparator()
    sim_matrix = comparator.calculate_similarity_matrix(networks_for_comparison)
    
    # Create table
    table_header = [html.Thead(html.Tr([html.Th("")] + [html.Th(name) for name in sim_matrix.columns]))]
    
    table_body = [
        html.Tbody([
            html.Tr([html.Th(row_name)] + [html.Td(f"{val:.2f}") for val in row])
            for row_name, row in sim_matrix.iterrows()
        ])
    ]
    
    return html.Table(
        table_header + table_body,
        style={'width': '100%', 'borderCollapse': 'collapse'}
    )


@callback(
    Output('comparison-metrics', 'children'),
    Input('networks-store', 'data'),
    prevent_initial_call=True
)
def update_comparison_metrics(networks_data):
    """Update comparison metrics display."""
    if not networks_data:
        return html.Div("No networks to display.")
    
    # Create metrics table
    metrics_list = []
    for name, data in networks_data.items():
        metrics_list.append({
            'Network': name,
            'Nodes': data['metrics'].get('nodes', 'N/A'),
            'Edges': data['metrics'].get('edges', 'N/A'),
            'Density': f"{data['metrics'].get('density', 0):.4f}",
            'Avg Degree': f"{data['metrics'].get('avg_degree', 0):.2f}"
        })
    
    df = pd.DataFrame(metrics_list)
    
    return html.Table(
        [
            html.Thead(html.Tr([html.Th(col) for col in df.columns]))
        ] + [
            html.Tbody([
                html.Tr([html.Td(df.iloc[i][col]) for col in df.columns])
                for i in range(len(df))
            ])
        ],
        style={'width': '100%', 'borderCollapse': 'collapse'}
    )


@callback(
    [Output('normalizing-target-network', 'options'),
     Output('normalizing-reference-network', 'options')],
    Input('networks-store', 'data'),
    prevent_initial_call=True
)
def update_normalizing_selectors(networks_data):
    """Update network selectors in normalizing tab."""
    if networks_data:
        options = [{'label': name, 'value': name} for name in networks_data.keys()]
        return options, options
    return [], []


@callback(
    [Output('normalizing-status', 'children'),
     Output('normalizing-results', 'children')],
    Input('find-targets-button', 'n_clicks'),
    State('normalizing-target-network', 'value'),
    State('normalizing-reference-network', 'value'),
    State('networks-store', 'data'),
    prevent_initial_call=True
)
def find_normalizing_targets(n_clicks, target_network, reference_network, networks_data):
    """Find normalizing targets between two networks."""
    if n_clicks == 0 or not target_network or not reference_network:
        raise PreventUpdate
    
    if target_network == reference_network:
        return "Please select different networks.", ""
    
    if not networks_data:
        return "No networks available.", ""
    
    # Get networks
    target_nd = NetworkData.from_dict(networks_data[target_network])
    reference_nd = NetworkData.from_dict(networks_data[reference_network])
    
    if target_nd.graph is None or reference_nd.targ_curve is None:
        return "Network data incomplete.", ""
    
    # Run normalizing target analysis
    finder = NormalizingTargetFinder(max_workers=4)
    
    # Prepare data
    target_data = {
        'graph': target_nd.graph,
        'targ_curve': target_nd.targ_curve
    }
    
    reference_curve = reference_nd.targ_curve
    
    try:
        results_df = finder.find_targets(
            target_data,
            reference_curve,
            strategy='degree',
            stop_threshold=0.05
        )
        
        # Get reference gene list for annotation
        reference_genes = list(reference_nd.graph.nodes()) if reference_nd.graph else []
        top_results = finder.get_top_targets(results_df, top_n=20, reference_genes=reference_genes)
        
        # Create results display
        if len(top_results) > 0:
            results_html = html.Div([
                html.H4(f"Top Normalizing Targets: {target_network} → {reference_network}", 
                       style={'color': '#2C3E50', 'marginBottom': '20px'}),
                
                # Summary
                html.Div([
                    html.Strong("Original Similarity: "),
                    html.Span(f"{results_df['similarity_score'].iloc[0]:.3f}" if len(results_df) > 0 else "N/A"),
                    html.Span(" | "),
                    html.Strong("Best Target Similarity: "),
                    html.Span(f"{top_results['similarity_score'].iloc[0]:.3f}" if len(top_results) > 0 else "N/A")
                ], style={'marginBottom': '15px', 'fontSize': '14px'}),
                
                # Table
                html.Div(
                    style={
                        'overflowX': 'auto',
                        'marginBottom': '20px'
                    },
                    children=[
                        html.Table(
                            [
                                html.Thead(html.Tr([
                                    html.Th('Rank'),
                                    html.Th('Gene'),
                                    html.Th('Similarity Score'),
                                    html.Th('Effect Size'),
                                    html.Th('Degree'),
                                    html.Th('In Reference')
                                ]))
                            ] + [
                                html.Tbody([
                                    html.Tr([
                                        html.Td(row['rank']),
                                        html.Td(row['gene']),
                                        html.Td(f"{row['similarity_score']:.3f}"),
                                        html.Td(f"{row['effect_size']:+.3f}"),
                                        html.Td(row['degree']),
                                        html.Td('✓' if row['in_reference'] else '✗')
                                    ])
                                    for _, row in top_results.iterrows()
                                ])
                            ],
                            style={'width': '100%', 'borderCollapse': 'collapse'}
                        )
                    ]
                ),
                
                # Export buttons
                html.Div([
                    html.Button(
                        '📥 Export to CSV',
                        id='export-normalizing-csv',
                        n_clicks=0,
                        style={
                            'backgroundColor': '#3498DB',
                            'color': 'white',
                            'border': 'none',
                            'padding': '10px 20px',
                            'borderRadius': '5px',
                            'cursor': 'pointer',
                            'marginRight': '10px'
                        }
                    ),
                ]),
                
                # Before/After plot for top target
                html.Div(
                    style={
                        'background': '#F8F9FA',
                        'padding': '15px',
                        'borderRadius': '8px',
                        'marginTop': '20px'
                    },
                    children=[
                        html.H4(f"Effect of Removing {top_results['gene'].iloc[0]}", 
                               style={'color': '#2C3E50'}),
                        dcc.Graph(
                            id='normalizing-before-after-plot',
                            figure=self._create_before_after_figure(
                                target_nd, reference_nd, top_results.iloc[0]
                            )
                        )
                    ]
                ) if len(top_results) > 0 else None,
                
                # Store results for export
                dcc.Store(id='normalizing-results-store', data=top_results.to_dict('records'))
            ])
        else:
            results_html = html.Div("No normalizing targets found.")
        
        status = html.Div(f"✓ Analysis complete. Found {len(top_results)} potential normalizing targets.", 
                         style={'color': '#27AE60', 'fontWeight': 'bold'})
        
        return status, results_html
        
    except Exception as e:
        logger.error(f"Error finding normalizing targets: {e}")
        return f"Error: {str(e)}", ""


# Export callbacks
@callback(
    Output('download-normalizing-csv', 'data'),
    Input('export-normalizing-csv', 'n_clicks'),
    State('normalizing-results-store', 'data'),
    prevent_initial_call=True
)
def export_normalizing_csv(n_clicks, results_data):
    """Export normalizing targets to CSV."""
    if n_clicks == 0 or not results_data:
        raise PreventUpdate
    
    import pandas as pd
    df = pd.DataFrame(results_data)
    return dcc.send_data_frame(df.to_csv, f"normalizing_targets_{time.time()}.csv", index=False)


@callback(
    Output('graph-network-selector', 'options'),
    Input('networks-store', 'data'),
    prevent_initial_call=True
)
def update_graph_selector(networks_data):
    """Update network selector in graphs tab."""
    if networks_data:
        return [{'label': name, 'value': name} for name in networks_data.keys()]
    return []


@callback(
    Output('network-graph-display', 'figure'),
    Input('graph-network-selector', 'value'),
    State('networks-store', 'data'),
    prevent_initial_call=True
)
def update_network_graph(selected_network, networks_data):
    """Update network graph display."""
    if not selected_network or not networks_data:
        return go.Figure()
    
    nd = NetworkData.from_dict(networks_data[selected_network])
    return self._create_network_figure(nd)


@callback(
    Output('graph-info', 'children'),
    Input('graph-network-selector', 'value'),
    State('networks-store', 'data'),
    prevent_initial_call=True
)
def update_graph_info(selected_network, networks_data):
    """Update graph info display."""
    if not selected_network or not networks_data:
        return ""
    
    nd = NetworkData.from_dict(networks_data[selected_network])
    return html.Div([
        html.Strong(f"{selected_network} Network: "),
        html.Span(f"{nd.metrics.get('nodes', 0)} nodes, {nd.metrics.get('edges', 0)} edges")
    ])


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _create_robustness_figure(nd: NetworkData) -> go.Figure:
    """Create a robustness curve figure for a single network."""
    fig = go.Figure()
    
    if nd.rand_curve:
        x, y = nd.rand_curve
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='lines+markers',
            name='Random Removal',
            line=dict(color='green', width=3),
            marker=dict(size=6),
            hovertemplate='Fraction Removed: %{x:.2f}<br>Integrity: %{y:.2f}<extra></extra>'
        ))
    
    if nd.targ_curve:
        x, y = nd.targ_curve
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='lines+markers',
            name='Targeted Removal',
            line=dict(color='red', width=3, dash='dash'),
            marker=dict(size=6),
            hovertemplate='Fraction Removed: %{x:.2f}<br>Integrity: %{y:.2f}<extra></extra>'
        ))
    
    fig.update_layout(
        title=f'Robustness Curves: {nd.factor_name}',
        xaxis_title='Fraction of Network Removed',
        yaxis_title='Network Integrity (LCC Size)',
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        plot_bgcolor='white',
        height=500,
        hovermode='closest'
    )
    
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    
    return fig


def _create_network_figure(nd: NetworkData) -> go.Figure:
    """Create an interactive network graph figure."""
    if nd.graph is None or nd.graph.number_of_nodes() == 0:
        return go.Figure()
    
    # Compute positions
    try:
        pos = nx.spring_layout(nd.graph, seed=42)
    except:
        pos = nx.random_layout(nd.graph)
    
    # Edges
    edge_x = []
    edge_y = []
    for edge in nd.graph.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        line=dict(width=0.5, color='grey'),
        hoverinfo='none',
        mode='lines'
    )
    
    # Nodes
    node_x = []
    node_y = []
    node_text = []
    node_degree = dict(nd.graph.degree())
    max_degree = max(node_degree.values()) if node_degree else 1
    node_size = [10 + 8 * (node_degree[node] / max_degree) for node in nd.graph.nodes()]
    node_color = [node_degree[node] for node in nd.graph.nodes()]
    
    for node in nd.graph.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        deg = node_degree[node]
        node_text.append(f'Node: {node}<br>Degree: {deg}')
    
    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers',
        hoverinfo='text',
        text=node_text,
        marker=dict(
            showscale=True,
            colorscale='YlGnBu',
            reversescale=True,
            color=node_color,
            size=node_size,
            colorbar=dict(
                thickness=15,
                title='Node Degree',
            ),
            line_width=1
        )
    )
    
    fig = go.Figure(data=[edge_trace, node_trace])
    
    fig.update_layout(
        title=f'Network Graph: {nd.factor_name}<br><sub>Size by degree, color by degree</sub>',
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='white',
        height=700
    )
    
    return fig


def _create_before_after_figure(target_nd: NetworkData, reference_nd: NetworkData, target_info: pd.Series) -> go.Figure:
    """Create before/after plot for normalizing target."""
    fig = go.Figure()
    
    # Original target curve
    if target_nd.targ_curve:
        x, y = target_nd.targ_curve
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='lines+markers',
            name=f'{target_nd.factor_name} (Original)',
            line=dict(color='red', width=3, dash='dash'),
            marker=dict(size=6)
        ))
    
    # Modified target curve (from target_info)
    if 'curve_x' in target_info and 'curve_y' in target_info:
        fig.add_trace(go.Scatter(
            x=target_info['curve_x'],
            y=target_info['curve_y'],
            mode='lines+markers',
            name=f'{target_nd.factor_name} - {target_info["gene"]}',
            line=dict(color='orange', width=3),
            marker=dict(size=6)
        ))
    
    # Reference curve
    if reference_nd.targ_curve:
        x, y = reference_nd.targ_curve
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='lines+markers',
            name=f'{reference_nd.factor_name} (Reference)',
            line=dict(color='green', width=3),
            marker=dict(size=6)
        ))
    
    fig.update_layout(
        title=f'Normalizing Effect: Removing {target_info["gene"]}',
        xaxis_title='Fraction of Network Removed',
        yaxis_title='Network Integrity (LCC Size)',
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        plot_bgcolor='white',
        height=500,
        hovermode='closest'
    )
    
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    
    return fig


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    logger.info("Starting NetFragile 2.0...")
    logger.info("Open http://localhost:8050 in your browser")
    app.run_server(debug=True, host='0.0.0.0', port=8050)
