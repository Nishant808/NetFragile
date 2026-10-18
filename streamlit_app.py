"""
NetFragile 2.0 - Streamlit App for Multi-Network Comparative Analysis
Upload EdgeR DEG files and find normalizing targets across conditions
"""

import streamlit as st
import networkx as nx
import plotly.graph_objects as go
import requests
import io
import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
import tempfile
import os
import time

# Import utilities
from utils.file_parser import EdgeRParser, parse_deg_file
from utils.network_builder import NetworkBuilder, build_network_from_genes
from utils.robustness import RobustnessAnalyzer
from utils.comparison import NetworkComparator
from utils.normalizing import NormalizingTargetFinder

# Page configuration
st.set_page_config(
    page_title="NetFragile 2.0 - Multi-Network Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for modern look
st.markdown("""
    <style>
    .main {
        padding: 0rem 1rem;
    }
    .stAlert {
        margin-top: 1rem;
    }
    .block-container {
        padding-top: 1rem;
        padding-bottom: 1rem;
    }
    h1 {
        color: #2C3E50;
        padding-bottom: 0.5rem;
        border-bottom: 3px solid #667eea;
    }
    h2 {
        color: #2C3E50;
        margin-top: 1.5rem;
    }
    h3 {
        color: #4A5568;
    }
    .metric-card {
        background: linear-gradient(135deg, rgba(102, 126, 234, 0.1) 0%, rgba(118, 75, 162, 0.1) 100%);
        padding: 1.5rem;
        border-radius: 12px;
        border-left: 4px solid #667eea;
        box-shadow: 0 4px 12px rgba(102, 126, 234, 0.2);
        margin: 1rem 0;
    }
    .gene-tag {
        display: inline-block;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 0.35rem 0.75rem;
        border-radius: 12px;
        margin: 0.25rem;
        font-size: 0.875rem;
        box-shadow: 0 2px 4px rgba(102, 126, 234, 0.3);
    }
    .pathway-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1.5rem;
        border-radius: 12px;
        color: white;
        margin: 1rem 0;
        box-shadow: 0 4px 12px rgba(102, 126, 234, 0.3);
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 2rem;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 12px 24px;
        border-radius: 8px 8px 0 0;
        font-weight: 600;
    }
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
    }
    .stButton>button {
        border-radius: 8px;
        font-weight: 600;
        transition: all 0.3s ease;
        border: none;
        box-shadow: 0 4px 8px rgba(0,0,0,0.15);
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(0,0,0,0.2);
    }
    .upload-box {
        border: 2px dashed #667eea;
        border-radius: 12px;
        padding: 2rem;
        text-align: center;
        background: rgba(102, 126, 234, 0.05);
        margin: 1rem 0;
    }
    .result-box {
        background: rgba(255, 255, 255, 0.8);
        padding: 1.5rem;
        border-radius: 12px;
        border: 1px solid #E2E8F0;
        margin: 1rem 0;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }
    .success-box {
        background: linear-gradient(135deg, rgba(39, 174, 96, 0.1) 0%, rgba(34, 153, 84, 0.1) 100%);
        padding: 1.5rem;
        border-radius: 12px;
        border-left: 4px solid #27AE60;
        margin: 1rem 0;
        box-shadow: 0 4px 12px rgba(39, 174, 96, 0.2);
    }
    .warning-box {
        background: linear-gradient(135deg, rgba(231, 76, 60, 0.1) 0%, rgba(192, 57, 43, 0.1) 100%);
        padding: 1.5rem;
        border-radius: 12px;
        border-left: 4px solid #E74C3C;
        margin: 1rem 0;
        box-shadow: 0 4px 12px rgba(231, 76, 60, 0.2);
    }
    .network-card {
        background: rgba(255, 255, 255, 0.9);
        padding: 1.5rem;
        border-radius: 12px;
        border: 1px solid #E2E8F0;
        margin: 1rem 0;
        box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        transition: transform 0.3s ease, box-shadow 0.3s ease;
    }
    .network-card:hover {
        transform: translateY(-4px);
        box-shadow: 0 8px 16px rgba(0,0,0,0.1);
    }
    .target-highlight {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 8px;
        display: inline-block;
        margin: 0.25rem;
        box-shadow: 0 2px 4px rgba(102, 126, 234, 0.3);
    }
    </style>
""", unsafe_allow_html=True)


class NetworkData:
    """Container for network data"""
    def __init__(self, factor_name: str):
        self.factor_name = factor_name
        self.gene_list = []
        self.graph = None
        self.degree_cent = None
        self.betweenness_cent = None
        self.rand_curve = None
        self.targ_curve = None
        self.metrics = {}
        self.metadata = {}


class AppState:
    """Global state for the Streamlit app"""
    def __init__(self):
        self.networks = {}
        self.uploaded_files = {}
        self.analysis_complete = False


# Initialize state
if 'state' not in st.session_state:
    st.session_state.state = AppState()


def parse_edgeR_file(uploaded_file, factor_name: str = None) -> Optional[Dict[str, Any]]:
    """Parse uploaded EdgeR DEG file"""
    try:
        # Read CSV
        df = pd.read_csv(uploaded_file)
        
        # Auto-detect gene column
        parser = EdgeRParser()
        gene_col = parser._find_gene_column(df)
        
        if gene_col is None:
            gene_col = df.columns[0]
        
        # Extract genes
        gene_list = df[gene_col].astype(str).tolist()
        gene_list = [g for g in gene_list if g and g != 'nan' and g != 'NaN']
        
        # Remove duplicates
        gene_list = list(dict.fromkeys(gene_list))
        
        return {
            'factor_name': factor_name or uploaded_file.name.replace('.csv', '').replace('.tsv', ''),
            'gene_list': gene_list,
            'gene_count': len(gene_list),
            'columns': list(df.columns),
            'file_name': uploaded_file.name
        }
    except Exception as e:
        st.error(f"Error parsing {uploaded_file.name}: {str(e)}")
        return None


def build_network_from_genes(gene_list: List[str], factor_name: str, 
                              required_score: int = 700, add_nodes: int = 200) -> Optional[NetworkData]:
    """Build network from gene list"""
    try:
        builder = NetworkBuilder()
        graph, metadata = builder.build_network(
            gene_list,
            required_score=required_score,
            add_nodes=add_nodes
        )
        
        if graph.number_of_nodes() == 0:
            st.warning(f"No network could be built for {factor_name}. Try lowering the confidence score.")
            return None
        
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
        
        # Calculate metrics
        nd.metrics = {
            'nodes': graph.number_of_nodes(),
            'edges': graph.number_of_edges(),
            'density': nx.density(graph) if graph.number_of_nodes() > 1 else 0,
            'avg_degree': sum(dict(graph.degree()).values()) / graph.number_of_nodes() if graph.number_of_nodes() > 0 else 0,
        }
        
        # Calculate centrality
        if graph.number_of_nodes() > 0:
            nd.degree_cent = nx.degree_centrality(graph)
            nd.betweenness_cent = nx.betweenness_centrality(graph)
        
        return nd
        
    except Exception as e:
        st.error(f"Error building network for {factor_name}: {str(e)}")
        return None


def create_robustness_figure(networks_data: Dict[str, NetworkData]) -> go.Figure:
    """Create overlay robustness figure"""
    fig = go.Figure()
    
    colors = ['#667eea', '#764ba2', '#f093fb', '#4facfe', '#43e97b', '#fa709a']
    
    for i, (name, nd) in enumerate(networks_data.items()):
        if nd.targ_curve:
            x, y = nd.targ_curve
            fig.add_trace(go.Scatter(
                x=x, y=y,
                mode='lines+markers',
                name=name,
                line=dict(color=colors[i % len(colors)], width=3),
                marker=dict(size=6),
                hovertemplate=f'{name}<br>Fraction Removed: %{{x:.2%}}<br>Integrity: %{{y:.2%}}<extra></extra>'
            ))
    
    fig.update_layout(
        title={
            'text': 'Comparative Robustness: Targeted Removal',
            'y': 0.95,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': dict(size=20, color='#2C3E50')
        },
        xaxis_title='Fraction of Network Removed',
        yaxis_title='Network Integrity (LCC Size)',
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        plot_bgcolor='rgba(255,255,255,0.9)',
        paper_bgcolor='rgba(255,255,255,0.9)',
        font=dict(size=14, family='Arial'),
        hovermode='closest',
        height=600,
        margin=dict(t=80)
    )
    
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    
    return fig


def create_network_figure(nd: NetworkData) -> go.Figure:
    """Create network visualization"""
    if nd.graph is None or nd.graph.number_of_nodes() == 0:
        return go.Figure()
    
    try:
        pos = nx.spring_layout(nd.graph, seed=42, k=0.15, iterations=50)
    except:
        pos = nx.random_layout(nd.graph)
    
    # Edges
    edge_x, edge_y = [], []
    for edge in nd.graph.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='rgba(102, 126, 234, 0.3)'),
        hoverinfo='none',
        mode='lines'
    )
    
    # Nodes
    node_x, node_y, node_text = [], [], []
    node_degree = dict(nd.graph.degree())
    max_degree = max(node_degree.values()) if node_degree else 1
    node_size = [10 + 12 * (node_degree[node] / max_degree) for node in nd.graph.nodes()]
    node_color = [node_degree[node] for node in nd.graph.nodes()]
    
    for node in nd.graph.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        deg = node_degree[node]
        node_text.append(f'<b>{node}</b><br>Degree: {deg}')
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        text=node_text,
        marker=dict(
            showscale=True,
            colorscale='YlOrRd',
            color=node_color,
            size=node_size,
            colorbar=dict(thickness=15, title='Node Degree'),
            line_width=1,
            line_color='DarkSlateGrey'
        )
    )
    
    fig = go.Figure(data=[edge_trace, node_trace])
    
    fig.update_layout(
        title=f'Network: {nd.factor_name}',
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=5, r=5, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='rgba(255,255,255,0.9)',
        paper_bgcolor='rgba(255,255,255,0.9)',
        height=700
    )
    
    return fig


def create_normalizing_figure(target_nd: NetworkData, reference_nd: NetworkData, 
                              target_info: pd.Series) -> go.Figure:
    """Create before/after figure for normalizing target"""
    fig = go.Figure()
    
    # Original target
    if target_nd.targ_curve:
        x, y = target_nd.targ_curve
        fig.add_trace(go.Scatter(
            x=x, y=y,
            mode='lines+markers',
            name=f'{target_nd.factor_name} (Original)',
            line=dict(color='#E74C3C', width=3, dash='dash'),
            marker=dict(size=6)
        ))
    
    # Modified target
    if 'curve_x' in target_info and 'curve_y' in target_info:
        fig.add_trace(go.Scatter(
            x=target_info['curve_x'],
            y=target_info['curve_y'],
            mode='lines+markers',
            name=f'{target_nd.factor_name} - {target_info["gene"]}',
            line=dict(color='#F39C12', width=3),
            marker=dict(size=6)
        ))
    
    # Reference
    if reference_nd.targ_curve:
        x, y = reference_nd.targ_curve
        fig.add_trace(go.Scatter(
            x=x, y=y,
            mode='lines+markers',
            name=f'{reference_nd.factor_name} (Reference)',
            line=dict(color='#27AE60', width=3),
            marker=dict(size=6)
        ))
    
    fig.update_layout(
        title=f'Normalizing Effect: Removing {target_info["gene"]}',
        xaxis_title='Fraction of Network Removed',
        yaxis_title='Network Integrity (LCC Size)',
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        plot_bgcolor='rgba(255,255,255,0.9)',
        paper_bgcolor='rgba(255,255,255,0.9)',
        height=500,
        hovermode='closest'
    )
    
    return fig


def main():
    """Main Streamlit app"""
    
    # Title with emoji
    st.markdown("""
        <div style='text-align: center; padding: 2rem 0;'>
            <h1 style='color: #2C3E50; font-size: 48px; font-weight: 700;'>
                🧬 NetFragile 2.0
            </h1>
            <h3 style='color: #4A5568; font-size: 24px; font-weight: 500; margin-top: 10px;'>
                Multi-Network Comparative Analysis
            </h3>
            <p style='color: #718096; font-size: 16px; margin-top: 10px;'>
                Upload EdgeR DEG files → Build Networks → Find Normalizing Targets
            </p>
        </div>
    """, unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.markdown("## ⚙️ Settings")
        
        # Species (locked to human)
        st.markdown("### Species")
        species_id = st.selectbox(
            "Species",
            options=[("Human (Homo sapiens)", 9606)],
            index=0,
            disabled=True,
            help="Currently only human (TaxID: 9606) is supported"
        )
        species_id = 9606  # Extract the ID
        
        # STRING settings
        st.markdown("### STRING API Settings")
        required_score = st.slider(
            "Confidence Score",
            min_value=100,
            max_value=900,
            value=700,
            step=100,
            help="Minimum STRING confidence score (higher = fewer, more reliable interactions)"
        )
        add_nodes = st.slider(
            "Additional Nodes",
            min_value=0,
            max_value=500,
            value=200,
            step=50,
            help="Number of additional nodes to include from STRING database"
        )
        
        # Filter settings
        st.markdown("### DEG Filtering")
        apply_pvalue = st.checkbox("Filter by p-value", value=True)
        pvalue_threshold = st.number_input(
            "p-value threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.01,
            disabled=not apply_pvalue
        )
        apply_logfc = st.checkbox("Filter by |logFC|", value=True)
        logfc_threshold = st.number_input(
            "|logFC| threshold",
            min_value=0.0,
            max_value=10.0,
            value=1.0,
            step=0.5,
            disabled=not apply_logfc
        )
        
        # Info
        st.markdown("---")
        st.markdown("""
            <div style='text-align: center; padding: 1rem;'>
                <p style='color: #718096; font-size: 12px; margin: 0;'>
                    NetFragile 2.0
                </p>
                <p style='color: #A0AEC0; font-size: 10px; margin: 5px 0 0 0;'>
                    Built with Streamlit
                </p>
            </div>
        """, unsafe_allow_html=True)
    
    # Main tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📤 Upload", 
        "🔍 Individual Networks", 
        "📊 Comparison", 
        "🎯 Normalizing Targets",
        "🌐 Network Graphs"
    ])
    
    # TAB 1: Upload Files
    with tab1:
        st.markdown("## 📤 Upload EdgeR DEG Files")
        st.markdown("""
            Upload your EdgeR differential expression output CSV/TSV files.
            Each file represents a condition (e.g., Metastatic, Normal, Tumor).
            The tool will automatically parse gene lists and build networks.
        """)
        
        # File uploader
        uploaded_files = st.file_uploader(
            "Upload DEG Files",
            type=['csv', 'tsv', 'txt'],
            accept_multiple_files=True,
            help="Upload one or more EdgeR output files"
        )
        
        if uploaded_files:
            st.success(f"✅ {len(uploaded_files)} file(s) uploaded")
            
            # Parse files
            parsed_data = []
            for uploaded_file in uploaded_files:
                factor_name = st.text_input(
                    f"Factor name for {uploaded_file.name}",
                    value=uploaded_file.name.replace('.csv', '').replace('.tsv', '').replace('.txt', ''),
                    key=f"factor_{uploaded_file.name}"
                )
                
                parsed = parse_edgeR_file(uploaded_file, factor_name)
                if parsed:
                    parsed_data.append(parsed)
            
            if parsed_data:
                st.markdown("### 📋 Parsed Files Summary")
                for data in parsed_data:
                    st.markdown(f"""
                        <div class='network-card'>
                            <h4 style='color: #667eea; margin-top: 0;'>✅ {data['factor_name']}</h4>
                            <p style='color: #4A5568; margin: 5px 0;'>
                                <strong>File:</strong> {data['file_name']}<br>
                                <strong>Genes:</strong> {data['gene_count']} unique genes
                            </p>
                        </div>
                    """, unsafe_allow_html=True)
                
                # Build networks button
                if st.button("🚀 Build Networks & Run Analysis", 
                             type="primary",
                             use_container_width=True):
                    with st.spinner("Building networks... This may take a few minutes..."):
                        st.session_state.state.networks = {}
                        st.session_state.state.uploaded_files = {}
                        
                        progress_bar = st.progress(0)
                        
                        for i, data in enumerate(parsed_data):
                            progress = (i + 1) / len(parsed_data)
                            progress_bar.progress(progress)
                            
                            # Apply filters
                            gene_list = data['gene_list']
                            if apply_pvalue or apply_logfc:
                                # Re-parse with filters
                                try:
                                    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
                                    df = pd.read_csv(uploaded_files[i])
                                    gene_col = parsed_data[i]['columns'][0]
                                    
                                    if apply_pvalue:
                                        pval_col = next((c for c in df.columns if 'pvalue' in c.lower() or 'p.val' in c.lower()), None)
                                        if pval_col:
                                            df = df[df[pval_col] < pvalue_threshold]
                                    
                                    if apply_logfc:
                                        logfc_col = next((c for c in df.columns if 'logfc' in c.lower() or 'logFC' in c.lower()), None)
                                        if logfc_col:
                                            df = df[abs(df[logfc_col]) > logfc_threshold]
                                    
                                    gene_list = df[gene_col].astype(str).tolist()
                                    gene_list = [g for g in gene_list if g and g != 'nan']
                                    gene_list = list(dict.fromkeys(gene_list))
                                except:
                                    pass
                            
                            # Build network
                            nd = build_network_from_genes(
                                gene_list,
                                data['factor_name'],
                                required_score=required_score,
                                add_nodes=add_nodes
                            )
                            
                            if nd:
                                st.session_state.state.networks[data['factor_name']] = nd
                                st.session_state.state.uploaded_files[data['factor_name']] = data
                        
                        st.session_state.state.analysis_complete = True
                        st.success("✅ Analysis complete! Networks built successfully.")
                        st.rerun()
        else:
            st.info("👆 Upload your EdgeR DEG files to begin")
    
    # TAB 2: Individual Networks
    with tab2:
        st.markdown("## 🔍 Individual Network Analysis")
        
        if not st.session_state.state.analysis_complete or not st.session_state.state.networks:
            st.info("Please upload and analyze files first in the Upload tab.")
        else:
            network_names = list(st.session_state.state.networks.keys())
            selected_network = st.selectbox(
                "Select Network",
                options=network_names,
                key="individual_network_select"
            )
            
            if selected_network:
                nd = st.session_state.state.networks[selected_network]
                
                # Metrics
                st.markdown("### 📊 Network Metrics")
                cols = st.columns(4)
                cols[0].metric("Nodes", nd.metrics.get('nodes', 0))
                cols[1].metric("Edges", nd.metrics.get('edges', 0))
                cols[2].metric("Density", f"{nd.metrics.get('density', 0):.4f}")
                cols[3].metric("Avg Degree", f"{nd.metrics.get('avg_degree', 0):.2f}")
                
                # Centrality
                st.markdown("### 🎯 Top Hubs & Bottlenecks")
                if nd.degree_cent:
                    top_hubs = sorted(nd.degree_cent.items(), key=lambda x: x[1], reverse=True)[:10]
                    hub_genes = [n for n, v in top_hubs]
                    st.markdown("**Top Hubs (by Degree Centrality):**")
                    st.write(", ".join(f"`{g}`" for g in hub_genes[:10]))
                
                if nd.betweenness_cent:
                    top_bottlenecks = sorted(nd.betweenness_cent.items(), key=lambda x: x[1], reverse=True)[:10]
                    bottleneck_genes = [n for n, v in top_bottlenecks]
                    st.markdown("**Top Bottlenecks (by Betweenness Centrality):**")
                    st.write(", ".join(f"`{g}`" for g in bottleneck_genes[:10]))
                
                # Robustness curves
                st.markdown("### 📈 Robustness Curves")
                if nd.rand_curve and nd.targ_curve:
                    fig = go.Figure()
                    
                    x, y = nd.rand_curve
                    fig.add_trace(go.Scatter(x=x, y=y, mode='lines+markers', 
                                             name='Random Removal',
                                             line=dict(color='green', width=3)))
                    
                    x, y = nd.targ_curve
                    fig.add_trace(go.Scatter(x=x, y=y, mode='lines+markers',
                                             name='Targeted Removal',
                                             line=dict(color='red', width=3, dash='dash')))
                    
                    fig.update_layout(
                        title=f'Robustness: {selected_network}',
                        xaxis_title='Fraction Removed',
                        yaxis_title='Network Integrity',
                        height=500
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                # Network graph
                st.markdown("### 🌐 Network Visualization")
                network_fig = create_network_figure(nd)
                st.plotly_chart(network_fig, use_container_width=True)
    
    # TAB 3: Comparison
    with tab3:
        st.markdown("## 📊 Network Comparison")
        
        if not st.session_state.state.analysis_complete or len(st.session_state.state.networks) < 2:
            st.info("Please upload and analyze at least 2 networks to compare.")
        else:
            # Overlay plot
            st.markdown("### 📈 Overlay Robustness Curves")
            curve_type = st.radio("Curve Type", ["Targeted Removal", "Random Removal"], horizontal=True)
            
            fig = create_robustness_figure(st.session_state.state.networks)
            st.plotly_chart(fig, use_container_width=True)
            
            # Similarity matrix
            st.markdown("### 🔗 Network Similarity Matrix")
            networks_data = {
                name: {
                    'targ_curve': nd.targ_curve,
                    'rand_curve': nd.rand_curve,
                    'graph': nd.graph
                }
                for name, nd in st.session_state.state.networks.items()
            }
            
            comparator = NetworkComparator()
            sim_matrix = comparator.calculate_similarity_matrix(networks_data)
            st.dataframe(sim_matrix, use_container_width=True)
            
            # Metrics comparison
            st.markdown("### 📊 Network Metrics Comparison")
            metrics_data = []
            for name, nd in st.session_state.state.networks.items():
                metrics_data.append({
                    'Network': name,
                    'Nodes': nd.metrics.get('nodes', 0),
                    'Edges': nd.metrics.get('edges', 0),
                    'Density': f"{nd.metrics.get('density', 0):.4f}",
                    'Avg Degree': f"{nd.metrics.get('avg_degree', 0):.2f}"
                })
            st.dataframe(pd.DataFrame(metrics_data), use_container_width=True)
    
    # TAB 4: Normalizing Targets
    with tab4:
        st.markdown("## 🎯 Normalizing Target Discovery")
        st.markdown("""
            Find genes in a target network whose removal makes its robustness profile
            resemble a reference network. These are potential drug targets that can
            "normalize" pathological network behavior.
        """)
        
        if not st.session_state.state.analysis_complete or len(st.session_state.state.networks) < 2:
            st.info("Please upload and analyze at least 2 networks to find normalizing targets.")
        else:
            network_names = list(st.session_state.state.networks.keys())
            
            col1, col2 = st.columns(2)
            with col1:
                target_network = st.selectbox("Target Network", options=network_names, key="target_net")
            with col2:
                reference_network = st.selectbox("Reference Network", 
                                                 options=[n for n in network_names if n != target_network],
                                                 key="reference_net")
            
            if target_network and reference_network:
                target_nd = st.session_state.state.networks[target_network]
                reference_nd = st.session_state.state.networks[reference_network]
                
                if st.button("🔍 Find Normalizing Targets", type="primary", use_container_width=True):
                    with st.spinner("Analyzing genes... This may take a few minutes..."):
                        finder = NormalizingTargetFinder(max_workers=4)
                        
                        target_data = {
                            'graph': target_nd.graph,
                            'targ_curve': target_nd.targ_curve
                        }
                        reference_curve = reference_nd.targ_curve
                        
                        results_df = finder.find_targets(
                            target_data,
                            reference_curve,
                            strategy='degree',
                            stop_threshold=0.05,
                            min_network_size=5
                        )
                        
                        if len(results_df) > 0:
                            # Get reference genes
                            reference_genes = list(reference_nd.graph.nodes())
                            top_results = finder.get_top_targets(
                                results_df, top_n=20, reference_genes=reference_genes
                            )
                            
                            st.success(f"✅ Found {len(top_results)} potential normalizing targets")
                            
                            # Original similarity
                            original_sim = results_df['similarity_score'].iloc[0] if len(results_df) > 0 else 0
                            best_sim = top_results['similarity_score'].iloc[0] if len(top_results) > 0 else 0
                            
                            st.markdown(f"""
                                <div class='success-box'>
                                    <h4>📊 Results Summary</h4>
                                    <p><strong>Original Similarity:</strong> {original_sim:.3f}</p>
                                    <p><strong>Best Target Similarity:</strong> {best_sim:.3f}</p>
                                    <p><strong>Improvement:</strong> {best_sim - original_sim:+.3f}</p>
                                </div>
                            """, unsafe_allow_html=True)
                            
                            # Top targets table
                            st.markdown("### 🏆 Top Normalizing Targets")
                            
                            # Format dataframe
                            display_df = top_results.copy()
                            display_df['similarity_score'] = display_df['similarity_score'].round(3)
                            display_df['effect_size'] = display_df['effect_size'].round(3)
                            display_df['in_reference'] = display_df['in_reference'].apply(
                                lambda x: "✅" if x else "❌"
                            )
                            
                            st.dataframe(
                                display_df[['rank', 'gene', 'similarity_score', 'effect_size', 'degree', 'in_reference']],
                                use_container_width=True,
                                hide_index=True
                            )
                            
                            # Before/After plot for top target
                            if len(top_results) > 0:
                                st.markdown("### 📈 Normalizing Effect Visualization")
                                top_gene = top_results.iloc[0]
                                before_after_fig = create_normalizing_figure(
                                    target_nd, reference_nd, top_gene
                                )
                                st.plotly_chart(before_after_fig, use_container_width=True)
                                
                                st.markdown(f"""
                                    <div class='target-highlight'>
                                        Top Target: {top_gene['gene']} | 
                                        Similarity Score: {top_gene['similarity_score']:.3f} | 
                                        In Reference: {"✅" if top_gene['in_reference'] else "❌"}
                                    </div>
                                """, unsafe_allow_html=True)
                            
                            # Export
                            st.markdown("### 💾 Export Results")
                            csv = top_results.to_csv(index=False)
                            st.download_button(
                                label="📥 Download CSV",
                                data=csv,
                                file_name=f"normalizing_targets_{target_network}_to_{reference_network}.csv",
                                mime="text/csv"
                            )
                        else:
                            st.warning("No normalizing targets found. Try with different networks.")
    
    # TAB 5: Network Graphs
    with tab5:
        st.markdown("## 🌐 Network Visualizations")
        
        if not st.session_state.state.analysis_complete or not st.session_state.state.networks:
            st.info("Please upload and analyze files first in the Upload tab.")
        else:
            network_names = list(st.session_state.state.networks.keys())
            selected_network = st.selectbox(
                "Select Network",
                options=network_names,
                key="graph_network_select"
            )
            
            if selected_network:
                nd = st.session_state.state.networks[selected_network]
                network_fig = create_network_figure(nd)
                st.plotly_chart(network_fig, use_container_width=True)
                
                # Network info
                st.markdown(f"""
                    <div class='network-card'>
                        <h4>{selected_network}</h4>
                        <p><strong>Nodes:</strong> {nd.metrics.get('nodes', 0)} | 
                           <strong>Edges:</strong> {nd.metrics.get('edges', 0)} | 
                           <strong>Density:</strong> {nd.metrics.get('density', 0):.4f}</p>
                    </div>
                """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
