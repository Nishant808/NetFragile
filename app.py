import networkx as nx
import plotly.graph_objects as go
import requests
import io
import pandas as pd
import random
import dash
from dash import dcc, html, Input, Output, State, callback
from dash.exceptions import PreventUpdate


class MycothiolRealWorldModel:
    def __init__(self):
        self.graph = None
        self.degree_cent = None
        self.betweenness_cent = None
        self.rand_res = None
        self.targ_res = None

    def fetch_real_data(self, species_id, seed_genes):
        """
        Fetches real PPI data from STRING Database API.
        """
        print(f"[Step 1] Fetching real data for species (TaxID: {species_id})...")
        print(f" Seeding network with genes: {', '.join(seed_genes)}")

        # STRING API Endpoint
        url = "https://string-db.org/api/tsv/network"

        params = {
            "identifiers": "%0d".join(seed_genes),  # Input genes
            "species": species_id,
            "add_nodes": 200,  # Expand network to find neighbors (build a cluster)
            "required_score": 700,  # High confidence interactions only (>0.7)
            "network_type": "functional"
        }

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()

            # Parse TSV data into a Graph
            data = response.text
            df = pd.read_csv(io.StringIO(data), sep="\t")

            self.graph = nx.Graph()
            # Build graph from edge list
            for _, row in df.iterrows():
                p1 = row['preferredName_A']
                p2 = row['preferredName_B']
                score = row['score']
                self.graph.add_edge(p1, p2, weight=score)

            # Keep only the largest connected component for stable analysis
            largest_cc = max(nx.connected_components(self.graph), key=len)
            self.graph = self.graph.subgraph(largest_cc).copy()

            print(
                f" [Success] Real Graph Built: {self.graph.number_of_nodes()} nodes, {self.graph.number_of_edges()} edges.")
            return True, f"Success: Graph built with {self.graph.number_of_nodes()} nodes and {self.graph.number_of_edges()} edges."

        except Exception as e:
            print(f" [Error] Could not fetch data: {e}")
            return False, f"Error: Could not fetch data - {str(e)}"

    def analyze_centrality(self):
        """
        Calculates centrality to identify Real-World Targets (Hubs & Bottlenecks).
        """
        if self.graph is None:
            return "Graph not built yet."

        print("\n[Analysis] Calculating Centrality Metrics...")

        # 1. Degree Centrality (Hubs)
        self.degree_cent = nx.degree_centrality(self.graph)

        # 2. Betweenness Centrality (Bottlenecks)
        self.betweenness_cent = nx.betweenness_centrality(self.graph)

        # Sort and find top targets
        top_hubs = sorted(self.degree_cent.items(), key=lambda x: x[1], reverse=True)[:5]
        top_bottlenecks = sorted(self.betweenness_cent.items(), key=lambda x: x[1], reverse=True)[:5]

        output = "\n >>> TOP DRUG TARGETS (Real Proteins):\n"
        output += f" Top 5 Hubs (high degree nodes, potential key players in the network): {[n for n, v in top_hubs]}\n"
        output += f" Top 5 Bottlenecks (high betweenness, control information flow): {[n for n, v in top_bottlenecks]}\n"
        output += "\nExplanation: Hubs are proteins with many interactions, making them good targets for broad disruption. \nBottlenecks are critical bridges in the network; targeting them can fragment the network efficiently."

        return output

    def _simulate_removal(self, node_order):
        """
        Simulates sequential removal of nodes and tracks LCC size.
        """
        G_temp = self.graph.copy()
        lcc_sizes = []
        n_removed = []
        initial_size = len(G_temp)

        for i, node in enumerate(node_order):
            if node in G_temp:
                G_temp.remove_node(node)

            if len(G_temp) > 0:
                # Calculate fraction of remaining network
                largest_cc = len(max(nx.connected_components(G_temp), key=len))
                lcc_sizes.append(largest_cc / initial_size)
            else:
                lcc_sizes.append(0)

            n_removed.append((i + 1) / initial_size)

            # Stop if network is effectively destroyed (<5% remaining)
            if lcc_sizes[-1] < 0.05:
                break

        return n_removed, lcc_sizes

    def robustness_test(self):
        """
        Compare Random Drug vs. Targeted Drug on the real network.
        """
        if self.graph is None:
            return "Graph not built yet."

        print("\n[Step 2] Simulating Drug Attacks on Network...")

        nodes = list(self.graph.nodes())

        # Strategy A: Random Attack
        random_order = nodes.copy()
        random.shuffle(random_order)
        rand_x, rand_y = self._simulate_removal(random_order)
        self.rand_res = (rand_x, rand_y)

        # Strategy B: Targeted Attack (Degree Centrality)
        degree_dict = dict(self.graph.degree())
        targeted_order = sorted(degree_dict, key=degree_dict.get, reverse=True)
        targ_x, targ_y = self._simulate_removal(targeted_order)
        self.targ_res = (targ_x, targ_y)

        output = "\nExplanation: The robustness test simulates drug effects. Random attacks mimic non-specific drugs, removing nodes randomly.\n Targeted attacks mimic specific drugs hitting high-degree hubs, disrupting the network faster."

        return output

    def get_robustness_figure(self):
        if self.rand_res is None or self.targ_res is None:
            return None

        rand_x, rand_y = self.rand_res
        targ_x, targ_y = self.targ_res

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=rand_x,
            y=rand_y,
            mode='lines+markers',
            name='Random Drug (Non-specific)',
            line=dict(color='green', width=3),
            marker=dict(size=6),
            hovertemplate='Fraction Removed: %{x:.2f}<br>Integrity: %{y:.2f}<extra></extra>'
        ))

        fig.add_trace(go.Scatter(
            x=targ_x,
            y=targ_y,
            mode='lines+markers',
            name='Targeted Drug (Anti-Hub)',
            line=dict(color='red', width=3, dash='dash'),
            marker=dict(size=6),
            hovertemplate='Fraction Removed: %{x:.2f}<br>Integrity: %{y:.2f}<extra></extra>'
        ))

        fig.update_layout(
            title={
                'text': f'Robustness of Sub-network<br>(Nodes: {self.graph.number_of_nodes()}, Edges: {self.graph.number_of_edges()})',
                'y': 0.95,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top',
                'font': dict(size=20)
            },
            xaxis_title='Fraction of Proteome Knocked Out',
            yaxis_title='Network Integrity (LCC Size)',
            legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01, bgcolor='rgba(255,255,255,0.8)'),
            plot_bgcolor='white',
            font=dict(size=14),
            hovermode='closest',
            height=600
        )

        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey', zeroline=True, zerolinewidth=2,
                         zerolinecolor='Black')
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey', zeroline=True, zerolinewidth=2,
                         zerolinecolor='Black')

        return fig

    def get_network_figure(self):
        if self.graph is None:
            return None

        print("\n[Visualization] Generating interactive network graph...")

        # Compute positions
        pos = nx.spring_layout(self.graph, seed=42)

        # Edges
        edge_x = []
        edge_y = []
        edge_hover = []
        for edge in self.graph.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            weight = self.graph.edges[edge]['weight']
            edge_hover.extend(
                [f'Edge: {edge[0]}-{edge[1]}<br>Weight: {weight:.2f}'] * 3)  # Repeat for each segment including None

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='grey'),
            hoverinfo='text',
            text=edge_hover,
            mode='lines'
        )

        # Nodes
        node_x = []
        node_y = []
        node_text = []
        node_degree = dict(self.graph.degree())
        max_degree = max(node_degree.values()) if node_degree else 1
        node_size = [15 + 10 * (node_degree[node] / max_degree) for node in self.graph.nodes()]
        node_color = [node_degree[node] for node in self.graph.nodes()]

        for node in self.graph.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            deg = node_degree[node]
            node_text.append(f'Node: {node}<br>Degree: {deg}')

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
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
                    title=dict(text='Node Degree', side='right'),
                    xanchor='left'
                ),
                line_width=2
            ),
            textfont=dict(size=8)
        )

        fig = go.Figure(data=[edge_trace, node_trace],
                        layout=go.Layout(
                            title={
                                'text': 'Interactive Network Graph<br>Hover over nodes/edges for details. Zoom and pan enabled.',
                                'y': 0.95,
                                'x': 0.5,
                                'xanchor': 'center',
                                'yanchor': 'top',
                                'font': dict(size=20)
                            },
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=20, l=5, r=5, t=40),
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            plot_bgcolor='white',
                            height=800
                        ))

        return fig


# Initialize the model
model = MycothiolRealWorldModel()

# Create Dash app
app = dash.Dash(__name__)

app.layout = html.Div(
    style={
        'fontFamily': 'Arial, sans-serif',
        'padding': '20px',
        'maxWidth': '1200px',
        'margin': 'auto',
        'background': 'linear-gradient(to bottom right, #e0eafc, #cfdef3)',
    },
    children=[
        html.Div(
            dcc.Markdown("""
                <style>
                    @keyframes fadeIn {
                        from { opacity: 0; transform: translateY(20px); }
                        to { opacity: 1; transform: translateY(0); }
                    }
                    .fade-in {
                        animation: fadeIn 0.8s ease-out;
                    }
                    .button {
                        transition: all 0.3s ease;
                        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                    }
                    .button:hover {
                        transform: scale(1.05);
                        box-shadow: 0 4px 8px rgba(0,0,0,0.2);
                    }
                    .graph-container {
                        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                        border-radius: 8px;
                        background: white;
                        padding: 10px;
                        margin-top: 30px;
                    }
                    h1 {
                        text-shadow: 1px 1px 2px rgba(0,0,0,0.1);
                    }
                    .input-container {
                        background: white;
                        border-radius: 10px;
                        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                        padding: 20px;
                    }
                </style>
            """, dangerously_allow_html=True)
        ),


        html.H1("NetFragile", style={'textAlign': 'center', 'color': '#2C3E50'}),
        html.H3("Simulate targeted drug attacks on biological networks using real STRING PPI data" , style={'textAlign': 'center', 'color': '#2C3E50'}),

        html.Div(
            className='input-container',
            children=[
                html.Label("Species ID (e.g., 83332 for M. tuberculosis H37Rv):",
                           style={'fontWeight': 'bold', 'fontSize': '16px'}),
                dcc.Input(id='species-id', type='number', value=83332,
                          style={'width': '100%', 'marginBottom': '15px', 'padding': '10px', 'borderRadius': '5px',
                                 'border': '1px solid #ddd'}),

                html.Label("Seed Genes (comma-separated, e.g., Rv0486, Rv1170, Rv2130c, Rv0819, Rv1082):",
                           style={'fontWeight': 'bold', 'fontSize': '16px'}),
                dcc.Input(id='seed-genes', type='text', value='Rv0486, Rv1170, Rv2130c, Rv0819, Rv1082',
                          style={'width': '100%', 'marginBottom': '15px', 'padding': '10px', 'borderRadius': '5px',
                                 'border': '1px solid #ddd'}),

                html.Button('Fetch Data and Analyze', id='submit-button', n_clicks=0,
                            className='button',
                            style={'backgroundColor': '#3498DB', 'color': 'white', 'border': 'none',
                                   'padding': '12px 24px', 'borderRadius': '5px', 'cursor': 'pointer', 'width': '100%',
                                   'fontSize': '16px'}),
            ]
        ),

        dcc.Loading(
            id='loading-output',
            type='circle',
            children=html.Div(id='output-container', className='fade-in',
                              style={'marginTop': '30px', 'background': 'white', 'padding': '20px',
                                     'borderRadius': '10px', 'boxShadow': '0 4px 8px rgba(0,0,0,0.1)'}),
        ),

        html.Div(
            className='graph-container fade-in',
            children=dcc.Graph(id='robustness-graph')
        ),

        html.Div(
            id='network-explanation',
            className='fade-in',
            style={'marginTop': '20px', 'backgroundColor': '#E9ECEF', 'padding': '15px', 'borderRadius': '8px',
                   'boxShadow': '0 4px 8px rgba(0,0,0,0.1)'},
            children="Explanation: This graph shows proteins as nodes (size and color by degree) and interactions as edges. Larger, darker nodes are hubs. Interactivity allows zooming, panning, and hovering for details."
        ),

        html.Div(
            className='graph-container fade-in',
            children=dcc.Graph(id='network-graph')
        ),
    ]
)


@callback(
    [Output('output-container', 'children'),
     Output('robustness-graph', 'figure'),
     Output('network-graph', 'figure')],
    [Input('submit-button', 'n_clicks')],
    [State('species-id', 'value'),
     State('seed-genes', 'value')]
)
def update_analysis(n_clicks, species_id, seed_genes_str):
    if n_clicks == 0:
        raise PreventUpdate

    if not species_id or not seed_genes_str:
        return "Please provide species ID and seed genes.", {}, {}

    seed_genes = [gene.strip() for gene in seed_genes_str.split(',')]

    success, fetch_message = model.fetch_real_data(species_id, seed_genes)

    if not success:
        return fetch_message, {}, {}

    centrality_output = model.analyze_centrality()
    robustness_output = model.robustness_test()

    output = html.Div([
        html.H3("Fetch Result", style={'color': '#2980B9'}),
        html.Pre(fetch_message, style={'backgroundColor': '#F8F9FA', 'padding': '10px', 'borderRadius': '5px'}),

        html.H3("Centrality Analysis", style={'color': '#2980B9', 'marginTop': '20px'}),
        html.Pre(centrality_output, style={'backgroundColor': '#F8F9FA', 'padding': '10px', 'borderRadius': '5px'}),

        html.H3("Robustness Test", style={'color': '#2980B9', 'marginTop': '20px'}),
        html.Pre(robustness_output, style={'backgroundColor': '#F8F9FA', 'padding': '10px', 'borderRadius': '5px'}),
    ])

    robustness_fig = model.get_robustness_figure() or {}
    network_fig = model.get_network_figure() or {}

    return output, robustness_fig, network_fig


if __name__ == "__main__":
    app.run(debug=True)