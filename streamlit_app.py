"""
NetFragile Streamlit App - Network Robustness Analyzer
Ingest CSV with gene sets and perform batch network analysis
"""

import streamlit as st
import networkx as nx
import plotly.graph_objects as go
import requests
import io
import pandas as pd
import random
from datetime import datetime
import json
from pathlib import Path
import tempfile
import os

# Page configuration
st.set_page_config(
    page_title="NetFragile - Network Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
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
        color: #1f77b4;
        padding-bottom: 0.5rem;
        border-bottom: 3px solid #1f77b4;
    }
    h2 {
        color: #2c3e50;
        margin-top: 1.5rem;
    }
    h3 {
        color: #34495e;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
    .gene-tag {
        display: inline-block;
        background-color: #e8f4f8;
        color: #2c3e50;
        padding: 0.25rem 0.5rem;
        border-radius: 0.25rem;
        margin: 0.25rem;
        font-size: 0.875rem;
        border: 1px solid #b8daff;
    }
    .pathway-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1.5rem;
        border-radius: 0.5rem;
        color: white;
        margin: 1rem 0;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 2rem;
    }
    </style>
""", unsafe_allow_html=True)


class NetworkAnalyzer:
    """Core network analysis functionality"""
    
    def __init__(self):
        self.graph = None
        self.degree_cent = None
        self.betweenness_cent = None
        self.closeness_cent = None
        self.rand_res = None
        self.targ_res = None
        self.nodes = None
        self.edges = None
    
    def fetch_string_data(self, species_id, seed_genes, required_score=700, add_nodes=200):
        """
        Fetches PPI data from STRING Database API.
        """
        if not seed_genes or len(seed_genes) == 0:
            return False, "No seed genes provided"
        
        # STRING API Endpoint
        url = "https://string-db.org/api/tsv/network"
        
        params = {
            "identifiers": "%0d".join(seed_genes),
            "species": species_id,
            "add_nodes": add_nodes,
            "required_score": required_score,
            "network_type": "functional"
        }
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            if response.text.strip() == "":
                return False, f"No data found for species {species_id} with genes: {', '.join(seed_genes[:5])}..."
            
            # Parse TSV data
            data = response.text
            df = pd.read_csv(io.StringIO(data), sep="\t")
            
            if df.empty:
                return False, "Empty response from STRING API"
            
            # Build graph from edge list
            self.graph = nx.Graph()
            for _, row in df.iterrows():
                p1 = row['preferredName_A']
                p2 = row['preferredName_B']
                score = row['score']
                self.graph.add_edge(p1, p2, weight=score)
            
            # Keep only the largest connected component
            if self.graph.number_of_nodes() > 0:
                largest_cc = max(nx.connected_components(self.graph), key=len)
                self.graph = self.graph.subgraph(largest_cc).copy()
            
            self.nodes = list(self.graph.nodes())
            self.edges = list(self.graph.edges())
            
            return True, f"Success: {self.graph.number_of_nodes()} nodes, {self.graph.number_of_edges()} edges"
            
        except requests.exceptions.Timeout:
            return False, "STRING API request timed out. Please try again."
        except requests.exceptions.RequestException as e:
            return False, f"STRING API error: {str(e)}"
        except Exception as e:
            return False, f"Error processing data: {str(e)}"
    
    def analyze_centrality(self):
        """Calculate centrality metrics"""
        if self.graph is None or self.graph.number_of_nodes() == 0:
            return {}
        
        results = {}
        
        # Degree Centrality (Hubs)
        self.degree_cent = nx.degree_centrality(self.graph)
        results['degree_centrality'] = dict(sorted(self.degree_cent.items(), key=lambda x: x[1], reverse=True))
        
        # Betweenness Centrality (Bottlenecks)
        self.betweenness_cent = nx.betweenness_centrality(self.graph)
        results['betweenness_centrality'] = dict(sorted(self.betweenness_cent.items(), key=lambda x: x[1], reverse=True))
        
        # Closeness Centrality
        self.closeness_cent = nx.closeness_centrality(self.graph)
        results['closeness_centrality'] = dict(sorted(self.closeness_cent.items(), key=lambda x: x[1], reverse=True))
        
        # Additional metrics
        results['node_count'] = self.graph.number_of_nodes()
        results['edge_count'] = self.graph.number_of_edges()
        results['density'] = nx.density(self.graph)
        results['average_clustering'] = nx.average_clustering_coefficient(self.graph)
        
        # Top targets
        top_n = min(10, len(self.degree_cent))
        results['top_hubs'] = [n for n, v in sorted(self.degree_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]]
        results['top_bottlenecks'] = [n for n, v in sorted(self.betweenness_cent.items(), key=lambda x: x[1], reverse=True)[:top_n]]
        
        return results
    
    def simulate_removal(self, node_order):
        """Simulates sequential removal of nodes and tracks LCC size."""
        if self.graph is None:
            return [], []
        
        G_temp = self.graph.copy()
        lcc_sizes = []
        n_removed = []
        initial_size = len(G_temp)
        
        for i, node in enumerate(node_order):
            if node in G_temp:
                G_temp.remove_node(node)
            
            if len(G_temp) > 0:
                largest_cc = len(max(nx.connected_components(G_temp), key=len))
                lcc_sizes.append(largest_cc / initial_size)
            else:
                lcc_sizes.append(0)
            
            n_removed.append((i + 1) / initial_size)
            
            # Stop if network is effectively destroyed
            if lcc_sizes[-1] < 0.05:
                break
        
        return n_removed, lcc_sizes
    
    def robustness_test(self):
        """Compare Random vs. Targeted attack"""
        if self.graph is None or self.graph.number_of_nodes() == 0:
            return None, None, {}
        
        nodes = list(self.graph.nodes())
        
        # Random Attack
        random_order = nodes.copy()
        random.shuffle(random_order)
        rand_x, rand_y = self.simulate_removal(random_order)
        self.rand_res = (rand_x, rand_y)
        
        # Targeted Attack (Degree Centrality)
        degree_dict = dict(self.graph.degree())
        targeted_order = sorted(degree_dict, key=degree_dict.get, reverse=True)
        targ_x, targ_y = self.simulate_removal(targeted_order)
        self.targ_res = (targ_x, targ_y)
        
        return rand_x, rand_y, targ_x, targ_y
    
    def get_robustness_figure(self, pathway_name=""):
        """Generate robustness plot"""
        if self.rand_res is None or self.targ_res is None:
            return None
        
        rand_x, rand_y = self.rand_res
        targ_x, targ_y = self.targ_res
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=rand_x, y=rand_y,
            mode='lines+markers',
            name='Random Removal',
            line=dict(color='green', width=3),
            marker=dict(size=6),
            hovertemplate='Fraction Removed: %{x:.2%}<br>Integrity: %{y:.2%}<extra></extra>'
        ))
        
        fig.add_trace(go.Scatter(
            x=targ_x, y=targ_y,
            mode='lines+markers',
            name='Targeted (Hubs)',
            line=dict(color='red', width=3, dash='dash'),
            marker=dict(size=6),
            hovertemplate='Fraction Removed: %{x:.2%}<br>Integrity: %{y:.2%}<extra></extra>'
        ))
        
        title = f'Network Robustness' + (f' - {pathway_name}' if pathway_name else '')
        
        fig.update_layout(
            title={
                'text': title,
                'y': 0.95,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top',
                'font': dict(size=18)
            },
            xaxis_title='Fraction of Nodes Removed',
            yaxis_title='Network Integrity (LCC Size)',
            legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01, bgcolor='rgba(255,255,255,0.8)'),
            plot_bgcolor='white',
            font=dict(size=14),
            hovermode='closest',
            height=500,
            margin=dict(t=80)
        )
        
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
        
        return fig
    
    def get_network_figure(self, pathway_name=""):
        """Generate interactive network visualization"""
        if self.graph is None or self.graph.number_of_nodes() == 0:
            return None
        
        # Compute positions
        pos = nx.spring_layout(self.graph, seed=42, k=0.15, iterations=50)
        
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
            edge_hover.extend([f'Edge: {edge[0]}-{edge[1]}<br>Weight: {weight:.0f}'] * 3)
        
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='grey', opacity=0.5),
            hoverinfo='text',
            text=edge_hover,
            mode='lines',
            showlegend=False
        )
        
        # Nodes
        node_x = []
        node_y = []
        node_text = []
        node_degree = dict(self.graph.degree())
        max_degree = max(node_degree.values()) if node_degree else 1
        node_size = [10 + 15 * (node_degree[node] / max_degree) for node in self.graph.nodes()]
        node_color = [node_degree[node] for node in self.graph.nodes()]
        
        for node in self.graph.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            deg = node_degree[node]
            node_text.append(f'<b>{node}</b><br>Degree: {deg}')
        
        # Highlight top hubs
        top_hubs = sorted(node_degree.items(), key=lambda x: x[1], reverse=True)[:5]
        top_hub_names = [n for n, d in top_hubs]
        node_colors = ['#ff7f0e' if node in top_hub_names else node_color[i] for i, node in enumerate(self.graph.nodes())]
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            text=node_text,
            marker=dict(
                showscale=True,
                colorscale='YlOrRd',
                reversescale=False,
                color=node_colors,
                size=node_size,
                colorbar=dict(
                    thickness=15,
                    title=dict(text='Node Degree', side='right'),
                ),
                line_width=1,
                line_color='DarkSlateGrey',
                opacity=0.8
            ),
            name='Nodes'
        )
        
        fig = go.Figure(data=[edge_trace, node_trace])
        
        title = f'Network Visualization' + (f' - {pathway_name}' if pathway_name else '')
        
        fig.update_layout(
            title={
                'text': title,
                'y': 0.95,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top',
                'font': dict(size=18)
            },
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=80),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            height=600
        )
        
        return fig
    
    def get_summary_stats(self):
        """Get summary statistics"""
        if self.graph is None:
            return {}
        
        return {
            'Nodes': self.graph.number_of_nodes(),
            'Edges': self.graph.number_of_edges(),
            'Density': f"{nx.density(self.graph):.4f}",
            'Avg Clustering': f"{nx.average_clustering_coefficient(self.graph):.4f}",
            'Is Connected': nx.is_connected(self.graph),
            'Components': nx.number_connected_components(self.graph),
            'Diameter': nx.diameter(self.graph) if nx.is_connected(self.graph) else 'N/A',
            'Avg Path Length': f"{nx.average_shortest_path_length(self.graph):.4f}" if nx.is_connected(self.graph) else 'N/A'
        }


def parse_csv_genes(csv_file, gene_column='leadingEdge'):
    """Parse CSV file and extract gene sets from specified column"""
    try:
        df = pd.read_csv(csv_file)
        
        if gene_column not in df.columns:
            st.error(f"Column '{gene_column}' not found in CSV. Available columns: {list(df.columns)}")
            return None
        
        # Extract pathways and their genes
        results = []
        for idx, row in df.iterrows():
            pathway = row.get('pathway', f'Pathway_{idx}')
            genes_str = str(row[gene_column])
            
            # Parse genes (semicolon-separated)
            genes = [g.strip() for g in genes_str.split(';') if g.strip()]
            
            if genes:
                results.append({
                    'pathway': pathway,
                    'genes': genes,
                    'gene_count': len(genes),
                    'pval': row.get('pval', 'N/A'),
                    'padj': row.get('padj', 'N/A'),
                    'NES': row.get('NES', 'N/A'),
                    'size': row.get('size', len(genes))
                })
        
        return pd.DataFrame(results)
        
    except Exception as e:
        st.error(f"Error parsing CSV: {str(e)}")
        return None


def main():
    """Main Streamlit app"""
    
    # Title
    st.title("🧬 NetFragile - Network Robustness Analyzer")
    st.markdown("""
        **Upload a CSV with gene sets and analyze their network robustness**
        
        This tool processes pathway/gene set data from GSEA or similar analyses,
        builds protein-protein interaction networks using STRING database,
        and evaluates network robustness through targeted and random attacks.
    """)
    
    # Sidebar
    with st.sidebar:
        st.header("⚙️ Settings")
        
        # Species selection
        species_id = st.selectbox(
            "Species (TaxID)",
            [83332, 9606, 10090, 10116],
            index=1,
            help="""
            Select species TaxID:
            - 9606: Homo sapiens (Human)
            - 10090: Mus musculus (Mouse)
            - 10116: Rattus norvegicus (Rat)
            - 83332: Mycobacterium tuberculosis
            """
        )
        
        species_names = {
            83332: "Mycobacterium tuberculosis",
            9606: "Homo sapiens (Human)",
            10090: "Mus musculus (Mouse)",
            10116: "Rattus norvegicus (Rat)"
        }
        
        if species_id in species_names:
            st.write(f"**Selected:** {species_names[species_id]}")
        
        st.divider()
        
        # STRING API settings
        st.subheader("STRING API Settings")
        required_score = st.slider(
            "Minimum Interaction Score",
            0, 900, 700, 50,
            help="Higher scores = higher confidence interactions"
        )
        
        add_nodes = st.slider(
            "Additional Nodes",
            0, 500, 200, 50,
            help="Number of additional nodes to include in network expansion"
        )
        
        st.divider()
        
        # Analysis settings
        st.subheader("Analysis Settings")
        max_pathways = st.number_input(
            "Max Pathways to Process",
            min_value=1, max_value=50, value=10,
            help="Limit number of pathways for batch processing"
        )
        
        st.divider()
        
        # Info
        st.markdown("""
        ### 📖 About
        
        **NetFragile** analyzes network robustness by:
        
        1. ✅ Building PPI networks from STRING database
        2. ✅ Identifying hubs and bottlenecks
        3. ✅ Simulating random vs. targeted attacks
        4. ✅ Visualizing network integrity curves
        
        **Citation:** Inspired by network robustness studies in systems biology
        
        ---
        
        **Version:** 2.0 (Streamlit Edition)
        """)
    
    # Main content
    tab1, tab2, tab3, tab4 = st.tabs(["📁 Upload & Analyze", "📊 Results", "💾 Export", "ℹ️ Help"])
    
    # Initialize session state
    if 'analyzer' not in st.session_state:
        st.session_state.analyzer = NetworkAnalyzer()
    if 'pathway_results' not in st.session_state:
        st.session_state.pathway_results = []
    if 'current_pathway' not in st.session_state:
        st.session_state.current_pathway = None
    if 'analysis_complete' not in st.session_state:
        st.session_state.analysis_complete = False
    if 'uploaded_file' not in st.session_state:
        st.session_state.uploaded_file = None
    
    with tab1:
        st.header("Upload CSV and Analyze")
        
        # File upload
        uploaded_file = st.file_uploader(
            "Upload CSV File",
            type=['csv'],
            help="Upload a CSV file containing pathway/gene set data (e.g., GSEA results)",
            key='csv_uploader'
        )
        
        if uploaded_file:
            st.session_state.uploaded_file = uploaded_file
            
            # Parse CSV
            with st.spinner("Parsing CSV file..."):
                pathway_df = parse_csv_genes(uploaded_file)
            
            if pathway_df is not None and len(pathway_df) > 0:
                st.success(f"✅ CSV parsed successfully! Found {len(pathway_df)} pathways")
                
                # Show preview
                with st.expander("📋 Preview Pathway Data"):
                    st.dataframe(pathway_df.head(), use_container_width=True)
                    st.write(f"Total pathways: {len(pathway_df)}")
                
                # Select pathways to analyze
                st.subheader("Select Pathways to Analyze")
                
                selected_pathways = st.multiselect(
                    "Choose pathways (or leave empty to analyze top pathways)",
                    pathway_df['pathway'].tolist(),
                    default=pathway_df['pathway'].head(min(5, len(pathway_df))).tolist(),
                    help="Select specific pathways to analyze, or leave empty to analyze the top pathways"
                )
                
                # If no selection, use top N
                if not selected_pathways:
                    selected_pathways = pathway_df['pathway'].head(min(max_pathways, len(pathway_df))).tolist()
                    st.info(f"Analyzing top {len(selected_pathways)} pathways")
                else:
                    # Limit to max_pathways
                    selected_pathways = selected_pathways[:max_pathways]
                    st.info(f"Analyzing {len(selected_pathways)} selected pathways (limited to {max_pathways})")
                
                # Filter dataframe
                analysis_df = pathway_df[pathway_df['pathway'].isin(selected_pathways)]
                
                # Start analysis button
                if st.button("🚀 Start Network Analysis", type="primary", use_container_width=True):
                    st.session_state.analysis_complete = False
                    st.session_state.pathway_results = []
                    st.session_state.current_pathway = None
                    
                    # Initialize progress bar
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    total = len(analysis_df)
                    
                    for i, (idx, row) in enumerate(analysis_df.iterrows()):
                        pathway_name = row['pathway']
                        genes = row['genes']
                        
                        # Update progress
                        progress = (i + 1) / total
                        progress_bar.progress(progress)
                        status_text.text(f"Processing {pathway_name} ({i+1}/{total})...")
                        
                        try:
                            # Analyze this pathway
                            analyzer = NetworkAnalyzer()
                            
                            success, message = analyzer.fetch_string_data(
                                species_id, 
                                genes, 
                                required_score=required_score,
                                add_nodes=add_nodes
                            )
                            
                            if success:
                                # Calculate centrality
                                centrality_results = analyzer.analyze_centrality()
                                
                                # Run robustness test
                                rand_x, rand_y, targ_x, targ_y = analyzer.robustness_test()
                                
                                # Get figures
                                robustness_fig = analyzer.get_robustness_figure(pathway_name)
                                network_fig = analyzer.get_network_figure(pathway_name)
                                
                                # Get summary stats
                                summary_stats = analyzer.get_summary_stats()
                                
                                # Store results
                                result = {
                                    'pathway': pathway_name,
                                    'genes': genes,
                                    'gene_count': len(genes),
                                    'network_nodes': summary_stats.get('Nodes', 0),
                                    'network_edges': summary_stats.get('Edges', 0),
                                    'message': message,
                                    'centrality': centrality_results,
                                    'robustness_fig': robustness_fig,
                                    'network_fig': network_fig,
                                    'summary_stats': summary_stats,
                                    'analyzer': analyzer,
                                    'pval': row.get('pval', 'N/A'),
                                    'padj': row.get('padj', 'N/A'),
                                    'NES': row.get('NES', 'N/A')
                                }
                                
                                st.session_state.pathway_results.append(result)
                            else:
                                st.warning(f"⚠️ {pathway_name}: {message}")
                                
                        except Exception as e:
                            st.error(f"❌ Error processing {pathway_name}: {str(e)}")
                    
                    st.session_state.analysis_complete = True
                    status_text.text(f"✅ Analysis complete! Processed {len(st.session_state.pathway_results)} pathways")
                    progress_bar.empty()
                    
                    # Show summary
                    st.balloons()
                    st.success(f"Analysis complete! Successfully processed {len(st.session_state.pathway_results)} pathways")
                    
                    # Quick summary
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        total_nodes = sum(r['network_nodes'] for r in st.session_state.pathway_results)
                        st.metric("Total Nodes", f"{total_nodes:,}")
                    with col2:
                        total_edges = sum(r['network_edges'] for r in st.session_state.pathway_results)
                        st.metric("Total Edges", f"{total_edges:,}")
                    with col3:
                        avg_nodes = total_nodes / len(st.session_state.pathway_results) if st.session_state.pathway_results else 0
                        st.metric("Avg Nodes/Pathway", f"{avg_nodes:.1f}")
    
    with tab2:
        st.header("📊 Analysis Results")
        
        if not st.session_state.analysis_complete:
            st.info("⬆️ Please upload a CSV and start analysis in the first tab")
        else:
            if not st.session_state.pathway_results:
                st.warning("No results to display. Please run analysis first.")
            else:
                # Sort by network size
                sorted_results = sorted(
                    st.session_state.pathway_results,
                    key=lambda x: x['network_nodes'],
                    reverse=True
                )
                
                # Summary table
                st.subheader("📈 Summary Table")
                
                summary_data = []
                for r in sorted_results:
                    summary_data.append({
                        'Pathway': r['pathway'],
                        'Genes': r['gene_count'],
                        'Nodes': r['network_nodes'],
                        'Edges': r['network_edges'],
                        'Density': r['summary_stats'].get('Density', 'N/A'),
                        'Top Hub': r['centrality'].get('top_hubs', [])[0] if r['centrality'].get('top_hubs') else 'N/A',
                        'p-value': r.get('pval', 'N/A'),
                        'NES': r.get('NES', 'N/A')
                    })
                
                summary_df = pd.DataFrame(summary_data)
                st.dataframe(summary_df, use_container_width=True, hide_index=True)
                
                st.divider()
                
                # Individual pathway analysis
                st.subheader("🔍 Detailed Pathway Analysis")
                
                # Pathway selector
                pathway_options = [r['pathway'] for r in sorted_results]
                selected_pathway = st.selectbox(
                    "Select a pathway to view details",
                    pathway_options,
                    key='pathway_selector'
                )
                
                # Find selected pathway
                selected_result = next(
                    (r for r in sorted_results if r['pathway'] == selected_pathway),
                    None
                )
                
                if selected_result:
                    st.session_state.current_pathway = selected_result
                    
                    # Display pathway info
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("Pathway", selected_result['pathway'])
                    with col2:
                        st.metric("Input Genes", selected_result['gene_count'])
                    with col3:
                        st.metric("Network Nodes", selected_result['network_nodes'])
                    with col4:
                        st.metric("Network Edges", selected_result['network_edges'])
                    
                    # Additional metrics
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("Density", selected_result['summary_stats'].get('Density', 'N/A'))
                    with col2:
                        st.metric("Avg Clustering", selected_result['summary_stats'].get('Avg Clustering', 'N/A'))
                    with col3:
                        st.metric("Is Connected", selected_result['summary_stats'].get('Is Connected', 'N/A'))
                    with col4:
                        st.metric("Components", selected_result['summary_stats'].get('Components', 'N/A'))
                    
                    st.divider()
                    
                    # Top hubs and bottlenecks
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.subheader("🎯 Top Hubs (Degree Centrality)")
                        top_hubs = selected_result['centrality'].get('top_hubs', [])
                        for i, hub in enumerate(top_hubs[:10], 1):
                            degree = selected_result['analyzer'].degree_cent.get(hub, 0)
                            st.markdown(f"**{i}. {hub}** (Degree: {degree:.4f})")
                    
                    with col2:
                        st.subheader("🔄 Top Bottlenecks (Betweenness)")
                        top_bottlenecks = selected_result['centrality'].get('top_bottlenecks', [])
                        for i, bn in enumerate(top_bottlenecks[:10], 1):
                            betweenness = selected_result['analyzer'].betweenness_cent.get(bn, 0)
                            st.markdown(f"**{i}. {bn}** (Betweenness: {betweenness:.4f})")
                    
                    st.divider()
                    
                    # Visualizations
                    st.subheader("📊 Visualizations")
                    
                    # Robustness curve
                    if selected_result['robustness_fig']:
                        st.plotly_chart(
                            selected_result['robustness_fig'],
                            use_container_width=True,
                            key=f"robustness_{selected_pathway}"
                        )
                    
                    # Network graph
                    if selected_result['network_fig']:
                        st.plotly_chart(
                            selected_result['network_fig'],
                            use_container_width=True,
                            key=f"network_{selected_pathway}"
                        )
                    
                    # Gene list
                    with st.expander("🧬 Input Genes"):
                        genes_str = ", ".join(selected_result['genes'][:50])
                        if len(selected_result['genes']) > 50:
                            genes_str += f" ... (+{len(selected_result['genes']) - 50} more)"
                        st.write(genes_str)
    
    with tab3:
        st.header("💾 Export Data")
        
        if not st.session_state.analysis_complete:
            st.info("⬆️ Please run analysis first")
        else:
            st.subheader("Export Analysis Results")
            
            # Export options
            export_format = st.radio(
                "Select export format",
                ["CSV", "JSON", "Excel"],
                horizontal=True
            )
            
            # Prepare export data
            export_data = []
            for r in st.session_state.pathway_results:
                row = {
                    'pathway': r['pathway'],
                    'input_genes': len(r['genes']),
                    'network_nodes': r['network_nodes'],
                    'network_edges': r['network_edges'],
                    'density': r['summary_stats'].get('Density', 'N/A'),
                    'avg_clustering': r['summary_stats'].get('Avg Clustering', 'N/A'),
                    'is_connected': r['summary_stats'].get('Is Connected', 'N/A'),
                    'components': r['summary_stats'].get('Components', 'N/A'),
                    'diameter': r['summary_stats'].get('Diameter', 'N/A'),
                    'avg_path_length': r['summary_stats'].get('Avg Path Length', 'N/A'),
                    'pval': r.get('pval', 'N/A'),
                    'padj': r.get('padj', 'N/A'),
                    'NES': r.get('NES', 'N/A'),
                    'top_hubs': ";".join(r['centrality'].get('top_hubs', [])[:5]),
                    'top_bottlenecks': ";".join(r['centrality'].get('top_bottlenecks', [])[:5]),
                    'message': r['message']
                }
                export_data.append(row)
            
            export_df = pd.DataFrame(export_data)
            
            # Export buttons
            col1, col2, col3 = st.columns(3)
            
            with col1:
                csv_data = export_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download CSV",
                    data=csv_data,
                    file_name=f"netfragile_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            
            with col2:
                json_data = export_df.to_json(orient='records', indent=2)
                st.download_button(
                    label="📥 Download JSON",
                    data=json_data,
                    file_name=f"netfragile_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                    mime="application/json",
                    use_container_width=True
                )
            
            with col3:
                try:
                    excel_data = io.BytesIO()
                    export_df.to_excel(excel_data, index=False, engine='openpyxl')
                    st.download_button(
                        label="📥 Download Excel",
                        data=excel_data.getvalue(),
                        file_name=f"netfragile_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True
                    )
                except ImportError:
                    st.warning("Excel export requires openpyxl. Install with: pip install openpyxl")
            
            st.divider()
            
            # Export individual pathway data
            st.subheader("Export Individual Pathway Data")
            
            if st.session_state.current_pathway:
                current = st.session_state.current_pathway
                
                # Centrality data
                cent_data = {
                    'node': [],
                    'degree_centrality': [],
                    'betweenness_centrality': [],
                    'closeness_centrality': []
                }
                
                for node in current['analyzer'].graph.nodes():
                    cent_data['node'].append(node)
                    cent_data['degree_centrality'].append(
                        current['analyzer'].degree_cent.get(node, 0)
                    )
                    cent_data['betweenness_centrality'].append(
                        current['analyzer'].betweenness_cent.get(node, 0)
                    )
                    cent_data['closeness_centrality'].append(
                        current['analyzer'].closeness_cent.get(node, 0)
                    )
                
                cent_df = pd.DataFrame(cent_data)
                
                col1, col2 = st.columns(2)
                
                with col1:
                    csv_cent = cent_df.to_csv(index=False)
                    st.download_button(
                        label=f"📥 Download {current['pathway']} Centrality Data (CSV)",
                        data=csv_cent,
                        file_name=f"{current['pathway']}_centrality.csv",
                        mime="text/csv",
                        use_container_width=True
                    )
                
                with col2:
                    json_cent = cent_df.to_json(orient='records', indent=2)
                    st.download_button(
                        label=f"📥 Download {current['pathway']} Centrality Data (JSON)",
                        data=json_cent,
                        file_name=f"{current['pathway']}_centrality.json",
                        mime="application/json",
                        use_container_width=True
                    )
    
    with tab4:
        st.header("ℹ️ Help & Documentation")
        
        st.markdown("""
        ## 📖 About NetFragile
        
        **NetFragile** is a network robustness analysis tool that helps you understand
        how biological networks respond to targeted disruptions.
        
        ### 🎯 Key Features
        
        - **CSV Ingestion**: Upload GSEA or similar pathway analysis results
        - **Network Construction**: Builds PPI networks from STRING database
        - **Centrality Analysis**: Identifies hubs and bottlenecks
        - **Robustness Testing**: Simulates random vs. targeted attacks
        - **Interactive Visualization**: Explore networks and robustness curves
        - **Data Export**: Download results in multiple formats
        
        ### 📥 Input CSV Format
        
        Your CSV should contain at least:
        - A `pathway` column with pathway names
        - A `leadingEdge` column with semicolon-separated gene lists
        
        Example columns that will be used:
        - `pathway`: Pathway name (e.g., "HALLMARK_ADIPOGENESIS")
        - `leadingEdge`: Genes in the pathway (semicolon-separated)
        - `pval`, `padj`, `NES`: Optional statistical values for display
        
        The app will extract genes from the `leadingEdge` column and build networks for each pathway.
        
        ### 🔬 How It Works
        
        1. **Upload**: Upload your CSV file with pathway data
        2. **Parse**: The app extracts gene lists from each pathway
        3. **Fetch**: For each pathway, fetches PPI data from STRING database
        4. **Analyze**: Calculates network metrics and simulates attacks
        5. **Visualize**: Generates interactive plots and tables
        6. **Export**: Download results for further analysis
        
        ### 📊 Understanding the Results
        
        **Robustness Curve**: Shows how the network maintains integrity as nodes are removed.
        - **Green line**: Random removal (non-specific drug effects)
        - **Red dashed line**: Targeted removal of hubs (specific drug effects)
        - Steeper decline = more vulnerable network
        
        **Network Graph**: Interactive visualization of the PPI network.
        - Node size/color = degree (number of connections)
        - Orange nodes = top hubs
        - Hover for details
        
        **Centrality Metrics**:
        - **Degree Centrality**: Number of connections (hubs)
        - **Betweenness Centrality**: Controls information flow (bottlenecks)
        - **Closeness Centrality**: Central position in network
        
        ### ⚙️ Settings
        
        **Species (TaxID)**: Select the organism for your analysis
        
        **STRING API Settings**:
        - **Minimum Interaction Score**: Higher = more confident interactions only
        - **Additional Nodes**: How many extra nodes to include in network expansion
        
        **Analysis Settings**:
        - **Max Pathways**: Limit batch processing to prevent long wait times
        
        ### 🐛 Troubleshooting
        
        **No data found**: The genes may not have interactions in STRING for the selected species.
        Try a different species or lower the interaction score threshold.
        
        **Slow processing**: STRING API has rate limits. For large gene sets, consider:
        - Reducing the number of pathways
        - Using a lower interaction score
        - Adding fewer additional nodes
        
        **Empty CSV**: Make sure your CSV has a `leadingEdge` column with semicolon-separated genes.
        
        ### 📚 References
        
        - STRING Database: https://string-db.org/
        - NetworkX: https://networkx.org/
        - GSEA: https://www.gsea-msigdb.org/gsea/index.jsp
        
        ---
        
        **Version**: 2.0 (Streamlit Edition)  
        **Built with**: Streamlit, NetworkX, Plotly, STRING API
        """)


if __name__ == "__main__":
    main()
