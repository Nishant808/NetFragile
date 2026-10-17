# NetFragile Streamlit App

A **Streamlit-based network robustness analyzer** that ingests CSV files with gene sets (e.g., from GSEA) and performs batch network analysis using STRING database.

## 🚀 Quick Start

### Option 1: Run Locally

```bash
# Clone or navigate to the repository
cd NetFragile

# Install dependencies
pip install -r streamlit_requirements.txt

# Run the Streamlit app
streamlit run streamlit_app.py
```

The app will open in your browser at `http://localhost:8501`

### Option 2: Deploy to Streamlit Cloud

1. **Fork this repository** to your GitHub account
2. **Go to [Streamlit Community Cloud](https://share.streamlit.io/)**
3. **Click "New App"** and select your repository
4. **Set the file path** to `streamlit_app.py`
5. **Deploy!**

### Option 3: Use Dev Container

The repository includes a `.devcontainer` configuration for VS Code:

1. Open the project in VS Code
2. Click "Reopen in Container" when prompted
3. The app will be available at port 8501

## 📥 Input CSV Format

The app expects a CSV file with pathway/gene set data, typically from GSEA analysis. Required columns:

- **`pathway`**: Pathway name (e.g., "HALLMARK_ADIPOGENESIS")
- **`leadingEdge`**: Semicolon-separated list of genes in the pathway

Optional columns (for display):
- `pval`: p-value
- `padj`: Adjusted p-value
- `NES`: Normalized Enrichment Score
- `size`: Gene set size

### Example CSV Structure:

```csv
pathway,pval,padj,NES,leadingEdge
HALLMARK_ADIPOGENESIS,1.94e-06,6.46e-06,-1.798,ITIH5;PTGER3;COL4A1;COL15A1;...
HALLMARK_ANGIOGENESIS,1.12e-06,3.97e-06,-2.187,FGFR1;FSTL1;VCAN;THBD;...
```

## ✨ Features

### 📁 CSV Ingestion
- Upload any CSV with gene sets
- Automatic parsing of semicolon-separated gene lists
- Preview data before analysis
- Select specific pathways or analyze all

### 🧬 Network Analysis
- **STRING Database Integration**: Fetches real PPI data
- **Network Construction**: Builds induced sub-networks
- **Centrality Analysis**: Identifies hubs and bottlenecks
  - Degree Centrality (hubs)
  - Betweenness Centrality (bottlenecks)
  - Closeness Centrality
- **Network Metrics**: Density, clustering, path length, etc.

### 🎯 Robustness Testing
- **Random Attack Simulation**: Non-specific drug effects
- **Targeted Attack Simulation**: Hub-focused drug effects
- **Interactive Robustness Curves**: Compare network vulnerability

### 📊 Visualizations
- **Robustness Plots**: Track network integrity during attacks
- **Interactive Network Graphs**: Zoom, pan, hover for details
- **Top Hub Highlighting**: Visual identification of key targets

### 💾 Data Export
- **Summary Tables**: CSV, JSON, Excel formats
- **Individual Pathway Data**: Centrality metrics per pathway
- **Timestamped Files**: Automatic naming with date/time

## 🎛️ Settings

### Species Selection
Choose from common model organisms:
- **9606**: Homo sapiens (Human)
- **10090**: Mus musculus (Mouse)
- **10116**: Rattus norvegicus (Rat)
- **83332**: Mycobacterium tuberculosis

### STRING API Settings
- **Minimum Interaction Score**: Filter by confidence (0-900)
- **Additional Nodes**: Expand network beyond seed genes

### Analysis Settings
- **Max Pathways**: Limit batch processing for performance

## 📊 Output Interpretation

### Robustness Curve
- **X-axis**: Fraction of nodes removed
- **Y-axis**: Network integrity (LCC size as % of original)
- **Green line**: Random removal (non-specific effects)
- **Red dashed line**: Targeted removal of hubs (specific effects)
- **Steeper decline**: More vulnerable network

### Network Graph
- **Node size/color**: Degree centrality (number of connections)
- **Orange nodes**: Top 5 hubs
- **Hover**: View gene name and degree

### Centrality Metrics
- **Degree Centrality**: Direct connections (hubs are good drug targets)
- **Betweenness Centrality**: Controls information flow (bottlenecks)
- **Closeness Centrality**: Central position in network

## 🔧 Requirements

### Python Packages
```
streamlit>=1.20.0
networkx>=3.0.0
pandas>=2.0.0
plotly>=5.0.0
requests>=2.28.0
openpyxl>=3.0.0  # For Excel export
```

### Dependencies
- Python 3.8+
- Internet connection (for STRING API)

## 📁 Project Structure

```
NetFragile/
├── streamlit_app.py          # Main Streamlit application
├── streamlit_requirements.txt # Dependencies for Streamlit app
├── app.py                    # Original Dash application
├── README.md                 # Original project documentation
└── .devcontainer/            # VS Code Dev Container config
```

## 🎯 Use Cases

### 1. Drug Target Identification
- Upload GSEA results
- Identify hubs in enriched pathways
- Prioritize targets based on network centrality

### 2. Pathway Vulnerability Assessment
- Compare robustness across pathways
- Find pathways most sensitive to targeted attacks
- Understand network-level effects of gene perturbations

### 3. Batch Network Analysis
- Process multiple gene sets at once
- Generate comparative reports
- Export data for downstream analysis

### 4. Educational Tool
- Visualize network robustness concepts
- Understand hub vs. random attack effects
- Explore PPI network properties

## 🐛 Troubleshooting

### Common Issues

**"No data found for species"**
- Check that your genes exist in STRING for the selected species
- Try lowering the interaction score threshold
- Verify gene names are correct (STRING uses specific identifiers)

**"STRING API request timed out"**
- STRING has rate limits; try again later
- Reduce the number of pathways being processed
- Use fewer additional nodes

**"Column 'leadingEdge' not found"**
- Your CSV must have a column with gene lists
- If using a different column name, modify the code or rename your column

### Performance Tips

1. **Limit pathways**: Use the "Max Pathways" setting to process fewer at once
2. **Lower confidence**: Reduce the interaction score to get more results
3. **Fewer nodes**: Reduce "Additional Nodes" for faster processing
4. **Small gene sets**: Pathways with 20-100 genes work best

## 📚 References

- **STRING Database**: [https://string-db.org/](https://string-db.org/)
- **NetworkX**: [https://networkx.org/](https://networkx.org/)
- **Streamlit**: [https://streamlit.io/](https://streamlit.io/)
- **Plotly**: [https://plotly.com/](https://plotly.com/)
- **GSEA**: [https://www.gsea-msigdb.org/](https://www.gsea-msigdb.org/)

## 📜 License

This project is open-source and available for academic and research use.

## 🙏 Acknowledgments

- STRING Database team for the PPI data API
- Streamlit team for the amazing web app framework
- NetworkX developers for the graph analysis library

---

**Version**: 2.0 (Streamlit Edition)  
**Maintainer**: Nishant808  
**Contact**: For issues or questions, please open a GitHub issue
