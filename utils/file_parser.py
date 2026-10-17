"""
EdgeR DEG File Parser for NetFragile 2.0

Handles parsing of EdgeR output CSV/TSV files to extract gene lists.
Supports various EdgeR output formats and filtering options.
"""

import pandas as pd
import numpy as np
from typing import List, Optional, Union
import logging

logger = logging.getLogger(__name__)


class EdgeRParser:
    """
    Parser for EdgeR differential expression output files.
    
    Supports:
    - CSV and TSV formats
    - Various column naming conventions
    - Filtering by p-value, adjusted p-value, logFC
    """
    
    # Common gene column names in EdgeR output
    GENE_COLUMN_NAMES = [
        'gene', 'Gene', 'GeneID', 'gene_id', 'GeneID',
        'symbol', 'Symbol', 'gene_name', 'GeneName',
        'ensembl', 'Ensembl', 'ensembl_id', 'EnsemblID'
    ]
    
    # Common significance columns
    PValue_COLUMN_NAMES = [
        'pvalue', 'p.value', 'PValue', 'P.Value',
        'p_val', 'p-val', 'PVal'
    ]
    
    ADJ_PVALUE_COLUMN_NAMES = [
        'padj', 'adj_p', 'adj.p', 'AdjP', 'Adj.P',
        'fdr', 'FDR', 'qvalue', 'QValue', 'q_value',
        'p.adjust', 'p.adjust'
    ]
    
    LOGFC_COLUMN_NAMES = [
        'logfc', 'logFC', 'log2fc', 'log2FC',
        'log_fc', 'log_fold_change', 'LogFoldChange'
    ]
    
    def __init__(self):
        """Initialize the EdgeR parser."""
        pass
    
    def parse_file(
        self,
        file_path: str,
        gene_column: Optional[str] = None,
        pvalue_threshold: Optional[float] = None,
        padj_threshold: Optional[float] = None,
        logfc_threshold: Optional[float] = None,
        min_genes: int = 10
    ) -> List[str]:
        """
        Parse an EdgeR output file and extract gene list.
        
        Args:
            file_path: Path to CSV/TSV file
            gene_column: Specific column name for genes (auto-detected if None)
            pvalue_threshold: Filter genes by p-value < threshold
            padj_threshold: Filter genes by adjusted p-value < threshold
            logfc_threshold: Filter genes by |logFC| > threshold
            min_genes: Minimum number of genes required (raises error if fewer)
            
        Returns:
            List of gene symbols
            
        Raises:
            ValueError: If file cannot be parsed or too few genes
        """
        # Read the file
        try:
            if file_path.endswith('.tsv') or file_path.endswith('.txt'):
                df = pd.read_csv(file_path, sep='\t')
            else:
                df = pd.read_csv(file_path)
        except Exception as e:
            logger.error(f"Failed to read file {file_path}: {e}")
            raise ValueError(f"Could not read file: {e}")
        
        if df.empty:
            logger.error(f"File {file_path} is empty")
            raise ValueError("File is empty or contains no data")
        
        # Find gene column
        gene_col = gene_column
        if gene_col is None:
            gene_col = self._find_gene_column(df)
        
        if gene_col is None:
            logger.error(f"Could not find gene column in file {file_path}")
            raise ValueError("Could not identify gene column. Please specify manually.")
        
        # Extract gene list
        gene_list = df[gene_col].astype(str).tolist()
        
        # Remove empty/NaN genes
        gene_list = [g for g in gene_list if g and g != 'nan' and g != 'NaN' and g != 'NA']
        
        # Apply filters if specified
        if pvalue_threshold is not None:
            pval_col = self._find_column(df, self.PValue_COLUMN_NAMES)
            if pval_col:
                mask = df[pval_col] < pvalue_threshold
                gene_list = [g for g, keep in zip(gene_list, mask) if keep]
            else:
                logger.warning(f"p-value column not found, skipping p-value filter")
        
        if padj_threshold is not None:
            padj_col = self._find_column(df, self.ADJ_PVALUE_COLUMN_NAMES)
            if padj_col:
                mask = df[padj_col] < padj_threshold
                gene_list = [g for g, keep in zip(gene_list, mask) if keep]
            else:
                logger.warning(f"Adjusted p-value column not found, skipping padj filter")
        
        if logfc_threshold is not None:
            logfc_col = self._find_column(df, self.LOGFC_COLUMN_NAMES)
            if logfc_col:
                mask = np.abs(df[logfc_col]) > logfc_threshold
                gene_list = [g for g, keep in zip(gene_list, mask) if keep]
            else:
                logger.warning(f"logFC column not found, skipping logFC filter")
        
        # Remove duplicates while preserving order
        seen = set()
        gene_list = [g for g in gene_list if not (g in seen or seen.add(g))]
        
        # Check minimum genes
        if len(gene_list) < min_genes:
            logger.error(f"Only {len(gene_list)} genes found (minimum {min_genes} required)")
            raise ValueError(f"Insufficient genes: {len(gene_list)} < {min_genes}")
        
        logger.info(f"Parsed {len(gene_list)} genes from {file_path}")
        return gene_list
    
    def _find_gene_column(self, df: pd.DataFrame) -> Optional[str]:
        """Auto-detect the gene column from DataFrame."""
        for col in self.GENE_COLUMN_NAMES:
            if col in df.columns:
                return col
        
        # Fallback: use first column
        return df.columns[0]
    
    def _find_column(self, df: pd.DataFrame, possible_names: List[str]) -> Optional[str]:
        """Find a column from a list of possible names."""
        for name in possible_names:
            if name in df.columns:
                return name
        return None
    
    def get_file_info(self, file_path: str) -> dict:
        """
        Get metadata about an EdgeR file without parsing fully.
        
        Returns:
            Dictionary with file info (columns, row count, etc.)
        """
        try:
            if file_path.endswith('.tsv') or file_path.endswith('.txt'):
                df = pd.read_csv(file_path, sep='\t', nrows=5)
            else:
                df = pd.read_csv(file_path, nrows=5)
            
            return {
                'file_path': file_path,
                'columns': list(df.columns),
                'row_count': len(df),
                'detected_gene_column': self._find_gene_column(df),
                'detected_pvalue_column': self._find_column(df, self.PValue_COLUMN_NAMES),
                'detected_padj_column': self._find_column(df, self.ADJ_PVALUE_COLUMN_NAMES),
                'detected_logfc_column': self._find_column(df, self.LOGFC_COLUMN_NAMES),
            }
        except Exception as e:
            return {'file_path': file_path, 'error': str(e)}


def parse_deg_file(
    file_path: str,
    gene_column: Optional[str] = None,
    pvalue_threshold: Optional[float] = None,
    padj_threshold: Optional[float] = None,
    logfc_threshold: Optional[float] = None,
    min_genes: int = 10
) -> List[str]:
    """
    Convenience function to parse a DEG file.
    
    Args:
        file_path: Path to the DEG file
        gene_column: Column name for genes (auto-detected if None)
        pvalue_threshold: Filter by p-value
        padj_threshold: Filter by adjusted p-value
        logfc_threshold: Filter by |logFC|
        min_genes: Minimum number of genes required
        
    Returns:
        List of gene symbols
    """
    parser = EdgeRParser()
    return parser.parse_file(
        file_path,
        gene_column=gene_column,
        pvalue_threshold=pvalue_threshold,
        padj_threshold=padj_threshold,
        logfc_threshold=logfc_threshold,
        min_genes=min_genes
    )
