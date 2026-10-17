"""
Normalizing Target Finder for NetFragile 2.0

Identifies genes in a target network whose removal makes its robustness
profile resemble a reference network. This is the core feature for finding
druggable targets that "normalize" pathological networks.
"""

import networkx as nx
import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Any, Optional
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import time

logger = logging.getLogger(__name__)


class NormalizingTargetFinder:
    """
    Finder for genes that normalize a target network toward a reference network.
    
    For each gene in the target network, simulates its removal and compares
    the resulting robustness curve to the reference network's curve.
    """
    
    def __init__(self, max_workers: int = 4):
        """
        Initialize the normalizing target finder.
        
        Args:
            max_workers: Maximum number of parallel workers for gene analysis
        """
        self.max_workers = max_workers
        self._rng = np.random.RandomState(42)
    
    def find_targets(
        self,
        target_network_data: Dict[str, Any],
        reference_curve: Tuple[List[float], List[float]],
        strategy: str = 'degree',
        stop_threshold: float = 0.05,
        num_points: int = 100,
        min_network_size: int = 10
    ) -> pd.DataFrame:
        """
        Find genes in target network that normalize it toward reference.
        
        Args:
            target_network_data: Dictionary containing:
                - 'graph': NetworkX graph of target network
                - 'targ_curve': Original targeted removal curve
            reference_curve: Tuple of (x, y) for reference network's targeted curve
            strategy: Removal strategy for robustness simulation
            stop_threshold: Stop simulation if LCC < threshold
            num_points: Number of points for curve interpolation
            min_network_size: Skip genes if resulting network is too small
            
        Returns:
            DataFrame with results sorted by similarity score
        """
        graph = target_network_data.get('graph')
        original_curve = target_network_data.get('targ_curve', ([], []))
        
        if graph is None:
            raise ValueError("Target network graph is required")
        
        if not reference_curve or len(reference_curve[0]) == 0:
            raise ValueError("Reference curve is required")
        
        # Calculate original similarity
        original_similarity = self._calculate_curve_similarity(
            original_curve, reference_curve, num_points
        )
        
        # Get all genes in target network
        all_genes = list(graph.nodes())
        logger.info(f"Testing {len(all_genes)} genes in target network")
        
        # Prepare results list
        results = []
        
        # Process genes in parallel batches
        batch_size = min(50, len(all_genes) // self.max_workers + 1)
        
        for i in range(0, len(all_genes), batch_size):
            batch = all_genes[i:i + batch_size]
            logger.info(f"Processing batch {i//batch_size + 1}/{(len(all_genes)-1)//batch_size + 1}")
            
            batch_results = self._process_batch(
                graph, batch, reference_curve, strategy, 
                stop_threshold, num_points, min_network_size
            )
            results.extend(batch_results)
        
        # Create DataFrame
        df = pd.DataFrame(results)
        
        if len(df) > 0:
            # Sort by similarity score (descending)
            df = df.sort_values('similarity_score', ascending=False)
            
            # Calculate effect size
            df['effect_size'] = df['similarity_score'] - original_similarity
            
            # Add rank
            df['rank'] = range(1, len(df) + 1)
        
        return df
    
    def _process_batch(
        self,
        graph: nx.Graph,
        gene_batch: List[str],
        reference_curve: Tuple[List[float], List[float]],
        strategy: str,
        stop_threshold: float,
        num_points: int,
        min_network_size: int
    ) -> List[Dict[str, Any]]:
        """
        Process a batch of genes in parallel.
        
        Args:
            graph: Target network graph
            gene_batch: List of genes to test
            reference_curve: Reference network's robustness curve
            strategy: Removal strategy
            stop_threshold: Stop threshold for simulation
            num_points: Number of points for interpolation
            min_network_size: Minimum network size after gene removal
            
        Returns:
            List of result dictionaries
        """
        results = []
        
        # Create partial function for parallel execution
        worker_func = partial(
            self._analyze_single_gene,
            graph=graph,
            reference_curve=reference_curve,
            strategy=strategy,
            stop_threshold=stop_threshold,
            num_points=num_points,
            min_network_size=min_network_size
        )
        
        # Use ThreadPoolExecutor for parallel processing
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(worker_func, gene): gene 
                for gene in gene_batch
            }
            
            for future in as_completed(futures):
                gene = futures[future]
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                except Exception as e:
                    logger.warning(f"Error processing gene {gene}: {e}")
        
        return results
    
    def _analyze_single_gene(
        self,
        gene: str,
        graph: nx.Graph,
        reference_curve: Tuple[List[float], List[float]],
        strategy: str,
        stop_threshold: float,
        num_points: int,
        min_network_size: int
    ) -> Optional[Dict[str, Any]]:
        """
        Analyze the effect of removing a single gene.
        
        Args:
            gene: Gene to remove
            graph: Target network graph
            reference_curve: Reference network's robustness curve
            strategy: Removal strategy
            stop_threshold: Stop threshold
            num_points: Number of points for interpolation
            min_network_size: Minimum network size
            
        Returns:
            Result dictionary or None if analysis fails
        """
        try:
            # Create modified graph
            G = graph.copy()
            
            if gene not in G:
                logger.debug(f"Gene {gene} not in network, skipping")
                return None
            
            G.remove_node(gene)
            
            # Check if network is still large enough
            if len(G) < min_network_size:
                logger.debug(f"Network too small after removing {gene}, skipping")
                return None
            
            # Simulate robustness on modified network
            from .robustness import RobustnessAnalyzer
            analyzer = RobustnessAnalyzer()
            x, y = analyzer._simulate_removal(G, strategy, stop_threshold)
            
            if not x or len(x) == 0:
                return None
            
            # Calculate similarity to reference curve
            similarity = self._calculate_curve_similarity(
                (x, y), reference_curve, num_points
            )
            
            # Get gene degree in original network
            degree = graph.degree(gene)
            
            return {
                'gene': gene,
                'similarity_score': similarity,
                'degree': degree,
                'network_size_after_removal': len(G),
                'curve_x': x,
                'curve_y': y
            }
            
        except Exception as e:
            logger.warning(f"Error analyzing gene {gene}: {e}")
            return None
    
    def _calculate_curve_similarity(
        self,
        curve1: Tuple[List[float], List[float]],
        curve2: Tuple[List[float], List[float]],
        num_points: int = 100
    ) -> float:
        """
        Calculate similarity between two robustness curves.
        
        Args:
            curve1: First curve (x, y)
            curve2: Second curve (x, y)
            num_points: Number of points for interpolation
            
        Returns:
            Pearson correlation coefficient (0-1)
        """
        x1, y1 = curve1
        x2, y2 = curve2
        
        if not x1 or not x2 or len(x1) == 0 or len(x2) == 0:
            return 0.0
        
        # Interpolate to common x-values
        common_x = np.linspace(0, 1, num_points)
        y1_interp = np.interp(common_x, x1, y1, left=1.0, right=0.0)
        y2_interp = np.interp(common_x, x2, y2, left=1.0, right=0.0)
        
        # Calculate Pearson correlation
        try:
            corr = np.corrcoef(y1_interp, y2_interp)[0, 1]
            # Handle NaN (perfect correlation or anti-correlation)
            if np.isnan(corr):
                return 0.0
            # Convert to 0-1 scale
            return (corr + 1) / 2
        except Exception:
            return 0.0
    
    def get_top_targets(
        self,
        results_df: pd.DataFrame,
        top_n: int = 10,
        include_shared: bool = True,
        reference_genes: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Get top normalizing targets with additional annotations.
        
        Args:
            results_df: DataFrame from find_targets()
            top_n: Number of top targets to return
            include_shared: Whether to include genes shared with reference
            reference_genes: List of genes in reference network (for annotation)
            
        Returns:
            DataFrame with top targets and annotations
        """
        if len(results_df) == 0:
            return pd.DataFrame()
        
        # Get top N
        top_results = results_df.head(top_n).copy()
        
        # Add annotation for whether gene is in reference network
        if reference_genes:
            reference_set = set(reference_genes)
            top_results['in_reference'] = top_results['gene'].apply(
                lambda x: x in reference_set
            )
        else:
            top_results['in_reference'] = False
        
        return top_results
