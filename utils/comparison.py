"""
Network Comparison Module for NetFragile 2.0

Handles comparative analysis between multiple networks,
including similarity metrics and visualization data.
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Any, Optional
import logging
from scipy.stats import pearsonr, spearmanr
from itertools import combinations

logger = logging.getLogger(__name__)


class NetworkComparator:
    """
    Comparator for analyzing multiple networks together.
    
    Computes similarity metrics between robustness curves
    and provides data for comparative visualizations.
    """
    
    def __init__(self):
        """Initialize the network comparator."""
        pass
    
    def calculate_similarity_matrix(
        self,
        networks_data: Dict[str, Dict[str, Any]]
    ) -> pd.DataFrame:
        """
        Calculate pairwise similarity between all networks.
        
        Args:
            networks_data: Dictionary mapping factor names to network data
                Each network data should contain:
                - 'targ_curve': Tuple of (x, y) for targeted removal
                - 'rand_curve': Tuple of (x, y) for random removal
                
        Returns:
            DataFrame with pairwise similarity scores
        """
        factor_names = list(networks_data.keys())
        n = len(factor_names)
        
        # Initialize matrix
        matrix = pd.DataFrame(
            np.ones((n, n)),
            index=factor_names,
            columns=factor_names
        )
        
        # Calculate pairwise similarities
        for i, j in combinations(range(n), 2):
            factor1 = factor_names[i]
            factor2 = factor_names[j]
            
            # Use targeted removal curves for comparison
            curve1 = networks_data[factor1].get('targ_curve', ([], []))
            curve2 = networks_data[factor2].get('targ_curve', ([], []))
            
            similarity = self._calculate_curve_similarity(curve1, curve2)
            matrix.loc[factor1, factor2] = similarity
            matrix.loc[factor2, factor1] = similarity
        
        return matrix
    
    def _calculate_curve_similarity(
        self,
        curve1: Tuple[List[float], List[float]],
        curve2: Tuple[List[float], List[float]]
    ) -> float:
        """
        Calculate similarity between two robustness curves.
        
        Uses Pearson correlation after interpolating to common x-values.
        
        Args:
            curve1: Tuple of (x_values, y_values)
            curve2: Tuple of (x_values, y_values)
            
        Returns:
            Similarity score (0-1, higher is more similar)
        """
        x1, y1 = curve1
        x2, y2 = curve2
        
        if not x1 or not x2:
            return 0.0
        
        # Interpolate both curves to common x-values
        common_x, y1_interp, y2_interp = self._interpolate_curves(x1, y1, x2, y2)
        
        if len(common_x) < 2:
            return 0.0
        
        # Calculate Pearson correlation
        try:
            corr, _ = pearsonr(y1_interp, y2_interp)
            # Convert to 0-1 scale (correlation can be -1 to 1)
            similarity = (corr + 1) / 2
            return float(similarity)
        except Exception as e:
            logger.warning(f"Could not calculate correlation: {e}")
            return 0.0
    
    def _interpolate_curves(
        self,
        x1: List[float],
        y1: List[float],
        x2: List[float],
        y2: List[float],
        num_points: int = 100
    ) -> Tuple[List[float], List[float], List[float]]:
        """
        Interpolate two curves to common x-values.
        
        Args:
            x1, y1: First curve
            x2, y2: Second curve
            num_points: Number of points for interpolation
            
        Returns:
            Tuple of (common_x, y1_interp, y2_interp)
        """
        # Create common x-values
        min_x = min(min(x1), min(x2)) if x1 and x2 else 0
        max_x = max(max(x1), max(x2)) if x1 and x2 else 1
        common_x = np.linspace(min_x, max_x, num_points)
        
        # Interpolate y-values
        y1_interp = np.interp(common_x, x1, y1, left=1.0, right=0.0)
        y2_interp = np.interp(common_x, x2, y2, left=1.0, right=0.0)
        
        return list(common_x), list(y1_interp), list(y2_interp)
    
    def get_overlay_plot_data(
        self,
        networks_data: Dict[str, Dict[str, Any]],
        curve_type: str = 'targ'
    ) -> Dict[str, Any]:
        """
        Prepare data for overlay robustness curve plot.
        
        Args:
            networks_data: Dictionary of network data
            curve_type: 'targ' for targeted or 'rand' for random
            
        Returns:
            Dictionary with plot-ready data for all networks
        """
        plot_data = {}
        
        for factor_name, data in networks_data.items():
            curve = data.get(f'{curve_type}_curve', ([], []))
            if curve and len(curve[0]) > 0:
                x, y = curve
                plot_data[factor_name] = {
                    'x': x,
                    'y': y,
                    'name': factor_name
                }
        
        return plot_data
    
    def calculate_network_metrics(
        self,
        networks_data: Dict[str, Dict[str, Any]]
    ) -> pd.DataFrame:
        """
        Calculate comparative metrics for all networks.
        
        Args:
            networks_data: Dictionary of network data
            
        Returns:
            DataFrame with metrics for each network
        """
        metrics_list = []
        
        for factor_name, data in networks_data.items():
            graph = data.get('graph')
            
            if graph:
                metrics = {
                    'Factor': factor_name,
                    'Nodes': graph.number_of_nodes(),
                    'Edges': graph.number_of_edges(),
                    'Density': nx.density(graph) if graph.number_of_nodes() > 1 else 0,
                    'Average Degree': sum(dict(graph.degree()).values()) / graph.number_of_nodes() if graph.number_of_nodes() > 0 else 0,
                }
                
                # Add robustness metrics
                targ_curve = data.get('targ_curve', ([], []))
                if targ_curve and len(targ_curve[0]) > 0:
                    x, y = targ_curve
                    # Area under curve (approximate)
                    auc = np.trapz(y, x)
                    metrics['Targeted AUC'] = auc
                    
                    # Percolation threshold (fraction removed when LCC drops below 50%)
                    for i, (x_val, y_val) in enumerate(zip(x, y)):
                        if y_val < 0.5:
                            metrics['Percolation Threshold'] = x_val
                            break
                    else:
                        metrics['Percolation Threshold'] = 1.0
                
                metrics_list.append(metrics)
        
        return pd.DataFrame(metrics_list)


def calculate_similarity_matrix(
    networks_data: Dict[str, Dict[str, Any]]
) -> pd.DataFrame:
    """
    Convenience function to calculate similarity matrix.
    
    Args:
        networks_data: Dictionary of network data
        
    Returns:
        DataFrame with pairwise similarity scores
    """
    comparator = NetworkComparator()
    return comparator.calculate_similarity_matrix(networks_data)
