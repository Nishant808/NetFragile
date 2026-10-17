"""
Robustness Analyzer for NetFragile 2.0

Handles simulation of network attacks (random and targeted removal)
and tracking of network integrity (Largest Connected Component size).
"""

import networkx as nx
import numpy as np
from typing import List, Tuple, Optional, Dict, Any
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import random

logger = logging.getLogger(__name__)


class RobustnessAnalyzer:
    """
    Analyzer for network robustness under node removal attacks.
    
    Simulates both random and targeted (degree-based) node removal
    and tracks the size of the largest connected component (LCC).
    """
    
    def __init__(self):
        """Initialize the robustness analyzer."""
        self._rng = np.random.RandomState(42)
    
    def analyze_network(
        self,
        graph: nx.Graph,
        removal_strategies: Optional[List[str]] = None,
        stop_threshold: float = 0.05
    ) -> Dict[str, Tuple[List[float], List[float]]]:
        """
        Analyze network robustness under different removal strategies.
        
        Args:
            graph: NetworkX graph to analyze
            removal_strategies: List of strategies to test
                Options: 'random', 'degree', 'betweenness'
                Default: ['random', 'degree']
            stop_threshold: Stop simulation if LCC < threshold * original_size
            
        Returns:
            Dictionary mapping strategy names to (x, y) curves
            where x = fraction of nodes removed, y = LCC size fraction
        """
        if removal_strategies is None:
            removal_strategies = ['random', 'degree']
        
        results = {}
        original_size = len(graph)
        
        if original_size == 0:
            logger.warning("Empty graph provided")
            for strategy in removal_strategies:
                results[strategy] = ([0.0], [1.0])
            return results
        
        for strategy in removal_strategies:
            logger.info(f"Simulating {strategy} removal on network with {original_size} nodes")
            x, y = self._simulate_removal(graph, strategy, stop_threshold)
            results[strategy] = (x, y)
        
        return results
    
    def _simulate_removal(
        self,
        graph: nx.Graph,
        strategy: str,
        stop_threshold: float = 0.05
    ) -> Tuple[List[float], List[float]]:
        """
        Simulate node removal using a specific strategy.
        
        Args:
            graph: NetworkX graph
            strategy: Removal strategy ('random', 'degree', 'betweenness')
            stop_threshold: Stop if LCC < threshold * original_size
            
        Returns:
            Tuple of (x_values, y_values) where:
            - x_values: Fraction of nodes removed
            - y_values: LCC size as fraction of original
        """
        # Make a copy to avoid modifying the original
        G = graph.copy()
        original_size = len(G)
        
        if original_size == 0:
            return [0.0], [1.0]
        
        # Determine node removal order
        if strategy == 'random':
            nodes = list(G.nodes())
            self._rng.shuffle(nodes)
        elif strategy == 'degree':
            degree_dict = dict(G.degree())
            nodes = sorted(degree_dict.keys(), key=lambda x: degree_dict[x], reverse=True)
        elif strategy == 'betweenness':
            betweenness = nx.betweenness_centrality(G)
            nodes = sorted(betweenness.keys(), key=lambda x: betweenness[x], reverse=True)
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
        
        # Initialize tracking lists
        x_values = []
        y_values = []
        
        # Simulate removal
        for i, node in enumerate(nodes):
            # Remove node if it still exists
            if node in G:
                G.remove_node(node)
            
            # Calculate current LCC size
            if len(G) > 0:
                largest_cc = max(nx.connected_components(G), key=len)
                lcc_size = len(largest_cc)
            else:
                lcc_size = 0
            
            # Store results
            fraction_removed = (i + 1) / original_size
            lcc_fraction = lcc_size / original_size
            
            x_values.append(fraction_removed)
            y_values.append(lcc_fraction)
            
            # Check stopping condition
            if lcc_fraction < stop_threshold:
                logger.debug(f"Stopping {strategy} simulation at {fraction_removed:.2f} removed")
                break
        
        return x_values, y_values
    
    def simulate_single_removal(
        self,
        graph: nx.Graph,
        gene_to_remove: str,
        strategy: str = 'degree',
        stop_threshold: float = 0.05
    ) -> Tuple[List[float], List[float]]:
        """
        Simulate robustness after removing a specific gene.
        
        This is used for normalizing target discovery.
        
        Args:
            graph: Original NetworkX graph
            gene_to_remove: Gene/node to remove before simulation
            strategy: Removal strategy for the simulation
            stop_threshold: Stop simulation threshold
            
        Returns:
            Tuple of (x_values, y_values) for the modified network
        """
        # Create modified graph
        G = graph.copy()
        
        if gene_to_remove in G:
            G.remove_node(gene_to_remove)
            logger.debug(f"Removed {gene_to_remove}, network now has {len(G)} nodes")
        else:
            logger.warning(f"Gene {gene_to_remove} not found in network")
        
        # Run simulation on modified graph
        return self._simulate_removal(G, strategy, stop_threshold)


def simulate_removal(
    graph: nx.Graph,
    strategy: str = 'degree',
    stop_threshold: float = 0.05
) -> Tuple[List[float], List[float]]:
    """
    Convenience function to simulate node removal on a graph.
    
    Args:
        graph: NetworkX graph
        strategy: Removal strategy ('random', 'degree', 'betweenness')
        stop_threshold: Stop if LCC < threshold * original_size
        
    Returns:
        Tuple of (x_values, y_values)
    """
    analyzer = RobustnessAnalyzer()
    result = analyzer.analyze_network(graph, [strategy], stop_threshold)
    return result[strategy]
