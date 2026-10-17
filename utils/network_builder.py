"""
Network Builder for NetFragile 2.0

Handles fetching PPI data from STRING database and constructing NetworkX graphs.
Includes caching for performance optimization.
"""

import networkx as nx
import pandas as pd
import requests
import io
from typing import List, Optional, Tuple, Dict, Any
import logging
from functools import lru_cache
import time

logger = logging.getLogger(__name__)

# Human species TaxID (locked as per requirements)
HUMAN_TAXID = 9606


class NetworkBuilder:
    """
    Builder for constructing biological networks from gene lists.
    
    Uses STRING database API to fetch protein-protein interactions.
    """
    
    STRING_API_URL = "https://string-db.org/api/tsv/network"
    
    def __init__(self, species_id: int = HUMAN_TAXID):
        """
        Initialize the network builder.
        
        Args:
            species_id: STRING species TaxID (default: 9606 for human)
        """
        self.species_id = species_id
        self._cache = {}
        self._last_api_call = 0
        self._api_call_count = 0
        
        # STRING API rate limiting: max 1 request per second
        self.API_DELAY = 1.0
    
    def build_network(
        self,
        gene_list: List[str],
        required_score: int = 700,
        add_nodes: int = 200,
        network_type: str = "functional"
    ) -> Tuple[nx.Graph, Dict[str, Any]]:
        """
        Build a network from a gene list using STRING API.
        
        Args:
            gene_list: List of gene symbols
            required_score: Minimum STRING confidence score (0-900)
            add_nodes: Number of additional nodes to include
            network_type: Type of network ("functional" or "physical")
            
        Returns:
            Tuple of (networkx.Graph, metadata_dict)
            
        Raises:
            ValueError: If network construction fails
        """
        if not gene_list:
            raise ValueError("Gene list cannot be empty")
        
        # Clean gene list
        gene_list = [str(g).strip() for g in gene_list if g and str(g).strip()]
        
        if len(gene_list) == 0:
            raise ValueError("No valid genes in list")
        
        # Create cache key
        cache_key = self._create_cache_key(gene_list, required_score, add_nodes, network_type)
        
        # Check cache
        if cache_key in self._cache:
            logger.info(f"Using cached network for {len(gene_list)} genes")
            return self._cache[cache_key]
        
        # Enforce rate limiting
        self._enforce_rate_limit()
        
        # Prepare STRING API parameters
        params = {
            "identifiers": "%0d".join(gene_list),
            "species": self.species_id,
            "required_score": required_score,
            "add_nodes": add_nodes,
            "network_type": network_type
        }
        
        logger.info(f"Fetching STRING data for {len(gene_list)} genes...")
        
        try:
            response = requests.get(self.STRING_API_URL, params=params)
            response.raise_for_status()
            
            data = response.text
            
            if not data.strip():
                logger.warning("STRING API returned empty response")
                raise ValueError("STRING API returned no data. Check gene symbols.")
            
            # Parse TSV data
            df = pd.read_csv(io.StringIO(data), sep="\t")
            
            # Build graph
            graph = nx.Graph()
            
            for _, row in df.iterrows():
                p1 = str(row['preferredName_A'])
                p2 = str(row['preferredName_B'])
                score = int(row['score'])
                
                # Only add if meets score requirement
                if score >= required_score:
                    graph.add_edge(p1, p2, weight=score)
            
            # Get largest connected component
            if graph.number_of_nodes() > 0:
                largest_cc = max(nx.connected_components(graph), key=len)
                graph = graph.subgraph(largest_cc).copy()
            
            # Store metadata
            metadata = {
                'original_gene_count': len(gene_list),
                'final_node_count': graph.number_of_nodes(),
                'final_edge_count': graph.number_of_edges(),
                'required_score': required_score,
                'add_nodes': add_nodes,
                'network_type': network_type,
                'species_id': self.species_id
            }
            
            # Cache result
            self._cache[cache_key] = (graph, metadata)
            
            logger.info(f"Built network: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
            
            return graph, metadata
            
        except requests.exceptions.RequestException as e:
            logger.error(f"STRING API request failed: {e}")
            raise ValueError(f"STRING API request failed: {e}")
        except Exception as e:
            logger.error(f"Error building network: {e}")
            raise ValueError(f"Error building network: {e}")
    
    def _create_cache_key(
        self,
        gene_list: List[str],
        required_score: int,
        add_nodes: int,
        network_type: str
    ) -> str:
        """Create a unique cache key for the given parameters."""
        gene_str = ",".join(sorted(gene_list))
        return f"{gene_str}|{required_score}|{add_nodes}|{network_type}"
    
    def _enforce_rate_limit(self):
        """Enforce STRING API rate limiting (max 1 request per second)."""
        elapsed = time.time() - self._last_api_call
        if elapsed < self.API_DELAY:
            sleep_time = self.API_DELAY - elapsed
            logger.debug(f"Sleeping for {sleep_time:.2f}s to respect rate limit")
            time.sleep(sleep_time)
        
        self._last_api_call = time.time()
        self._api_call_count += 1
    
    def clear_cache(self):
        """Clear the network cache."""
        self._cache = {}
        logger.info("Network cache cleared")
    
    def get_cache_stats(self) -> dict:
        """Get cache statistics."""
        return {
            'cached_networks': len(self._cache),
            'api_calls': self._api_call_count
        }


def build_network_from_genes(
    gene_list: List[str],
    species_id: int = HUMAN_TAXID,
    required_score: int = 700,
    add_nodes: int = 200,
    network_type: str = "functional"
) -> Tuple[nx.Graph, Dict[str, Any]]:
    """
    Convenience function to build a network from a gene list.
    
    Args:
        gene_list: List of gene symbols
        species_id: STRING species TaxID
        required_score: Minimum STRING confidence score
        add_nodes: Number of additional nodes
        network_type: Type of network
        
    Returns:
        Tuple of (networkx.Graph, metadata_dict)
    """
    builder = NetworkBuilder(species_id)
    return builder.build_network(
        gene_list,
        required_score=required_score,
        add_nodes=add_nodes,
        network_type=network_type
    )
