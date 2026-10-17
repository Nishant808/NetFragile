# NetFragile Utilities Package
# Multi-network comparative analysis for biological networks

__version__ = "2.0.0"
__author__ = "NetFragile Team"

from .file_parser import EdgeRParser, parse_deg_file
from .network_builder import NetworkBuilder, build_network_from_genes
from .robustness import RobustnessAnalyzer, simulate_removal
from .comparison import NetworkComparator, calculate_similarity_matrix
from .normalizing import NormalizingTargetFinder

__all__ = [
    "EdgeRParser",
    "parse_deg_file",
    "NetworkBuilder",
    "build_network_from_genes",
    "RobustnessAnalyzer",
    "simulate_removal",
    "NetworkComparator",
    "calculate_similarity_matrix",
    "NormalizingTargetFinder",
]
