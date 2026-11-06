"""Patient-level association estimates for multi-patient single-cell studies."""

from .association import associate
from .result import AssociationResult
from . import datasets

__version__ = "1.1.0"
__all__ = ["associate", "AssociationResult", "datasets", "__version__"]
