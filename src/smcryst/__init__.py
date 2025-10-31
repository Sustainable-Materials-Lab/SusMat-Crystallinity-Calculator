"""
SusMat Crystallinity Calculator

A collection of Python scripts for determining the crystallinity of cellulosic samples
using WAXS (Wide Angle X-ray Scattering) data processed through TOPAS Academic.

This package provides tools for:
- Automated crystallinity calculation using the Ruland method
- Support for multiple cellulose phases (Cellulose I, Cellulose II)
- Local and remote (iRODS) data processing
- Automatic XYE format detection
- SVG and PNG plot generation

Main modules:
- cell_cryst: Local file processing with CLI interface
- cell_cryst_remote: Remote iRODS data processing
- data_importer: Centralized data parsing and import functionality
"""

__version__ = "0.3.0"
__author__ = "Samuel Eyley"
__email__ = "samuel.eyley@kuleuven.be"

# Import main classes and functions for easy access
from .data_importer import DataImporter
from .cell_cryst import cli as local_cli
from .cell_cryst_remote import cli as remote_cli

# Define what gets imported with "from smcryst import *"
__all__ = [
    "DataImporter",
    "local_cli", 
    "remote_cli",
    "__version__",
    "__author__",
    "__email__"
]