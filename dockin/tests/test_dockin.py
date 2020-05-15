"""
Unit and regression test for the dockin package.
"""

# Import package, test suite, and other packages as needed
import dockin
import pytest
import sys

def test_dockin_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "dockin" in sys.modules
