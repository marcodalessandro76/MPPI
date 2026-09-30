"""
Common fixtures of the MPPI test suite.

The tests use the reference data of the tutorials in sphinx_source/tutorials/Reference_data.
Tests that need the QuantumESPRESSO or Yambo executables are marked with ``requires_qe`` or
``requires_yambo`` and are skipped automatically if the executables are not found in the PATH.
"""
import os, shutil
import matplotlib
matplotlib.use('Agg')
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TUTORIALS = os.path.join(ROOT,'sphinx_source','tutorials')

def pytest_configure(config):
    config.addinivalue_line('markers','requires_qe: the test needs the pw.x executable')
    config.addinivalue_line('markers','requires_yambo: the test needs the yambo executable')

def pytest_collection_modifyitems(config, items):
    skip_qe = pytest.mark.skip(reason='pw.x not found in PATH')
    skip_yambo = pytest.mark.skip(reason='yambo not found in PATH')
    for item in items:
        if 'requires_qe' in item.keywords and shutil.which('pw.x') is None:
            item.add_marker(skip_qe)
        if 'requires_yambo' in item.keywords and shutil.which('yambo') is None:
            item.add_marker(skip_yambo)

@pytest.fixture
def ref_dir():
    """Folder with the reference data of the tutorials"""
    return os.path.join(TUTORIALS,'Reference_data')

@pytest.fixture
def io_dir():
    """Folder with the input files of the tutorials"""
    return os.path.join(TUTORIALS,'IO_files')
