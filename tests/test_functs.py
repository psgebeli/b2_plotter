import pytest as pt
import argparse as ap
import root_pandas as rp
from b2_plotter.Plotter import Plotter, parse_cmd, construct_dfs, get_fom
from unittest.mock import patch

def test_parse_cmd():

    # Use argparse.Namespace to simulate the parsed command line arguments
    parsed_args = ap.Namespace(input = 'path/to/MC', prefix = 'xic_prefix_name')

    # Patch the argparse.ArgumentParser to return the simulated parsed_args
    with patch('ap.ArgumentParser.parse_args', return_value = parsed_args):
        result = parse_cmd()

    # Check if the returned result matches the expected result
    assert result == parsed_args

def test_construct_dfs():
