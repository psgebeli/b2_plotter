
# Preamble
import pytest as pt
from b2_plotter.Plotter import Plotter, parse_cmd, construct_dfs, get_fom
import uproot as up
import os
import argparse as ap
from unittest.mock import patch
import pandas as pd
import matplotlib.pyplot as plt 

# Define a testfile and test columns
mixed = 'mc/xipipi_miprompt_700fb.root'
ccbar = 'mc/xipipi_ccprompt_700fb.root'

mycols= ['xipipi_xic_M', 'xipipi_xi_significanceOfDistance', 'xipipi_lambda_p_protonID', 'xipipi_xi_M', 'xipipi_xic_isSignal']

# Define a test cut
xicmassrangeloose = '2.3 < xipipi_xic_M < 2.65'

# Create dataframes
with up.open(mixed) as mixed:
    tree = mixed['xic_tree']
    df_mixed = tree.arrays(filter_name = mycols + ['xi03pi_xic_isSignal'], library = "pd")

with up.open(ccbar) as ccbar:
    tree = ccbar['xic_tree']
    df_ccbar = tree.arrays(filter_name = mycols + ['xi03pi_xic_isSignal'], library = "pd")

plotter = Plotter(isSigvar = 'xic_isSignal', mcdfs={'mixed': df_mixed}, signaldf = df_ccbar)

# CLASS TESTING
# -------------
def test_constructor():
    assert isinstance(plotter, Plotter)

def test_plot():
    for var in mycols[:-1]:
        plotter.plot(var, cuts = xicmassrangeloose).savefig(f'plot_{var}.png')
        assert os.path.isfile(f'plot_{var}.png')

def test_plotFom():
    for var in mycols[:-1]:
        fom, cut = plotter.plotFom(var, cuts = xicmassrangeloose, massvar = 'xipipi_xic_M', signalregion = (2.46, 2.475))
        fom.savefig(f'fom_{var}.png')
        assert os.path.isfile(f'fom_{var}.png')
        assert isinstance(cut, float)

def test_plotStep():
    for var in mycols[:-1]:
        plotter.plotStep(var, cuts = xicmassrangeloose).savefig(f'step_{var}.png')
        assert os.path.isfile(f'step_{var}.png')

def test_getpurity():
    assert isinstance(plotter.get_purity(xicmassrangeloose, 'xipipi_xic_M', (2.46, 2.475)), float)

def test_getsigeff():
    assert isinstance(plotter.get_sigeff(xicmassrangeloose, 'xipipi_xic_M', (2.46, 2.475)), float)

def test_errors():
    with pt.raises(TypeError):

        # Test isSigvar type errors
        plotter1 = Plotter(isSigvar = 5, mcdfs={'ccbar': df_mixed}, signaldf = df_mixed)
        plotter2 = Plotter(isSigvar = xic_M, mcdfs={'ccbar': df_mixed}, signaldf = df_mixed)

        # Test mcdfs 
        plotter3 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs = 5, signaldf = df_mixed)
        plotter4 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs = 'hello', signaldf = df_mixed)
        plotter5 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs={5 : df_mixed}, signaldf= df_mixed)
        plotter7 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs={'label': 5}, signaldf=df_mixed)
        plotter8 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs={'label': 'hello'}, signaldf=df_mixed)

        # Test signaldf
        plotter9 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs={'ccbar' : df_mixed}, signaldf = 5)
        plotter10 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs={'ccbar' : df_mixed}, signaldf = 'hello')
# -------------

# FUNCTION TESTING
# ----------------
        
def test_parse_cmd():

    # Use argparse.Namespace to simulate the parsed command line arguments
    parsed_args = ap.Namespace(input = 'path/to/MC', prefix = 'xic_prefix_name')

    # Patch the argparse.ArgumentParser to return the simulated parsed_args
    with patch('ap.ArgumentParser.parse_args', return_value = parsed_args):
        result = parse_cmd()

    # Check if the returned result matches the expected result
    assert result == parsed_args

def test_construct_dfs():

    mcdfs = construct_dfs('mc/', mycols = mycols, prefix = 'xipipi_xic')

    assert isinstance(mcdfs, dict)
    for df in mcdfs:
        assert isinstance(df, pd.DataFrame)

def test_get_fom():

    lessfom, lesscut, greaterfom, greatercut = get_fom(cuts = xicmassrangeloose, var = 'xipipi_xi_significanceOfDistance', prefix = 'xipipi_xic', plotter = plotter)

    assert isinstance(lessfom, plt.Figure)
    assert isinstance(greaterfom, plt.Figure)

    assert isinstance(lesscut, float)
    assert isinstance(greatercut, float)
    