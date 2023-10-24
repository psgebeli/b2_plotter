import pytest as pt
from b2_plotter.Plotter import Plotter
import root_pandas as rp
import os

ccbar = '/belle2work/psgebeli/samples/gmc/mc15rib/xipipi/ccbar.root'

mycols= ['xipipi_xic_M', 'xipipi_xi_significanceOfDistance', 'xipipi_lambda_p_protonID', 'xipipi_xi_M', 'xipipi_xic_isSignal']

xicmassrangeloose = '2.3 < xipipi_xic_M < 2.65'
df_ccbar = rp.read_root(ccbar, key='xic_tree', columns = mycols)

plotter = Plotter(isSigvar='xic_isSignal', mcdfs={'ccbar': df_ccbar}, signaldf = df_ccbar)

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
        plotter1 = Plotter(isSigvar=5, mcdfs={'ccbar': df_ccbar}, signaldf = df_ccbar)
        plotter2 = Plotter(isSigvar=xic_M, mcdfs={'ccbar': df_ccbar}, signaldf = df_ccbar)

        # Test mcdfs 
        plotter3 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs = 5, signaldf = df_ccbar)
        plotter4 = Plotter(isSigvar = 'xipipi_xic_isSignal', mcdfs = 'hello', signaldf = df_ccbar)
        plotter5 = Plotter(isSigvar='xipipi_xic_isSignal', mcdfs={5 : df_ccbar}, signaldf= df_ccbar)
        plotter7 = Plotter(isSigvar='xipipi_xic_isSignal', mcdfs={'label': 5}, signaldf=df_ccbar)
        plotter8 = Plotter(isSigvar='xipipi_xic_isSignal', mcdfs={'label': 'hello'}, signaldf=df_ccbar)

        # Test signaldf
        plotter9 = Plotter(isSigvar='xipipi_xic_isSignal', mcdfs={'ccbar' : df_ccbar}, signaldf = 5)
        plotter10 = Plotter(isSigvar='xipipi_xic_isSignal', mcdfs={'ccbar' : df_ccbar}, signaldf = 'hello')
        