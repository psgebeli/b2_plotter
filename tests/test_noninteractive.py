import pytest 
from b2_plotter.Plotter import Plotter
import root_pandas as rp
import os

ccbar = '/belle2work/psgebeli/samples/gmc/mc15rib/xipipi/ccbar.root'
mycols= ['xic_M', 'xic_significanceOfDistance','xi_significanceOfDistance', 
         'lambda0_p_protonID', 'xi_M', 'xic_mcFlightTime', 'xic_chiProb', 'xic_isSignal']
xicmassrangeloose = '2.3 < xic_M < 2.65'
df_ccbar = rp.read_root(ccbar, key='xic_tree', columns = mycols)

plotter = Plotter(isSigvar='xic_isSignal', mcdfs={'ccbar': df_ccbar}, signaldf = df_ccbar, interactive = False)

def test_constructor():
    assert isinstance(plotter, Plotter)

def test_plot():
    for var in mycols[:-3]:
        plotter.plot(var, cuts = xicmassrangeloose)
        assert os.path.isfile(f'plot_{var}.png')

def test_plotFom():
    for var in mycols[:-3]:
        optimal, fommax = plotter.plotFom(var, massvar = 'xic_M', signalregion = (2.46, 2.475))
        assert os.path.isfile(f'fom_{var}.png')
        assert isinstance(optimal, float) and isinstance(fommax, float)

def test_plotStep():
    for var in mycols[:-3]:
        plotter.plotStep(var, cuts = xicmassrangeloose)
        assert os.path.isfile(f'step_{var}.png')

def test_getpurity():
    assert isinstance(plotter.get_purity(xicmassrangeloose, 'xic_M', (2.46, 2.475)), float)

def test_getsigeff():
    assert isinstance(plotter.get_sigeff(xicmassrangeloose, 'xic_M', (2.46, 2.475)), float)

def test_errors():
    with pytest.raises(TypeError):

        # Test isSigvar type errors
        plotter1 = Plotter(isSigvar=5, mcdfs={'ccbar': df_ccbar}, signaldf = df_ccbar, interactive = False)
        plotter2 = Plotter(isSigvar=xic_M, mcdfs={'ccbar': df_ccbar}, signaldf = df_ccbar, interactive = False)

        # Test mcdfs 
        plotter3 = Plotter(isSigvar = 'xic_M', mcdfs = 5, signaldf = df_ccbar, interactive = False)
        plotter4 = Plotter(isSigvar = 'xic_M', mcdfs = 'hello', signaldf = df_ccbar, interactive = False)
        plotter5 = Plotter(isSigvar='xic_M', mcdfs={5 : df_ccbar}, signaldf= df_ccbar, interactive=False)
        plotter6 = Plotter(isSigvar='xic_M', mcdfs={'hello' : df_ccbar}, signaldf= df_ccbar, interactive=False)
        plotter7 = Plotter(isSigvar='xic_M', mcdfs={'label': 5}, signaldf=df_ccbar, interactive=False)
        plotter8 = Plotter(isSigvar='xic_M', mcdfs={'label': 'hello'}, signaldf=df_ccbar, interactive=False)

        # Test signaldf
        plotter9 = Plotter(isSigvar='xic_M', mcdfs={'ccbar' : df_ccbar}, signaldf = 5, interactive=True)
        plotter10 = Plotter(isSigvar='xic_M', mcdfs={'ccbar' : df_ccbar}, signaldf = 'hello', interactive=True)
        
        # Test interactive
        plotter11 = Plotter(isSigvar='xic_M', mcdfs={'ccbar' : df_ccbar}, signaldf = df_ccbar, interactive = 5)
        plotter12 = Plotter(isSigvar='xic_M', mcdfs={'ccbar' : df_ccbar}, signaldf = df_ccbar, interactive = 'hello')