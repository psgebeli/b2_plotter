import pytest 
from b2_plotter.Plotter import Plotter
import root_pandas as rp
import os

ccbar = '/belle2work/psgebeli/samples/gmc/mc15rib/xipipi/ccbar.root'
mycols= ['xic_M', 'xic_significanceOfDistance','xi_significanceOfDistance', 
         'lambda0_p_protonID', 'xi_M', 'xic_mcFlightTime', 'xic_chiProb', 'xic_isSignal']
xicmassrangeloose = '2.3 < xic_M < 2.65'
df_ccbar = rp.read_root(ccbar, key='xic_tree', columns = mycols)

def test_constructor():

    

    assert Plotter(isSigvar='xic_isSignal', mcdfs={'ccbar': df_ccbar}, signaldf = df_ccbar, interactive = False)

def test_plot():
    
    for var in mycols[:-3]:
        Plotter.plot(var, cuts = xicmassrangeloose)
        assert os.path.isfile(f'xic_{var}.png')