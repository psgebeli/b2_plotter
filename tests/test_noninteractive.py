import pytest 
from b2_plotter.Plotter import Plotter
import root_pandas as rp

def test_constructor():

    ccbar = '/belle2work/psgebeli/samples/gmc/mc15rib/xipipi/ccbar.root'
    mycols= ['xic_M','xic_isSignal', 'xic_significanceOfDistance','xi_significanceOfDistance', 
         'lambda0_p_protonID', 'xi_M', 'xic_mcFlightTime', 'xic_chiProb']
    
    df_ccbar = rp.read_root(ccbar, key='xic_tree', columns = mycols)

    assert Plotter(isSigvar='hi', mcdfs={'hi': 5, 'hii' : 6}, signaldf = df_ccbar, )

