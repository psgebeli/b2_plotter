
# Preamble
import numpy
import matplotlib.pyplot as plt 
import uproot as up 
import pandas as pd
import argparse
import os
import csv

class Plotter():

    def __init__(self, isSigvar: str, mcdfs: dict, signaldf: pd.DataFrame, massvar: str, signalregion: tuple,
                 datadf: pd.DataFrame = None):
        
        '''
        Initialize a plotter object upon constructor call.

        :param isSigvar: name of isSignal variable 
        :type isSigvar: str
        :param mcdfs: Monte carlo dataframes constructed with root_pandas
        :type mcdfs: dict (key: label, value: df)
        :param signaldf: Monte carlo dataframe to be treated as signal
        :type signaldf: pandas dataframe
        :param massvar: Name of primary mass variable
        :type massvar: str
        :param signalregion: Signal region of primary mass variable
        :type signalregion: tuple

        :raise TypeError: If any parameters dont match expected type
        '''
        
        # Error checking. Check type of each parameter. If it doesnt match expectations, raise a typeerror
        if isinstance(mcdfs, dict):
            for label, df in mcdfs.items():
                if not isinstance(label, str):
                    raise TypeError(f'The key associated with the value {df} is not a str.')
                if not isinstance(df, pd.DataFrame):
                    raise TypeError(f'The value associated with the key "{label}" is not a pandas DataFrame')
            self.mcdfs = mcdfs
        else:
            raise TypeError('Mcdfs is not a dictionary')
        
        if isinstance(signaldf, pd.DataFrame):
            self.signaldf = signaldf 
        else:
            raise TypeError('Signal df is not a dataframe')
        
        if datadf is not None:
            if isinstance(datadf, pd.DataFrame):
                self.datadf = datadf
            else:
                raise TypeError('Datadf is not a pandas DataFrame.')
        else:
            self.datadf = None
        
        if isinstance(isSigvar, str):
            self.isSigvar = isSigvar
        else:
            raise TypeError('isSigvar is not a string.')
        
        if isinstance(massvar, str):
            self.massvar = massvar
        else:
            raise TypeError('massvar is not a string.')
        
        if isinstance(signalregion, tuple):
            self.signalregion = signalregion
        else:
            raise TypeError('signalregion is not a tuple.')
        
        
    def plotMC(self, var, cuts, myrange = (), nbins = 100, isLog = False, xlabel = '', scale = 1, bgscale = 1, 
               color = ['b', '#ffa500', 'g', 'r', 'c', 'y', '#a52a2a', 'm' ]):

        '''Create a matplotlib stacked histogram of a variable over a certain range using only Monte Carlo.

        :param var: The variable to be cut
        :type var: str
        :param cuts: All cuts to be applied to the dataframes before plotting
        :type cuts: str
        :param myrange: Range on x-axis
        :type myrange: tuple 
        :param nbins: Number of bins 
        :type nbins: int 
        :param isLog: Whether or not the plot should be on a logarithmic scale 
        :type isLog: bool
        :param xlabel: Label on x-axis 
        :type xlabel: str (usually raw str)
        :param scale: Factor by which to scale the signal
        :type scale: Float
        :param bgscale: Factor by which to scale the background
        :type bgscale: Float
        :param color: List of colors to apply to each stack of the histogram
        :param color: List'''

        # Set up matplotlib plot 
        ax = plt.subplot()

        # Set up empty dict of MC numpy arrays
        mcnps = {}

        # For each key/value pair in mcdfs, create a new entry in mcnps of 'label' : numpy array.
        # Also create an entry for the signal.
        for label, df in self.mcdfs.items():
            mcnps[label] = df.query(cuts + f'and {self.isSigvar} != 1')[var].to_numpy()
        mcnps['signal'] = self.signaldf.query(cuts + f'and {self.isSigvar} == 1')[var].to_numpy()
        
        # Set up empty dict of weights 
        wnps = {}

        # For each label : np array in mcnps, create a weight of the corresponding scale times the length of the array
        for label, np in mcnps.items():
            if label != 'signal':
                wnps[label] = [bgscale] * len(np)
            else:
                wnps['signal'] = [scale] * len(np)

        if myrange == ():
            # Calculate the dynamic range for the variable based on the data within the specified cuts
            all_mc = numpy.concatenate(list(mcnps.values()))
            myrange = (numpy.min(all_mc), numpy.max(all_mc))
        
        # Create stacked matplotlib histogram
        ax.hist(list(mcnps.values()), bins = nbins, range = myrange,
                label = list(mcnps.keys()),
                weights = list(wnps.values()),
                stacked = True,
                color = color)
        
        # Plot features 
        plt.yscale('log') if isLog else plt.yscale('linear')
        plt.xlim(myrange)
        plt.ylabel('Number of Events')
        plt.xlabel(var) if xlabel == '' else plt.xlabel(xlabel)
        plt.legend()

        return plt

    def plotData(self, var, cuts, myrange = (), nbins = 100, isLog = False, xlabel = '', scale = 1, 
             bgscale = 1, color = ['b', '#ffa500', 'g', 'r', 'c', 'y', '#a52a2a', 'm' ], addBlinding = True):

        '''Create a matplotlib stacked histogram of a variable over a certain range.
        If datadf is provided to constructor, data will be stacked on top of MC.

        :param var: The variable to be cut
        :type var: str
        :param cuts: All cuts to be applied to the dataframes before plotting
        :type cuts: str
        :param myrange: Range on x-axis
        :type myrange: tuple 
        :param nbins: Number of bins 
        :type nbins: int 
        :param isLog: Whether or not the plot should be on a logarithmic scale 
        :type isLog: bool
        :param xlabel: Label on x-axis 
        :type xlabel: str (usually raw str)
        :param scale: Factor by which to scale the signal
        :type scale: Float
        :param bgscale: Factor by which to scale the background
        :type bgscale: Float
        :param color: List of colors to apply to each stack of the histogram
        :param color: List
        :param addBlinding: Add blinding to signal region?
        :param addBlinding: bool'''

        # Set up matplotlib plot 
        ax = plt.subplot()

        # Set up empty dict of MC numpy arrays
        mcnps = {}

        # For each key/value pair in mcdfs, create a new entry in mcnps of 'label' : numpy array.
        # Also create an entry for the signal.
        for label, df in self.mcdfs.items():
            mcnps[label] = df.query(cuts + f'and {self.isSigvar} != 1')[var].to_numpy()
        mcnps['signal'] = self.signaldf.query(cuts + f'and {self.isSigvar} == 1')[var].to_numpy()

        # Create data arrays for outside of signal region if blinding is enabled and mass is being plotted, otherwise make one array
        if addBlinding and var == self.massvar:
            npdata_less = self.datadf.query(f'{cuts} and {self.massvar} < {self.signalregion[0]}')[self.massvar].to_numpy()
            npdata_greater = self.datadf.query(f'{cuts} and {self.massvar} > {self.signalregion[1]}')[self.massvar].to_numpy()
        else:
            npdata = self.datadf.query(cuts)[var].to_numpy()
        
        # Set up empty dict of weights 
        wnps = {}

        # For each label : np array in mcnps, create a weight of the corresponding scale times the length of the array
        for label, np in mcnps.items():
            if label != 'signal':
                wnps[label] = [bgscale] * len(np)
            else:
                wnps['signal'] = [scale] * len(np)

        if myrange == ():
            # Calculate the dynamic range for the variable based on the data within the specified cuts
            all_data = numpy.concatenate(list(mcnps.values()))
            myrange = (numpy.min(all_data), numpy.max(all_data))
        
        # Plot data/mc using conditional for blinding
        
        ax.hist(list(mcnps.values()), bins = nbins, range = myrange,
                label = list(mcnps.keys()),
                weights = list(wnps.values()),
                stacked = True, 
                color = color)

        if addBlinding:
            ydata_less, bin_edges_less = numpy.histogram(npdata_less, bins=nbins, range=myrange)
            ydata_greater, bin_edges_greater = numpy.histogram(npdata_greater, bins=nbins, range=myrange)
            
            bin_centers_less = 0.5 * (bin_edges_less[1:] + bin_edges_less[:-1])
            bin_centers_greater = 0.5 * (bin_edges_greater[1:] + bin_edges_greater[:-1])
            
            ax.errorbar(bin_centers_less, ydata_less, yerr=ydata_less**0.5, fmt='ko', label="Data")
            ax.errorbar(bin_centers_greater, ydata_greater, yerr=ydata_greater**0.5, fmt='ko', label="Data")
        else:
            ydata, bin_edges = numpy.histogram(npdata, bins=nbins, range=myrange)
            bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
            ax.errorbar(bin_centers, ydata, yerr=ydata**0.5, fmt='ko', label="Data")

            
        # Plot features 
        plt.yscale('log') if isLog else plt.yscale('linear')
        plt.xlim(myrange)
        plt.ylabel('Number of Events')
        plt.xlabel(var) if xlabel == '' else plt.xlabel(xlabel)
        plt.legend()

        return plt
    


    def plotFom(self, var, cuts, myrange = (), isGreaterThan = True, nbins = 100, xlabel = ''):

        '''Function to plot the figure of merit for cuts on a particular variable,
        where FOM = sqrt[signalevents/(signalevents + bkgevents)]. The maximum
        of the FOM curve is the cut which removes the most background while keeping
        the most signal. Purity (how much of signal region is signal) and signal
        efficiency (what % of signal is removed) are also included. 

        :param var: The variable to be cut
        :type var: str
        :param cuts: Cuts to be applied before the FOM is generated
        :type cuts: str
        :param myrange: The range over which cuts should be applied
        :type myrange: tuple 
        :param isGreaterThan: Expresses whether to apply testcuts where var > value, or greater than cuts
        :type isGreaterThan: bool
        :param nbins: The number of bins 
        :type nbins: int 
        :param xlabel: Label for the x-axis
        :type xlabel: str'''

        # Create a background dataframe as the concatenation of all of the individual monte carlo dataframes
        df_bkg = pd.concat(self.mcdfs)

        # Store the total signal and background as numpy arrays
        np_bkg = df_bkg.query(f'{cuts} and {self.signalregion[0]} < {self.massvar} < {self.signalregion[1]} and {self.isSigvar} != 1')[var].to_numpy()
        np_sig = self.signaldf.query(f'{cuts} and {self.signalregion[0]} < {self.massvar} < {self.signalregion[1]} and {self.isSigvar} == 1')[var].to_numpy()


        # Store the total amount of sig events in the signal region by the size of the numpy array
        total_sig = np_sig.size

        if myrange == () and isGreaterThan:
            # Calculate the dynamic range for the variable based on the data within the specified cuts
            myrange = (numpy.min(np_sig), numpy.max(np_sig))
        elif myrange == () and not isGreaterThan:
            myrange = (numpy.min(np_sig) + (numpy.max(np_sig) / 10), numpy.max(np_sig))

        # Initialize empty lists for testcuts, number of background/sig events after some global cut has been applied, and the figure of merit values
        testcuts, globalsig, globalbkg, fom = [], [], [], []

        # Define some interval to iterate over based on the range and the number of bins.
        interval = (myrange[1] - myrange[0]) / nbins

        # For each bin,
        for bin in range(0, nbins):

            # Define a test cut as (interval * bin) + x_min and append it to testcuts
            testcut = interval * bin + myrange[0]
            testcuts.append(testcut)

            # If the paramater isGreaterThan is True, the global cut is var > value. Otherwise, its var < value.
            # The global cuts is a string of this cut as well as constraining the mass to the signal region.
            if isGreaterThan:
                globalcuts = f'{cuts} and {self.signalregion[0]} < {self.massvar} < {self.signalregion[1]} and {var} > {testcut}'
            else:
                globalcuts = f'{cuts} and {self.signalregion[0]} < {self.massvar} < {self.signalregion[1]} and {var} < {testcut}'
            
            # Append the size of the array after the constraints to the globalsig/globalbkg lists respectively
            globalsig.append(self.signaldf.query(f'{globalcuts} and {self.isSigvar} == 1')[var].to_numpy().size)
            globalbkg.append(df_bkg.query(f'{globalcuts} and {self.isSigvar} != 1')[var].to_numpy().size)

            # Calculate the figure of merit for this bin and append it to fom list
            fom.append(globalsig[bin] / numpy.sqrt(globalsig[bin] + globalbkg[bin]))

        
        # Setup the figure of merit plot
        fig, ax = plt.subplots()

        # Twin the x-axis twice to make 2 independent y-axes and make some extra space for them.
        axes = [ax, ax.twinx(), ax.twinx()]
        fig.subplots_adjust(right=0.75)
        fig.subplots_adjust(right=0.75)

        # Move the last y-axis spine over to the right by 20% of the width of the axes
        axes[-1].spines['right'].set_position(('axes', 1.2))

        # To make the border of the right-most axis visible, we need to turn the frame
        # on. This hides the other plots, however, so we need to turn its fill off.
        axes[-1].set_frame_on(True)
        axes[-1].patch.set_visible(False)

        # Initialize empty lists for signal efficiency and purity and append values for each bin to the lists.
        sigeff = []
        purity = []
        for bin in range(0, (nbins - 1)):
            sigeff.append(globalsig[bin]/total_sig)
            purity.append(globalsig[bin]/(globalbkg[bin]+globalsig[bin]))
        
        # Append the signal efficiency and purity of the final bin again so the curves flatten out.
        sigeff.append(sigeff[nbins - 2])
        purity.append(purity[nbins - 2])

        # Plot the curves on their respective axes and label them.
        axes[0].plot(testcuts, fom, color='Red')
        axes[0].set_ylabel('Figure of merit', color='Red')
        axes[1].plot(testcuts, sigeff, color='Blue')
        axes[1].set_ylabel('Signal efficiency', color='Blue')
        axes[2].plot(testcuts, purity, color='Green')
        axes[2].set_ylabel('Purity', color='Green')

        # Label the x-axis according to the input parameter and add grid lines to the axis plot
        if xlabel == '' and isGreaterThan:
            axes[0].set_xlabel(f'{var} > ...')
        elif xlabel == '' and not isGreaterThan:
            axes[0].set_xlabel(f'{var} < ...')
        else:
            axes[0].set_xlabel(xlabel)
        ax.grid()
        
        # Convert fom to a numpy array for easier manipulation
        fom = numpy.array(fom)

        # Find the index of the maximum value in the fom array
        max_fom_index = numpy.argmax(fom)

        # Get the corresponding test cut value at the maximum FOM
        optimal_cut = testcuts[max_fom_index]

        return plt, optimal_cut

    def plotStep(self, var, cuts, myrange = (), nbins = 100, xlabel = ''):

        '''Function to plot an unstacked step histogram for a variable, which is useful in cases where you do not 
        need to see the individual background types and/or the signal is hidden underneath a sea of background.
        
        :param var: The variable to be cut
        :type var: str
        :param cuts: All cuts to be applied to the dataframes before plotting
        :type cuts: str
        :param myrange: Range on x-axis
        :type myrange: tuple 
        :param nbins: Number of bins 
        :type nbins: int 
        :param xlabel: Label on x-axis 
        :type xlabel: str (usually raw str)'''

        # Setup plot
        ax = plt.subplot()

        df_bkg = pd.concat(self.mcdfs)

        # Define bkg/true numpy arrays
        npbkg = df_bkg.query(f'{cuts} and {self.isSigvar} != 1')[var].to_numpy()
        npsig = self.signaldf.query(f'{cuts} and {self.isSigvar} == 1')[var].to_numpy()

        if myrange == ():
            # Calculate the dynamic range for the variable based on the data within the specified cuts
            myrange = (numpy.min(npbkg), numpy.max(npbkg))

        # Create the histogram
        ax.hist([npbkg, npsig], bins = nbins, range = myrange, label =  ['bkg', 'signal'], histtype = 'step', stacked = False)

        # Set plot features 
        plt.yscale('log')
        plt.xlim(myrange)
        plt.xlabel(var) if xlabel == '' else plt.xlabel(xlabel)
        
        # Create a legend and show plot
        plt.legend()
        
        return plt

    def get_purity(self, cuts):
        
        '''Function to return the purity, % of signal in signal region
        
        :param cuts: All cuts to be applied to the dataframes before plotting
        :type cuts: str'''

        df_bkg = pd.concat(self.mcdfs)

        npsig = self.signaldf.query(f'{self.signalregion[0]} < {self.massvar} < {self.signalregion[1]} and {cuts} and {self.isSigvar} == 1')[self.massvar].to_numpy()
        npbkg = df_bkg.query(f'{self.signalregion[0]} < {self.massvar} < {self.signalregion[1]} and {cuts} and {self.isSigvar} != 1')[self.massvar].to_numpy()

        sig_events, bkg_events = len(npsig), len(npbkg)
        total_events = sig_events + bkg_events

        return sig_events / total_events * 100
    
    def get_sigeff(self, cuts):

        '''Function to return the sigeff, % of signal lost from applying cuts
        
        :param cuts: All cuts to be applied to the dataframes before plotting
        :type cuts: str'''

        sig_before = len(self.signaldf.query(f'{self.signalregion[0]} < {self.massvar} < {self.signalregion[1]} and {self.isSigvar} == 1'))
        sig_after = len(self.signaldf.query(f'{self.signalregion[0]} < {self.massvar} < {self.signalregion[1]} and {cuts} and {self.isSigvar} == 1'))

        return sig_after / sig_before * 100


# ----------------------------------------------------------------------------------------------------------------------------

# Hard coded columns
cols = ['xi03pi_xi_significanceOfDistance', 'xi03pi_lambda_significanceOfDistance',
          'xi03pi_xi_pi0_M','xi03pi_lambda_p_protonID_noSVD', 'xi03pi_xi_M', 'xi03pi_lambda_M',
          'xi03pi_gamma1_beamBackgroundSuppression', 'xi03pi_gamma2_beamBackgroundSuppression',
          'xi03pi_gamma1_fakePhotonSuppression', 'xi03pi_gamma2_fakePhotonSuppression',
          'xi03pi_xic_M']
    
# Frequently used vars 
xicmassrangeloose = '2.3 < xi03pi_xic_M < 2.65'
xicmassrangetight = '2.46 < xi03pi_xic_M < 2.475'

# Vars that could potentially be useful 
potentially_useful_vars = cols[:-1]



def main():
    
    # Parse cmd line arguments and store them in variables
    args = parse_cmd()
    mcpath, prefix = args.input, args.prefix

    # Call construct_dfs with these columns and store return value
    mcdfs = construct_dfs(mcpath, cols, prefix)
    
    # Construct a plotter object 
    plotter = Plotter(isSigvar = f'{prefix}_isSignal', mcdfs = mcdfs, signaldf = pd.concat(mcdfs.values()))

    # Initialize cuts
    cuts = xicmassrangeloose

    # Initialize csv file to store cuts 
    with open('cuts.csv', 'a') as file:

        # Create a writer object from csv library with columns variable, lower_bound, upper_bound
        writer = csv.DictWriter(file, fieldnames = ['variable', 'lower_bound', 'upper_bound'])

        # Write a header to the csv
        writer.writeheader()

        # For each variable in a predefined column slice
        for var in potentially_useful_vars:

            # Save the return values of get fom
            less_fom, less_cut = get_fom(cuts, var, prefix, plotter)[0]
            greater_fom, greater_cut = get_fom(cuts, var, prefix, plotter)[1]

            # Save the plots as pngs 
            less_fom.savefig(f'{var}_lessfom.png')
            greater_fom.savefig(f'{var}_greaterfom.png')

            # Close the plots
            less_fom.close()
            greater_fom.close()

            # Get optimal cuts from FOM and write them to the csv
            writer.writerow({'variable' : var, 'lower_bound' : greater_cut, 'upper_bound' : less_cut})

        


# Read in args from cmd line
def parse_cmd():
    
    # Create an argument parser from argparse with a usage statement
    parser = argparse.ArgumentParser(usage = 'python3 Plotter.py -i path/to/MC [-d path/to/data] -p xic_prefix_name')

    # Search the command line for arguments following these flags and provide help statement for 
    # python3 Plotter.py --help 
    parser.add_argument('-i', '--input', help = 'Relative path to directory containing all MC root files', type = str)
    parser.add_argument('-p', '--prefix', help = 'Prefix of Xic+ variables', type = str)

    # Return the parsed arguments
    return parser.parse_args()

# Construct dataframes
def construct_dfs(mcpath, mycols, prefix):
    
    # Initialize an empty dictionary to hold monte carlo dataframes
    mcdfs = {}

    # For each file in the provided MC path
    for mcfile in os.listdir(mcpath):

        # If it ends with .root
        if mcfile.endswith('.root'):

            # Create the full path including the file's name
            path = os.path.join(mcpath, mcfile)

            # Constract a pandas dataframe with that file 
            with up.open(path) as file:
                tree = file['xic_tree']
                df = tree.arrays(filter_name = mycols + [f'{prefix}_isSignal'], library = "pd")

            # Create a pair in the mcdfs dictionary of filename : df
            mcdfs[mcfile] = df

    # Return the dict of mc dfs and the single data df 
    return mcdfs         


def get_fom(cuts, var, prefix, plotter):
    return plotter.plotFom(var = var, massvar = f'{prefix}_M', signalregion = (2.46, 2.475), cuts = cuts, isGreaterThan = False), plotter.plotFom(var = var, massvar = f'{prefix}_M', signalregion = (2.46, 2.475), cuts = cuts)

if __name__ == '__main__':
    main()
