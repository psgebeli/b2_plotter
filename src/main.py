
# Preamble
import numpy as np 
import matplotlib.pyplot as plt 
import root_pandas 
import pandas as pd 

class Plotter():

    def __init__(self, isSigvar: str, mcdfs: dict, signaldf: pd.DataFrame, datadf: pd.DataFrame = None):
        
        '''
        Initialize a plotter object upon constructor call.

        :param isSigvar: name of isSignal variable 
        :type isSigvar: str
        :param mcdfs: Monte carlo dataframes constructed with root_pandas
        :type mcdfs: dict (key: label, value: df)
        :param signaldf: Monte carlo dataframe to be treated as signal
        :type signaldf: pandas dataframe
        :param datadf: Data dataframe constructed with root_pandas
        :type datadf: pandas dataframe
        

        :raise TypeError: If any parameters dont match expected type
        '''
        
        # Error checking. Check type of each parameter. If it doesnt match expectations, raise a typeerror
        if isinstance(mcdfs, dict):
            for label, df in mcdfs.items():
                if not isinstance(df, pd.DataFrame):
                    raise TypeError(f'The value associated with the key "{label}" is not a pandas DataFrame')
            self.mcdfs = mcdfs
        else:
            raise TypeError('Mcdfs is not a dictionary')
        
        if isinstance(signaldf, pd.DataFrame):
            self.signaldf = signaldf 
        else:
            raise ValueError('Signal df is not a dataframe')
        
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
        
    def plot(self, var, cuts, myrange, nbins = 100, isLog = False, xlabel = '', scale = 1, bgscale = 1):

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
        :param label: Label on x-axis 
        :type label: str (usually raw str)'''

        # Set up matplotlib plot 
        ax = plt.subplot()

        # Set up empty dict of numpy arrays
        nps = {}

        for label, df in self.mcdfs.items():
            nps[label] = df.query(cuts + f'and {self.isSigvar} != 1')[var].to_numpy()
        nps['signal'] = self.signaldf.query(cuts + f'and {self.isSigvar} == 1')[var].to_numpy()

        if self.datadf is not None:
            nps['data'] = self.datadf.query(cuts)[var].to_numpy()
        
        # Set up empty dict of weights 
        wnps = {}

        for label, np in nps.items():
            if label != 'signal':
                wnps[label] = [bgscale] * len(np)
            else:
                wnps['signal'] = [scale] * len(np)
        
        # Create stacked matplotlib histogram
        ax.hist(list(nps.values()), bins = nbins, range = myrange,
                label = list(nps.keys()),
                weights = list(wnps.values()),
                stacked = True)
        
        # Plot features 

        plt.yscale('log') if isLog else plt.yscale('linear')

        plt.xlim(myrange)
        plt.ylabel('Number of Events')

        plt.xlabel(var) if xlabel == '' else plt.xlabel(xlabel)

        plt.legend()
        plt.show()

    def plotFom(self, var, massvar, myrange, signalregion, isGreaterThan = True, nbins = 100, xlabel = ''):

        '''Function to plot the figure of merit for cuts on a particular variable,
        where FOM = sqrt[signalevents/(signalevents + bkgevents)]. The maximum
        of the FOM curve is the cut which removes the most background while keeping
        the most signal. Purity (how much of signal region is signal) and signal
        efficiency (what % of signal is removed) are also included. 

        :param var: The variable to be cut
        :type var: str
        :param massvar: The name of the variable for your particle's mass
        :type massvar: str
        :param myrange: The range over which cuts should be applied
        :type myrange: tuple 
        :param signalregion: The signal region for the mass of your particular
        particle -- used to calculate purity and signal efficiency.
        :type signalregion: tuple 
        :param isGreaterThan: Expresses whether to apply testcuts where var > value, or greater than cuts
        :type isGreaterThan: bool
        :param nbins: The number of bins 
        :type nbins: int 
        :param xlabel: Label for the x-axis
        :type xlabel: str'''

        # Create a background dataframe as the concatenation of all of the individual monte carlo dataframes
        df_bkg = pd.concat(self.mcdfs)

        # Store the total amount of background/sig events in the signal region by the size of their numpy array
        total_bkg = df_bkg.query(f'{signalregion[0]} < {massvar} < {signalregion[1]} and {self.isSigvar} != 1')[var].to_numpy().size
        total_sig = self.signaldf.query(f'{signalregion[0]} < {massvar} < {signalregion[1]} and {self.isSigvar} == 1')[var].to_numpy().size

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
                globalcuts = f'{signalregion[0]} < {massvar} < {signalregion[1]} and {var} > {testcut}'
            else:
                globalcuts = f'{signalregion[0]} < {massvar} < {signalregion[1]} and {var} < {testcut}'
            
            # Append the size of the array after the constraints to the globalsig/globalbkg lists respectively
            globalsig.append(self.signaldf.query(f'{globalcuts} and {self.isSigvar} == 1')[var].to_numpy().size)
            globalbkg.append(df_bkg.query(f'{globalcuts} and {self.isSigvar} != 1')[var].to_numpy().size)

            # Calculate the figure of merit for this bin and append it to fom list
            fom.append(globalsig[bin] / np.sqrt(globalsig[bin] + globalbkg[bin]))

        
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
        axes[0].set_xlabel(xlabel)
        ax.grid()
        
        # Show the plot
        plt.show()
        
        # Convert fom to a numpy array for easier manipulation
        fom = np.array(fom)

        # Find the index of the maximum value in the fom array
        max_fom_index = np.argmax(fom)

        # Get the corresponding test cut value at the maximum FOM
        optimal_cut = testcuts[max_fom_index]

        # Print the optimal cut and corresponding maximum FOM
        print("Optimal Cut:", optimal_cut)
        print("Maximum Figure of Merit:", fom[max_fom_index])

    def plotStep(self, var, cuts, myrange, nbins = 100, xlabel = ''):

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

        # Create the histogram
        ax.hist([npbkg, npsig], bins = nbins, range = myrange, label =  ['bkg', 'signal'], histtype = 'step', stacked = False)

        # Set plot features 
        plt.yscale('log')
        plt.xlim(myrange)
        plt.xlabel(var) if xlabel == '' else plt.xlabel(xlabel)
        
        # Create a legend and show both of the plots
        plt.legend()
        plt.show()
        plt.show()

    def plotBias(self, compVar, lifetimeVar, massvar, myrange, signalregion, step, refLifetime, xlabel = ''):

        '''Function to plot a variable, var1, binned in terms of another variable, var2. This shows the correlation between the two variables
        and if any constraints on one variable will bias the other. The most common use case of this is to determine if constraints to variables
        introduce a bias to the particle lifetime.

        :param compVar: The variable that may introduce bias
        :type compVar: str
        :param lifetimeVar: The name of lifetime variable, usually mcFlightTime
        :type lifetimeVar: str
        :param massvar: The name of the variable for your particle's mass
        :type massvar: str
        :param myrange: Range on x-axis
        :type myrange: tuple
        :param signalregion: The signal region for the mass of your particular particle.
        :type signalregion: tuple 
        :param step: The x-axis width between bins
        :type step: float
        :param refLifetime: Average lifetime of particle in picoseconds
        :type refLifetime: float
        :param xlabel: Label on x-axis 
        :type xlabel: str (usually raw str)'''

        # Setup plot
        ax = plt.subplot()

        # Partition the bins
        mybins = np.arange(myrange[0], myrange[1], step)
        probs = np.arange(myrange[0], myrange[1] + step, step)

        # Create empty list for the mean, error bars, and widths
        means, errs, sigmas = [], [], []

        for i in range(0, len(probs) - 1):

            # Create a signal numpy array
            npsig = self.signaldf.query(f'{signalregion[0]} < {massvar} < {signalregion[1]} and {self.isSigvar} == 1 and {probs[i]} < {compVar} < {probs[i+1]}')[lifetimeVar].to_numpy() * 1000

            # Append the mean, error bar, and width
            means.append(np.average(npsig))
            errs.append(np.std(npsig) / np.sqrt(npsig.size - 1))
            sigmas.append((refLifetime - np.average(npsig)) * np.sqrt(npsig.size - 1) / np.std(npsig))

        # Create error bar
        ax.errorbar(mybins, means, yerr = errs, fmt = 'kp')

        # Perform constant chi^2 fit 
        p2, res2, _, _, _ = np.polyfit(mybins, means, deg = 0, full = True)
        p2_fn = np.poly1d(p2)
        p2bins = np.polyval(p2, mybins)
        cs2 = np.sum(np.square(p2bins - means) / np.square(errs))
        plt.plot(mybins, p2_fn(mybins), color = 'red', label = f'constant fit $\chi^2$/dof = {round(cs2, 5)}, dof = {len(means) - 1}')

        # Perform linear chi^2 fit
        p1, res, _, _, _ = np.polyfit(mybins, means, deg = 1, full = True)
        p1_fn = np.poly1d(p1)
        p1bins = np.polyval(p1, mybins)
        cs1 = np.sum(np.square(p1bins - means) / np.square(errs))
        plt.plot(mybins, p1_fn(mybins), color = 'red', label = f'linear fit $\chi^2$/dof = {round(cs1, 5)}, dof = {len(means) - 2}')

        print(cs1)
        print(cs2)

        # Set plot features
        plt.ylim(refLifetime - (refLifetime / 2), refLifetime + (2 * refLifetime))
        plt.legend()
        plt.xlabel(compVar) if xlabel == '' else plt.xlabel(xlabel)
        plt.ylabel('Lifetime (ps)')
        plt.show()

    def get_purity(self, cuts, massVar, myrange, signalregion):
        
        '''Function to return the purity, % of signal in signal region
        
        :param cuts: All cuts to be applied to the dataframes before plotting
        :type cuts: str
        :param massvar: The name of the variable for your particle's mass
        :type massvar: str
        :param myrange: The range over which cuts should be applied
        :type myrange: tuple 
        :param signalregion: The signal region for the mass of your particular particle.
        :type signalregion: tuple'''

        df_bkg = pd.concat(self.mcdfs)

        npsig = self.signaldf.query(f'{signalregion[0]} < {massVar} < {signalregion[1]} and {cuts} and {self.isSigvar} == 1').to_numpy()
        npbkg = df_bkg.query(f'{signalregion[0]} < {massVar} < {signalregion[1]} and {cuts} and {self.isSigvar} != 1')[massVar].to_numpy()

        sig_events, bkg_events = len(npsig), len(npbkg)
        total_events = sig_events + bkg_events

        return sig_events / total_events * 100
    
    def get_sigeff(self, cuts, massVar, myrange, signalregion):

        '''Function to return the sigeff, % of signal lost from applying cuts
        
        :param cuts: All cuts to be applied to the dataframes before plotting
        :type cuts: str
        :param massvar: The name of the variable for your particle's mass
        :type massvar: str
        :param myrange: The range over which cuts should be applied
        :type myrange: tuple 
        :param signalregion: The signal region for the mass of your particular particle.
        :type signalregion: tuple'''

        sig_before = len(self.signaldf.query(f'{signalregion[0]} < {massVar} < {signalregion[1]} and {self.isSigvar} == 1'))
        sig_after = len(self.signaldf.query(f'{signalregion[0]} < {massVar} < {signalregion[1]} and {cuts} and {self.isSigvar} == 1'))

        return sig_after / sig_before * 100
        



        



