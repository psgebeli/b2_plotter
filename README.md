# `b2_plotter`

## Description
The `b2_plotter` package consists of the class Plotter(), whereby several plots based on numpy arrays and pandas dataframes can be created easily from interpretation of Belle II ROOT mDST files. 

## Releases

### Version 3.0.8 
- Bugfix data plotting

### Version 3.0.7
- Add signalregion and massvar to constructor call and remove it as a parameter from plotFom, getpurity, getsigeff
- Split up plot into two functions, plotMC and plotData
    - Add blinding functionality to plotData function

### Version 3.0.6
- Add blinding to invariant mass distributions through introduction of addBlinding parameter to initialization

### Version 3.0.5
- Update example.ipynb

### Version 3.0.4
- Bugfix missing call to colors parameter

### Version 3.0.3
- Make colors() for plot an argument with default for all MC types, so that
it can be overridden if needed

### Version 3.0.2
- Bugfix additional root_pandas call


### Version 3.0.1
- Remove deprecated root_pandas library, replace with uproot
- Update main() with a working example usage for xi03pi
- Add unit tests

### Version 2.0.9 
- 2.0.8 bugfix (missing comma)

### Version 2.0.8
- Remove scaled signal feature -- will work on this on a later release
- Specify the colors of plot() histograms since the signal was defaulting to grey

### Version 2.0.7 
- Bugfix to version 2.0.6, where scaled signal was not properly showing

### Version 2.0.6
- Change functions to return plots instead of displaying them, and remove interactive feature
- Add a scaled signal histogram to the plot function
- Update dynamic range calculation in plotFom for isGreaterThan = False cuts to prevent division by zero errors

### Version 2.0.5
- Add `cuts` argument to `plotFom`

### Version 2.0.4
- Last push didnt include changes, so include bugfix
- Add dist/, .pytest_cache, and .ipynb_checkpoints to git ignore 

### Version 2.0.3
- Bugfix call to numpy.histogram() (was called np.histogram)

### Version 2.0.2
- Implement unit tests for constructor TypeError raises

### Version 2.0.1
- Implement data plotting in plot()
- Implement unit tests for functions in non-interactive mode

### Version 1.0.9
- Set myrange = () in plotFom to enable the dynamic range calculation from v.1.0.5 if the user does not provide a range.

### Version 1.0.8
- Implement dynamic range calculation in plotStep 
- Bugfix labels not showing up on plotFom for myrange = () (default)
- Bugfix get_sigeff (missing period)

### Version 1.0.7
- Bugfix (remove () from numpy.size, as size is a property, not a method)

### Version 1.0.6
- Remove deprecated var total_bkg from plotFom()

### Version 1.0.5 
- Create example files for interactive session
- Change massvar name to be consistent (massVar in get_purity but massvar in FOM, etc)
- Implement dynamic range calculation in plotFom

### Version 1.0.4
- Bugfix main and implement a for loop to plot multiple variables
- Create usage details in README.md

### Version 1.0.3 
- Change readme such that recent versions are at the top 
- Remove `plotBias` function, as it is very niche
- Remove unused parameter `myrange` in get_sigeff and get_purity
- Add main() function 
- Add interactivity boolean in constructor call to decide whether or not to save plots to a .png or show them directly

### Version 1.0.2 
- Edit README.md
- Remove tests/ directory (unit tests are not particularly useful for this package)
- Rename directories so imports are more intuitive

### Version 1.0.1
- Add unit tests
- Switch backend to hatchling
- Move metadeta to pyproject.toml and clean it up
- Remove brackets from LICENSE.txt
- Create tests/ directory

### Version 1.0.0
Initial launch.


