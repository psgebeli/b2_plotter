# `b2_plotter`

## Description
The `b2_plotter` package consists of the class Plotter(), whereby several plots based on numpy arrays and pandas dataframes can be created easily. 

## Releases

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


