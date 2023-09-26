# `b2_plotter`

## Description
The `b2_plotter` package consists of the class Plotter(), whereby several plots based on numpy arrays and pandas dataframes can be created easily. 

## Releases

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


