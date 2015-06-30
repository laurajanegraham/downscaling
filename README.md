# Downscaling species occurrence records

This code takes species records at the grid cell level and downscales to the patch level. To run, the following input data are required:
- `grids` shapefile (this is in the same resolution as the species data you wish to downscale)
- `species.info.csv` (at a minimum should contain the value for alpha for each species - 1/mean dispersal distance)
- A folder for each species containing:
  - `habitat` shapefile (this will be patches of suitable habitat for the species - I created mine by filtering Land Cover Map 2007 data by habitat type and minimum required patch size and dissolving boundaries)
  - `survey.csv` for each species (presence absence for each grid cell)


Edit `runDownscaling.R` following the comments and then run. The `runDownscaling()` function will then be available. This takes as input the species of interest and the number of reps (how many time you want to run the downscaling functions). The output of this function is a patch by downscaling method matrix showing the probability a patch is predicted to be occupied for each downscaling method.  
