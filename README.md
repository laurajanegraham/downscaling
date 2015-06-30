# Downscaling species occurrence records

This code takes species records at the grid cell level and downscales to the patch level. To run, the following input data are required:
- `grids` shapefile (this is in the same resolution as the species data you wish to downscale)
- `habitat` shapefile for each species (this will be patches of suitable habitat for the species)
- `survey.csv` for each species (presence absence for each grid cell)
- `species.info.csv` (at a minimum should contain the value for alpha for each species - 1/mean dispersal distance)

Edit `runDownscaling.R` following the comments and then run. The `runDownscaling()` function will then be available. This takes as input the species of interest and the number of reps (how many time you want to run the downscaling functions). The output of this function is a patch by downscaling method matrix showing the probability a patch is predicted to be occupied for each downscaling method.  
