# Run downscaling. This code calls on the ds.factors.R and ds.run.R scripts.
# Species data are downscaled as follows:
# 1. Load habitat data
# 2. Load species information
# 3. Gets the factors required for downscaling
# 4. Downscales the results to the patch
runDownscaling <- function(sp.name, reps) {
  # load required packages (install them with install.packages() if not done already)
  require(rgeos)
  require(rgdal)
  
  # load the other scripts
  source("ds.run.R")
  source("ds.factors.R")
  
  # landscape data (specific to species habitat requirements)
  # this needs to be a shapefile - placed into the species folder in the data folder
  habitat <- readOGR(paste0("data/",sp.name), "current", verbose=FALSE) 
  
  # species details (minimum area, dispersal etc.)
  sp.info <- read.csv("data/species.info.csv")
  
  # the grid - this must be at the same resolution as the species records. Here it
  # is 2km square cells
  grid.poly <- readOGR("data", "grids")
  
  # survey data - must be a matrix with years on cols and grid refs on rows. Data
  # are 0 and 1
  sp.rec <- read.csv(paste("data", sp.name, "survey.csv", sep="/"), stringsAsFactors = FALSE)
  
  # assign value of alpha for the species - will be used in ds.factors() (this is
  # 1/dispersal distance)
  alpha <- sp.info$alpha[sp.info$species==sp.name]
  
  # change the ref to grid.ref to the name of your grid.ref column and same with survey
  survey.col <- data.frame(grid.ref = sp.rec$grid.ref, survey = sp.rec$survey, 
                           stringsAsFactors = FALSE)
  
  if(sum(survey.col$survey)>0) {
    grids <- survey.col[survey.col$survey==1,1]
    survey.grid <- subset(grid.poly, grid_ref %in% grids)
    ds.fac <- ds.factors(habitat, survey.grid, alpha)
  }
  
  occ <- length(survey.grid)/length(grid.poly)
  
  ds.res <- list()
  for(r in 1:reps) {
    # DOWNSCALE
    ds.survey <- ds.run(ds.fac, survey.grid, sp.name, occ = occ, grid.poly = grid.poly)
    rownames(ds.survey) <- ds.survey$lcmID
    ds.survey <- ds.survey[,-1]
    ds.res[[r]] <- ds.survey        
  }
  
  # This gives back a data frame with the probability of occupancy for each patch.
  # The rownames are the IDs from the habitat shapefile and so the results can be
  # matched back using this
  ds.res <- Reduce("+", ds.res)/reps
  return(ds.res)
}
