# ds.factors function: this function creates the factors to weight by: proportion of patch 
# falling within each grid, distance to nearest grid cell not including the grid(s) it lies in, 
# and patch area. 
ds.factors <- function(lcm, grid, alpha) {
  # get area of lcm in each grid cell
  # this section replaces the gIntersection function
  lcm.grid <- lcm[grid,]
  int <- gIntersects(lcm.grid, grid, byid=TRUE)
  vec <- vector(mode="list", length=dim(int)[2])
  for (i in seq(along=vec)) vec[[i]] <- try(gIntersection(lcm.grid[i,], 
                                                          grid[int[,i],], byid=TRUE), silent=TRUE)
  # this section deals with the (very rare) polygons which match exactly with
  # grid lines. It adds a marginal buffer to the shape so that it no longer 
  # matches exactly
  for(i in 1:length(vec)){
    if(class(vec[[i]])=="try-error"){
      shp <- gBuffer(lcm.grid[i,], width=0.00000001) # smallest buffer that works
      vec[[i]] <- try(gIntersection(shp,grid[int[,i],], byid=TRUE))
    }
  }
  setAs("SpatialCollections", "SpatialPolygons", function(from) from@polyobj)
  vecP = lapply(vec,function(x){as(x,"SpatialPolygons")})
  lcm.by.grid <- do.call("rbind", vecP)
  ##
  
  lcm.by.grid.df <- data.frame(lcm.by.grid.id=names(lcm.by.grid))
  lcm.by.grid.df$lcm.by.grid.id <- as.character(lcm.by.grid.df$lcm.by.grid.id) # convert the name to character
  splitid <- strsplit(lcm.by.grid.df$lcm.by.grid.id, " ", fixed=TRUE) # split the names
  splitid <- do.call("rbind", splitid) # rbind those back together
  colnames(splitid) <- c("lcmID", "gridID") # now you have the lcm ID and the grid ID as separate variables in a dataframe 
  lcm.by.grid.df <- data.frame(lcm.by.grid.df, splitid) # make that into a dataframe
  lcm.by.grid.area <- data.frame(gArea(lcm.by.grid, byid=T)) # calculate area lcm falling in each grid cell
  lcm.by.grid.area <- lcm.by.grid.area/10000 # convert to ha
  lcm.by.grid.area$lcm.by.grid.id <- rownames(lcm.by.grid.area) # get field to merge on
  colnames(lcm.by.grid.area) <- c("lcm.by.grid.area", "lcm.by.grid.id")
  lcm.by.grid.df <- merge(lcm.by.grid.df, lcm.by.grid.area) # merge
  
  #get total area of lcm
  lcm.area <- data.frame(gArea(lcm, byid=TRUE))
  lcm.area <- lcm.area/10000 # convert units to hectares
  lcm.area$lcmID <- row.names(lcm.area)
  colnames(lcm.area) <- c("total.area", "lcmID")
  lcm.by.grid.df <- merge(lcm.by.grid.df, lcm.area)
  
  # create proportion of each lcm patch falling into grid cells
  lcm.by.grid.df$prop <- lcm.by.grid.df$lcm.by.grid.area / lcm.by.grid.df$total.area
  
  # get the distance to the closest patch within a DIFFERENT occupied grid cell
  dist.df <- data.frame(gDistance(lcm.by.grid, lcm.by.grid, byid=TRUE))
  dist.df$ID <- as.character(row.names(dist.df)) # convert the name to character
  splitid <- strsplit(dist.df$ID, " ", fixed=TRUE) # split the names
  splitid <- do.call("rbind", splitid) # rbind those back together
  colnames(splitid) <- c("lcmID", "gridID") # now you have the lcm ID and the grid ID as separate variables in a dataframe 
  dist.df <- data.frame(dist.df, splitid) # make that into a dataframe
  
  dist <- data.frame(lcmID = 0, dist = 0)
  for(i in 1:(nrow(dist.df))){
    dist[i, 1] <- as.numeric(as.character(dist.df[i, (ncol(dist.df) - 1)]))
    dist[i, 2] <- if(min(dist.df[dist.df[,ncol(dist.df)]!=dist.df[i,ncol(dist.df)],i]) == 0) 0.001 else 
      min(dist.df[dist.df[,ncol(dist.df)]!=dist.df[i,ncol(dist.df)],i]) 
    # replacing 0 with 0.001 so that it will work in the inverse (avoids div by zero issues)
    # JUSTTIFICATION to put in  write up - reason for leaving in the (almost) zeros is 
    # that if a patch straddles two occupied grid cells, it is probably occupied. 
  }

  # create different weightings
  lcm.info <- unique(merge(lcm.by.grid.df, dist))
  lcm.info$IDW <- 1 / lcm.info$dist # turn this into distance weighting based on inverse distance weighting
  lcm.info$ifm.dist <- exp(-alpha*(lcm.info$dist/1000)) # distance measure based on IFM connectivity (without area)
  lcm.info$ifm.dist.area <- lcm.info$ifm.dist*log(lcm.info$total.area) # distance measure based on IFM connectivity (with area)
  lcm.info$str.dist <- max(lcm.info$dist) - lcm.info$dist # distance measure that forms a straight line relationship
  lcm.info$log.IDW <- log(lcm.info$IDW) # log of the inverse distance measure
  lcm.info$log.area <- log(lcm.info$total.area)
  
  # scale all weightings to 0 to 1 to create probability.
  lcm.info$st.total.area <- (lcm.info$total.area - min(lcm.info$total.area))/(max(lcm.info$total.area) - min(lcm.info$total.area))
  lcm.info$st.IDW <- (lcm.info$IDW - min(lcm.info$IDW))/(max(lcm.info$IDW) - min(lcm.info$IDW))
  lcm.info$st.ifm.dist <- (lcm.info$ifm.dist - min(lcm.info$ifm.dist))/(max(lcm.info$ifm.dist) - min(lcm.info$ifm.dist))
  lcm.info$st.ifm.dist.area <- (lcm.info$ifm.dist.area - min(lcm.info$ifm.dist.area))/(max(lcm.info$ifm.dist.area) - min(lcm.info$ifm.dist.area))
  lcm.info$st.str.dist <- (lcm.info$str.dist - min(lcm.info$str.dist))/(max(lcm.info$str.dist) - min(lcm.info$str.dist))
  lcm.info$st.log.IDW <- (lcm.info$log.IDW - min(lcm.info$log.IDW))/(max(lcm.info$log.IDW) - min(lcm.info$log.IDW))
  lcm.info$st.log.area <- (lcm.info$log.area - min(lcm.info$log.area))/(max(lcm.info$log.area) - min(lcm.info$log.area))
  
  return(lcm.info)  
}
