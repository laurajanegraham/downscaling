ds.run <- function(lcm.info, grid, species, occ, grid.poly, occ.sensitivity = 0, uncertainty.pct = 5) { 

        # get occupancy level for sensitivity analysis
        occ <- ifelse(occ.sensitivity + occ < 0, 0, occ.sensitivity + occ) 
        
        if(occ==0) {
                samp <- data.frame(lcmID=lcm.info$lcmID, rand=0, dist=0, log.dist=0, 
                                    str.dist=0, ifm.dist=0, area=0, log.area=0, 
                                    logidw.logarea=0, str.dist.logarea=0, ifm.dist.area=0)
        } else {
                grids <- unique(lcm.info$gridID) # get list of grids for looping
                i <- 1
                d.samp <- list() # list to store IDs from distance weighted sampling
                a.samp <- list() # list to store IDs from area weighted sampling
                r.samp <- list() # list to store IDs from random sampling
                logd.samp <- list()
                loga.samp <- list()
                ifm.samp <- list()
                ifm.area.samp <- list()
                str.dist.samp <- list()
                logidw.logarea.samp <- list()
                str.dist.logarea.samp <- list()
                
                for(gridID in grids) { # loop through each grid sampling for each of the weightings
                        grid.info <- lcm.info[lcm.info$gridID==gridID,] # get lcm patches contained in the grid
                        if(ceiling(occ*nrow(grid.info))==nrow(grid.info)){
                                d.samp[[i]] <- grid.info$lcmID
                                a.samp[[i]] <- grid.info$lcmID
                                r.samp[[i]] <- grid.info$lcmID
                                logd.samp[[i]] <- grid.info$lcmID
                                loga.samp[[i]] <- grid.info$lcmID
                                ifm.samp[[i]] <- grid.info$lcmID
                                ifm.area.samp[[i]] <- grid.info$lcmID
                                str.dist.samp[[i]] <- grid.info$lcmID
                                logidw.logarea.samp[[i]] <- grid.info$lcmID
                                str.dist.logarea.samp[[i]] <- grid.info$lcmID
                        } else {
                                d.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*grid.info$st.IDW))
                                a.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*grid.info$st.total.area))
                                r.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop))
                                logd.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*grid.info$st.log.IDW))
                                loga.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*grid.info$st.log.area))
                                ifm.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*grid.info$st.ifm.dist))
                                ifm.area.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*grid.info$st.ifm.dist.area))
                                str.dist.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*grid.info$st.str.dist))
                                logidw.logarea.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*((grid.info$st.log.IDW+grid.info$st.log.area)/2)))
                                str.dist.logarea.samp[[i]] <- array(sample(grid.info$lcmID, ceiling(occ*nrow(grid.info)), prob=grid.info$prop*((grid.info$st.str.dist+grid.info$st.log.area)/2)))
                        }
                        i <- i + 1
                }
                
                d.samp <- unique(unlist(d.samp))
                a.samp <- unique(unlist(a.samp))
                r.samp <- unique(unlist(r.samp))
                logd.samp <- unique(unlist(logd.samp))
                loga.samp <- unique(unlist(loga.samp))
                ifm.samp <- unique(unlist(ifm.samp))
                ifm.area.samp <- unique(unlist(ifm.samp))
                str.dist.samp <- unique(unlist(str.dist.samp))
                logidw.logarea.samp <- unique(unlist(logidw.logarea.samp))
                str.dist.logarea.samp <- unique(unlist(str.dist.logarea.samp))
                
                samp <- data.frame(unique(lcm.info$lcmID)) # create a data frame to populate with the results
                colnames(samp) <- "lcmID"
                
                # update the column for distance, area and random sampling in the output file
                # random
                samp <- transform(samp, rand = ifelse(lcmID %in% r.samp, 1, 0)) 
                
                # distance based measures
                samp <- transform(samp, dist = ifelse(lcmID %in% d.samp, 1, 0))
                samp <- transform(samp, log.dist = ifelse(lcmID %in% logd.samp, 1, 0))
                samp <- transform(samp, str.dist = ifelse(lcmID %in% str.dist.samp, 1, 0))
                samp <- transform(samp, ifm.dist = ifelse(lcmID %in% ifm.samp, 1, 0))
                
                # area based measures
                samp <- transform(samp, area = ifelse(lcmID %in% a.samp, 1, 0))
                samp <- transform(samp, log.area = ifelse(lcmID %in% loga.samp, 1, 0))
                
                # area and distance based measures 
                samp <- transform(samp, logidw.logarea = ifelse(lcmID %in% logidw.logarea.samp, 1, 0))
                samp <- transform(samp, str.dist.logarea = ifelse(lcmID %in% str.dist.logarea.samp, 1, 0))
                samp <- transform(samp, ifm.dist.area = ifelse(lcmID %in% ifm.area.samp, 1, 0))
        }
        # return the results of the downscaling 
        return(samp)
}
