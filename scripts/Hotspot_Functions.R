# Spatial Analysis functions 
# Initial function provided by Allen Zhang
# Modified for use in TMAs

# Currently configured to work with default QuPath outputs
# Needs some minor additional modification to accept other 
# formats / marker sets

#' calculate Getis-Ord Gi* hotspots
#' 
#' @param x a spatstat ppp for hotspot calculation
#'
#' @export
get_localG <- function(x, eps=c(10, 10), 
                       neighbourhood_size=4, 
                       type = "pp", 
                       window_option = "none", 
                       nperms = 100, 
                       p_thres = 0.05
                       ) 
  {
  
  if (eps[1] != eps[2]) {
    ## Relies on pp.rescale, which does not accept separate x and y conversions
    stop("No longer valid.")
  }
  
  ## Save the original (possibly polygonal) window. as.mask does a bad job coercing
  if (type == "phenoDat") {
    pp <- x$pp
  } else if (type == "pp") {
    pp <- x
  } else {
    stop("Unrecognized type.")
  }
  
  origwin <- pp$window
  
  ## rasterize pp
  pp.pix <- spatstat::pixellate(pp, eps=eps)
  mat <- pp.pix$v
  dims <- dim(mat)
  
  ## Don't process if we don't have any hits
  if (pp$n == 0) {
    ## Use 0 since we do 2-sided p-value test later
    return(matrix(0, nrow=dims[1], ncol=dims[2]))
  }
  
  xyw <- setNames(data.frame(mat, stringsAsFactors=FALSE), c("y", "x", "npoints"))
  xyw$x <- pp.pix$xcol[as.numeric(as.character(xyw$x))] / eps[1]
  xyw$y <- pp.pix$yrow[as.numeric(as.character(xyw$y))] / eps[2]
  
  ## Points inside the window (complex for polygonal windows, due to bug in as.mask)
  if (window_option == "grid") {
    inside <- !is.na(xyw$npoints)
    #inside <- xyw$npoints > 0
  } else if (origwin$type == "polygonal") {
    spp <- sp::SpatialPoints(data.frame(x=xyw$x, y=xyw$y))
    pp.rescale <- spatstat::rescale(pp, s = eps[1])
    newwin <- pp.rescale$window
    # require(maptools)
    inside <- as.vector(polycontains(as(newwin, "SpatialPolygons"), spp, byid = TRUE))
  } else if (origwin$type == "mask") {
    # not sure why I need to move conversion to spp... Not enough pixel resolution I think
    spp <- sp::SpatialPoints(data.frame(x=xyw$x, y=xyw$y) * eps[1])
    #pp.rescale <- spatstat::rescale(pp, s = eps[1])
    newwin <- pp$window
    inside <- spatstat::inside.owin(x = spp@coords[,"x"], y = spp@coords[,"y"], w = newwin)
  } else {
    inside <- rep(TRUE, nrow(xyw))
  }
  
  xyw_filtered <- xyw[inside,]
  xyw_filtered$npoints[is.na(xyw_filtered$npoints)] <- 0
  
  G_vals <- rep(NA, nrow(xyw))
  
  xy <- with(xyw_filtered, cbind(x, y))
  nb <- spdep::dnearneigh(xy, 0, neighbourhood_size)
  
  ## Binary coding style
  G <- spdep::localG(xyw_filtered$npoints, spdep::nb2listw(spdep::include.self(nb), style="B", zero.policy = TRUE), zero.policy = TRUE)
  if (nperms > 0) {
    permutations <- lapply(1:nperms, function(x) sample(1:nrow(xyw_filtered)))
    permGs <- lapply(seq_along(permutations), function(i) {
      perm <- permutations[[i]]
      pG <- spdep::localG(xyw_filtered$npoints[perm], spdep::nb2listw(spdep::include.self(nb), style="B", zero.policy = TRUE), zero.policy = TRUE)
      return(pG)
    })
    
    pMat <- do.call('cbind', lapply(permGs, function(x) as.numeric(x)))
    pMatSorted <- t(apply(pMat, 1, grr::sort2))
    
    sig_col <- pracma::ceil(nperms * (1-p_thres))
    crit_values <- pMatSorted[, sig_col]
    
    binG <- as.numeric(G >= crit_values)
    
  } else {
    ## taken from Yinyin's lab
    sq.count <- length(G)
    
    if (p_thres != 0.05) {
      warning("Calculating with p = 0.05, cause critical values not determined for others ...")
    }
    
    if (sq.count >= 1 & sq.count < 50){ g <- 1.645}
    if (sq.count >= 50 & sq.count < 100){ g <- 3.083}
    if (sq.count >= 100 & sq.count < 1000){ g <- 3.289}
    if (sq.count >= 1000){ g <- 3.886}
    
    binG <- as.numeric(G > g)
  }
  
  G_vals[inside] <- binG
  G_mat <- matrix(G_vals, nrow=dims[1], ncol=dims[2])
  return(G_mat)
}

#' Polygon contains helper functions
#' 
polycontains <- function(sp, spp, byid = TRUE) {
  rgeos::gContains(sp, spp, byid = byid)
}


## Following uses the default kernel smoother in spatstat to
## mask out non-tumor areas (from tumor / stroma dual
## classification) ... assumes 'tissue' mark with 
## tumor class.

#' Mask spatstat ppp to retain only tumor class cell areas
#' 
#' @param pp a ppp object created by spatstat, with marks for tissue type
#'
#' @export
remask_tissue <- function(pp, 
                          tissue_type='Tumor',
                          cutoff=1.05) 
{
  message('re-masking tissue')
  smooth_pp <- Smooth(pp, sigma = bw.ppl, adjust=1.5)
  nmask <- smooth_pp$tissue$v
  nmask[is.na(nmask)] <- 1
  nmask <- nmask > cutoff
  pp.sub <- spatstat::subset.ppp(pp, tissue==tissue_type)
  nwin <- spatstat::owin(xrange = smooth_pp$tissue$xrange, 
                         yrange = smooth_pp$tissue$yrange, 
                         mask = nmask)
  pp.t <- spatstat::ppp(pp.sub$x, pp.sub$y,
                        window = nwin)
  marks(pp.t) <- spatstat::marks(pp.sub)[ , !names(spatstat::marks(pp.sub)) %in% 'tissue']
  return(pp.t)
}

# remask_tissue <- function(pp, 
#                           tissue_type='Tumor',
#                           sigma=NULL) 
#   {
#   
#   smooth_pp <- Smooth(pp)
#   nmask <- smooth_pp$tissue$v
#   # nmask[is.na(nmask)] <- 1
#   nmask <- !is.na(nmask)
#   pp.sub <- spatstat::subset.ppp(pp, tissue==tissue_type)
#   nwin <- spatstat::owin(xrange = smooth_pp$tissue$xrange, 
#                                   yrange = smooth_pp$tissue$yrange, 
#                                   mask = nmask)
#   pp.t <- spatstat::ppp(pp.sub$x, pp.sub$y,
#                         window = nwin)
#   marks(pp.t) <- spatstat::marks(pp.sub)[ , !names(spatstat::marks(pp.sub)) %in% 'tissue']
#   return(pp.t)
# }


# helper function from Allen
create_hotspot_object <- function(G_mat) {
  mat <- G_mat
  obj <- list(mat=mat)
  return(obj)
}

# Another helper: intersect hotspots between two objects
intersect_hotspots <- function(obj1, obj2) {
  # Assumes the same distribution of NAs -- this should be true.!
  eidx <- which(!is.na(obj1))
  eidx2 <- which(!is.na(obj2))
  if(! identical(eidx, eidx2)) stop('are you comparing the same tissue?')
  idx1 <- obj1 == 1
  idx2 <- obj2 == 1
  
  mat <- obj1
  mat[eidx] <- 0
  mat[idx1 & idx2] <- 1
  
  return(mat)
}

#' Compute fractional overlap statistics
#' 
#' @param objc Cancer hotspots
#' @param objl Lymphocyte hotspots
#' @param objcl Cancer-immune hotspots
#'
#' @export
compute_fractional_overlap <- function(objc, objl, objcl) {
  
  matc <- objc
  matl <- objl
  matcl <- objcl
  
  idx1 <- matc == 1
  idx2 <- matl == 1
  idx3 <- matcl == 1
  
  f_c <- length(which(idx3)) / length(which(idx1))
  f_i <- length(which(idx3)) / length(which(idx2))
  f_ci <- length(which(idx3)) / length(which(!is.na(matc)))
  
  return(list(c=f_c, i=f_i, ci=f_ci, 
              nc=length(which(idx1)), 
              ni=length(which(idx2)),
              nci = length(which(idx3)),
              ntot = length(which(!is.na(matc))))
         )
}


#' Run hotspot pipeline. 
#' Could generalize to more flexible markers/thresholds: 
#' For now just run with default params for QuPath/CD155 project
#' 
#' @param pp a point process to run hotspot analysis on
#' 
#' @export
run_hotspot <- function(pp) {
  
  if(sum(pp$marks$tissue=='Tumor') < 100 ) {
    message('Fewer than 100 tumor cells... Skipping.')
    return(NULL)
  }
  
  # thresholds set by inspection in QuPath
  pp.t <- remask_tissue(pp)
  pp.pd1 <- subset.ppp(pp.t, pd1 > 1)
  pp.cd8 <- subset.ppp(pp.t, cd8 > 1)
  pp.pdl1 <- subset.ppp(pp.t, pdl1 > 0.5)
  pp.cd155 <- subset.ppp(pp.t, cd155 > 1)
  
  if(pp.pd1$n < 5 | pp.pdl1$n < 5 | pp.cd155$n < 5) {
    message('One or more phenotypes has fewer than 5 cells... Skipping.')
    return(NULL)
  }
  
  gpd1 <- get_localG(pp.pd1, nperms = 100)
  gcd8 <- get_localG(pp.cd8, nperms = 100)
  gpdl1 <- get_localG(pp.pdl1, nperms = 100)
  gcd155 <- get_localG(pp.cd155, nperms = 100)
  
  l_out  <- list(g_pd1 = gpd1,
                 g_cd8 = gcd8, 
                 g_pdl1 = gpdl1, 
                 g_cd155 = gcd155)
  return(l_out)
}

#' generate point process from QuPath annotations (spatstat::ppp)
#' specific to current QuPath output format
#' 
#' @param qupath_detections qupath detection output
#' 
#' @export
generate_ppp <- function(qupath_detections) 
{
  # verbose for debugging
  message('## Reading ', qupath_detections)
  
  ff <- data.table::fread(qupath_detections)
  
  x <- ff$`Centroid X \xb5m`
  y <- ff$`Centroid Y \xb5m`
  
  # use convex hull to clean tissue boundaries
  ch <- chull( cbind(x,y) )
  coords <- cbind(x,y)[c(ch, ch[1]), ] 
  # reverse to clockwise traversal
  coords <- coords[nrow(coords):1, ]
  
  cd155 <- ff$`Cell: CD155 mean`
  pd1 <- ff$`Cell: PD-1 mean`
  cd8 <- ff$`Cell: CD8 mean`
  pdl1 <- ff$`Cell: PD-L1 mean`
  tissue <- ff$`Name`
  tissue[tissue!='Tumor'] <- 'Other'
  
  pwin <- spatstat::owin(poly = list(x=coords[,1], y=coords[,2]))
  pwin <- spatstat::expand.owin(pwin, distance=15)
  pp <- ppp(x=x, y=y, pwin)
  marks(pp)  <- data.frame(cd155, pd1, cd8, pdl1, tissue)
  
  return(pp)
  
}

#' Compute F stats from output of Getis-Ord Gi* wrapper. 
#' Just set up for CD155 project right now.
#' TODO: generalize to any marker sets
#' 
#' @param gmats list of Gi* matrices
#' 
#' @export
getFcs <- function(gmats) {
  
  f_pdl1_cd8 <- compute_fractional_overlap(gmats$g_pdl1, 
                                           gmats$g_cd8,
                                           intersect_hotspots(gmats$g_cd8, gmats$g_pdl1)
                                           )
  f_cd155_cd8 <- compute_fractional_overlap(gmats$g_cd155, 
                                            gmats$g_cd8, 
                                            intersect_hotspots(gmats$g_cd8, gmats$g_cd155)
                                            )
  f_pdl1_pd1 <- compute_fractional_overlap(gmats$g_pdl1, 
                                           gmats$g_pd1, 
                                           intersect_hotspots(gmats$g_pd1, gmats$g_pdl1)
                                           )
  f_cd155_pd1 <- compute_fractional_overlap(gmats$g_cd155, 
                                            gmats$g_pd1, 
                                            intersect_hotspots(gmats$g_pd1, gmats$g_cd155)
                                            )
  return( list(f_pdl1_cd8=f_pdl1_cd8, 
               f_cd155_cd8=f_cd155_cd8, 
               f_pdl1_pd1=f_pdl1_pd1, 
               f_cd155_pd1=f_cd155_pd1) 
  )
}
