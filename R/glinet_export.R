# Remove edges between home ranges
remove_edges_between_hr <- function(net, perc){
  net <- add_hr(net, perc)
  field_name <- paste("vol", perc*100, sep="")
  new_vlist <- list()
  new_indices <- rep(0, length(net$vlist))
  for (i in 1:length(net$vlist)){
    if (net$vlist[[i]]$hr[field_name] == 1){
      new_vlist[[length(new_vlist)+1]] <- net$vlist[[i]]
      new_indices[i] <- length(new_vlist)
    }
  }
  for (i in 1:length(new_vlist)){
    for (j in 1:length(new_vlist[[i]]$neibs)){
      new_vlist[[i]]$neibs[j] <- new_indices[new_vlist[[i]]$neibs[j]]
    }
    new_vlist[[i]]$neibs <- new_vlist[[i]]$neibs[!new_vlist[[i]]$neibs==0]
  }
  net$vlist <- new_vlist
  net$elist <- edge_list(new_vlist)
  return(net)
}

# Convert network edges to spatial lines
edges2lines <- function(net){
  ids <- 1:length(net$elist)
  ll <- lapply(ids, function(i){
    e <- net$elist[[i]]
    l <- Lines(list(Line(rbind(net$vlist[[e$from]]$coords, net$vlist[[e$to]]$coords))), ID=i)
  })
  return(SpatialLines(LinesList = ll))
}

# Convert network vertices to spatial points
verts2points <- function(net){
  coords <- do.call(rbind, lapply(net$vlist, function(.) .$coords))
  w <- get_vertex_weights(net)
  sp <- SpatialPointsDataFrame(coords, data=data.frame(w=w))
  return(sp)
}

# Convert a network into xy matrices, each corresponding to one "branch" of the network
net2segments <- function(net){

  segments <- list()

  # initialize processed indicators
  free_visits <- sapply(1:length(net$vlist), function(i) {
    ifelse(length(net$vlist[[i]]$neibs) %in% c(1,2), 1, length(net$vlist[[i]]$neibs))
  })
  vert_degrees <- sapply(net$vlist, function(.) length(.$neibs))

  continue <- TRUE
  while (continue){

    # find free end vertex
    processed <- match(TRUE, vert_degrees != 2 & free_visits > 0)

    # if no free end vertex
    if (is.na(processed)){
      continue <- FALSE

      # else
    } else {
      segment <- matrix(net$vlist[[processed]]$coords, nrow = 1)
      free_visits[processed] <- free_visits[processed] - 1
      previous <- -1

      step_forward <- TRUE
      while (step_forward){

        # find a free neighbor
        neib_i <- match(TRUE, free_visits[net$vlist[[processed]]$neibs] > 0 & net$vlist[[processed]]$neibs != previous)

        # if no free neighbor
        if (is.na(neib_i)){
          step_forward <- FALSE
          if (nrow(segment) > 1){
            segments[[length(segments) + 1]] <- segment
          }

          # else
        } else {
          previous <- processed
          processed <- net$vlist[[processed]]$neibs[neib_i]
          segment <- rbind(segment, net$vlist[[processed]]$coords)
          free_visits[processed] <- free_visits[processed] - 1

          # if it's end vertex
          if (vert_degrees[processed] != 2){
            step_forward <- FALSE
            if (nrow(segment) > 1){
              segments[[length(segments) + 1]] <- segment
            }
          }
        }
      }
    }
  }
  return(segments)
}

# Convert network to segment lines (from end to confluence etc...)
# If dissolve == TRUE, touching segments will form multipart spatial lines
net2segmLines <- function(net, dissolve = TRUE){
  s <- net2segments(net)
  if (dissolve) {
    s <- dissolve_segments(s)
    l <- msegments2sl(s)
  } else {
    l <- segments2sl(s)
  }
  return(l)
}

#' Convert the network to Spatial Points
#'
#' This function converts the glinet network into SpatialPointsDataFrame with weights in a field "w" and
#' eventual home range identifiers in appropriately named fields.
#'
#' @param net The glinet object.
#' @return The SpatialPointsDataFrame object.
#' @author Vojta Barták
#' @export
glinet_to_sp <- function(net){
  pp = verts2points(net)
  for (hr in names(net$vlist[[1]]$hr)){
    pp@data[,hr] <- sapply(net$vlist, function(.) .$hr[hr])
  }
  return(pp)
}

#' Rasterize the network
#'
#' This function converts the glinet network into raster, with raster values representing either the
#' network geometry or the vertex weights as well.
#'
#' @param net The glinet object.
#' @param cellsize The cellsize of the output raster.
#' @param use_vertex_weights If TRUE (default), the raster values represent vertext weights. If FALSE, the
#' raster values represent binary information about the network location.
#' @return A raster::RasterLayer object.
#' @author Vojta Barták
#' @export
rasterize_glinet <- function(net, cellsize, use_vertex_weights = TRUE){
  coords <- do.call(rbind, lapply(net$vlist, function(v) v$coords))
  if (use_vertex_weights){
    w <- get_vertex_weights(net)
  } else {
    w <- 1
  }
  pp <- SpatialPointsDataFrame(coords, data.frame(w=w))
  template <- raster(ext = extent(bbox(pp)), resolution = cellsize)
  r <- raster::rasterize(pp, template, field = "w", fun = max)
  return(r)
}

#' Extract home ranges
#'
#' This function extract home ranges from the glinet network in the form of spatial lines
#' each representing one home range.
#'
#' @param net The glinet object.
#' @param perc The volume percentage of the underlying utilization distribution that defines home ranges.
#' @return A SpatialLines (sp) object.
#' @author Vojta Barták
#' @export
glinet_to_hr_lines <- function(net, perc=.95){
  net <- remove_edges_between_hr(net, perc)
  if (length(net$elist) > 0){
    return(net2segmLines(net))
  } else {
    return(NA)
  }
}

#' Determine the number of home ranges
#'
#' This function determines the number of home ranges from the linear network with
#' vertex weights representing utilization distribution.
#' @param net The glinet object.
#' @param perc The volume percentage of the underlying utilization distribution that defines home ranges.
#' @param split If TRUE (default), then home ranges longer than 2*hr_size are splitted.
#' @param hr_size Length of home range size used for splitting longer home ranges. Not used when split
#' is FALSE.
#' @return The number of home ranges in the network (integer).
#' @author Vojta Barták
#' @export
hr_count <- function(net, perc=.95, split=TRUE, hr_size=2000){
  l <- net2hrLines(net, perc)
  if (!is(l, "SpatialLines")){
    return(0)
  } else {
    add <- 0
    if (split){
      for (line in l@lines){
        len <- sum(sapply(line@Lines, function(.) length_line(.@coords)))
        add <- add + max(len%/%hr_size - 1, 0)
      }
    }
    return(length(l@lines)+add)
  }
}
