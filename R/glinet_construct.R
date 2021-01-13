library(sp)

# Cut a segment of length dx from the beginning of the input line xy. Return a new line.
cut_segment <- function(xy, dx){
  al <- 0
  i <- 1
  while ((al + dist_point_to_point(xy[i,], xy[i+1,])) < dx){
    al <- al + dist_point_to_point(xy[i,], xy[i+1,])
    i <- i + 1
  }
  rest <- dx - al
  x <- xy[i,1] + rest*(xy[i+1,1]-xy[i,1])/dist_point_to_point(xy[i,], xy[i+1,])
  y <- xy[i,2] + rest*(xy[i+1,2]-xy[i,2])/dist_point_to_point(xy[i,], xy[i+1,])
  return(rbind(matrix(c(x,y), ncol=2), xy[(i+1):dim(xy)[1],]))
}

# Replace the input line xy by a new line made of segments of approximate length dx.
split_line <- function(xy, dx){
  segment <- matrix(c(xy[1,1], xy[1,2]), ncol = 2)
  ll <- length_line(xy)
  new_dx <- ll / (ll %/% dx)
  ns <- ll / new_dx
  if (ns > 1){
    for (i in 1:(ns-1)){
      xy <- cut_segment(xy, new_dx)
      segment <- rbind(segment, xy[1,])
    }
  }
  rbind(segment, xy[dim(xy)[1],])
}

# Determine whether a vertex with coordinates xy is already present in the list of vertices.
position_in_vertlist <- function(xy, vertlist){
  position <- 0
  if (length(vertlist) > 0){
    for (i in 1:length(vertlist)){
      if (min(vertlist[[i]]$coords==xy) == 1){
        position <- i
        break
      }
    }
  }
  return(position)
}

# Create a list of vertices and their neighbors from the input polyline (SpatialLines).
vert_list <- function(lines, dx){
  verts <- list()
  for (line in lines@lines){
    xy <- split_line(line@Lines[[1]]@coords, dx)
    start <- 2

    # prvni vrchol
    pos <- position_in_vertlist(xy[1,], verts)
    if (pos == 0){
      verts[[length(verts)+1]] <- list(coords=xy[1,], neibs=c(length(verts)+2), w=0, hr=c())
    }
    else {
      verts[[pos]]$neibs[length(verts[[pos]]$neibs)+1] <- length(verts)+1
      verts[[length(verts)+1]] <- list(coords=xy[2,], neibs=c(pos, length(verts)+2), w=0, hr=c())
      start <- 3
    }

    # mezivrcholy
    if (dim(xy)[1] > start){
      for (i in start:(dim(xy)[1]-1)){
        verts[[length(verts)+1]] <- list(coords=xy[i,], neibs=c(length(verts), length(verts)+2), w=0, hr=c())
      }
    }

    # posledni vrchol
    if (dim(xy)[1] >= start){
      pos <- position_in_vertlist(xy[dim(xy)[1],], verts)
      if (pos == 0){
        verts[[length(verts)+1]] <- list(coords=xy[dim(xy)[1],], neibs=c(length(verts)), w=0, hr=c())
      }
      else {
        verts[[length(verts)]]$neibs[2] <- pos
        verts[[pos]]$neibs[length(verts[[pos]]$neibs)+1] <- length(verts)
      }
    } else {
      verts[[length(verts)]]$neibs <- verts[[length(verts)]]$neibs[-length(verts[[length(verts)]]$neibs)]
    }
  }
  return(verts)
}

# Determine whether two edges are same. If complete_test is TRUE, even weights and
# orientation are tested.
is_same_edge <- function(e1, e2, complete_test = FALSE){
  same_xy <- FALSE
  same_other <- TRUE
  if ((e1$from == e2$from & e1$to == e2$to) | (e1$from == e2$to & e1$to == e2$from))
    same_xy <- TRUE
  if (complete_test == TRUE){
    if (e1$oriented != e2$oriented) same_other <- FALSE
    if (e1$w != e2$w) same_other <- FALSE
  }
  return(same_xy & same_other)
}

# Determine whether an edge e is already in the edge_list
is_in_edgelist <- function(e, edge_list){
  is_in <- FALSE
  for (edge in edge_list){
    if (is_same_edge(e, edge)) is_in <- TRUE
  }
  return(is_in)
}

# Create a list of edges from the list of vertices vlist
edge_list <- function(vlist){
  elist <- list()
  for (i in 1:length(vlist)){
    for (n in vlist[[i]]$neibs){
      if (i == n) print("AH!")
      e <- list(from=i, to=n, l=dist_point_to_point(vlist[[i]]$coords, vlist[[n]]$coords), w=NA, oriented=FALSE)
      if (!is_in_edgelist(e, elist))
        elist[[length(elist)+1]] <- e
    }
  }
  return(elist)
}

# glinet constructor
new_glinet <- function(lines, dx=double()){
  stopifnot(class(lines)[1]=="SpatialLines")
  stopifnot(is.double(dx))
  vl <- vert_list(lines, dx)
  el <- edge_list(vl)
  #structure(list(vlist=vl, elist=el), class = "glinet")
  structure(list(dx=dx, vlist=vl, elist=el), class = "glinet")
}

#' Creating the network
#'
#' Creates a discretized network from SpatialLines lines as a list of 2: vlist, elist
#'
#' @param lines The input SpatialLines or SpatialLinesDataFrame object
#' @param dx The edge-segment length (double)
#' @return A glinet object, which is a linear network represented as a list of two lists,
#' vlist - list of vertices, and elist - list of edges.
#' @details The glinet object consists of two lists: list of vertices (vlist) and list of edges (elist).
#' Vertices and edges in these lists have the following structure:
#'
#' Vertex: list:
#' - $coords: numeric vector of coordinates
#' - $neibs: integer vector of neighbors (each neighbor coded by its index in vlist)
#' - $w: vertex weight (numeric)
#' - $hr: named vector of home range identifiers (0/1)
#'
#' Edge: list:
#' - $from: index of first vertex (integer)
#' - $to: index of second  vertex (integer),
#' - $w: edge weight (numeric),
#' - $oriented: indicates whether the edge is oriented (boolean)
#' - $l: edge length (numeric)
#'
#' The weights and orientation of edges are intended for more general use in the future and are not utilized
#' in this version of the package.
#' @export
glinet <- function(lines, dx=double()){
  dx = as.double(dx)
  lines = sp::as.SpatialLines.SLDF(lines)
  new_glinet(lines, dx)
}


#
# Net: list: $vlist: list of vertices
#            $elist: list of edges
