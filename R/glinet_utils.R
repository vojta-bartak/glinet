#' Get the vertex weights
#'
#' This function extracts the vertex weights as a numeric vector.
#'
#' @param net glinet object
#' @return numeric vector of vertex weights
#' @author Vojta Barták
#' @export
get_vertex_weights <- function(net){
  sapply(net$vlist, function(v) v$w)
}

# Get the edge lengths as a numeric vector
get_edge_lengths <- function(net){
  sapply(net$elist, function(e) e$l)
}

# Get vertex number
get_vertex <- function(x, y, net, prec = 0.01){
  i <- 0
  for (j in 1:length(net$vlist)){
    v <- net$vlist[[j]]
    if (max(abs(c(v$coords[1] - x, v$coords[2] - y))) < prec)
      i <- j
  }
  return(i)
}

#' Plot glinet network
#'
#' This function plots the glinet network.
#'
#' @param x glinet object to plotted
#' @param col_edge Color of edges
#' @param col_vert Color of vertices
#' @param size_vert Size of vertices.
#' @param size_factor If size_vert not set, the size of vertices is determined by
#' their weights and multiplied by size_factor. The default is 2.
#' @export
plot.glinet <- function(x, y=0, col_edge="black", col_vert="red", size_vert=double(), size_factor=2){
  lines <- edges2lines(x)
  points <- verts2points(x)
  plot(lines, col = col_edge)
  if (length(size_vert) == 0) size_vert = points$w*size_factor/max(points$w)
  plot(points, add=TRUE, col=col_vert, pch=16, cex=size_vert)
}
