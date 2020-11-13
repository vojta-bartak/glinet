# Find the position of a point xy in the network net, i.e. the two nearest neighboring vertices
find_in_net <- function(xy, net){
  mindist <- Inf
  edge_index <- NA
  i <- 1
  for (e in net$elist){
    d <- dist_point_to_segment(xy, rbind(net$vlist[[e$from]]$coords, net$vlist[[e$to]]$coords))
    if (!is.na(d) & d < mindist){
      mindist <- d
      edge_index <- i
    }
    i <- i + 1
  }
  return(list(edge_index, mindist))
}

#' Initialize the network
#'
#' Initializes the network net by a set of points xy and their corresponding weights w.
#'
#' @param net A glinet object.
#' @param xy A table-like object (data.frame or matrix) with two columns for x and y coordinates,
#' representing the set of points (observations or events) whose density should be estimated on the
#' network.
#' @param w A numeric vector (of the same length as the number of points in xy) of weights. Each weight
#' will be distributed to the nearest two vertices of the network.
#' @return An initialized glinet object.
#' @export
init_glinet <- function(net, xy, w){
  for (i in 1:nrow(xy)){
    eimd <- find_in_net(xy[i,], net)
    ei <- eimd[[1]]
    if (!is.na(ei)){
      part1 <- dist_point_to_point(xy[i,], net$vlist[[net$elist[[ei]]$to]]$coords)/net$elist[[ei]]$l
      if (part1 > 1) print(paste("Too big! Mindist", eimd[[2]], "Point", xy[1,1], xy[1,2], "Start", net$elist[[ei]]$from, "End", net$elist[[ei]]$to, "Edge", ei))
      flush.console()
      #part2 <- dist_point_to_point(xy[i,], net$vlist[[net$elist[[ei]]$from]]$coords)/net$elist[[ei]]$l
      part2 <- 1 - part1
      net$vlist[[net$elist[[ei]]$from]]$w <- net$vlist[[net$elist[[ei]]$from]]$w + part1*w[i]
      net$vlist[[net$elist[[ei]]$to]]$w <- net$vlist[[net$elist[[ei]]$to]]$w + part2*w[i]
    }
  }
  return(net)
}

# Estimate the Anderson-Morley bound B for the network
estimate_B <- function(net){
  max(sapply(net$elist, function(e){
    length(net$vlist[[e$from]]$neibs) + length(net$vlist[[e$to]]$neibs)
  }))
}

# Update the network in the heat simulation
update_net <- function(net, alpha){
  new_net <- net
  for (i in 1:length(net$vlist)){
    deg <- length(net$vlist[[i]]$neibs)
    new_net$vlist[[i]]$w <- net$vlist[[i]]$w*(1-alpha*deg) + alpha*sum(
      sapply(net$vlist[[i]]$neibs, function(neib){
        net$vlist[[neib]]$w
      })
    )
  }
  return(new_net)
}

#' Simulate the heat equation
#'
#' This function simulates the heat equation on a linear network discretized into segments of same lengths.
#' It thus effectively estimates a density of point events (with weights) located on the network vertices.
#'
#' @param net A glinet object.
#' @param sigma ...
#' @param beta ...
#' @return Returns a new glinet object with vertex weights updated by heat propagation.
#' @author Vojta Barták
#' @references ...
#' @export
heat_glinet <- function(net, sigma, beta=.5){
  dx = net$dx
  B <- estimate_B(net)
  bound1 <- dx*dx/(beta*B)
  bound2 <- sigma*dx/3
  dt <- min(bound1, bound2)
  n <- floor(sigma*sigma/dt) + 1
  alpha <- beta*dt/(dx*dx)
  for (i in 1:n){
    net <- update_net(net, alpha)
  }
  return(net)
}

#' Add home range identifiers
#'
#' This function adds home range identifiers to glinet network vertices.
#'
#' @param net The glinet object.
#' @param perc Numeric vector of volume percentages to extract home ranges from the underlying utilization
#' distribution
#' @return The updated glinet object.
#' @author Vojta Barták
#' @export
add_hr <- function(net, perc=c(.5, .9, .95)){
  for (p in perc){
    den <- get_vertex_weights(net)
    z <- sort(den[!is.na(den)], decreasing=TRUE)
    y <- cumsum(as.numeric(z))
    i <- sum(y <= p * y[length(y)])
    field_name <- paste("vol", p*100, sep="")
    hr <- as.numeric(den >= z[i])
    for (i in 1:length(net$vlist)){
      hr_i <- match(field_name, names(net$vlist[[i]]$hr))
      if (!is.na(hr_i)) net$vlist[[i]]$hr <- net$vlist[[i]]$hr[-hr_i]
      net$vlist[[i]]$hr <- c(net$vlist[[i]]$hr, hr[i])
      names(net$vlist[[i]]$hr)[length(net$vlist[[i]]$hr)] <- field_name
    }
  }
  return(net)
}
