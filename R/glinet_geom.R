

# Determine Euclidian distance of two points p1, p2
dist_point_to_point <- function(p1, p2){
  sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2)
}

# Determine orthogonal projection of point p to a line defined by two points p1, p2
project_point_to_line <- function(p, p1, p2){
  a <- p1[2] - p2[2]
  b <- p2[1] - p1[1]
  c <- p1[1]*p2[2] - p2[1]*p1[2]
  d <- b*p[1] - a*p[2]
  #print(paste("new point", class(p)))
  #flush.console()
  return(solve(matrix(c(a, -b, b, a), ncol=2), c(-c, -d)))
}

# Determine orthogonal projection of point p to a segment segm
# If nearest=TRUE (default) the distance to the nearest point of a segment is returned,
# otherwise only orthogonal projection of the point is returned, or NA when the projection
# is outside the segment.
project_point_to_segment <- function(p, segm, nearest = TRUE){
  q <- project_point_to_line(p, segm[1,], segm[2,])
  d1 <- dist_point_to_point(segm[1,], segm[2,])
  d2 <- dist_point_to_point(q, segm[1,])
  d3 <- dist_point_to_point(q, segm[2,])
  if (max(d2, d3) > d1){
    if (nearest) {
      if (d2 < d3) return(segm[1,])
      else return(segm[2,])
    } else return(NA)
  } else {
    return(q)
  }
}

# Determine Euclidian distance between point p and segment line segm
# p: vector of coordinates, segm: matrix of coordinates
dist_point_to_segment <- function(p, segm){
  q <- project_point_to_segment(p, segm)
  if (!is.na(q[1])){
    return(dist_point_to_point(p, q))
  } else {
    return(NA)
  }
}

# Identifies whether two segments (xy matrices) touch or not (by ending points)
segm_touch <- function(segm1, segm2, precision = 0.001){
  (abs(segm1[1,1] - segm2[1,1]) < precision) & (abs(segm1[1,2] - segm2[1,2]) < precision) |
    (abs(segm1[1,1] - segm2[nrow(segm2),1]) < precision) & (abs(segm1[1,2] - segm2[nrow(segm2),2]) < precision) |
    (abs(segm2[1,1] - segm1[nrow(segm1),1]) < precision) & (abs(segm2[1,2] - segm1[nrow(segm1),2]) < precision) |
    (abs(segm2[nrow(segm2),1] - segm1[nrow(segm1),1]) < precision) & (abs(segm2[nrow(segm2),2] - segm1[nrow(segm1),2]) < precision)
}

# Identifies whether a segment is touching a multi-segment (list of segements) by ending points of segments
multi_segm_touch <- function(segm, segmList, precision = 0.001){
  max(sapply(segmList, function(s) segm_touch(segm, s, precision))) == 1
}

# Reorganize segments (xy matrices) to group those touching (by ending points)
dissolve_segments <- function(segments, precision=0.1){
  new_segments <- list()
  processed <- rep(FALSE, length(segments))
  continue <- TRUE
  while (continue){
    segment_i <- match(FALSE, processed)
    if (is.na(segment_i)){
      continue <- FALSE
    } else {
      new_segments[[length(new_segments) + 1]] <- list()
      todo <- list(segment_i)
      processed[segment_i] <- TRUE
      while (length(todo) > 0){
        i <- todo[[1]]
        todo <- todo[-1]
        new_segments[[length(new_segments)]][[length(new_segments[[length(new_segments)]]) + 1]] <- segments[[i]]
        for (j in 1:length(segments)){
          if (!processed[j] & segm_touch(segments[[i]], segments[[j]], precision)){
            todo[[length(todo) + 1]] <- j
            processed[j] <- TRUE
          }
        }
      }
    }
  }
  return(new_segments)
}

# Convert segments (list of xy matrices) to spatial lines
segments2sl <- function(segments){
  lineslist <- lapply(1:length(segments), function(i) Lines(list(Line(segments[[i]])), i))
  sl <- SpatialLines(lineslist)
  return(sl)
}

# Convert multipart segments (list of lists of xy matrices) to spatial lines
msegments2sl <- function(segments){
  lineslist <- lapply(1:length(segments), function(i) Lines(lapply(segments[[i]], Line), i))
  sl <- SpatialLines(lineslist)
  return(sl)
}

# Total length of a line given as a matrix of XY coordinates
length_line <- function(xy){
  sum(sapply(2:dim(xy)[1], function(i) sqrt((xy[i-1, 1]-xy[i, 1])^2 + (xy[i-1, 2]-xy[i, 2])^2)))
}

# Distance from point to line ------------------------------------------------------------------------------------------------
# line is object of class Lines (package sp), point is matrix of xy coords
dist_to_line <- function(point, line){
  min(sapply(line@Lines, function(l){
    min(sapply(1:(nrow(l@coords)-1), function(i){
      dist_point_to_segment(point, l@coords[i:(i+1),])
    }))
  }))
}

# Distance to nearest line, returns vector of distances of the same length as input SpatialPoint layer -----------------------
dist_to_nearest_hr <- function(points, lines){
  sapply(1:nrow(points@coords), function(i) {
    min(sapply(lines@lines, function(l) {
      dist_to_line(points@coords[i,], l)
    }))
  })
}
