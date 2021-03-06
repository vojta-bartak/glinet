% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glinet_utils.R
\name{plot.glinet}
\alias{plot.glinet}
\title{Plot glinet network}
\usage{
\method{plot}{glinet}(
  x,
  y = 0,
  col_edge = "black",
  col_vert = "red",
  size_vert = double(),
  size_factor = 2
)
}
\arguments{
\item{x}{glinet object to plotted}

\item{col_edge}{Color of edges}

\item{col_vert}{Color of vertices}

\item{size_vert}{Size of vertices.}

\item{size_factor}{If size_vert not set, the size of vertices is determined by
their weights and multiplied by size_factor. The default is 2.}
}
\description{
This function plots the glinet network.
}
