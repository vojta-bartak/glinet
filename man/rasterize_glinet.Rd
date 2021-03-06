% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glinet_export.R
\name{rasterize_glinet}
\alias{rasterize_glinet}
\title{Rasterize the network}
\usage{
rasterize_glinet(net, cellsize, use_vertex_weights = TRUE)
}
\arguments{
\item{net}{The glinet object.}

\item{cellsize}{The cellsize of the output raster.}

\item{use_vertex_weights}{If TRUE (default), the raster values represent vertext weights. If FALSE, the
raster values represent binary information about the network location.}
}
\value{
A raster::RasterLayer object.
}
\description{
This function converts the glinet network into raster, with raster values representing either the
network geometry or the vertex weights as well.
}
\author{
Vojta Barták
}
