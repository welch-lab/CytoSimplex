% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ternary.R
\name{plotTernary.simMat}
\alias{plotTernary.simMat}
\title{Create quaternary simplex plot with pre-calculated similarity matrix}
\usage{
\method{plotTernary}{simMat}(
  x,
  title = NULL,
  veloMat = NULL,
  nGrid = 10,
  radius = 0.1,
  dotSize = NULL,
  dotShuffle = NULL,
  colorArg = NULL,
  labelColors = c("#3B4992FF", "#EE0000FF", "#008B45FF"),
  vertexLabelSize = NULL,
  vertexLabelDrift = 0.03,
  axisBreak = 5,
  axisTextShow = TRUE,
  axisTextSize = NULL,
  axisTextDrift = 0.01,
  gridLineAlpha = 0.6,
  arrowLinewidth = NULL,
  arrowAngle = 20,
  arrowLen = 0.2,
  titleSize = NULL,
  equilateral = TRUE,
  margin = 0.1,
  interactive = FALSE,
  ...
)
}
\arguments{
\item{x}{simMat object, n cells by 4 vertices, each row summing to 1.}

\item{title}{Title text of the plot. Default \code{NULL}.}

\item{veloMat}{Aggregated velocity matrix. Output of \code{aggrVeloGraph}.}

\item{nGrid}{Number of grids along the bottom side of the equilateral
triangle. Default \code{10}.}

\item{radius}{Arrow length of unit velocity. Lower this when arrows point
outside of the coordinate. Default \code{0.1}.}

\item{dotSize}{Dot aesthetics passed to \code{geom_point}. Default \code{0.6}
when not interactive, \code{4} when interactive.}

\item{dotShuffle}{Whether to shuffle the order of dots being added to the
plot, useful when categorical colors are used and mixing of categories is
expected. Default \code{NULL} does shuffle when \code{dotColorBy} given is
categorical and does not otherwise.}

\item{colorArg}{A "colorArg" object, internally prepared by
\code{\link{plotQuaternary.default}}. Default \code{NULL}.}

\item{labelColors}{Colors of the axis lines and vertex labels.
Default \code{c("#3B4992FF", "#EE0000FF", "#008B45FF")} (blue, red and green)}

\item{vertexLabelSize}{Size of vertex labels. Default \code{6} when not
interactive, \code{16} when interactive.}

\item{vertexLabelDrift}{Position adjustment of the vertex labels, only
applied to non-interactive view. Default \code{0.03}.}

\item{axisBreak}{Number of breaks to be labeled along axis. Default
\code{5}.}

\item{axisTextShow}{Logical, whether to show axis text. Default \code{TRUE}.}

\item{axisTextSize}{Size of text along each axis break. Default \code{4} for
non-interactive view, \code{12} for interactive view.}

\item{axisTextDrift}{Position adjustment of the axis text, only applied to
non-interactive view. Default \code{0.01}.}

\item{gridLineAlpha}{Transparency of background grid lines. Default
\code{0.6}.}

\item{arrowLinewidth}{Line width of the velocity arrows. Default \code{0.25}
for non-interactive view, \code{2} for interactive view.}

\item{arrowAngle}{Controls the angle of the arrowhead, only applied to
non-interactive view. Default \code{20}.}

\item{arrowLen}{Control length in centimetre from arrow tip to arrow tail,
only applied to non-interactive view. Default \code{0.2}.}

\item{titleSize}{Size of title text. Default \code{14} for non-interactive
view, \code{20} for interactive view.}

\item{equilateral}{Logical, whether to always display the triangle as
equilateral. Default \code{TRUE}.}

\item{margin}{Margin allowed around of the triangle plotting region when
\code{equilateral = TRUE}}

\item{interactive}{Logical. Whether to display plotly interactive view.
Default \code{FALSE}.}

\item{...}{Not used}
}
\description{
Create quaternary simplex plot with pre-calculated similarity matrix
}
