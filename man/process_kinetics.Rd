% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_kinetics.R
\name{process_kinetics}
\alias{process_kinetics}
\title{Process Kinetics}
\usage{
process_kinetics(
  x,
  y = NULL,
  data = NULL,
  x0 = 0,
  method = c("monoexponential", "sigmoidal", "half_time", "peak_slope"),
  ...
)
}
\arguments{
\item{x}{A numeric vector specifying the predictor variable for \code{y} if \code{y} is
defined. Otherwise, \code{x} is assumed to define the response variable.}

\item{y}{(\emph{Optional}). A numeric vector specifying the response variable. If \code{y}
is not defined, the response variable is assumed to be defined by \code{x}, with
\code{idx = seq_along(x)} as the predictor variable.}

\item{data}{(\emph{Optional}). A dataframe containing at least the response variable
(\code{x}), or the predictor and response variables (\code{x} and \code{y}).}

\item{x0}{A numeric scalar indicating the value of the predictor variable \code{x}
or \code{idx} representing the start of the kinetics event.}

\item{method}{Indicates how to process the kinetics.
\describe{
\item{\code{method = "monoexponential"}}{...}
\item{\code{method = "sigmoidal"}}{...}
\item{\code{method = "half_time"}}{...}
\item{\code{method = "peak_slope"}}{...}
}}

\item{...}{Additional arguments. Used to define fixed parameters which will
not be optimised by the kinetics methods. e.g. \code{A = 10} will define
\code{SSmonoexp(x, A = 10, B, TD, tau)}}
}
\value{
A list \code{L} of class \code{mNIRS.kinetics} with components \code{L$...}:
\item{\code{model}}{The model object.}
\item{\code{data}}{A dataframe of the input and fitted model data.}
\item{\code{fitted}}{A vector of the fitted values returned by the model.}
\item{\code{coefs}}{A dataframe of the model coefficients, including manually
fixed parameters.}
\item{\code{fit_criteria}}{A dataframe of the model fit criteria
(\code{AIC}, \code{BIC}, \code{R2}, \code{RMSE}, \code{RSE}, \code{MAE}, \code{MAPE}).}
}
\description{
Fit mNIRS kinetics vector data with a parametric or non-parametric curve fitting
model.
}
\details{
\code{method \%in\% c("monoexponential", "sigmoidal")} use \link[stats:nls]{nls()}
for nonlinear (weighted) least-squares estimates.
}
\seealso{
\code{\link[stats:nls]{stats::nls()}}, \code{\link[stats:SSasymp]{stats::SSasymp()}}, \code{\link[stats:SSlogis]{stats::SSlogis()}},
}
