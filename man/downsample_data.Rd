% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downsample_data.R
\name{downsample_data}
\alias{downsample_data}
\title{Downsample Data}
\usage{
downsample_data(
  data,
  sample_column = NULL,
  sample_rate = NULL,
  downsample_rate = NULL,
  downsample_time = NULL,
  verbose = TRUE
)
}
\arguments{
\item{data}{A dataframe.}

\item{sample_column}{\emph{(Optional)}. A character scalar indicating the name of
the time or sample data column. Must match exactly.}

\item{sample_rate}{\emph{(Optional)}. A numeric scalar for the sample rate in Hz.}

\item{downsample_rate}{\emph{(Optional)}. A numeric scalar indicating the desired
output sample rate (in Hz) to convert the dataframe.}

\item{downsample_time}{\emph{(Optional)}. A numeric scalar indicating the desired
sample time (in seconds) to convert the dataframe.}

\item{verbose}{A logical. \code{TRUE} (\emph{default}) will return warnings and
messages which can be used for data error checking. \code{FALSE} will silence these
messages. Errors will always be returned.}
}
\value{
A \link[tibble:tibble-package]{tibble} of class \code{mNIRS.data} with
metadata available with \code{attributes()}.
}
\description{
Condense the number of samples of a dataframe using time-weighted averaging.
}
\details{
\code{sample_column} and \code{sample_rate} will be taken from metadata for an
mNIRS dataframe, if not defined explicitly.

If not present in metadata, \code{sample_column} must be defined explicitly.
\code{sample_rate} will be estimated based on the mean difference between values
in the \code{sample_column}. If \code{sample_column} contains integer sample numbers,
then \code{sample_rate} will be incorrectly estimated  to be 1 Hz, and should be
defined explicitly.

\emph{TODO include upsample interpolation?}
}
