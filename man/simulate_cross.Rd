% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_cross.R
\name{simulate_cross}
\alias{simulate_cross}
\title{Simulate Biparental Crosses}
\usage{
simulate_cross(n.ind, h1, h2, cm.map, prob1 = NULL, prob2 = NULL, seed = NULL)
}
\arguments{
\item{n.ind}{Integer: The number of offspring to generate.}

\item{h1}{Matrix: Homology group for parent 1.}

\item{h2}{Matrix: Homology group for parent 2.}

\item{cm.map}{Vector: Genetic map positions in centiMorgans.}

\item{prob1}{Vector: Probability distribution for parent 1 bivalent configuration.}

\item{prob2}{Vector: Probability distribution for parent 2 bivalent configuration.}

\item{seed}{Numeric: Optional seed for reproducibility.}
}
\value{
A list containing:
\describe{
\item{offspring}{An array of offspring formed by combining gametes from both parents.}
\item{rf.calc}{Recombination fractions calculated for the cross.}
\item{ph1}{The homology group of parent 1.}
\item{ph2}{The homology group of parent 2.}
}
}
\description{
This function simulates offspring from a biparental cross by combining gametes
from two parents, each represented by their homology groups.
}
