% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapping_functions.R
\name{imf_haldane_cM}
\alias{imf_haldane_cM}
\title{Inverse Haldane Mapping Function (returns cM)}
\usage{
imf_haldane_cM(r)
}
\arguments{
\item{r}{Numeric vector. Recombination frequencies (\eqn{0 \leq r \leq 0.5}).}
}
\value{
Numeric vector of genetic distances in centiMorgans (cM).
}
\description{
Converts recombination frequency (\eqn{r}) into genetic distance in centiMorgans (cM)
according to the Haldane mapping function, assuming no crossover interference.
}
\details{
The Haldane mapping function assumes no crossover interference and relates
recombination fraction (\eqn{r}) to map distance (\eqn{d}) by:
\deqn{d = -\frac{1}{2} \log(1 - 2r)}{d = -0.5 * log(1 - 2r)}
in Morgans. Multiplying by 100 gives centiMorgans.
}
\examples{
imf_haldane_cM(0.1)   # 11.06 cM
imf_haldane_cM(c(0.01, 0.1, 0.2))  # multiple recombination frequencies

}
