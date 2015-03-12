\name{BOD2}
\alias{BOD2}
\non_function{}
\title{Biochemical Oxygen Demand}
\description{
  The \code{BOD2} data frame has 8 rows and 2 columns giving the
  biochemical oxygen demand versus time in an evaluation of water
  quality.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Time}{
      A numeric vector giving the time of the measurement (days).
    }
    \item{demand}{
      A numeric vector giving the biochemical oxygen demand (mg/l).
    }
  }
}
\source{
  Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
    Applications}, Wiley, (Appendix A4.1).

  Originally from Marske (1967), M.Sc. Thesis, University of Wisconsin - Madison.
}
\examples{
\dontrun{
str(BOD2)
# simplest form of fitting a first-order model to these data
(fm1 <- nls(demand ~ A*(1-exp(-exp(lrc)*Time)), data = BOD2,
   start = c(A = 2.2, lrc = log(0.25))))
coef(fm1)
# using the plinear algorithm
(fm2 <- nls(demand ~ (1-exp(-exp(lrc)*Time)), data = BOD2,
   start = c(lrc = log(0.25)), algorithm = "plinear", trace = TRUE))
# using a self-starting model
(fm3 <- nls(demand ~ SSasympOrig(Time, A, lrc), data = BOD2))
plotfit(fm3, xlab = "Time (days)",
        ylab = "Biochemical Oxygen Demand")

}
}
\keyword{datasets}
