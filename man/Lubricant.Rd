\name{Lubricant}
\alias{Lubricant}
\non_function{}
\title{Viscosity of lubricants}
\description{
    The \code{Lubricant} data frame has 53 rows and 3 columns on the
    viscosity of a lubricant at different pressures and temperatures.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{pressure}{
	    a numeric vector of pressures (stokes).
	}
	\item{viscos}{
	    a numeric vector of observed log(kinematic viscosity).
	}
	\item{tempC}{
	    a numeric vector of temperatures (degrees Celsius).
	}
    }
}
\source{
    Bates and Watts (1998), \emph{Nonlinear Regression Analysis and Its
	Applications}, Wiley (Appendix A1.5).
}
\examples{
\dontrun{
str(Lubricant)
tempf <- as.factor(Lubricant$tempC)
levels(tempf) <- paste(c(0.0,25.0,37.8,98.9),"C")
xyplot(viscos ~ pressure, Lubricant, groups = tempf,
     xlab = "Pressure (atm)", type = c("g", "p"),
     ylab = "kinematic viscosity (stokes)",
     scales = list(y = list(log = 2)),
     auto.key = list(space = "top", columns = 4))
rm(tempf)
}}
\keyword{datasets}
