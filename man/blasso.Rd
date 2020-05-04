\name{blassocox}
\alias{blassocox_f}
\title{blasso placeholder}
\usage{
	blassocox(x)
}
\arguments{
	\item{x}{a vector of values}
}
\description{
	R interface via C to Cox DBDC algorithm in Fortran
}
\examples{
	data(ex)
	blassocox(x=ex_X[1:5,1:5], y=ex_Y[1:5,], kits=ex_K[1:5,1:5], costs=ex_c[1:5])
}
