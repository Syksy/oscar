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
	blassocox(x=ex_X, y=ex_Y, kits=ex_K, costs=ex_c)
}
