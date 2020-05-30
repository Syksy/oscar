\name{casso}
\title{casso placeholder}
\usage{
	casso(x, y, k, w)
}
\arguments{
	\item{x}{a vector of values}
}
\description{
	R interface via C to Cox DBDC algorithm in Fortran
}
\examples{
	data(ex)
	casso(x=ex_X, y=ex_Y, k=ex_K, w=ex_c)
}
