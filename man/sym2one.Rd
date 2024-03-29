% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_sym2one.R
\name{sym2one}
\alias{sym2one}
\title{Switching between representations of the derivative matrices at knots of splines}
\usage{
sym2one(S, supp = NULL, inv = F)
}
\arguments{
\item{S}{\code{(m+2) x (k+1)} numeric matrix, the derivatives in one of the two representations;}

\item{supp}{\code{(Nsupp x 2)} or \code{NULL} matrix; row-wise the endpoint indices of the interval support; If it 
is equal to \code{NULL} (which is also the default), then the full support is assumed.}

\item{inv}{logical, if \code{F} (default), then the functiona assumes that the input is 
in the symmetric format and transforms it to the left-to-right format, if \code{T}, then 
the inverse transformation is applied.}
}
\value{
a matrix that is respective transformation of the input.
}
\description{
A technical but useful transformation of the matrix of derivatives form the one-sided 
to symmetric representations, or a reverse one. It allows for switching between the standard representation of the matrix 
of the derivatives for \code{Splinets} which is symmetric around the central knot(s) to the one-sided that yields 
the RHS limits at the knots, which is more convenient for computations.
}
\details{
The transformation essentially changes only the last column in \code{S}, i.e. the highest (discontinuous) derivatives so that 
the one-sided representation yields the right-hand-side limit. 
It is expected that the number of rows in \code{S} is the same as the total size of the support
as indicated by \code{supp}, i.e. if \code{supp!=NULL}, then \code{sum(supp[,2]-supp[,1]+1)=m+2}. 
If this is true all derivative submatrices of the components in \code{S} will be reversed.
However, this condition is not checked, which may lead to switch of 
the representations only for parts of the matrix \code{S}.
}
\section{References}{

Liu, X., Nassar, H., Podg\eqn{\mbox{\'o}}{o}rski, K. (2019) Splinets -- efficient orthonormalization of the B-splines.  ArXiv, abs/1910.07341.

???? (2020) \code{Splinets} -- Local Taylor expansion driven representation of splines, their support sets and orthogonal bases. (Preprint)
}

\examples{
#-----------------------------------------------------#
#-------Representations of derivatives at knots-------#
#-----------------------------------------------------#
n=10; k=3; xi=seq(0,1,by=1/(n+1)) #the even number of equally spaced knots 
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) #construction of a spline
a=spl@der[[1]]
b=sym2one(a)
aa=sym2one(b,inv=TRUE) # matrix 'aa' is the same as 'a'

n=11; xi2=seq(0,1,by=1/(n+1)) #the odd number of knots case
S2=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl2=construct(xi2,k,S2) #construction of a spline
a2=spl2@der[[1]]
b2=sym2one(a2)
aa2=sym2one(b2, inv=TRUE) # matrix 'aa2' is the same as 'a2'

#-----------------------------------------------------#
#--------------More complex support sets--------------#
#-----------------------------------------------------#
#Zero order splines, non-equidistant case, support with three components
n=43; xi=seq(0,1,by=1/(n+1)); k=3; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1;
support=list(matrix(c(2,14,17,30,32,43),ncol=2,byrow = TRUE))
#Third order splines
ssp=new("Splinets",knots=xi,supp=support,smorder=k) #with partial support

m=sum(ssp@supp[[1]][,2]-ssp@supp[[1]][,1]+1) #the total number of knots in the support
ssp@der=list(matrix(rnorm(m*(k+1)),ncol=(k+1)))  #the derivative matrix at random
IS=is.splinets(ssp) 
IS$robject@der
IS$robject@supp
b=sym2one(IS$robject@der[[1]],IS$robject@supp[[1]]) #the RHS limits at the knots
a=sym2one(b,IS$robject@supp[[1]],inv=TRUE) #is the same as the SLOT supp in IS@robject
}
\seealso{
\code{\link{Splinets-class}} for the description of the \code{Splinets}-class; 
\code{\link{is.splinets}} for diagnostic of \code{Splinets}-objects;
}
