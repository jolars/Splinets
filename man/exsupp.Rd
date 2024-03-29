% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_exsupp.R
\name{exsupp}
\alias{exsupp}
\title{Correcting support sets and reshaping derivative matrix}
\usage{
exsupp(S, supp = NULL, epsilon = 1e-07)
}
\arguments{
\item{S}{\code{(m+2)x(k+1)} matrix, values of derivatives at the knots in the input support set 
which has the cardinality \code{m+2}; The matrix is assumed to be in the symmetric around center form with respect to provided 
support component.}

\item{supp}{\code{NULL} or \code{Nsupp x2} matrix of integers, of the indices of the endpoints of 
the input support intervals, where \code{Nsupp} is the number of the intervals; If the parameter is \code{NULL},
than the full support is assumed.}

\item{epsilon}{small positive number, threshold value of the norm of rows of \code{S}; If the norm
of a row of \code{S} is less than \code{epsilon}, then it will be viewed as a neglible and the knot is excluded from 
the support set.}
}
\value{
The list of two elements: \code{exsupp$rS} is the reduced derivative matrix from which the neglible rows, if any, have been removed
 and \code{exsupp$rsupp} the corresponding reduced support;
The output matrix is in the symmetric around the center form.
}
\description{
The function is adjusting for a potential reduction in the support sets due to negligibly small values of rows
in the derivative matrix. If the derivative matrix has a row equal to zero (or smaller than a neglible positive value) in the one-sided representation
of it (see the references and \code{\link{sym2one}}), then the corresponding knot should be removed
from the support set. The function can be used to reduce the neglible support components from a \code{Splinets}-object.
}
\details{
This function typically would be applied to an element in the list given by SLOT 
\code{der} of a \code{Splinets}-object. This reduces support sets by regions of negligible values
of a corresponding spline and its derivatives.
}
\section{References}{

Liu, X., Nassar, H., Podg\eqn{\mbox{\'o}}{o}rski, K. (2019) Splinets -- efficient orthonormalization of the B-splines.  ArXiv, abs/1910.07341.

???? (2020) \code{Splinets} -- Local Taylor expansion driven representation of splines, their support sets and orthogonal bases. (Preprint)
}

\examples{
#----------------------------------------------------#
#---Correcting support sets in a derivative matrix---#
#----------------------------------------------------#
n=20; k=3; xi=seq(0,1,by=1/(n+1)) #an even number of equally spaced knots 
set.seed(5)
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
spl=construct(xi,k,S) #this spline will be used below to construct a 'sparse' spline
is.splinets(spl) #verification
plot(spl)


xxi=seq(0,20,by=1/(n+1)) #large set of knots for construction of a sparse spline 
nn=length(xxi)-2
spspl=new('Splinets',knots=xxi,smorder=k) #generic object from the 'Splinets'-class
spspl@der[[1]]=matrix(0,ncol=(k+1),nrow=(nn+2)) #starting with zeros everywhere

spspl@der[[1]][1:(n+2),]=sym2one(spl@der[[1]]) #assigning local spline to a sparse spline at 
spspl@der[[1]][nn+3-(1:(n+2)),]=spspl@der[[1]][(n+2):1,] #the beginning and the same at the end
spspl@der[[1]]=sym2one(spspl@der[[1]],inv=TRUE) 
                                #at this point the object does not account for the sparsity

is.splinets(spspl) #a sparse spline on 421 knots with a non-zero terms at the first 22 
                   #and at the last 22 knots, the actual support set is not yet reported
plot(spspl)
plot(spspl,xlim=c(0,1)) #the local part of the sparse spline

exsupp(spspl@der[[1]]) #the actual support of the spline given the sparse derivative matrix

#Expanding the previous spline by building a slightly more complex support set
spspl@der[[1]][(n+1)+(1:(n+2)),]=sym2one(spl@der[[1]]) #double the first component of the
                                              #support because these are tangent supports 
spspl@der[[1]][(2*n+3)+(1:(n+2)),]=sym2one(spl@der[[1]]) #tdetect a single component of
                                                #the support with no internal knots removed
is.splinets(spspl)
plot(spspl)


es=exsupp(spspl@der[[1]])
es[[2]]   #the new support made of three components with the two first ones
          #separated by an interval with no knots in it 
          

spspl@der[[1]]=es[[1]]         #defining the spline on the evaluated actual support
spspl@supp[[1]]=es[[2]]
#Example with reduction of not a full support. 

xi1=seq(0,14/(n+1),by=1/(n+1)); n1=13; #the odd number of equally spaced knots 
S1=matrix(rnorm((n1+2)*(k+1)),ncol=(k+1))
spl1=construct(xi1,k,S1) #construction of a local spline
xi2=seq(16/(n+1),42/(n+1),by=1/(n+1)); n2=25; #the odd number of equally spaced knots 

S2=matrix(rnorm((n2+2)*(k+1)),ncol=(k+1))
spl2=construct(xi2,k,S2) #construction of a local spline

spspl@der[[1]][1:15,]=sym2one(spl1@der[[1]])
spspl@der[[1]][16,]=rep(0,k+1)
spspl@der[[1]][17:43,]=sym2one(spl2@der[[1]])
spspl@der[[1]][1:43,]=sym2one(spspl@der[[1]][1:43,],inv=TRUE)

is.splinets(spspl) #three intervals in the support are repported 

exsupp(spspl@der[[1]],spspl@supp[[1]])
}
\seealso{
\code{\link{Splinets-class}} for the description of the \code{Splinets}-class;
 \code{\link{sym2one}} the function
switching between the representations of a derivative matrix over a general support set;
\code{\link{lincomb}} evaluating a linear transformation of splines in a \code{Splinets}-object; 
\code{\link{is.splinets}} a diagnostic tool for \code{Splinets}-objects;
}
