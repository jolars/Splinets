% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_evspline.R
\name{evspline}
\alias{evspline}
\title{Evaluating splines at given arguments.}
\usage{
evspline(object, sID = NULL, x = NULL, N = 250)
}
\arguments{
\item{object}{\code{Splinets} object;}

\item{sID}{vector of integers, the indicies specifying splines in the \code{Splinets} list 
to be evaluated; If \code{sID=NULL}, then all elements in the list are evaluated. Default 
is \code{NULL}}

\item{x}{vector, the arguments at which the splines are evaluated; If \code{x} is
\code{NULL}, then the splines are evaluated over a regular grid. Default is \code{x=NULL}.}

\item{N}{integer, the average number of points per an interval between to consequiteve knots, default \code{N = 250};}
}
\value{
An \code{length(x) x length(sID+1)} matrix. The argument values, in the first column, 
then, columnwise, values of the splines.
}
\description{
For a given \code{Splinets} object \code{S} and a vector of arguments \code{t}, 
the function returns the matrix of values for the splines in \code{S}. The evaluations are done 
through the Taylor expansion on the \eqn{i}th interval 
\eqn{[\xi_i,\xi_{i+1}]}{[\xi_i,\xi_[i+1]]}
\deqn{S(t)=\sum_{j=0}^{k} s_{i j} \frac{(t-\xi_{i})^j}{j!}.}{S(t) = \sum(l=1:k) (t-\xi_i)^l * 1/l! * s_il.}
For the zero order splines which are discontinuous at the knots, at the LHS knots the value is taken as the RHS-limit, and at the RHS knots as the LHS-limit. 
The value at the central knot for the zero order and an odd number of knots case is zero.
}
\section{References}{

Liu, X., Nassar, H., Podg\eqn{\mbox{\'o}}{o}rski, K. (2019) Splinets -- efficient orthonormalization of the B-splines.  ArXiv, abs/1910.07341.

???? (2020) \code{Splinets} -- Local Taylor expansion driven representation of splines, their support sets and orthogonal bases. (Preprint)
}

\examples{
#---------------------------------------------#
#-- Example piecewise polynomial vs. spline --#
#---------------------------------------------#
n=20; k=3; xi=sort(runif(n+2))
sp=new("Splinets",knots=xi) 

#Randomly assigning the derivatives -- a very 'wild' function.
S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1))
sp@supp=list(t(c(1,n+2))); sp@smorder=k; sp@der[[1]]=S

y = evspline(sp)
plot(y,type = 'l',col='red')

#A correct spline object
nsp=is.splinets(sp)
sp2=nsp$robject 
y = evspline(sp2)
lines(y,type='l')

#---------------------------------------------#
#-- Example piecewise polynomial vs. spline --#
#---------------------------------------------#
#Gathering three 'Splinets' objects using three different
#method to correct the derivative matrix
n=17; k=4; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1

S=matrix(rnorm((n+2)*(k+1)),ncol=(k+1)) # generate a random matrix S

spl=construct(xi,k,S) #constructing the first correct spline
spl=gather(spl,construct(xi,k,S,mthd='CRFC')) #the second and the first ones
spl=gather(spl,construct(xi,k,S,mthd='CRLC')) #the third is added
y = evspline(spl, sID= 1)
plot(y,type = 'l',col='red')

y = evspline(spl, sID = c(1,3))
plot(y[,1:2],type = 'l',col='red')
points(y[,c(1,3)],type = 'l',col='blue')

#sID = NULL
y = evspline(spl)
plot(y[,1:2],type = 'l',col='red',ylim=range(y[,2:4]))
points(y[,c(1,3)],type = 'l',col='blue')
points(y[,c(1,4)],type = 'l',col='green')

#---------------------------------------------#
#--- Example with different support ranges ---#
#---------------------------------------------#
n=25; k=3; xi=sort(runif(n+2)); xi[1]=0; xi[n+2]=1
#Defining support ranges for three splines
supp=matrix(c(2,12,4,20,6,25),byrow=TRUE,ncol=2)
#Initial random matrices of the derivative for each spline
SS1=matrix(rnorm((supp[1,2]-supp[1,1]+1)*(k+1)),ncol=(k+1)) 
SS2=matrix(rnorm((supp[2,2]-supp[2,1]+1)*(k+1)),ncol=(k+1)) 
SS3=matrix(rnorm((supp[3,2]-supp[3,1]+1)*(k+1)),ncol=(k+1)) 
spl=construct(xi,k,SS1,supp[1,]) #constructing the first correct spline
nspl=construct(xi,k,SS2,supp[2,],'CRFC')
spl=gather(spl,nspl) #the second and the first ones
nspl=construct(xi,k,SS3,supp[3,],'CRLC')
spl=gather(spl,nspl) #the third is added

y = evspline(spl, sID= 1)
plot(y,type = 'l',col='red')

y = evspline(spl, sID = c(1,3))
plot(y[,1:2],type = 'l',col='red')
points(y[,c(1,3)],type = 'l',col='blue')

#sID = NULL -- all splines evaluated
y = evspline(spl)
plot(y[,c(1,3)],type = 'l',col='red',ylim=c(-1,1))
points(y[,1:2],type = 'l',col='blue')
points(y[,c(1,4)],type = 'l',col='green')

}
\seealso{
\code{\link{is.splinets}} for diagnostic of \code{Splinets}-objects;
\code{\link{plot,Splinets-method}} for plotting \code{Splinets}-objects;
}
