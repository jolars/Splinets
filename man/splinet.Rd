% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_splinet.R
\name{splinet}
\alias{splinet}
\title{B-splines and their orthogonalization}
\usage{
splinet(knots = NULL, smorder = 3, type = "spnt", Bsplines = NULL, norm = F)
}
\arguments{
\item{knots}{\code{n+2} vector, the knots (presented in the increasing order); It is not needed, when
\code{Bsplines} argumment if the latter is not \code{NULL}, in which the case the knots from this \code{Splinets}-object are inherited.}

\item{smorder}{integer,  the order of the splines, default is \code{smorder=3}; Again it is inherited from the
\code{Bsplines} argumment if the latter is not \code{NULL}.}

\item{type}{string, the type of the basis; The following choices are available 
\itemize{
  \item \code{'bs'} for the unorthogonalized B-splines,
  \item \code{'spnt'} for the orthogonal splinet (the default),
  \item \code{'gsob'} for the Gramm-Schmidt (one-sided) O-splines,
  \item \code{'twob'} for the two-sided O-splines.
 }}

\item{Bsplines}{\code{Splinet}-object, the basis of B-splines (if not \code{NULL}); The B-splines in the \code{Splinet}-object. 
When this argument isnot \code{NULL} the first two arguments  
are not needed as they will be inherited from  \code{Bsplines}.}

\item{norm}{logical, a flag to indicate if the output B-splines will be normalized;}
}
\value{
Either a list \code{list("bs"=Bsplines)} made of a single \code{Splinet}-object \code{Bsplines} 
if \code{type=='bs'}, which represents the B-splines (the B-splines are normalized or not, which depends
on the \code{norm}-flag) or a list of two \code{Splinets}-objects: \code{list("bs"=Bsplines,"os"=Splinet)}, 
where \code{Bsplines} are either computed (in the input \code{Bspline= NULL}) or taken from the input \code{Bspline}
(this output will be normalized or not depending on the \code{norm}-flag),
\code{Splinet} is the B-spline orthognalization determined by the input argument \code{type}.
}
\description{
The B-splines are either given in the input or generated inside the routine. Then, given 
the B-splines and the argument \code{type}, the routine additionally generates a \code{Splinets}-object
representing an orthonormal spline basis obtained from the 
orthonormalization of the B-splines. The orthonormal spline bases are obtained by one of the methods:
the Gram-Schmidt method, the two-sided method, and/or the splinet algorithm, which is the default method.
All spline bases are kept in the format of \code{Splinets}-objects.
}
\details{
The B-spline basis, if not given in 
the input, is computed 
from  the following recurrent (with respect to the smoothness order of the B-splines) formula
\deqn{
B_{l,k}^{\boldsymbol \xi }(x)=
\frac{x- {\xi_{l}}
 }{
{\xi_{l+k}}-{\xi_{l}}
}
B_{l,k-1}^{\boldsymbol \xi}(x)
+
 \frac{{\xi_{l+1+k}}-x }{ {\xi_{l+1+k}}-{\xi_{l+1}}}
 B_{l+1,k-1}^{\boldsymbol \xi}(x), l=0,\dots, n-k.
}{
 B_lk(x)=(x-\xi_l)/(\xi_{l+k}-\xi_l) * B_{lk-1}(x)
 +
 (\xi_{l+1+k}-x)/(\xi_{l+1+k}-\xi_{l+1}) * B_{l+1k-1}(x), l=0,\dots, n-k
 } 
 The dyadic algorithm that is implemented takes into account efficiencies due to the equally space knots 
(exhibited in the Toeplitz form of the Gram matrix) only if the problem is fully dyadic, i.e. if the number of 
the internal knots is \eqn{smorder*2^N-1}, for some integer \eqn{N}. To utilize this efficiency it may be advantageous 
for a large number of equally spaced knots to choose them so that their number follows the fully dyadic form.
Additional advantage of the dyadic form is the complete symmetry at all levels of the support.
}
\section{References}{

Liu, X., Nassar, H., Podg\eqn{\mbox{\'o}}{o}rski, K. (2019) Splinets -- efficient orthonormalization of the B-splines.  ArXiv, abs/1910.07341.

???? (2020) \code{Splinets} -- Local Taylor expansion driven representation of splines, their support sets and orthogonal bases. (Preprint)
}

\examples{
#--------------------------------------#
#----Splinet, equally spaced knots-----#
#--------------------------------------#
k=2 # order
n_knots = 12 # number of knots
xi = seq(0, 1, length.out = n_knots)

so = splinet(xi, k)

plot(so$bs) #Plotting B-splines
plot(so$os) #Plotting Splinet

#Verifying the orthogonalization
gm = gramian(so$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm)))

#An example of the dyadic structure with equally spaced knots
k=3
N=3
n_knots=2^N*k-1 #the number of internal knots for the dyadic case

xi = seq(0, 1, length.out = n_knots+2)

so = splinet(xi) 

plot(so$bs,type="simple",vknots=FALSE,lwd=3) #Plotting B-splines in a single simple plot
plot(so$os,type="simple",vknots=FALSE,lwd=3) 

plot(so$os,lwd=3,mrgn=2) #Plotting the splinet on the dyadic net of support intervals

so=splinet(xi, Bsplines=so$bs, type='gsob') #Obtaining the Gram-Schmidt orthogonalization
plot(so$os,type="simple",vknots=FALSE)      #Without computing B-splines again

so=splinet(xi, Bsplines=so$bs, type='twob') #Obtaining the symmetrize orthogonalization
plot(so$os,type="simple",vknots=FALSE)     

#-------------------------------------#
#---Splinet, unequally spaced knots---#
#-------------------------------------#
n_knots=30

xi = c(0, sort(runif(n_knots)), 1)

sone = splinet(xi, k)

plot(sone$bs, type='dyadic') #Plotting B-splines
plot(sone$os) #Plotting Splinet

#Verifying the orthogonalization
gm = gramian(sone$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm)))

#------------------------------------------#
#---Dyadic splinet, equally spaced knots---#
#------------------------------------------#
k = 2 # order
N = 3 # support level
n_so = k*(2^N-1) # number of splines in a dyadic structure with N and k
n_knots = n_so + k + 1 # number of knots
xi = seq(0, 1, length.out = n_knots)

sodyeq = splinet(xi, k)

plot(sodyeq$bs) #Plotting B-splines
plot(sodyeq$os) #Plotting Splinet

#Verifying the orthogonalization
gm = gramian(sodyeq$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm)))

#--------------------------------------------#
#---Dyadic splinet, unequally spaced knots---#
#--------------------------------------------#
xi = c(0, sort(runif(n_knots)), 1)

sody = splinet(xi, k)


plot(sody$bs) #Plotting B-splines
plot(sody$os) #Plotting Splinet

#Verifying the orthogonalization
gm = gramian(sody$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm)))

#-----------------------------------------#
#---Bspline basis, equally spaced knots---#
#-----------------------------------------#
n = 15
xi = seq(0,1,length.out = n+2)
order = 2

bs = splinet(xi, order, type = 'bs')

plot(bs$bs)

#---------------------------------------------#
#---Bspline basis, non-equally spaced knots---#
#---------------------------------------------#
n = 6
xi = c(0,sort(runif(n)),1)
order = 3

so = splinet(xi, order, type = 'bs') #unnormalized version
plot(so$bs) 

so1 = splinet(type='bs',Bsplines=so$bs,norm=TRUE) #normalized version
plot(so1$bs)

#-------------------------------------------------#
#---Gram-Schmidt osplines, equally spaced knots---#
#-------------------------------------------------#
so = splinet(xi, order,  type = 'gsob')

plot(so$bs)
plot(so$os)

#Using the previously generated B-splines and normalizing them
so1 = splinet(Bsplines=so$bs, type = "gsob",norm=TRUE) 

plot(so1$bs) #normalized B-splines
plot(so1$os) #the one sided osplines

gm = gramian(so1$os) #evaluation of the inner products
diag(gm)
sum(gm - diag(diag(gm))) #verification of the orghonoalization of the matrix

#-----------------------------------------------------#
#---Gram-Schmidt osplines, non-equally spaced knots---#
#-----------------------------------------------------#

so = splinet(Bsplines=sody$bs, type = 'gsob') #previously genereted Bsplines

plot(so$bs)
plot(so$os)

gm = gramian(so$os)
diag(gm)
sum(gm - diag(diag(gm)))

#---------------------------------------------#
#---Twosided osplines, equally spaced knots---#
#---------------------------------------------#

so = splinet(Bsplines=bs$bs, type = 'twob')
plot(so$os)

gm = gramian(so$os) #verification of the orthogonality
diag(gm)
sum(gm - diag(diag(gm)))

#-------------------------------------------------#
#---Twosided osplines, non equally spaced knots---#
#-------------------------------------------------#

so = splinet(Bsplines=sody$bs, type = 'twob')
plot(so$os)

gm = gramian(so$os) #verification of the orthogonality
diag(gm)
sum(gm - diag(diag(gm)))


}
\seealso{
\code{\link{project}} for projecting into the functional spaces spanned by the spline bases; 
\code{\link{lincomb}} for evaluation of a linear combination of splines;
\code{\link{seq2dyad}} for building the dyadic structure for a splinet of a given smoothness order;
\code{\link{plot,Splinets-method}} for visualisation of splinets; 
Splinets'-class defintion
}
