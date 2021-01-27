#' @title Correcting support sets and reshaping derivative matrix 
#'
#'@description The function is adjusting for a potential reduction in the support sets due to negligibly small values of rows
#'in the derivative matrix. If the derivative matrix has a row equal to zero (or smaller than a neglible positive value) in the one-sided representation
#'of it (see the references and \code{\link{sym2one}}), then the corresponding knot should be removed
#'from the support set. The function can be used to reduce the neglible support components from a \code{Splinets}-object.  
#'
#' @param S  \code{(m+2)x(k+1)} matrix, values of derivatives at the knots in the input support set 
#' which has the cardinality \code{m+2}; The matrix is assumed to be in the symmetric around center form with respect to provided 
#' support component.
#' @param supp \code{NULL} or \code{Nsupp x2} matrix of integers, of the indices of the endpoints of 
#' the input support intervals, where \code{Nsupp} is the number of the intervals; If the parameter is \code{NULL},
#' than the full support is assumed.
#' @param epsilon small positive number, threshold value of the norm of rows of \code{S}; If the norm
#' of a row of \code{S} is less than \code{epsilon}, then it will be viewed as a neglible and the knot is excluded from 
#' the support set.  
#'
#' @return The list of two elements: \code{exsupp$rS} is the reduced derivative matrix from which the neglible rows, if any, have been removed
#'  and \code{exsupp$rsupp} the corresponding reduced support;
#' The output matrix is in the symmetric around the center form.
#' 
#' @details This function typically would be applied to an element in the list given by SLOT 
#' \code{der} of a \code{Splinets}-object. This reduces support sets by regions of negligible values
#' of a corresponding spline and its derivatives.  
#' @export
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{Splinets-class}} for the description of the \code{Splinets}-class;
#'  \code{\link{sym2one}} the function
#' switching between the representations of a derivative matrix over a general support set;
#' \code{\link{lincomb}} evaluating a linear transformation of splines in a \code{Splinets}-object; 
#' \code{\link{is.splinets}} a diagnostic tool for \code{Splinets}-objects;
#' @example R/Examples/ExExsupp.R
#' 
#' 
exsupp = function(S, supp=NULL , epsilon = 1e-7){
  m=dim(S)[1]-2
  old_supp=1:(m+2)
  k=dim(S)[2]-1 #The order of the spline
  SS=sym2one(S)
  if(!is.null(supp)){
    Nsupp=dim(supp)[1] #the number of support components
    dd=supp[,2]-supp[,1] #sizes of support intervals
    if(sum(supp[,2]-supp[,1]+1)!=(m+2)){
      stop("The support set is not compatible with the dimension of the derivative matrix.")
    }
    st = 1 # evaluating all the indices in the support set
    for(j in 1:Nsupp){ 
      en = st + dd[j]
      old_supp[st:en]=supp[j,1]:supp[j,2]
      st = en+1
    }
    SS=sym2one(S,supp) #Transforming from the symmetric form
  }
  #two subsequent zero rows in 1:k, zero at k+1 (highest derivative) of the first one, 
  #and non-zero highest derivative in the second one 
  #indicates the locations of the beginning of the new intervals
  SSS=rbind(rep(0,k+1),SS,c(rep(0,k),1)) #adding vectors of zeros at the beginning and zeros-and-one at the 
  #end to treat the beginning and the end the same as the inner points. 
  #Evaluating the sum of squares in the matrix of derivatives
  if(k>1){
    sqsums=(k+1)*rowMeans((SSS)^2) #this is a standard function in R 
    sqsums1=k*rowMeans((SSS[,1:k])^2) # and, as usually and stupidly, it does not work for one dimension. 
  }else{
    if(k==1){
      sqsums=(SSS[,1])^2+(SSS[,2])^2 
      sqsums1=(SSS[,1])^2
    }else{ #the zero order case
      sqsums=(SSS[,1])^2 #Sum of squares of values of the function at the knots
      sqsums1=rep(0,m+4)
    }
  }
  
  #Detecting the endings: the first term detects two consequitive rows with the zeroes except the last value in the second row
  zerosL=((sqsums[1:(m+2)]+sqsums1[2:(m+3)])<epsilon) * ((SSS[2:(m+3),(k+1)])^2>epsilon) #the LHS 
  #The two rows, with the second being all zeros and the first one not.
  zerosR = (sqsums[1:(m+2)]>epsilon) * (sqsums[2:(m+3)] < epsilon) #the RHS
  if(sum(zerosL)!=sum(zerosR)){stop('Something wrong, check if the matrix is a valid spline matrix, or try to change the threshold parameter.')}
  rsupp=cbind(old_supp[zerosL==1],old_supp[zerosR==1])
  
  #Removing the second from two consequitive zero rows
  zeros=(sqsums[1:(m+2)]+sqsums[2:(m+3)])<epsilon
  rS = SS[!zeros, ,drop=F] # reduced der matrix
  rS = sym2one(rS,rsupp,inv=T) #Transforming back to the symmetric form
  return(list(rS = rS, rsupp = rsupp))
}