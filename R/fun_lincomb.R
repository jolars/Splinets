#' @title Linear transformation of splines.
#'
#' @description The following linear combination of the splines \eqn{S_j}  is computed  
#' \deqn{R_i=\sum_{j=0}^{d} a_{i j} S_j,\, i=1,\dots, l}
#'
#' @param object \code{Splinets} object containing \code{d} splines;
#' @param A \code{l x d} matrix; 
#' @param reduced logical; If \code{TRUE} (default), then the linear combination will be 
#' calculated through the der matrix based on actual support sets (recommended in the sparse function case), 
#' if \code{FALSE}, then the 
#' full support computations are used (can be faster for a lower dimension or non-sparse cases).
#' @param SuppContr logical; If \code{TRUE} (default), the true support is extracted, otherwise, full range 
#' is reported as the support. Applies only to the case when \code{reduced=FALSE}.   
#'
#' @return A \code{Splinet}-object that contains \code{l} splines obtained by linear combinations of  
#' using coefficients in rows of \code{A}. The  \code{SLOT type} of the output splinet objects is \code{sp}.
#' 
#' 
#' @export
#' @inheritSection Splinets-class References
#' @seealso \code{\link{exsupp}} for extracting the correct support; 
#' \code{\link{construct}} for building a valid splines; 
#' \code{\link{rspline}} for random generation of splines;
#' @example R/Examples/ExLincomb.R
#' 
#' 

lincomb = function(object, A, reduced = TRUE, SuppContr = TRUE){
  k = object@smorder
  S = object@der
  xi = object@knots
  supp = object@supp
  n = length(xi)-2
  if(dim(A)[2] != length(S)){
    stop("The number of columns of input A matrix is different from the number of splines in the input object")
  }
  
  if(reduced){
    # do linear combination using reduced der matrix
    temp = lincomb_supp(S, supp, A, n)
  } else{
    # do linear combination using full der matrix
    temp = lincomb_full(S, supp, k, n, A, SuppContr)
  }
  
  object@type = "sp"
  object@der = temp$S
  object@supp = temp$supp
  
  return(object)
}