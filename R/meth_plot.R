#' @title Plotting splines 
#'
#' @description The method provides graphical visualization of a \code{Splinets}-class object. 
#' @param object \code{Splinets} object;
#' @param x vector, specifying where the splines will be evaluated for the plots.
#' @param sID vector, specifying indices of splines in the splinet object to be plotted;
#' @param vknots logic, indicates if add auxiliary lines to highlight the positions of knots; The default is \code{TRUE}. 
#' @param type string, controls the layout of graphs; The following options are available
#' #' \itemize{
#'   \item \code{"stnd"} \code{object@type="dspnt"} or \code{="spnt"}, then it is ploted over 
#'   the dyadic net of supports, other bases are on a single plot with information about 
#'   the basis,
#'   \item \code{"simple"} all the objects are ploted in a single plot with minimal information,
#'   \item \code{"dyadic"} if not \code{object@type="sp"}, then the plot is over the dyadic net of supports.
#' } 
#' @param mrgn number, specyfying the margin size in the dydadic structure plot;
#' @param lwd positive integer, the line width;
#' @param ... also other graphical parameters can be passed;
#' @return A plot visualizing a Splinet object. The entire set of splines will 
#' be displayed in a plot. 
#' @details Standard method of plotting splines in a \code{Splinet}-object. 
#' It plots a single graph with all splines in the object except if \code{type} field of the
#' object represents a splinet. In the latter case, the default is the (dyadic) net plot of 
#' the basis. The string argument \code{type} can overide this to produce a single plot.
#' Most of the standard R graphical parameters can be passed by this function. 
#' 
#' @export
#' @inheritSection Splinets-class References
#' 
#' @seealso \code{\link{evspline}} for manually evaluating splines in a \code{Splinet}-object;
#' \code{\link{Splinets-class}} for the definition of the \code{Splinet}-class;
#' \code{\link{lines,Splinets-method}} for adding graphs to exisiting plots;
#'
#' @example R/Examples/ExPlot.R 
#' 
#' 


setMethod(
  "plot","Splinets",
  function(object, x = NULL, sID = NULL, vknots = TRUE, type = "stnd",mrgn=2, lwd=2, ...){
    if(object@type == "sp" | type=="simple" ){
      plot.spline(object, x, sID, vknots, mrgn=mrgn, ...) #if x is NULL then the default sampling of 'evaluate' function is used for x
    } else if(object@type == "bs" & type== "stnd"){
      plot.basis(object)
    } else if(object@type %in% c("gsob", "twob") & type=="stnd"){
      plot.obasis(object)
    } else if(object@type %in% c("dspnt", "spnt") | (object@type == "bs" & type== "dyadic") ){
      plot.splinet(object,lwd,mrgn)
    }
  })

# The above method allowed to use existing generics 'plot' to create a new method. 
# The following that created a new generics stopped working for R version 4.03
# setGeneric("plot.Splinets",
#            function(object, x = NULL, sID = NULL, vknots = TRUE, type = "stnd", mrgn=2, ...){
#              standardGeneric("plot.Splinets")
#            }
# )
# 
# setMethod(
#   "plot.Splinets","Splinets",
#   function(object, x = NULL, sID = NULL, vknots = TRUE, type = "stnd",mrgn=2, lwd=2, ...){
#     if(object@type == "sp" | type=="simple" ){
#       plot.spline(object, x, sID, vknots, ...)
#     } else if(object@type == "bs" & type== "stnd"){
#       plot.basis(object)
#     } else if(object@type %in% c("gsob", "tob") & type=="stdn"){
#       plot.obasis(object)
#     } else if(object@type %in% c("dspnt", "ndspnt")){
#       plot.splinet(object,lwd,mrgn)
#     }
# })
