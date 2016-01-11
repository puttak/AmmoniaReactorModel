pr.effectiveness.factor <- function(x){return(1)}
effectiveness.factor <- function(stream, reaction, catalyst){return(1)}

catalystParams <- c("intrPorosity", "porosity", "pelletD")
setClass("Catalyst",
         representation = representation(
             geometry = "numeric"
         ),
         validity = function(object){
             if ( !element.match(names(object@geometry), catalystParams) )
                 stop("Wrong names in geometry")
             if ( any(object@geometry < 0))
                 stop("Negative parameters of catalyst geometry")
             if (object@geometry[["intrPorosity"]] > 1)
                 stop("IntrPorosity is greater than 0")
             if (object@geometry[["porosity"]] > 1)
                 stop("porosity is greater than 0")
             return(TRUE)
         }
)

setMethod("initialize",
          "Catalyst",
          definition = function(.Object, geometry){
                .Object@geometry = geometry
                validObject(.Object)
                return(.Object)
          }
)