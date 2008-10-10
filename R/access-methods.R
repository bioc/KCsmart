setMethod("[", "probeAnnotation", function(x, i, j = "missing", drop = "missing"){
    x@chromosome <- x@chromosome[i]
    x@maploc <- x@maploc[i]
    x@name <- x@name[i]
    x
})

setMethod("[[", signature("KCData"), function(x, i, j = "missing"){
   x@data[[i]]
})

setReplaceMethod("[[", signature("KCData"), function(x, i, j = "missing", value){
   x@data[[i]] <- value
   x
})