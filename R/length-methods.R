setMethod("length", signature("KCData"), function(x){
    length(x@data)
})