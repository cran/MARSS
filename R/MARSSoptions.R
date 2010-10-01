MARSSoptions <- function(..., method="kem") {
    if(!(method %in% allowed.methods)) stop("Stopped in MARSSoptions: method not allowed.",call.=FALSE)
    if ( missing(...) ) return(alldefaults[[method]])
    current <- alldefaults
    temp <- list(...)
    if (is.null(names(temp)) && is.null(names(unlist(temp)))) {
        arg <- unlist(temp)
    if (!all(arg %in% names(current[[method]]))) 
      stop(paste("Stopped in MARSSoptions: ", arg[!(arg %in% names(current[[method]]))],"is not an arg in alldefaults.  \nSee ?MARSSoptions for proper formatting."),call.=FALSE)
      switch(mode(arg),
            character = return(alldefaults[[method]][arg]),
            stop("invalid argument: ", sQuote(arg)))
    }
    n <- names(temp)
    # A series of error-checks
    if (is.null(n)) stop("Stopped in MARSSoptions: options must be listed by the format name=value. See ?MARSSoptions.", call.=FALSE)
    if (!all(n %in% names(current[[method]]))) 
      stop(paste("Stopped in MARSSoptions: ", n[!(n %in% names(current[[method]]))],"is not an arg in alldefaults.  \nSee ?MARSSoptions for proper formatting."),call.=FALSE)
    #set up the names in all defaults which must be list
    list.names=c()
    for(el.name in names(current[[method]])){
    if(is.list(alldefaults[[method]][[el.name]])) list.names=c(list.names,el.name)
    }
    for(el.name in n){
    if(el.name %in% list.names && !is.list(temp[[el.name]]))
      stop(paste("Stopped in MARSSoptions: ",el.name, " must be passed in as a list."), call.=FALSE)
    if(!(el.name %in% list.names) && is.list(temp[[el.name]]))
      stop(paste("Stopped in MARSSoptions: ",el.name, " cannot be a list."), call.=FALSE)
    }
    
    changed <- current[[method]][n]
    current[[method]][n[!(n %in% list.names)]] <- temp[n[!(n %in% list.names)]]
    for(el.name in n[n %in% list.names]){
       temp.n = names(temp[[el.name]])
       if (is.null(temp.n)) 
          stop("Stopped in MARSSoptions: options must be listed by the format name=value. See ?MARSSoptions.", call.=FALSE)        
       if (!all(temp.n %in% names(current[[method]][[el.name]]))) 
          stop(paste("Stopped in MARSSoptions: ", temp.n[!(temp.n %in% names(current[[method]][[el.name]]))]," arg is not in alldefaults$",el.name,".",sep=""),call.=FALSE)
       current[[method]][[el.name]][temp.n] <- temp[[el.name]][temp.n]
    }

    if (sys.parent() == 0) 
        env <- asNamespace("MARSS") 
    else 
        env <- parent.frame()
    unlockBinding("alldefaults",asNamespace("MARSS"))
    assign("alldefaults", current, envir = env)
    lockBinding("alldefaults",asNamespace("MARSS"))
    invisible(current)
}
