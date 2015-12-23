novelsetup <- setup.model(DRera=TRUE, treatSL=TRUE, treatnovel=TRUE)
values <- set.values()
mergedvalues <- append(append(values[[1]], values[[2]]), append(values[[3]], values[[4]]))
tallynames <- colnames(equilib()$log)[-(1:(length(dssetup$statenames)+1))]
elementnames <- set.novelvalues()$elementnames

drout <- screendrout()


# limits output to those simulations with DR fraction of incidence in targetepi range

evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DS", DST="DSTall") # can also specify ids and idr to run just a subset of drout
evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DR", DST="DSTall", idr=1:10) # can also specify ids and idr to run just a subset of drout
evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DS", DST="DSTnone", idr=1:10) # can also specify ids and idr to run just a subset of drout
evaltrp(genericvalues = mergedvalues, drsetup = drsetup, drout=drout, targetpt="DR", DST="DSTnone", idr=1:10) # can also specify ids and idr to run just a subset of drout using ids and idr
