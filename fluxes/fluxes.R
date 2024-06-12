#library(openxlsx)
#library(tidyverse)
library(fluxweb)

# This script relies on products of "attributes.R" and "foodwebs.R". 
# Make sure to run those first.

fluxes = vector(mode = "list", 
                length=length(att))
fluxxes = data.frame(Unit_quarter = names(att),
                     Total.flux = numeric(length(att)),
                     Herb.flux = numeric(length(att)),
                     Detr.flux = numeric(length(att)),
                     Micr.flux = numeric(length(att)),
                     Pred.flux = numeric(length(att)))
for (i in 1:length(att)) {
  
  fluxes[[i]] <- fluxing(# the interaction matrix
    # columns are expected diet composition
    # of each consumer, given 
    # general preference for resource types,
    # predator-prey mass ratios, 
    # prey defences, 
    # vertical stratification
    # and relative availability (biomass)
    mat = web[[i]],
    biomasses = NULL, 
    # population level losses
    losses = att[[i]]$pop.Loss.Jh,
    # assimilation efficiency of node as resource
    efficiencies = att[[i]]$efficiency,
    # the matrix already accounts for relative biomass
    bioms.prefs = F,
    # losses are already of the population
    bioms.losses = F,
    ef.level = "prey")
  
  fluxxes[i,"Total.flux"] = sum(fluxes[[i]])
  fluxxes[i,"Herb.flux"] = sum(fluxes[[i]]["plants",])
  fluxxes[i,"Detr.flux"] = sum(fluxes[[i]]["detritus",])
  fluxxes[i,"Micr.flux"] = sum(fluxes[[i]]["microbes",])
  fluxxes[i,"Pred.flux"] = sum(fluxes[[i]][4:nrow(fluxes[[i]]),])

}

#write.csv(fluxxes, "flux_aggregates.csv", row.names = FALSE)



