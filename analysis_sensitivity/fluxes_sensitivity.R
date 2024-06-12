library(openxlsx)
library(tidyverse)
library(fluxweb)

# This script relies on products of "attributes_sensitivity.R" and "foodwebs_sensitivity.R". 
# Make sure to run those first.


ffluxes = vector(mode = "list", 
                 length=length(attt))
for  (i in 1:length(attt)) {
  ffluxes[[i]] = vector(mode = "list", 
                   length=length(att))
}
ffluxxes = vector(mode = "list", 
                     length=length(attt))

for (i in 1: length(attt)) { 
  ffluxxes[[i]] = data.frame(Unit_quarter = names(att),
                             Total.flux = numeric(length(att)),
                             Herb.flux = numeric(length(att)),
                             Detr.flux = numeric(length(att)),
                             Micr.flux = numeric(length(att)),
                             Pred.flux = numeric(length(att)))
}

for (h in 1:length(attt)) { 
for (i in 1:length(att)) {
  
  ffluxes[[h]][[i]] <- fluxing(# the interaction matrix
    # columns are expected diet composition
    # of each consumer, given 
    # general preference for resource types,
    # predator-prey mass ratios, 
    # prey defences, 
    # vertical stratification
    # and relative availability (biomass)
    mat = webb[[h]][[i]],
    biomasses = NULL, 
    # population level losses
    losses = attt[[h]][[i]]$pop.Loss.Jh,
    # assimilation efficiency of node as resource
    efficiencies = attt[[h]][[i]]$efficiency,
    # the matrix already accounts for relative biomass
    bioms.prefs = F,
    # losses are already of the population
    bioms.losses = F,
    ef.level = "prey")
  
  ffluxxes[[h]][i,"Total.flux"] = sum(ffluxes[[h]][[i]])
  ffluxxes[[h]][i,"Herb.flux"] =  sum(ffluxes[[h]][[i]]["plants",])
  ffluxxes[[h]][i,"Detr.flux"] =  sum(ffluxes[[h]][[i]]["detritus",])
  ffluxxes[[h]][i,"Micr.flux"] =  sum(ffluxes[[h]][[i]]["microbes",])
  ffluxxes[[h]][i,"Pred.flux"] =  sum(ffluxes[[h]][[i]][4:nrow(ffluxes[[h]][[i]]),])

}
  ############################## Show loop progress ##############################
  cat('\014')
  cat(paste0(round((h/250) * 100), '% completed'))
  Sys.sleep(.05)
  if (h == 250) cat(': Done')
  ################################################################################
}
