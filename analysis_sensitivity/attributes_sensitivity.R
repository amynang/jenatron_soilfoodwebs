library(openxlsx)
library(tidyverse)

nem.mes.mac.taxonomy = read.csv("datafiles/Potapov2022.csv")


nematodes = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "Nematodes"]
meso = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "mesofauna"]
macro = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "macrofauna"]
# indexing for the metabolic loss regressions
insecta      = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class == "Insecta"]
araneae      = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Araneae"]
chilopoda    = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class == "Chilopoda"]
progoneata   = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class %in% c("Pauropoda","Symphyla","Diplopoda")]
mesostigmata = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Mesostigmata"]
oribatida    = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Sarcoptiformes"]
prostigmata  = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Trombidiformes"]

# per mesocosm abundance data
nem.mes.mac.dens = read_csv("datafiles/jenatron_soil_fauna_abundance.csv")
# meso and macro fauna bodymass
mes.mac.mean.masses = read_csv("datafiles/jenatron_soil_fauna_bodymass.csv")
across.unit.avg = mes.mac.mean.masses %>% 
  group_by(taxon_name) %>% 
  summarise(Bodymass.mg = mean(Bodymass.mg)) %>% 
  ungroup()
# nematode bodymass
ecophys = read.csv("datafiles/jenatron_nematode_bodymass.csv")

# split density dataframe to mesocosm-specific dataframes
att = nem.mes.mac.dens %>% 
  add_column(.before = "Achromadora",
             plants = 1, # dummy values so that basals are part of the foodweb
             detritus = 1,
             microbes = 1) %>% 
  pivot_longer(2:length(names(.)),
               names_to = "taxon_name",
               values_to = "density")%>% 
  filter(density>0) %>% 
  # finally, create a list: each element is a mesocosm
  split(., with(.,Unit_quarter))



attt = vector(mode = "list", 250)
set.seed(321)
# with Poisson samples
#for (j in 1:nrow(attt[[i]])) {
#  attt[[i]][j,"density"] = rpois(1, round(attt[[i]][j,"density"]))
#}

# a list of lists; 250 versions of 96 mesocosm specific dataframes
for (i in 1:length(attt)) {
  attt[[i]] = att
  for (j in 1:length(att)) {
    attt[[i]][[j]] = attt[[i]][[j]] %>% 
      mutate(.after = taxon_name,
             Agility = nem.mes.mac.taxonomy$Agility[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
             PhysicalProtection = nem.mes.mac.taxonomy$PhysicalProtection[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
             Metabolites = nem.mes.mac.taxonomy$Metabolites[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
             above = nem.mes.mac.taxonomy$above[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
             epi = nem.mes.mac.taxonomy$epi[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
             hemi = nem.mes.mac.taxonomy$hemi[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)],
             eu = nem.mes.mac.taxonomy$eu[match(.$taxon_name, nem.mes.mac.taxonomy$taxon_name)]) %>% 
      mutate(Bodymass.mg = case_when(.$taxon_name %in% nematodes ~ ecophys$AvgMass[match(.$taxon_name, ecophys$taxon_name)],
                                     .$taxon_name %in% c(meso,macro) ~ mes.mac.mean.masses$Bodymass.mg[match(paste(.$Unit_quarter, .$taxon_name), 
                                                                                                             paste(mes.mac.mean.masses$Unit_quarter, mes.mac.mean.masses$taxon_name))],
                                     TRUE ~ 1)) %>% 
      # mean imputation for taxa  that were counted in a mesocosm but not measured (11 cases)
      mutate(Bodymass.mg = ifelse(is.na(Bodymass.mg), across.unit.avg$Bodymass.mg[match(.$taxon_name, 
                                                                                        across.unit.avg$taxon_name)], 
                                  Bodymass.mg)) %>% 
      # note to future me: relax length(density) does the right thing https://stackoverflow.com/a/63534386  
      # with Poisson samples
      #mutate(density = rpois(length(density), round(density))) %>% 
      # with uniform samples from [n-0.2*n, n+0.2*n]
      mutate(density = runif(length(density),
                             density -
                               .2*density,
                             density +
                               .2*density)) %>% 
      mutate(ind.Loss.Jh = case_when(# based on Ehnes 2011 10.1111/j.1461-0248.2011.01660.x
        # group specific coefficients (phylogenetic model)
        .$taxon_name %in% araneae ~ exp(24.581475 + .5652537*log(Bodymass.mg) - .7093476*(1/(8.62*1e-5*(20+273.15)))),
        .$taxon_name %in% insecta ~ exp(21.972050 + .7588950*log(Bodymass.mg) - .6574038*(1/(8.62*1e-5*(20+273.15)))),
        .$taxon_name %in% chilopoda ~ exp(28.252911 + .5580991*log(Bodymass.mg) - .8030069*(1/(8.62*1e-5*(20+273.15)))),
        .$taxon_name %in% progoneata ~ exp(22.347024 + .5713411*log(Bodymass.mg) - .6700449*(1/(8.62*1e-5*(20+273.15)))),
        .$taxon_name %in% oribatida ~ exp(22.022770 + .6793706*log(Bodymass.mg) - .7060855*(1/(8.62*1e-5*(20+273.15)))),
        .$taxon_name %in% prostigmata ~ exp(10.281495 + .6599399*log(Bodymass.mg) - .4125318*(1/(8.62*1e-5*(20+273.15)))),
        .$taxon_name %in% mesostigmata ~ exp(9.6740230 + .6904864*log(Bodymass.mg) - .3792541*(1/(8.62*1e-5*(20+273.15)))),
        # The general relationship (linear model)
        TRUE ~ exp(23.055335 + .6950710*log(Bodymass.mg) - .6864200*(1/(8.62*1e-5*(20+273.15))))),
        ind.Loss.Jh = runif(length(ind.Loss.Jh), 
                            ind.Loss.Jh - .2*ind.Loss.Jh, 
                            ind.Loss.Jh + .2*ind.Loss.Jh),
        Biomass.mg = density*Bodymass.mg,
        pop.Loss.Jh = if_else(.$taxon_name %in% c("plants","detritus","microbes"), 0, # not necessary, basal losses are ignored
                              density*ind.Loss.Jh)) %>% 
      #https://www.ecologycenter.us/species-richness/the-importance-of-transfer-efficiencies-in-determining-energy-pathways.html
      mutate(efficiency = case_when(# based on Lang 2017 10.1111/oik.04419
        .$taxon_name == "detritus" ~ .158,
        .$taxon_name == "plants" ~ .545,
        TRUE ~ .906))
  }
  ############################## Show loop progress ##############################
  cat('\014')
  cat(paste0(round((i/250) * 100), '% completed'))
  Sys.sleep(.05)
  if (i == 250) cat(': Done')
  ################################################################################ 
}

# with uniform samples from [n-0,25*n, n+0,25*n]
for (i in 1:length(attt)) {
  attt[[i]] = as.data.frame(att)
  for (j in 1:nrow(attt[[i]])) {
    attt[[i]][j,"density"] = runif(1,
                                   round(attt[[i]][j,"density"]) -
                                     .25*round(attt[[i]][j,"density"]),
                                   round(attt[[i]][j,"density"]) +
                                     .25*round(attt[[i]][j,"density"]))
  }
}

longg = vector(mode = "list", 250)
for (i in 1:length(attt)) {
  longg[[i]] = do.call(rbind.data.frame, attt[[i]])
}

