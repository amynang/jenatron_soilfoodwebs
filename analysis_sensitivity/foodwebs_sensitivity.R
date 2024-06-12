library(tidyverse)

# This script relies on products of "attributes_sensitivity.R". Make sure to run that one first. 

mat = matrix(0,
             nrow = length(nem.mes.mac.taxonomy$taxon_name)+3,
             ncol = length(nem.mes.mac.taxonomy$taxon_name)+3,
             dimnames = list(c("plants","detritus","microbes",
                               nem.mes.mac.taxonomy$taxon_name),
                             c("plants","detritus","microbes",
                               nem.mes.mac.taxonomy$taxon_name)))

nematodes = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "Nematodes"]
meso = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "mesofauna"]
macro = nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$group == "macrofauna"]

# Nematodes
mat["plants",    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$Abbreviation == "Ne-H"]] = 1
mat["microbes", nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$Abbreviation == "Ne-B"]] = 1
mat["microbes",    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$Abbreviation == "Ne-F"]] = 1
mat[c("microbes",
      nematodes),
    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$Abbreviation == "Ne-O"]] = c(1/2,
                                                                                      rep((1/2)/length(nematodes),
                                                                                          length(nematodes)))
mat[nematodes,    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$Abbreviation == "Ne-P"]] = rep(1/length(nematodes),
                                                                                                      length(nematodes))
# Gamasids
mat[c(nematodes,
      meso,macro),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$Abbreviation == "Me"]] = rep(1/length(c(nematodes,meso,macro)),  
                                                                                                        length(c(nematodes,meso,macro)))
# Oribatids
# fungivores Maraun 2023
mat["microbes",    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name %in% c("Astegistidae", 
                                                                                          "Autognethidae",
                                                                                          "Protoribatidae")]] = 1
# fungivores and decomp Maraun 2023
mat[c("microbes",
      "detritus"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name %in% c("Euphthiracaridae",
                                                                                             "Oribatulidae")]] = .5
# fungi & pred
mat[c("microbes",
      nematodes),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Oppiidae"]] = c(.5,
                                                                                                         rep(.5/length(c(nematodes)),  
                                                                                                             length(c(nematodes))))
mat[c("microbes",
      nematodes),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Hypochthoniidae"]] = c(.5,
                                                                                                                rep(.5/length(c(nematodes)), 
                                                                                                                    length(c(nematodes))))
mat[c("detritus","microbes",
      nematodes),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Scheloribatidae"]] = c(1/3,1/3,
                                                                                                                rep((1/3)/length(c(nematodes)), 
                                                                                                                    length(c(nematodes))))
mat[c("detritus",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$Abbreviation == "Ori" &
                                                        !(nem.mes.mac.taxonomy$taxon_name %in% c("Astegistidae",
                                                                                                 "Autognethidae",
                                                                                                 "Protoribatidae",
                                                                                                 "Euphthiracaridae",
                                                                                                 "Oribatulidae",
                                                                                                 "Oppiidae","Scheloribatidae",
                                                                                                 "Hypochthoniidae"))]] = .5
# Astigmatids
mat[c("detritus",
      "microbes",
      nematodes),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Acaridae"]] = c(.4,.4,
                                                                                                         rep(.2/length(nematodes),
                                                                                                             length(nematodes)))

# Prostigmatids
mat[c("microbes",
      nematodes),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name %in% c("Eupodidae",	
                                                                                            "Tydaeidae")]] = c(.5,
                                                                                                               rep(.5/length(nematodes),
                                                                                                                   length(nematodes)))
mat["plants",    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Penthaleidae"]] = 1

#10.1657/1523-0430(2006)38[292:DOPRMA]2.0.CO;2  #10.11646/zootaxa.4176.1.1
mat[c(nematodes,
      meso,
      macro),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name %in%  c("Rhagidiidae",	    
                                                                                         "Stigmaeidae")]] = rep(1/length(c(nematodes,meso,macro)), 
                                                                                                                length(c(nematodes,meso,macro)))
mat[c("plants",
      "microbes",
      nematodes, 
      meso,
      macro),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$Abbreviation == "Prost" &           
                                                   !(nem.mes.mac.taxonomy$taxon_name %in% c("Rhagidiidae",	
                                                                                            "Stigmaeidae",
                                                                                            "Penthaleidae",
                                                                                            "Eupodidae",	
                                                                                            "Tydaeidae"))]] =c(.2,.4,
                                                                                                               rep((.4)/length(c(nematodes,meso,macro)),
                                                                                                                   length(c(nematodes,meso,macro))))
# Collembola

mat[c("plants",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% c("Entomobryidae", 
                                                                                         "Katiannidae", 
                                                                                         "Sminthurididae")]] = c(1/3,2/3)
mat[c("plants",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Isotomidae_large"]] = c(1/3,2/3)

mat[c("detritus",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name %in% c("Isotomidae_medium", 
                                                                                             "Isotomidae_small")]] = c(1/3,2/3)
mat[c("detritus",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% c("Neelidae", 
                                                                                         "Arrhopalitidae")]] = c(1/3,2/3)

mat["microbes",    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% "Hypogastruridae"]] = 1

mat[c("plants",                                             
      "detritus",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% c("Onychiuridae")]] = 1/3

mat[c("detritus",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% c("Tullbergiidae")]] = c(1/3,2/3)


# Pauropoda
mat[c("plants",
      "detritus",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Pauropoda"]] = 1/3
# Protura
mat["microbes",    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Protura"]] = 1
# Symphyla
mat[c("plants",
      "detritus",
      nematodes, 
      meso,
      macro),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Symphyla"]] = c(1/3,1/3,rep((1/3)/length(c(nematodes,meso,macro)),
                                                                                                                 length(c(nematodes,meso,macro))))
# Spiders
mat[c(meso,macro),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Araneae"]] = rep(1/length(c(meso,macro)), 
                                                                                                      length(c(meso,macro)))
# Centipedes
mat[c(meso,macro),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class == "Chilopoda"]] = rep(1/length(c(meso,macro)), 
                                                                                                        length(c(meso,macro)))
# Diplopoda
mat[c("detritus",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class == "Diplopoda"]] = c(.75,.25)
# Diplura
mat[c("plants",
      "detritus",
      "microbes",
      nematodes),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% "Campodeidae"]] = c(.25,.25,.4, rep((.1)/length(c(nematodes)),
                                                                                                                          length(c(nematodes))))  
# Gastropoda
mat[c("plants",
      "detritus",
      "microbes"),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$class == "Gastropoda"]] = c(.2,.4,.4)

# Coleoptera
mat["plants",    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% "Curculionidae"]] = 1

mat[c("plants",
      "detritus",
      "microbes",
      meso,macro),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% "Elateridae"]] = c(1/4,1/4,1/4,rep((1/4)/length(c(meso,macro)), 
                                                                                                                          length(c(meso,macro))))
mat[c("microbes",
      meso,macro),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% "Staphylinidae"]] = c(.2,rep(.8/length(c(meso,macro)), 
                                                                                                                    length(c(meso,macro))))
# Diptera
mat[c("detritus",
      "microbes",
      nematodes),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$family %in% "Chironomidae"]] = c(1/3,1/3,rep((1/3)/length(nematodes),
                                                                                                                       length(nematodes)))
mat[c("plants",
      "detritus",
      "microbes",                                                                  
      nematodes),    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$taxon_name == "Nematocera_larva"]] = c(1/4,1/4,1/4,rep((1/4)/length(nematodes),
                                                                                                                                 length(nematodes)))
# Hemiptera (Sernoryncha)
mat["plants",    nem.mes.mac.taxonomy$taxon_name[nem.mes.mac.taxonomy$order == "Hemiptera"]] = 1


# now we split to mesocosm-specific matrices 

webb = vector(mode = "list", 
             length=length(attt))

for (h in 1:length(webb)) { 
  for (i in 1:length(attt[[h]])) {
    webb[[h]][[i]] = mat[att[[i]]$taxon_name, 
                         att[[i]]$taxon_name]
  }
}

# Li et al. 2023
# re-estimated for terrestrial invertebrates 
P_ij = function(pred_mass.mg, prey_mass.mg) {
  
  m_i = log10(pred_mass.mg/1e3)
  m_j = log10(prey_mass.mg/1e3)
  

  
  beta_mu0       =  .35
  beta_mu1       = -.2235
  beta_sigma_sq0 = -.1167
  beta_sigma_sq1 = -.0106
  beta_theta0    =  .3742
  beta_theta1    =  .0105
  
  # log10 optimal size ratio of predator i
  mu_i = beta_mu0 + beta_mu1*m_i
  
  # feeding range of predator i
  squared.sigma_i = exp(beta_sigma_sq0 + beta_sigma_sq1*m_i)
  
  # feeding probability of i, if prey j matches the optimal size ratio
  theta_i =  (exp(beta_theta0 + beta_theta1*m_i)) / 
    (1 + exp(beta_theta0 + beta_theta1*m_i))
  
  p_ij = theta_i*exp( - (m_i - m_j - mu_i)^2/(2*squared.sigma_i))
  
  return(p_ij)
}

set.seed(321)
for (h in 1:length(webb)) { 
for (i in 1:length(att)) {
  n = nrow(att[[i]]) - 3
  bodymat = matrix(0,n,n)
  
  fauna = intersect(colnames(webb[[h]][[i]]), c(nematodes,meso,macro))
  # vertical stratification similarity
  ver = attt[[h]][[i]] %>% 
    filter(taxon_name %in% fauna) %>% 
    select(above,epi,hemi,eu) %>% 
    # dissimilarity
    vegan::vegdist(method = "bray") %>% 
    as.matrix() 
  # similarity
  vertical = 1- ver
  
  # which cells in the interaction matrix are non zero
  ind = which(webb[[h]][[i]][fauna,
                             fauna] != 0, arr.ind = T)
  # we will use this to standardize animal predation in omnivores
  std = colSums(mat[c(nematodes,meso,macro), 
                    c(nematodes,meso,macro)])[fauna]
  
  # vector of bodymasses
  # this needs to be foodweb-specific and same order as in the web[[i]]
  meanmass = attt[[h]][[i]]$Bodymass.mg[attt[[h]][[i]]$taxon_name %in% fauna]
  
  for (j in 1:nrow(ind)) { # for every predator-prey pair
    # we calculate prey suitability based on Li et al. 2023
    prop = P_ij(meanmass[ind[j,][2]], # predator
                meanmass[ind[j,][1]]) # prey
    bodymat[ind[j,][1],ind[j,][2]] = prop
  }
  
  # now we multiply by three vectors that modify this relationship based on prey
  # agility, physical protection or metabolites 
  # and finally, a matrix of vertical stratification similarity
  props = bodymat*
    attt[[h]][[i]]$Agility[attt[[h]][[i]]$taxon_name %in% fauna]*
    attt[[h]][[i]]$PhysicalProtection[attt[[h]][[i]]$taxon_name %in% fauna]*
    attt[[h]][[i]]$Metabolites[attt[[h]][[i]]$taxon_name %in% fauna]*
    vertical
  
  # regulating cannibalism 
  # (prevents "Error in fluxing : model chosen is unable to determine fluxes accoringly to data")
  diag(props) = 0*diag(props)
  
  # now we make predation on different taxa also depend on their relative biomass
  rel.biomass = as.vector(vegan::decostand(attt[[h]][[i]]$Biomass.mg[attt[[h]][[i]]$taxon_name %in% fauna],"total", 2))
  props = props*rel.biomass
  
  # first standardization: animal preferences sum to 1
  props = vegan::decostand(props,"total", 2)
  # second, we have them sum to the complement of whatever else they eat
  # so if we expect that an omnivore eats 50% plants and 50% animals, animal preferences sum to .5
  # we do this by multiplying every column in the matrix with the std vector
  props = props*rep(std, each=ncol(props))
  webb[[h]][[i]][fauna,
                 fauna] = props
  
}
  ############################## Show loop progress ##############################
  cat('\014')
  cat(paste0(round((h/250) * 100), '% completed'))
  Sys.sleep(.05)
  if (h == 250) cat(': Done')
  ################################################################################
}



