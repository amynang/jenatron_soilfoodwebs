library(brms)
library(emmeans)
library(modelr)
library(tidybayes)

# create an empty list
datalist = vector(mode = "list", 
                  length=length(attt))
# mesocosm information
d.1 = read.xlsx("datafiles/PlotList_Jenatron.xlsx",
                sheet = "Tabelle4")[,-(8:9)]
d.1.1 = d.1 %>% mutate(.after = subunit,
                       Uq = paste0(unit,"_",subunit)) %>% 
  
  mutate(.after = soil.history,
         soil_history = case_when(grepl("B_", soil.history) ~ 0,
                                  grepl("E_", soil.history) ~ 1),
         plant_history = case_when(plant.history == "+" ~ 1,
                                   plant.history == "-" ~ 0)) %>% 
  mutate(.after = soil.history,
         treatment = case_when(soil_history == 0 &
                                 plant_history == 0 ~ "00",
                               soil_history == 0 &
                                 plant_history == 1 ~ "01",
                               soil_history == 1 &
                                 plant_history == 0 ~ "10",
                               soil_history == 1 &
                                 plant_history == 1 ~ "11") %>% as.factor()) %>% 
  mutate(.before = unit,
         Unit_quarter = str_c("U",unit,"q",subunit))
# datalist will be populated by 250 versions of the dataset.

for (i in 1:length(attt)) {
  datalist[[i]] = d.1.1 %>% 
    select(hall.block, unit, Unit_quarter, treatment, plant.diversity) %>% 
    right_join(., ffluxxes[[i]], by = join_by(Unit_quarter)) %>% 
    mutate(hall.block = as.factor(hall.block),
           unit= as.factor(unit),
           plant.diversity.sc = scale(plant.diversity))
}

# brm_multiple() compiles the model once and then fits it 
# to each of the 250 versions of the data contained in the datalist
total.m.2 = brm_multiple(bf(log(Total.flux) ~ 
                     treatment
                   + (1|hall.block/unit)),
                family = student(),
                chains = 4,
                cores = 4,
                iter = 8000,
                control = list(adapt_delta = 0.999),
                backend = "cmdstanr",
                seed = 321,
                data = datalist,
                # we want the results of 250 models
                combine = FALSE)

saveRDS(total.m.2, "totalm2sens.rds")
rm(total.m.2)
gc()

pred.m.2 = brm_multiple(bf(log(Pred.flux) ~ 
                              treatment
                            + (1|hall.block/unit)),
                         family = student(),
                         chains = 4,
                         cores = 4,
                         iter = 8000,
                         control = list(adapt_delta = 0.999),
                         backend = "cmdstanr",
                         seed = 321,
                         data = datalist,
                         # we want the results of 250 models
                         combine = FALSE)

saveRDS(pred.m.2, "predm2sens.rds")
rm(pred.m.2)
gc()

herb.m.2 = brm_multiple(bf(log(Herb.flux) ~ 
                             treatment
                           + (1|hall.block/unit)),
                        family = student(),
                        chains = 4,
                        cores = 4,
                        iter = 8000,
                        control = list(adapt_delta = 0.999),
                        backend = "cmdstanr",
                        seed = 321,
                        data = datalist,
                        # we want the results of 250 models
                        combine = FALSE)

saveRDS(herb.m.2, "herbm2sens.rds")
rm(herb.m.2)
gc()



micr.m.2 = brm_multiple(bf(log(Micr.flux) ~ 
                             treatment
                           + (1|hall.block/unit)),
                        family = gaussian(),
                        chains = 4,
                        cores = 4,
                        iter = 8000,
                        control = list(adapt_delta = 0.999),
                        backend = "cmdstanr",
                        seed = 321,
                        data = datalist,
                        # we want the results of 250 models
                        combine = FALSE)

saveRDS(micr.m.2, "micrm2sens.rds")
rm(micr.m.2)
gc()

detr.m.2 = brm_multiple(bf(log(Detr.flux) ~ 
                             treatment
                           + (1|hall.block/unit)),
                        family = student(),
                        chains = 4,
                        cores = 4,
                        iter = 8000,
                        control = list(adapt_delta = 0.999),
                        backend = "cmdstanr",
                        seed = 321,
                        data = datalist,
                        # we want the results of 250 models
                        combine = FALSE)

saveRDS(detr.m.2, "detrm2sens.rds")
rm(detr.m.2)
gc()


total.m.2 = readRDS("totalm2sens.rds")

#summary(total.m.2);pp_check(total.m.2, ndraws = 100)
#ci <- emmeans(total.m.2, specs = pairwise ~ treatment)
#summary(ci, point.est = mean)

# we want to examine the sensitivity of the mean difference
# between the two soil treatments to biomass and metabolic loss uncertainty
# we extract from each model the mean estimate, the lower and upper limits of
# the 95% credible interval
total.soil.contrast = data.frame(mean.est = numeric(length(total.m.2)),
                                 hdi95.lower = numeric(length(total.m.2)),
                                 hdi95.upper = numeric(length(total.m.2)))

for (i in 1:length(total.m.2)) {
  total.soil.contrast[i,] = total.m.2[[i]] %>%
    ref_grid() %>%
    add_grouping("soil",
                 "treatment",
                 factor(c("0", "0", "1", "1"),
                        levels = c("1", "0"))) %>%
    emmeans(~ soil) %>%
    contrast(method = "pairwise") %>%
    gather_emmeans_draws() %>%
    mean_hdi(.width = 0.95) %>% .[,2:4]
}

# we plot their distributions
s1 = total.soil.contrast %>%
  pivot_longer(1:3,
               names_to = "statistic",
               values_to = "value") %>%
  ggplot(aes(x = value,
             # conditionally color dots above/below the zero threshold
             fill = (value > 0),
             color = (value > 0))) +
  geom_dots(side = "topright",
            orientation = "y",
            dotsize = 1) +
  facet_wrap(~factor(statistic,
                     levels=c('hdi95.lower',
                              'mean.est',
                              'hdi95.upper')),
             scales = 'free_x') +
  scale_color_manual(values = c("#b3cf99","#B22222")) +
  scale_fill_manual(values = c("#b3cf99","#B22222")) +
  ggtitle("Total flux - mean difference of soil treatments\nmain text: 0.71 [0.44, 0.98]") +
  theme_bw() +
  theme(legend.position = "none")

rm(total.m.2)
gc()

pred.m.2 = readRDS("predm2sens.rds")

pred.soil.contrast = data.frame(mean.est = numeric(length(pred.m.2)),
                                 hdi95.lower = numeric(length(pred.m.2)),
                                 hdi95.upper = numeric(length(pred.m.2)))

for (i in 1:length(pred.m.2)) {
  pred.soil.contrast[i,] = pred.m.2[[i]] %>%
    ref_grid() %>%
    add_grouping("soil",
                 "treatment",
                 factor(c("0", "0", "1", "1"),
                        levels = c("1", "0"))) %>%
    emmeans(~ soil) %>%
    contrast(method = "pairwise") %>%
    gather_emmeans_draws() %>%
    mean_hdi(.width = 0.95) %>% .[,2:4]
}

# we plot their distributions
s2 = pred.soil.contrast %>%
  pivot_longer(1:3,
               names_to = "statistic",
               values_to = "value") %>%
  ggplot(aes(x = value,
             # conditionally color dots above/below the zero threshold
             fill = (value > 0),
             color = (value > 0))) +
  geom_dots(side = "topright",
            orientation = "y",
            dotsize = 1) +
  facet_wrap(~factor(statistic,
                     levels=c('hdi95.lower',
                              'mean.est',
                              'hdi95.upper')),
             scales = 'free_x') +
  scale_color_manual(values = c("#b3cf99","#B22222")) +
  scale_fill_manual(values = c("#b3cf99","#B22222")) +
  ggtitle("Predation flux - mean difference of soil treatments\nmain text: 0.6 [0.3, 0.93]") +
  theme_bw() +
  theme(legend.position = "none")


rm(pred.m.2)
gc()
herb.m.2 = readRDS("herbm2sens.rds")

herb.soil.contrast = data.frame(mean.est = numeric(length(herb.m.2)),
                                hdi95.lower = numeric(length(herb.m.2)),
                                hdi95.upper = numeric(length(herb.m.2)))

for (i in 1:length(herb.m.2)) {
  herb.soil.contrast[i,] = herb.m.2[[i]] %>%
    ref_grid() %>%
    add_grouping("soil",
                 "treatment",
                 factor(c("0", "0", "1", "1"),
                        levels = c("1", "0"))) %>%
    emmeans(~ soil) %>%
    contrast(method = "pairwise") %>%
    gather_emmeans_draws() %>%
    mean_hdi(.width = 0.95) %>% .[,2:4]
}

# we plot their distributions
s3 = herb.soil.contrast %>%
  pivot_longer(1:3,
               names_to = "statistic",
               values_to = "value") %>%
  ggplot(aes(x = value,
             # conditionally color dots above/below the zero threshold
             fill = (value > 0),
             color = (value > 0))) +
  geom_dots(side = "topright",
            orientation = "y",
            dotsize = 1) +
  facet_wrap(~factor(statistic,
                     levels=c('hdi95.lower',
                              'mean.est',
                              'hdi95.upper')),
             scales = 'free_x') +
  scale_color_manual(values = c("#b3cf99","#B22222")) +
  scale_fill_manual(values = c("#b3cf99","#B22222")) +
  ggtitle("Herbivory flux - mean difference of soil treatments\nmain text: 0.96 [0.53, 1.39]") +
  theme_bw() +
  theme(legend.position = "none")



rm(herb.m.2)
gc()
micr.m.2 = readRDS("micrm2sens.rds")

micr.soil.contrast = data.frame(mean.est = numeric(length(micr.m.2)),
                                hdi95.lower = numeric(length(micr.m.2)),
                                hdi95.upper = numeric(length(micr.m.2)))

for (i in 1:length(micr.m.2)) {
  micr.soil.contrast[i,] = micr.m.2[[i]] %>%
    ref_grid() %>%
    add_grouping("soil",
                 "treatment",
                 factor(c("0", "0", "1", "1"),
                        levels = c("1", "0"))) %>%
    emmeans(~ soil) %>%
    contrast(method = "pairwise") %>%
    gather_emmeans_draws() %>%
    mean_hdi(.width = 0.95) %>% .[,2:4]
}

# we plot their distributions
s4 = micr.soil.contrast %>%
  pivot_longer(1:3,
               names_to = "statistic",
               values_to = "value") %>%
  ggplot(aes(x = value,
             # conditionally color dots above/below the zero threshold
             fill = (value > 0),
             color = (value > 0))) +
  geom_dots(side = "topright",
            orientation = "y",
            dotsize = 1) +
  facet_wrap(~factor(statistic,
                     levels=c('hdi95.lower',
                              'mean.est',
                              'hdi95.upper')),
             scales = 'free_x') +
  scale_color_manual(values = c("#b3cf99","#B22222")) +
  scale_fill_manual(values = c("#b3cf99","#B22222")) +
  ggtitle("Microbivory flux - mean difference of soil treatments\nmain text: 0.59 [0.4, 0.78]") +
  theme_bw() +
  theme(legend.position = "none")


rm(micr.m.2)
gc()
detr.m.2 = readRDS("detrm2sens.rds")

detr.soil.contrast1 = data.frame(mean.est = numeric(length(detr.m.2)),
                                hdi95.lower = numeric(length(detr.m.2)),
                                hdi95.upper = numeric(length(detr.m.2)))
detr.soil.contrast2 = data.frame(mean.est = numeric(length(detr.m.2)),
                                 hdi90.lower = numeric(length(detr.m.2)),
                                 hdi90.upper = numeric(length(detr.m.2)))
for (i in 1:length(detr.m.2)) {
  contr = detr.m.2[[i]] %>%
    ref_grid() %>%
    add_grouping("soil",
                 "treatment",
                 factor(c("0", "0", "1", "2"),
                        levels = c("1", "0", "2"))) %>%
    emmeans(~ soil) %>%
    contrast(method = "pairwise") %>%
    gather_emmeans_draws() %>%
    mean_hdi(.width = c(.9,.95))
  
  detr.soil.contrast1[i,] = contr[5,2:4]
  detr.soil.contrast2[i,] = contr[3,2:4]
  
}

# we plot their distributions
s5 = detr.soil.contrast1 %>%
  pivot_longer(1:3,
               names_to = "statistic",
               values_to = "value") %>%
  ggplot(aes(x = value,
             # conditionally color dots above/below the zero threshold
             fill = (value > 0),
             color = (value > 0))) +
  geom_dots(side = "topright",
            orientation = "y",
            dotsize = 1) +
  facet_wrap(~factor(statistic,
                     levels=c('hdi95.lower',
                              'mean.est',
                              'hdi95.upper')),
             scales = 'free_x') +
  scale_color_manual(values = c("#b3cf99","#B22222")) +
  scale_fill_manual(values = c("#b3cf99","#B22222")) +
  ggtitle("Detritivory flux - mean difference of 'soil but no plant history' to 'no soil history'\nmain text: -0.83 [-1.22, -0.45]") +
  theme_bw() +
  theme(legend.position = "none")

s6 = detr.soil.contrast2 %>%
  pivot_longer(1:3,
               names_to = "statistic",
               values_to = "value") %>%
  ggplot(aes(x = value,
             # conditionally color dots above/below the zero threshold
             fill = (value > 0),
             color = (value > 0))) +
  geom_dots(side = "topright",
            orientation = "y",
            dotsize = 1) +
  facet_wrap(~factor(statistic,
                     levels=c('hdi90.lower',
                              'mean.est',
                              'hdi90.upper')),
             scales = 'free_x') +
  scale_color_manual(values = c("#b3cf99","#B22222")) +
  scale_fill_manual(values = c("#b3cf99","#B22222")) +
  ggtitle("Detritivory flux - mean difference of plant treatments in communities with soil history\nmain text: -0.44 [-0.87, -0.01]") +
  theme_bw() +
  theme(legend.position = "none")


library(patchwork)

s1/s2/s3

s1
s2
s3
s4
s5
s6

commm = vector(mode = "list", length(longg))
for (i in 1:length(longg)) {
  commm[[i]] = longg[[i]] %>% 
    filter(!(taxon_name %in% c("plants","detritus","microbes"))) %>% 
    group_by(Unit_quarter) %>% 
    summarise(Biomass = sum(Biomass.mg),
              Diversity = exp(vegan::diversity(Biomass.mg, MARGIN = 2)),
              CWM_Bodymass.mg = weighted.mean(Bodymass.mg, density))
}

datalist2 = vector(mode = "list", length(longg))
for (i in 1:length(longg)) {
datalist2[[i]] = right_join(datalist[[i]],commm[[i]])
}


biom.m.2 = brm_multiple(bf(log(Biomass) ~ 
                              treatment
                            + (1|hall.block/unit)),
                         family = student(),
                         chains = 4,
                         cores = 4,
                         iter = 8000,
                         control = list(adapt_delta = 0.99),
                         backend = "cmdstanr",
                         seed = 321,
                         data = datalist2,
                         # we want the results of 250 models
                         combine = FALSE)

saveRDS(biom.m.2, "biomm2sens.rds")
rm(biom.m.2)
gc()


cwm.m.2 = brm_multiple(bf(log(CWM_Bodymass.mg) ~ 
                             treatment
                           + (1|hall.block/unit)),
                        family = student(),
                        chains = 4,
                        cores = 4,
                        iter = 8000,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        seed = 321,
                        data = datalist2,
                        # we want the results of 250 models
                        combine = FALSE)

saveRDS(cwm.m.2, "cwmm2sens.rds")
rm(cwm.m.2)
gc()




biom.m.2 = readRDS("biomm2sens.rds")

biom.soil.contrast = data.frame(mean.est = numeric(length(biom.m.2)),
                                hdi95.lower = numeric(length(biom.m.2)),
                                hdi95.upper = numeric(length(biom.m.2)))

for (i in 1:length(biom.m.2)) {
  biom.soil.contrast[i,] = biom.m.2[[i]] %>%
    ref_grid() %>%
    add_grouping("soil",
                 "treatment",
                 factor(c("0", "0", "1", "1"),
                        levels = c("1", "0"))) %>%
    emmeans(~ soil) %>%
    contrast(method = "pairwise") %>%
    gather_emmeans_draws() %>%
    mean_hdi(.width = 0.95) %>% .[,2:4]
}

# we plot their distributions
s7 = biom.soil.contrast %>%
  pivot_longer(1:3,
               names_to = "statistic",
               values_to = "value") %>%
  ggplot(aes(x = value,
             # conditionally color dots above/below the zero threshold
             fill = (value > 0),
             color = (value > 0))) +
  geom_dots(side = "topright",
            orientation = "y",
            dotsize = 1) +
  facet_wrap(~factor(statistic,
                     levels=c('hdi95.lower',
                              'mean.est',
                              'hdi95.upper')),
             scales = 'free_x') +
  scale_color_manual(values = c("#b3cf99","#B22222")) +
  scale_fill_manual(values = c("#b3cf99","#B22222")) +
  ggtitle("Fauna community biomass - mean difference of soil treatments\nmain text: 0.41 [0.125, 0.705]") +
  theme_bw() +
  theme(legend.position = "none")






cwm.m.2 = readRDS("cwmm2sens.rds")

cwm.soil.contrast = data.frame(mean.est = numeric(length(cwm.m.2)),
                                hdi95.lower = numeric(length(cwm.m.2)),
                                hdi95.upper = numeric(length(cwm.m.2)))

for (i in 1:length(cwm.m.2)) {
  cwm.soil.contrast[i,] = cwm.m.2[[i]] %>%
    ref_grid() %>%
    add_grouping("soil",
                 "treatment",
                 factor(c("0", "0", "1", "1"),
                        levels = c("1", "0"))) %>%
    emmeans(~ soil) %>%
    contrast(method = "pairwise") %>%
    gather_emmeans_draws() %>%
    mean_hdi(.width = 0.95) %>% .[,2:4]
}

# we plot their distributions
s8 = cwm.soil.contrast %>%
  pivot_longer(1:3,
               names_to = "statistic",
               values_to = "value") %>%
  ggplot(aes(x = value,
             # conditionally color dots above/below the zero threshold
             fill = (value > 0),
             color = (value > 0))) +
  geom_dots(side = "topright",
            orientation = "y",
            dotsize = 1) +
  facet_wrap(~factor(statistic,
                     levels=c('hdi95.lower',
                              'mean.est',
                              'hdi95.upper')),
             scales = 'free_x') +
  scale_color_manual(values = c("#b3cf99","#B22222")) +
  scale_fill_manual(values = c("#b3cf99","#B22222")) +
  ggtitle("CWM Body-mass - mean difference of soil treatments\nmain text: -0.64 [-0.94, -0.36]") +
  theme_bw() +
  theme(legend.position = "none")

s7
s8
