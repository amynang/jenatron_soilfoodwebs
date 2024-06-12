library(tidyverse)
library(brms)
library(emmeans)
library(modelr)
library(tidybayes)

# flux aggregates per mesocosm
fluxxes = read.csv("datafiles/flux_aggregates.csv")
# community data
long = read_csv("datafiles/jenatron_longformat_community.csv")

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

fluxdata = d.1.1 %>% 
  select(hall.block, unit, Unit_quarter, treatment, plant.diversity) %>% 
  right_join(., fluxxes, by = join_by(Unit_quarter)) %>% 
  mutate(hall.block = as.factor(hall.block),
         unit= as.factor(unit),
         plant.diversity.sc = scale(plant.diversity))

comm = long %>% 
  filter(!(taxon_name %in% c("plants","detritus","microbes"))) %>% 
  group_by(Unit_quarter) %>% 
  summarise(Biomass = sum(Biomass.mg),
            Diversity = exp(vegan::diversity(Biomass.mg, MARGIN = 2)),
            CWM_Bodymass.mg = weighted.mean(Bodymass.mg, density))

commm = right_join(fluxdata,comm)


####### Total flux ######

total.m.1 = brm(bf(log(Total.flux) ~ 
                     treatment*plant.diversity.sc 
                   + (1|hall.block/unit)),
                family = student(),
                chains = 4,
                cores = 4,
                iter = 8000,
                control = list(adapt_delta = 0.99),
                backend = "cmdstanr",
                seed = 321,
                data = fluxdata,
                file = "totalm1")

total.m.2 = brm(bf(log(Total.flux) ~ 
                     treatment
                   + (1|hall.block/unit)),
                family = student(),
                chains = 4,
                cores = 4,
                iter = 8000,
                control = list(adapt_delta = 0.999),
                backend = "cmdstanr",
                seed = 321,
                data = fluxdata,
                file = "totalm2")

summary(total.m.1);pp_check(total.m.1, ndraws = 100)
summary(total.m.2);pp_check(total.m.2, ndraws = 100)

# slopes for each treatment
tot.emt = emtrends(total.m.1, "treatment", var = "plant.diversity.sc")
summary(tot.emt, point.est = mean, level = .95)
# pairwise comparisons
tot.pairs = pairs(tot.emt)
summary(tot.pairs, point.est = mean, level = .95)

(loo1 = loo(total.m.1))   
(loo2 = loo(total.m.2))
# compare both models
loo_compare(loo1, loo2)

ci <- emmeans(total.m.2, specs = pairwise ~ treatment)
summary(ci, point.est = mean)
plot(conditional_effects(total.m.2, effects = "treatment"),
     points = TRUE)

c_eff <- conditional_effects(total.m.1)
fig1a <- plot(c_eff, plot = FALSE, points = FALSE)[[3]] + 
  aes(linetype=treatment) +
  geom_point(
    aes(x = plant.diversity.sc, y = log(Total.flux),
        colour = treatment, shape = treatment), 
    data = fluxdata, 
    position = position_jitter(.05),
    size = 2,
    # This tells it to ignore the ymin and ymax settings used elsewhere
    inherit.aes = FALSE
  ) +
  scale_x_continuous(breaks = c(-1.0386217, -0.2832605, 0.4721008, 2.7381844),
                     labels = c('1', '2', '3', '6')) +
  scale_fill_manual(values = c('#cb992e40','#cb992e40', '#41291440', '#41291440'), 
                    name = "history treatment",
                    labels = c("soil (–), plant (–)",
                               "soil (–), plant (+)",
                               "soil (+), plant (–)", 
                               "soil (+), plant (+)")) +
  scale_linetype_manual(values = c("dashed","solid","dashed","solid"), 
                        name = "history treatment",
                        labels = c("soil (–), plant (–)",
                                   "soil (–), plant (+)",
                                   "soil (+), plant (–)", 
                                   "soil (+), plant (+)")) +
  scale_shape_manual(values = c(19, 15, 19, 15), 
                     name = "history treatment",
                     labels = c("soil (–), plant (–)",
                                "soil (–), plant (+)",
                                "soil (+), plant (–)", 
                                "soil (+), plant (+)")) +
  labs(
    y = "Community level energy flux log(J/h)",
    x = "Plant richness") +
  scale_color_manual(values = c('#cb992e','#cb992e', '#412914', '#412914'), 
                     name = "history treatment",
                     labels = c("soil (–), plant (–)",
                                "soil (–), plant (+)",
                                "soil (+), plant (–)", 
                                "soil (+), plant (+)")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))
fig1a

fig1b = fluxdata %>%
  select(hall.block, 
         unit,
         treatment, 
         Total.flux) %>% 
  data_grid(hall.block, unit, treatment) %>%
  add_epred_draws(total.m.2, dpar = c("mu", 
                                      "sigma"),
                  re_formula = NA) %>%
  ggplot(aes(x = treatment)) +
  stat_pointinterval(aes(y = .epred),
                     point_interval = "mean_qi",
                     .width = c(.95),
                     position = position_nudge(x = -.15)) +
  geom_point(aes(y = log(Total.flux),
                 fill = treatment,
                 shape = treatment), 
             data = fluxdata, 
             alpha = .7,
             size = 2, 
             position = position_jitter(.05)) +
  scale_shape_manual(values=c(21, 22, 21, 22))+
  scale_fill_manual(values=c('#cb992e','#cb992e', '#412914', '#412914')) +
  scale_x_discrete(labels = c('--', '-+', '+-', '++')) +
  annotate('text', x=c(.85,1.85,2.85,3.85), y=5.5, label=c('a','a','b','b'),
           color = "darkred", fontface = 'bold') +
  labs(
    y = "Community level energy flux log(J/h)",
    x = "History treatment") +
  theme_classic() +
  theme(axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))



################################## Pred flux ###################################

pred.m.1 = brm(bf(log(Pred.flux) ~ 
                    treatment*plant.diversity.sc 
                  + (1|hall.block/unit)),
               family = student(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = fluxdata,
               file = "predm1")

pred.m.2 = brm(bf(log(Pred.flux) ~ 
                    treatment
                  + (1|hall.block/unit)),
               family = student(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = fluxdata,
               file = "predm2")

summary(pred.m.1);pp_check(pred.m.1, ndraws = 100)
summary(pred.m.2);pp_check(pred.m.2, ndraws = 100)

# slopes for each treatment
pred.emt = emtrends(pred.m.1, "treatment", var = "plant.diversity.sc")
summary(pred.emt, point.est = mean, level = .95)
# pairwise comparisons
pred.pairs = pairs(pred.emt)
summary(pred.pairs, point.est = mean, level = .95)

(loo1 = loo(pred.m.1))   
(loo2 = loo(pred.m.2))
# compare both models
loo_compare(loo1, loo2)

summary(pred.m.2)
ci <- emmeans(pred.m.2, specs = pairwise ~ treatment)
summary(ci, point.est = mean)


fig1c = fluxdata %>%
  select(hall.block, 
         unit,
         treatment, 
         Pred.flux) %>% 
  data_grid(hall.block, unit, treatment) %>%
  add_epred_draws(pred.m.2, dpar = c("mu", 
                                     "sigma"),
                  re_formula = NA) %>%
  ggplot(aes(x = treatment)) +
  stat_pointinterval(aes(y = .epred),
                     point_interval = "mean_qi",
                     .width = c(.95),
                     position = position_nudge(x = -.15)) +
  geom_point(aes(y = log(Pred.flux),
                 fill = treatment,
                 shape = treatment), 
             data = fluxdata, 
             alpha = .7,
             size = 2, 
             position = position_jitter(.05)) +
  scale_shape_manual(values=c(21, 22, 21, 22))+
  scale_fill_manual(values=c('#cb992e','#cb992e', '#412914', '#412914')) +
  scale_x_discrete(labels = c('--', '-+', '+-', '++')) +
  annotate('text', x=c(.85,1.85,2.85,3.85), y=4.5, label=c('a','a','b','b'),
           color = "darkred", fontface = 'bold') +
  labs(
    y = "Predatory flux log(J/h)",
    x = "History treatment") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))


################################## Herb flux ###################################

herb.m.1 = brm(bf(log(Herb.flux) ~ 
                    treatment*plant.diversity.sc 
                  + (1|hall.block/unit)),
               family = student(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = fluxdata,
               file = "herbm1")

herb.m.2 = brm(bf(log(Herb.flux) ~ 
                    treatment 
                  + (1|hall.block/unit)),
               family = student(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = fluxdata,
               file = "herbm2")

summary(herb.m.1);pp_check(herb.m.1, ndraws = 100)
summary(herb.m.2);pp_check(herb.m.2, ndraws = 100)

# slopes for each treatment
herb.emt = emtrends(herb.m.1, "treatment", var = "plant.diversity.sc")
summary(herb.emt, point.est = mean, level = .95)
# pairwise comparisons
herb.pairs = pairs(herb.emt)
summary(herb.pairs, point.est = mean, level = .95)

(loo1 = loo(herb.m.1))   
(loo2 = loo(herb.m.2))
# compare both models
loo_compare(loo1, loo2)

summary(herb.m.2)
ci <- emmeans(herb.m.2, specs = pairwise ~ treatment)
summary(ci, point.est = mean)



fig1e = fluxdata %>%
  select(hall.block, 
         treatment, 
         unit,
         Herb.flux) %>% 
  data_grid(hall.block, unit, treatment) %>%
  add_epred_draws(herb.m.2, dpar = c("mu", 
                                     "sigma"),
                  re_formula = NA) %>%
  ggplot(aes(x = treatment)) +
  stat_pointinterval(aes(y = .epred),
                     point_interval = "mean_qi",
                     .width = c(.95),
                     position = position_nudge(x = -.15)) +
  geom_point(aes(y = log(Herb.flux),
                 fill = treatment,
                 shape = treatment), 
             data = fluxdata, 
             alpha = .7,
             size = 2, 
             position = position_jitter(.05)) +
  scale_shape_manual(values=c(21, 22, 21, 22))+
  scale_fill_manual(values=c('#cb992e','#cb992e', '#412914', '#412914')) +
  scale_x_discrete(labels = c('--', '-+', '+-', '++')) +
  annotate('text', x=c(.85,1.85,2.85,3.85), y=6.3, label=c('a','a','b','b'),
           color = "darkred", fontface = 'bold') +
  labs(
    y = "Herbivory flux log(J/h)",
    x = "History treatment") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))


################################## Micr flux ###################################

micr.m.1 = brm(bf(log(Micr.flux) ~ 
                    treatment*plant.diversity.sc 
                  + (1|hall.block/unit)),
               family = gaussian(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = fluxdata,
               file = "micrm1")

micr.m.2 = brm(bf(log(Micr.flux) ~ 
                    treatment 
                  + (1|hall.block/unit)),
               family = gaussian(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = fluxdata,
               file = "micrm2")

summary(micr.m.1);pp_check(micr.m.1, ndraws = 100)
summary(micr.m.2);pp_check(micr.m.2, ndraws = 100)

# slopes for each treatment
micr.emt = emtrends(micr.m.1, "treatment", var = "plant.diversity.sc")
summary(micr.emt, point.est = mean, level = .95)
# pairwise comparisons
micr.pairs = pairs(micr.emt)
summary(micr.pairs, point.est = mean, level = .95)

(loo1 = loo(micr.m.1))   
(loo2 = loo(micr.m.2))
# compare both models
loo_compare(loo1, loo2)

summary(micr.m.2)
ci <- emmeans(micr.m.2, specs = pairwise ~ treatment)
summary(ci, point.est = mean)

bayestestR::p_direction(ci)


fig1d = fluxdata %>%
  select(hall.block, 
         unit,
         treatment, 
         Micr.flux) %>% 
  data_grid(hall.block, unit, treatment) %>%
  add_epred_draws(micr.m.2, dpar = c("mu", 
                                     "sigma"),
                  re_formula = NA) %>%
  ggplot(aes(x = treatment)) +
  stat_pointinterval(aes(y = .epred),
                     point_interval = "mean_qi",
                     .width = c(.95),
                     position = position_nudge(x = -.15)) +
  geom_point(aes(y = log(Micr.flux),
                 fill = treatment,
                 shape = treatment), 
             data = fluxdata, 
             alpha = .7,
             size = 2, 
             position = position_jitter(.05)) +
  scale_shape_manual(values=c(21, 22, 21, 22))+
  scale_fill_manual(values=c('#cb992e','#cb992e', '#412914', '#412914')) +
  scale_x_discrete(labels = c('--', '-+', '+-', '++')) +
  annotate('text', x=c(.85,1.85,2.85,3.85), y=4, label=c('a','a','b','b'),
           color = "darkred", fontface = 'bold') +
  labs(
    y = "Microbivory flux log(J/h)",
    x = "History treatment") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))



################################## Detr flux ###################################

detr.m.1 = brm(bf(log(Detr.flux) ~ 
                    treatment*plant.diversity.sc 
                  + (1|hall.block/unit)),
               family = student(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = fluxdata,
               file = "detrm1")

detr.m.2 = brm(bf(log(Detr.flux) ~ 
                    treatment 
                  + (1|hall.block/unit)),
               family = student(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = fluxdata,
               file = "detrm2")

summary(detr.m.1);pp_check(detr.m.1, ndraws = 100)
summary(detr.m.2);pp_check(detr.m.2, ndraws = 100)

# slopes for each treatment
detr.emt = emtrends(detr.m.1, "treatment", var = "plant.diversity.sc")
summary(detr.emt, point.est = mean, level = .95)
# pairwise comparisons
detr.pairs = pairs(detr.emt)
summary(detr.pairs, point.est = mean, level = .95)

(loo1 = loo(detr.m.1))   
(loo2 = loo(detr.m.2))
# compare both models
loo_compare(loo1, loo2)

summary(detr.m.2)
ci <- emmeans(detr.m.2, specs = pairwise ~ treatment)
summary(ci, point.est = mean)
plot(conditional_effects(detr.m.2, effects = "treatment"),
     points = TRUE)

fig1f = fluxdata %>%
  select(hall.block, 
         unit,
         treatment, 
         Detr.flux) %>% 
  data_grid(hall.block, unit, treatment) %>%
  add_epred_draws(detr.m.2, dpar = c("mu", 
                                     "sigma"),
                  re_formula = NA) %>%
  ggplot(aes(x = treatment)) +
  stat_pointinterval(aes(y = .epred),
                     point_interval = "mean_qi",
                     .width = c(.95),
                     position = position_nudge(x = -.15)) +
  geom_point(aes(y = log(Detr.flux),
                 fill = treatment,
                 shape = treatment), 
             data = fluxdata, 
             alpha = .7,
             size = 2, 
             position = position_jitter(.05)) +
  scale_shape_manual(values=c(21, 22, 21, 22))+
  scale_fill_manual(values=c('#cb992e','#cb992e', '#412914', '#412914')) +
  scale_x_discrete(labels = c('--', '-+', '+-', '++')) +
  annotate('text', x=c(.85,1.85,2.85,3.85), y=3.5, label=c('a','a','b','ab'),
           color = "darkred", fontface = 'bold') +
  labs(
    y = "Detritivory flux log(J/h)",
    x = "History treatment") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))


library(patchwork)

(fig1a + fig1c + fig1e)/(fig1b + fig1d + fig1f) + 
  plot_annotation(tag_levels = list(c('a','c','e',
                                      'b','d','f'))) +
  plot_layout(widths = c(2  , 1, 1,
                         1.5, 1, 1))

(fig1a / fig1b) + 
  plot_annotation(tag_levels = list(c('a','b')))
(fig1c + fig1e)/(fig1d + fig1f) + 
  plot_annotation(tag_levels = list(c('c','e',
                                      'd','f')))

ggsave("Figure1ab.png",
       bg = "white",
       scale = 1.6,
       width = 70,
       height = 90,
       units = "mm",
       dpi = 600)
ggsave("Figure1cdef.png",
       bg = "white",
       scale = 1.6,
       width = 110,
       height = 90,
       units = "mm",
       dpi = 600)
################################## Biomass #####################################

biom.m.1 = brm(bf(log(Biomass) ~ 
                    treatment*plant.diversity.sc 
                  + (1|hall.block/unit)),
               family = student(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = commm,
               file = "biomm1")

biom.m.2 = brm(bf(log(Biomass) ~ 
                    treatment
                  + (1|hall.block/unit)),
               family = student(),
               chains = 4,
               cores = 4,
               iter = 8000,
               control = list(adapt_delta = 0.99),
               backend = "cmdstanr",
               seed = 321,
               data = commm,
               file = "biomm2")

summary(biom.m.1);pp_check(biom.m.1, ndraws = 100)
summary(biom.m.2);pp_check(biom.m.2, ndraws = 100)

# slopes for each treatment
biom.emt = emtrends(biom.m.1, "treatment", var = "plant.diversity.sc")
summary(biom.emt, point.est = mean, level = .95)
# pairwise comparisons
biom.pairs = pairs(biom.emt)
summary(biom.pairs, point.est = mean, level = .95)

(loo1 = loo(biom.m.1))   
(loo2 = loo(biom.m.2))
# compare both models
loo_compare(loo1, loo2)

summary(biom.m.2)
ci <- emmeans(biom.m.2, specs = pairwise ~ treatment)
summary(ci, point.est = mean)
summary(ci, point.est = mean, level = .9)

fig2c = commm %>%
  select(hall.block,
         unit,
         treatment, 
         Biomass) %>% 
  data_grid(hall.block, unit, treatment) %>%
  add_epred_draws(biom.m.2, dpar = c("mu", 
                                     "sigma"),
                  re_formula = NA) %>%
  ggplot(aes(x = treatment)) +
  stat_pointinterval(aes(y = .epred),
                     point_interval = "mean_qi",
                     .width = c(.95),
                     position = position_nudge(x = -.15)) +
  geom_point(aes(y = log(Biomass),
                 fill = treatment,
                 shape = treatment), 
             data = commm, 
             alpha = .7,
             size = 2, 
             position = position_jitter(.05)) +
  scale_shape_manual(values=c(21, 22, 21, 22))+
  scale_fill_manual(values=c('#cb992e','#cb992e', '#412914', '#412914')) +
  scale_x_discrete(labels = c('--', '-+', '+-', '++')) +
  annotate('text', x=c(.85,1.85,2.85,3.85), y=9, label=c('a','a','b','b'),
           color = "darkred", fontface = 'bold') +
  labs(
    y = "Community Biomass log(mg)",
    x = "History treatment") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))

################################# CWM Bodymass #################################

cwm.m.1 = brm(bf(log(CWM_Bodymass.mg) ~ 
                   treatment*plant.diversity.sc 
                 + (1|hall.block/unit)),
              family = student(),
              chains = 4,
              cores = 4,
              iter = 8000,
              control = list(adapt_delta = 0.99),
              backend = "cmdstanr",
              seed = 321,
              data = commm,
              file = "cwmm1")

cwm.m.2 = brm(bf(log(CWM_Bodymass.mg) ~ 
                   treatment
                 + (1|hall.block/unit)),
              family = student(),
              chains = 4,
              cores = 4,
              iter = 8000,
              control = list(adapt_delta = 0.999),
              backend = "cmdstanr",
              seed = 321,
              data = commm,
              file = "cwmm2")

summary(cwm.m.1);pp_check(cwm.m.1, ndraws = 100)
summary(cwm.m.2);pp_check(cwm.m.2, ndraws = 100)

# slopes for each treatment
cwm.emt = emtrends(cwm.m.1, "treatment", var = "plant.diversity.sc")
summary(cwm.emt, point.est = mean, level = .95)
# pairwise comparisons
cwm.pairs = pairs(cwm.emt)
summary(cwm.pairs, point.est = mean, level = .95)

(loo1 = loo(cwm.m.1))   
(loo2 = loo(cwm.m.2))
# compare both models
loo_compare(loo1, loo2)

summary(cwm.m.2)
ci <- emmeans(cwm.m.2, specs = pairwise ~ treatment)
summary(ci, point.est = mean)

fig2b = commm %>%
  select(hall.block, 
         unit,
         treatment, 
         CWM_Bodymass.mg) %>% 
  data_grid(hall.block, unit, treatment) %>%
  add_epred_draws(cwm.m.2, dpar = c("mu", 
                                    "sigma"),
                  re_formula = NA) %>%
  ggplot(aes(x = treatment)) +
  stat_pointinterval(aes(y = .epred),
                     point_interval = "mean_qi",
                     .width = c(.95),
                     position = position_nudge(x = -.15)) +
  geom_point(aes(y = log(CWM_Bodymass.mg),
                 fill = treatment,
                 shape = treatment), 
             data = commm, 
             alpha = .7,
             size = 2, 
             position = position_jitter(.05)) +
  scale_shape_manual(values=c(21, 22, 21, 22))+
  scale_fill_manual(values=c('#cb992e','#cb992e', '#412914', '#412914')) +
  scale_x_discrete(labels = c('--', '-+', '+-', '++')) +
  annotate('text', x=c(.85,1.85,2.85,3.85), y=-3.8, label=c('a','a','b','b'),
           color = "darkred", fontface = 'bold') +
  labs(
    y = "CWM Bodymass log(mg)",
    x = "History treatment") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))


################################# Shannon Entropy ##############################

div.m.1 = brm(bf(Diversity ~ 
                   treatment*plant.diversity.sc 
                 + (1|hall.block/unit)),
              family = gaussian(),
              chains = 4,
              cores = 4,
              iter = 8000,
              control = list(adapt_delta = 0.99),
              backend = "cmdstanr",
              seed = 321,
              data = commm,
              file = "divm1")

div.m.2 = brm(bf(Diversity ~ 
                   treatment 
                 + (1|hall.block/unit)),
              family = gaussian(),
              chains = 4,
              cores = 4,
              iter = 8000,
              control = list(adapt_delta = 0.99),
              backend = "cmdstanr",
              seed = 321,
              data = commm,
              file = "divm2")

summary(div.m.1);pp_check(div.m.1, ndraws = 100)
summary(div.m.2);pp_check(div.m.2, ndraws = 100)

# slopes for each treatment
div.emt = emtrends(div.m.1, "treatment", var = "plant.diversity.sc")
summary(div.emt, point.est = mean, level = .95)
# pairwise comparisons
div.pairs = pairs(div.emt)
summary(div.pairs, point.est = mean, level = .95)

(loo1 = loo(div.m.1))   
(loo2 = loo(div.m.2))
# compare both models
loo_compare(loo1, loo2)

summary(div.m.2)
ci <- emmeans(div.m.2, specs = pairwise ~ treatment)
summary(ci, point.est = mean)

fig2a = commm %>%
  select(hall.block, 
         treatment, 
         unit,
         Diversity) %>% 
  data_grid(hall.block, unit, treatment) %>%
  add_epred_draws(div.m.2, dpar = c("mu", 
                                    "sigma"),
                  re_formula = NA) %>%
  ggplot(aes(x = treatment)) +
  stat_pointinterval(aes(y = .epred),
                     point_interval = "mean_qi",
                     .width = c(.95),
                     position = position_nudge(x = -.15)) +
  geom_point(aes(y = (Diversity),
                 fill = treatment,
                 shape = treatment), 
             data = commm, 
             alpha = .7,
             size = 2, 
             position = position_jitter(.05)) +
  scale_shape_manual(values=c(21, 22, 21, 22))+
  scale_fill_manual(values=c('#cb992e','#cb992e', '#412914', '#412914')) +
  scale_x_discrete(labels = c('--', '-+', '+-', '++')) +
  labs(
    y = "Diversity",
    x = "History treatment") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.text.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.0)))


(fig2a + fig2b + fig2c) + 
  plot_annotation(tag_levels = list(c('a','b','c')))

ggsave("Figure2abc.png",
       bg = "white",
       scale = 1.6,
       width = 180,
       height = 60,
       units = "mm",
       dpi = 600)
################################ Something #####################################

m.tot.biom.cwm = brm(bf(log(Total.flux) ~ 
                          scale(logBiomass)*scale(logCWM)
                        + (1|hall.block/unit)),
                     family = student(),
                     chains = 4,
                     cores = 4,
                     iter = 8000,
                     control = list(adapt_delta = 0.999),
                     backend = "cmdstanr",
                     seed = 321,
                     data = commm %>% 
                       mutate(logBiomass = log(Biomass),
                              logCWM = log(CWM_Bodymass.mg)),
                     file = "fluxbiomcwm")
pp_check(m.tot.biom.cwm,ndraws = 100)
summary(m.tot.biom.cwm)
plot(conditional_effects(m.tot.biom.cwm),
     points = TRUE)

c_eff <- conditional_effects(m.tot.biom.cwm)
fig3a <- plot(c_eff, plot = FALSE, points = FALSE)[[1]] + 
  #aes(colour = "")
  geom_line(color = "black",size =1.1) +
  geom_point(
    aes(x = logBiomass, y = log(Total.flux),
        colour = treatment, shape = treatment), 
    data = commm %>% 
      mutate(logBiomass = log(Biomass),
             logCWM = log(CWM_Bodymass.mg)), 
    size = 2,
    # This tells it to ignore the ymin and ymax settings used elsewhere
    inherit.aes = FALSE
  ) +
  scale_shape_manual(values = c(19, 15, 19, 15), 
                     name = "history treatment",
                     labels = c("soil (–), plant (–)",
                                "soil (–), plant (+)",
                                "soil (+), plant (–)", 
                                "soil (+), plant (+)")) +
  labs(
    y = "Community level energy flux log(J/h)",
    x = "Plant richness") +
  scale_color_manual(values = c('#cb992e','#cb992e', '#412914', '#412914'), 
                     name = "history treatment",
                     labels = c("soil (–), plant (–)",
                                "soil (–), plant (+)",
                                "soil (+), plant (–)", 
                                "soil (+), plant (+)")) +
  labs(
    y = "Community level energy flux log(J/h)",
    x = "Biomass log(mg)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.title.y = element_text(size = rel(1.0)))
fig3a

fig3b <- plot(c_eff, plot = FALSE, points = FALSE)[[2]] + 
  #aes(colour = "")
  geom_line(color = "black",size =1.1) +
  geom_point(
    aes(x = logCWM, y = log(Total.flux),
        colour = treatment, shape = treatment), 
    data = commm %>% 
      mutate(logBiomass = log(Biomass),
             logCWM = log(CWM_Bodymass.mg)), 
    size = 2,
    # This tells it to ignore the ymin and ymax settings used elsewhere
    inherit.aes = FALSE
  ) +
  scale_shape_manual(values = c(19, 15, 19, 15), 
                     name = "history treatment",
                     labels = c("soil (–), plant (–)",
                                "soil (–), plant (+)",
                                "soil (+), plant (–)", 
                                "soil (+), plant (+)")) +
  labs(
    y = "Community level energy flux log(J/h)",
    x = "Plant richness") +
  scale_color_manual(values = c('#cb992e','#cb992e', '#412914', '#412914'), 
                     name = "history treatment",
                     labels = c("soil (–), plant (–)",
                                "soil (–), plant (+)",
                                "soil (+), plant (–)", 
                                "soil (+), plant (+)")) +
  labs(
    y = "Community level energy flux log(J/h)",
    x = "CWM Body-mass log(mg)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = rel(1.0)),
        axis.title.y = element_text(size = rel(1.0)))
fig3b

(fig3b + fig3a) + 
  plot_annotation(tag_levels = list(c('a','b')))

ggsave("Figure3ab.png",
       bg = "white",
       scale = 1.6,
       width = 180*(2/3),
       height = 60,
       units = "mm",
       dpi = 600)


detr.m.2 %>%
  
  emmeans(~ treatment) %>%
  gather_emmeans_draws() %>%
  median_qi()

soil = factor(c("0", "0", "1", "1"), levels = c("0", "1"))
detr.m.2 %>%
  ref_grid() %>% 
  add_grouping("soil",
               "treatment",
               factor(c("0", "0", "1", "1"), 
                      levels = c("0", "1"))) %>% 
  emmeans(~ soil) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  mean_qi()


pred.m.2 %>%
  ref_grid() %>% 
  add_grouping("soil",
               "treatment",
               factor(c("0", "0", "1", "1"), 
                      levels = c("0", "1"))) %>% 
  emmeans(~ soil) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  mean_qi()

total.m.2 %>%
  ref_grid() %>% 
  add_grouping("soil",
               "treatment",
               factor(c("0", "0", "1", "1"), 
                      levels = c("0", "1"))) %>% 
  emmeans(~ soil) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  mean_qi()

total.m.2 %>%
  ref_grid() %>% 
  add_grouping("soil",
               "treatment",
               factor(c("0", "0", "1", "1"), 
                      levels = c("1", "0"))) %>% 
  emmeans(~ soil) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws() %>%
  mean_qi(.width = 0.95)
