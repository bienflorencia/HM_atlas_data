# Flo Grattarola
# 2024-10-04
# Code based on https://doserlab.com/files/spoccupancy-web/articles/spacetimemodelshtml
# For a more basic start, check: https://doserlab.com/files/spoccupancy-web/articles/modelfitting

library(janitor)
library(patchwork)
library(sf)
sf_use_s2(FALSE)
library(spOccupancy)
library(tidyverse)

atlas_birds_CZ <- readRDS('data/Birds_Atlas_Czechia_beast_data.rds')
atlas_birds_CZ %>% 
  filter(start_year==2014 & end_year==2017 & cell_grouping ==1) %>%
  group_by(verbatim_name) %>% count() %>% arrange((n)) %>% print(n=100)

cell1grid <- sf::st_read('data/Birds_Atlas_Czechia_grid.gpkg', layer='cell1grid')
cell1grid_coords <- bind_cols(st_drop_geometry(cell1grid) %>% select(cell_label,area),
                              st_centroid(cell1grid) %>% st_coordinates())
cell1grid_effort <- left_join(cell1grid %>% select(cell_label), 
                              atlas_birds_CZ %>% 
                                filter(start_year==2014 & 
                                         end_year==2017 & 
                                         cell_grouping ==1) %>% 
                                distinct(cell_label, effort, area))

acce_grids <- readRDS('data/acce_grids.rds')

V.det.covs <- st_join(acce_grids, cell1grid_effort, left = T, join = st_equals) %>% 
  select(cell_label=cell_label.x, acce=acc_50k, effort)

################
# to choose species, check the verbatim_name
atlas_birds_CZ %>%
  filter(start_year==2014 & end_year==2017 & cell_grouping ==1) %>%
  group_by(verbatim_name) %>% count() %>% arrange((n)) %>%
  filter(grepl('mega', verbatim_name))

sp_list <- c('anas_strepera', 'ciconia_ciconia',
             'tyto_alba', 'haliaeetus_albicilla',
             'upupa_epops', 'merops_apiaster',
             'jynx_torquilla', 'dendrocopos_minor',
             'larus_ridibundus', 'charadrius_dubius',
             'tetrastes_bonasia', 'streptopelia_turtur',
             'pyrrhula_pyrrhula', 'oriolus_oriolus',
             'luscinia_megarhynchos')

# run for all species
for(sp in sp_list) {
  
  species_data <- atlas_birds_CZ %>% 
    filter(start_year==2014 & end_year==2017 & cell_grouping ==1) %>% 
    select(cell_label, verbatim_name, effort, area) %>% 
    mutate(period=2, presence=1) %>% 
    pivot_wider(names_from = verbatim_name,
                values_from = presence,
                values_fill = c(list(presence=0))) %>% 
    janitor::clean_names() %>% 
    select(presence=all_of(sp),   
           cell_label, effort, area) %>% 
    left_join(., acce_grids) %>% 
    arrange(cell_label)
  
  coords.species <- right_join(cell1grid %>% select(cell_label, area), 
                               species_data) %>%
    st_centroid(cell1grid) %>% st_coordinates()
  
  
  ## data for analysis
  y.species <- data.frame(y=species_data$presence)
  # occ.covs <- data.frame(area = matrix(species_data$area))
  det.covs <- data.frame(effort = matrix(species_data$effort),
                         acce = matrix(species_data$acc_50k))
  coords <- data.frame(X = matrix(coords.species[,1]), 
                       Y = matrix(coords.species[,2]))
  
  data.species <- list(y = y.species,
                       # occ.covs = occ.covs,
                       det.covs = det.covs,
                       coords = coords)
  
  str(data.species)
  
  ## spatial component
  # Pair-wise distances between all sites
  dist.species <- dist(data.species$coords)
  
  # Exponential covariance model
  cov.model <- "exponential"
  
  # Specify list of inits
  species.inits <- list(alpha = 0, 
                        beta = 0, 
                        z = data.species$y,  
                        sigma.sq = 2, 
                        phi = 3 / mean(dist.species), 
                        w = rep(0, nrow(data.species$y)))
  
  batch.length <- 100
  n.batch <- 1000
  n.burn <- 5000
  n.thin <- 20
  n.chains <- 3
  
  species.tuning <- list(phi = 1)
  # accept.rate = 0.43 by default, so we do not specify it.
  
  min.dist <- min(dist.species)
  max.dist <- max(dist.species)
  species.sp.priors <- list(beta.normal = list(mean = 0, var = 2.72),
                            alpha.normal = list(mean = 0, var = 2.72),
                            sigma.sq.ig = c(2, 1),
                            phi.unif = c(3/max.dist, 3/min.dist))
  
  n.omp.threads <- 1
  verbose <- TRUE
  n.report <- 100 # Report progress at every 100th batch.
  
  ######  ######
  # FIT MODEL
  
  out.sp.species <- spPGOcc(occ.formula = ~ 1, 
                            det.formula = ~ scale(effort) + scale(acce),  
                            data = data.species, 
                            # inits = talba.inits, 
                            n.batch = n.batch, 
                            batch.length = batch.length, 
                            priors = species.sp.priors, 
                            cov.model = cov.model, 
                            NNGP = TRUE, 
                            n.neighbors = 5,
                            tuning = species.tuning, 
                            n.report = n.report, 
                            n.burn = n.burn, 
                            n.thin = n.thin, 
                            n.chains = n.chains)
  
  summary(out.sp.species)
  # plot(out.sp.species, 'beta', density = FALSE)
  # plot(out.sp.species, 'alpha', density = FALSE)
  # plot(out.sp.species, 'theta', density = FALSE)
  
  ######  ######
  # Predictive check
  
  ppc.species.out <- ppcOcc(out.sp.species, fit.stat = 'freeman-tukey', group = 1)
  summary(ppc.species.out)
  
  # ppc.species.df <- data.frame(fit = ppc.species.out$fit.y,
  #                              fit.rep = ppc.species.out$fit.y.rep,
  #                              color = 'lightskyblue1')
  # ppc.species.df$color[ppc.species.df$fit.rep > ppc.species.df$fit] <- 'lightsalmon'
  # plot(ppc.species.df$fit, ppc.species.df$fit.rep, bg = ppc.species.df$color, pch = 21,
  #      ylab = 'Fit', xlab = 'True')
  # lines(ppc.species.df$fit, ppc.species.df$fit, col = 'black')

  ######  ######
  # Predictions
  
  # area.pred <- cell1grid_coords$area
  # X.0 <- cbind(1, area.pred) # with area
  # X.0 <- matrix(1, nrow(species_data), 1)
  X.0 <- matrix(1, nrow(cell1grid_coords), 1)
  
  coords.0 <- as.matrix(cell1grid_coords[, c('X', 'Y')]) # some unsampled grids
  # coords.0 <- coords.species
  
  acce.pred <- (V.det.covs$acce - mean(V.det.covs$acce)) / sd(V.det.covs$acce)
  effort.pred <- (V.det.covs$effort - mean(V.det.covs$effort, na.rm=T)) / sd(V.det.covs$effort, na.rm=T)
  V.0 <- cbind(1,effort.pred, acce.pred)
  
  # prediction: occupancy
  out.sp.species.pred.occ <- predict(out.sp.species, X.0, coords.0, 
                                     verbose = FALSE, 
                                     type='occupancy')
  
  # prediction: detection
  out.sp.species.pred.det <- predict(out.sp.species, V.0, coords.0, 
                                     verbose = FALSE,
                                     type = 'detection')
  
  # Produce a species distribution map
  dat.sp <- data.frame(x = cell1grid_coords$X,
                       y = cell1grid_coords$Y, 
                       mean.psi = apply(out.sp.species.pred.occ$psi.0.samples, 2, mean),
                       sd.psi = apply(out.sp.species.pred.occ$psi.0.samples, 2, sd),
                       mean.w = apply(out.sp.species.pred.occ$w.0.samples, 2, mean),
                       sd.w = apply(out.sp.species.pred.occ$w.0.samples, 2, sd),
                       mean.p = apply(out.sp.species.pred.det$p.0.samples, 2, mean, na.rm = TRUE),
                       sd.p = apply(out.sp.species.pred.det$p.0.samples, 2, sd, na.rm = TRUE))
  
  dat.sp.sf <- st_as_sf(dat.sp, coords = c('x', 'y'), crs=4326)
  
  sp_name <- str_to_sentence(str_replace(sp, '_', ' '))
  raw.occ <- round(sum(species_data$presence) / nrow(species_data),2)
  alpha.det <- round(plogis(mean(out.sp.species$alpha.samples)),3)
  beta.occ <- round(plogis(mean(out.sp.species$beta.samples)),3)
  
  plot_psi <- ggplot() + 
    geom_sf(data = st_join(cell1grid,dat.sp.sf), aes(fill = mean.psi)) +
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(title= str_glue('{sp_name} \u03A8 (occupancy / \u03B2 = {beta.occ})')) +
    theme_bw()
  
  plot_p <- ggplot() + 
    geom_sf(data = st_join(cell1grid,dat.sp.sf), aes(fill = mean.p)) +
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(title= str_glue('{sp_name} p (detection / \u03B1 = {alpha.det})' )) +
    theme_bw()
  
  plot_y <- ggplot() + 
    geom_sf(data = left_join(cell1grid_effort, species_data), 
            aes(fill = presence)) +
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(title= str_glue('{sp_name} presence/absence (occ % = {raw.occ})' )) +
    theme_bw()
  
  plot_effort <- ggplot() + 
    geom_sf(data = cell1grid_effort, aes(fill = log(effort))) +
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(title= 'Sampling effort') +
    theme_bw()
  
  plot_acce <- ggplot() + 
    geom_sf(data = acce_grids, aes(fill = log(acc_50k))) +
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(title= 'Accessibility to cities') +
    theme_bw()
  
  (plot_y | plot_effort | plot_acce)/ (plot_psi | plot_p) 
  ggsave(filename = str_glue('docs/figs/{sp}_occ.png'), width = 18, height = 10)
  
}
