################################################
######## R Code for Johnston et al. 2025, ######
###### Diverse historical fire disturbance #####
### and successional dynamics in Douglas-fir ###
## forests of the western Oregon Cascades, USA ##
#################################################

#All code created by James Johnston.  Contact jamesjoh@uoregon.edu with comments or suggestions for improvement (there is lots of room for improvement and I probably will update and streamline some of this code, so drop me a line for those updates).
#Some of this code is computationally intensive.  Running everything from beginning to end may take an hour—longer if you specify the same number of simulations as I used in the paper.
#This code should run end to end just fine. But it is probably best practice to restart your session after every section to ensure there's no unforeseen conflicts (there is at least one package conflict with the gstat and automap libraries which is noted below).  There is redundancy in the code within sections, i.e., the same files are opened multiple times and sometimes manipulated in the same way to facilitate running individual sections independently without creating a hopelessly cluttered workspace.   
#The graphical outputs are optimized for my MacBook Pro and some of the graphical workflow is a little idiosyncratic.  You will want to customize the graphical output for a good fit with your workflow.
#This code consists of the following sections:

#Create maps
#Evaluate environmental space
#Create fire histories
#Summarize data, calculate fire history metrics, and bootstrap MFRIs
#Create fire occurrence file
#Create ignition density file
#spaMM model for fire occurrence
#Evaluate watershed cohorts
#Evaluate large fires

###############################
######### Create maps #########
###############################

#Here we're going to make the maps in Figure 1.  This code will download rasters, plot rasters, and combine maps.  It may take 5-15 minutes to run.

#Clear workspace and set directory
rm(list=ls())
setwd("/YOURPATH")

#Load libraries for this section
library(tidyverse)
library(geosphere)
library(elevatr)
library(sf)
library(terra)
library(tidyterra)
library(ggspatial)
library(ggnewscale)
library(ggpubr)

#Read site centroids and environmental variables files
site_coords <- read.csv("site_centroids.csv")
env_vars <- read.csv("env_data.csv") %>%
  select(site, huc8_name) %>%
  group_by(site) %>%
  count(huc8_name) %>%
  slice(which.max(n))
site_coords <- site_coords %>% left_join(env_vars)

#Define path to GeoPackages
gpkg_study <- "wfi_layers.gpkg"
gpkg_region <- "wfi_layers_region.gpkg"

#Define a study area extent and download a raster
wfi_study_extent <- data.frame(x = c(-123, -121), y = c(46, 43))
prj_dd <- "EPSG:4326"
wfi_dem <- get_elev_raster(wfi_study_extent, prj = prj_dd, z = 10)
wfi_rast <- rast(wfi_dem)
wfi_rast_df <- as.data.frame(wfi_rast, xy = TRUE)
names(wfi_rast_df)[3] <- "elev_m"

#Calculate slope, aspect, and hillshade, and convert to dataframe.
slope <- terrain(wfi_rast, "slope", unit = "radians")
aspect <- terrain(wfi_rast, "aspect", unit = "radians")
wfi_hill <- shade(slope, aspect, angle = 45, direction = 300, normalize = TRUE)
wfi_hill_df <- as.data.frame(wfi_hill, xy = TRUE)

#Read spatial layers for the study area map
wfi_layers <- "wfi_layers.gpkg"
study_area_nfs <- st_read(wfi_layers, layer = "national_forests", quiet = TRUE)
ugbs <- st_read(wfi_layers, layer = "urban_growth_boundaries", quiet = TRUE)
rivers <- st_read(wfi_layers, layer = "large_rivers", quiet = TRUE)
lakes_res <- st_read(wfi_layers, layer = "lakes_reservoirs", quiet = TRUE)
large_fires <- st_read(wfi_layers, layer = "large_fires", quiet = TRUE)

#Define a region extent and download a raster
region_extent <- data.frame(x = c(-124.53, -120.75), y = c(48.6, 40.5))
region_dem <- get_elev_raster(region_extent, prj = prj_dd, z = 10)
region_rast <- rast(region_dem)
region_rast_df <- as.data.frame(region_rast, xy = TRUE)
names(region_rast_df)[3] <- "elev_m"

#Calculate slope, aspect, and hillshade, and convert to dataframe 
slope <- terrain(region_rast, "slope", unit = "radians")
aspect <- terrain(region_rast, "aspect", unit = "radians")
region_hill <- shade(slope, aspect, angle = 45, direction = 300, normalize = TRUE)
region_hill_df <- as.data.frame(region_hill, xy = TRUE)

#Read spatial layers for the region map
wfi_layers_region <- "wfi_layers_region.gpkg"
states <- st_read(wfi_layers_region, layer = "us_states_no_coast", quiet = TRUE)
cities <- st_read(wfi_layers_region, layer = "major_cities", quiet = TRUE)
ocean <- st_read(wfi_layers_region, layer = "ocean", quiet = TRUE)
wor_cascades <- st_read(wfi_layers_region, layer = "western_cascades", quiet = TRUE)
large_fires_region <- st_read(wfi_layers_region, layer = "large_fires", quiet = TRUE)

#Subset cities
eugene <- cities %>% filter(TOWN == "Eugene")
portland <- cities %>% filter(TOWN == "Portland")
seattle <- cities %>% filter(TOWN == "Seattle")

#Plot the study area map
study_map <- ggplot() +
  geom_raster(data = wfi_hill_df,
              aes(x, y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_raster(data = wfi_rast_df,
              aes(x, y, fill = elev_m),
              alpha = .2) +
  scale_fill_gradientn(colours = c("white", "grey100", "cornsilk2", "darkolivegreen4", "darkolivegreen", "darkolivegreen", "darkolivegreen", "whitesmoke", "white", "white", "white"), breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200), guide = "none") +
  geom_sf(data=study_area_nfs, color="black", fill="darkolivegreen", linewidth=.2, alpha=.2) +
  geom_sf(data = ugbs, color="black", fill="grey70", linewidth=.1, alpha=.5) +
  geom_sf(data = large_fires, color="darkred", fill="darkred", linewidth=.05, alpha=.25) +
  geom_sf(data=lakes_res, fill="royalblue4", color="royalblue4", alpha=.85) +
  geom_sf(data=rivers, color="royalblue4", linewidth=.2, alpha=.85) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude, color=huc8_name), shape=22, fill="black", size=1.75, stroke=1.5) +
  scale_color_manual(name = "", values=c("Lower Columbia-Sandy" = "green2", "Clackamas"="orangered", "North Santiam"="cyan2", "South Santiam"="darkorchid2", "Mckenzie"="goldenrod1", "Middle Fork Willamette"="violetred1"), guide = "none") +
  geom_text(data=site_coords, aes(x=longitude, y=latitude, label=site), color="black", size=3, hjust = 0.0, nudge_x = 0.03, vjust = 0.0, nudge_y = 0.0025, family = "Arial", fontface = "bold") +
  annotate("text", x=-122.637, y=45.5, label= "Portland", size=4, family = "Arial") +
  scale_x_continuous("", breaks = c(-123, -122, -121), limits=c(-123, -121)) +
  scale_y_continuous(name="", breaks = c(43, 44, 45, 46), limits=c(43, 46), position = 'right') +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey10", "white"),
    text_family = "Arial",
    line_width = 1.5,
    height = unit(0.225, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(.9, "cm"), 
    width = unit(.9, "cm"),
    pad_x = unit(0.0, "cm"), pad_y = unit(0.15, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "grey80"),
      line_col = "grey10",
      text_family = "Arial"
    )) +
  guides(fill = guide_colorsteps(barwidth = 25,
                                 barheight = 1,
                                 title.position = "right")) +
  labs(fill = "m") +
  coord_sf(xlim = c(-122.75, -121.555), ylim = c(45.54, 43.62)) +
  theme(legend.position="none", axis.text=element_text(size=10), panel.border = element_rect(color="black", fill = NA, linewidth=.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin = margin(t=.218, r=.5, b=.764, l=0, unit = "cm"))

#Plot the region map
region_map <- ggplot() +
  geom_raster(data = region_hill_df,
              aes(x, y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_raster(data = region_rast_df,
              aes(x, y, fill = elev_m),
              alpha = .1) +
  scale_fill_gradientn(colours = c("white", "white", "white", "white", "white", "white", "white", "white", "white", "white", "cornsilk2", "darkolivegreen4", "darkolivegreen", "khaki3", "whitesmoke", "whitesmoke", "grey100", "grey100", "grey100"), breaks = c(-4000, -3500, -3000, -2500, -2000, -1500, -1000, -500, 0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500), guide = "none") +
  geom_sf(data = wor_cascades, color="grey5", fill="olivedrab4", linewidth=.15, alpha=.5) +
  geom_sf(data = ocean, color=NA, fill="steelblue", alpha=1) +
  geom_sf(data = states, color="black", fill=NA, linewidth=.2, alpha=1) +
  geom_sf(data = large_fires, color="darkred", fill="darkred", linewidth=.025, alpha=.2) +
  #geom_point(data=site_coords, aes(x=longitude, y=latitude, color=huc8_name), shape=22, fill="black", size=.5, stroke=.15) +
  #scale_color_manual(name = "", values=c("Middle Fork Willamette"="violetred1", "Clackamas"="orangered", "Lower Columbia-Sandy" = "green2", "North Santiam"="cyan2", "Mckenzie"="goldenrod1", "South Santiam"="darkorchid2"), guide = "none") +
  geom_rect(aes(xmin = -122.835, xmax = -121.521, ymin = 43.525, ymax = 45.6), fill = "transparent", color = "black", linewidth = .9) +
  #geom_sf_label(data=cities, aes(label = TOWN)) +
  annotate("text", x=-121.92, y=47.6, label= "Seattle", size=3.5, family = "Arial") +
  annotate("text", x=-122.33, y=45.5, label= "Portland", size=3.5, family = "Arial") +
  annotate("text", x=-123.3, y=44.05, label= "Eugene", size=3.5, family = "Arial") + 
  annotate("text", x=-121.22, y=44.025, label= "Bend", size=3.5, family = "Arial") + 
  annotate("text", x=-122.25, y=46.72, label= "Washington", size=5, family = "Arial") +
  annotate("text", x=-121.74, y=43.1, label= "Oregon", size=5, family = "Arial") +
  annotate("text", x=-122.4, y=41.4, label= "California", size=5, family = "Arial") +
  scale_x_continuous("", breaks = c(-124, -122)) +
  scale_y_continuous("", breaks = c(42, 44, 46, 48)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey10", "white"),
    text_family = "Arial",
    line_width = 1.5,
    height = unit(0.225, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(.9, "cm"), 
    width = unit(.9, "cm"),
    pad_x = unit(0.0, "cm"), pad_y = unit(0.15, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "grey80"),
      line_col = "grey10",
      text_family = "Arial"
    )) +
  guides(fill = guide_colorsteps(barwidth = 25,
                                 barheight = 1,
                                 title.position = "right")) +
  labs(fill = "m") +
  coord_sf(xlim = c(-124.53, -120.75), y = c(48.2, 41)) +
  theme(legend.position="none", axis.text=element_text(size=10), panel.border = element_rect(color="black", fill = NA, linewidth=.5), plot.margin = margin(t=.65, r=0, b=.4, l=0, unit = "cm"))

#Combine both graphs and print... it may take 5-10 minutes to print (plotting rasters is time intensive).
quartz(width=8, height=8.25)
combine_maps <- ggpubr::ggarrange(region_map, study_map, labels = c("", ""), ncol = 2, nrow = 1, align = "v", widths = c(1, 1.19), heights = c(1, 1.19))
combine_maps
ggsave(path="/YOURPATH", "combine_maps.png", width = 8, height = 8.25, units = "in", dpi=500)



######################################
#### Evaluate environmental space ####
######################################

#Here we're going to compare the environmental space of the present study's data collection sites to the western United States and previously published fire history studies in the western US.

#Clear workspace and set working directory
rm(list=ls())
setwd("/YOURPATH")

#Load libraries for this section
library(tidyverse)
library(PMCMRplus)

#Read environmental space data.  This data consists of net primary productivity, 30-year normal mean monthly temperature, 30-year normal maximum monthly temperature, total annual precipitation, and net primary productivity (NPP) for approximately 932,000 systematically located points in forested regions of the western United States (11 western states; "wstates_grid"), as well as the same data for the location of 36 study sites located for the present study ("wsf"), and 1,666 NAFSN fire history study sites found across the western United States (see Margollis et al. 2022; "nafsn").  This is a large file. 
points <- read.csv("points_space.csv", header = TRUE)
nrow(points)
table(points$study)
tail(points)

#Count -9999 PPT observations... -9999 indicates no value for NPP. temp, or precip.  
names(points)
nrow(filter(points, latitude == -9999))
nrow(filter(points, longitude == -9999))
nrow(filter(points, prism_tmax_c == -9999))
nrow(filter(points, prism_ppt_mm == -9999))
nrow(filter(points, modis_npp == -9999))

#There are missing values for npp, possibly because these points fall on waterbodies or other areas where these values can't be calculated.  We'll remove them
nrow(points)
points <- points[rowSums(points == -9999)==0, ,drop = FALSE]
nrow(points)

#Summarize value of points by type and parameter
points_summary <- points %>%
  group_by(study) %>%
  summarize(
    mean_PPT = mean(prism_ppt_mm), 
    mean_tmean = mean(prism_tmean_c),
    mean_Tmax = mean(prism_tmax_c),
    mean_npp = mean(modis_npp))
points_summary

#Graph all the variables
#First melt the data into long format
head(points)
points_melt <- reshape::melt(points, id.vars='study', measure.vars=c('prism_tmean_c', 'prism_tmax_c', 'prism_ppt_mm', 'modis_npp'))
head(points_melt)

#Boxplots of all variables faceted by variables
all_variable_box <- ggplot(points_melt) +
  geom_boxplot(aes(x=study, y=value, fill=study)) +
  #scale_fill_manual(name = "", values = c('All NF Land' = 'goldenrod4', 'WFI Samples' = 'cornflowerblue')) +
  scale_x_discrete("") +
  scale_y_continuous("Value") +
  facet_wrap(~variable, scales = "free", ) +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(legend.position="none")
#all_variable_box

#Environmental space graphic pt. 1

#Plot Panel A of Figure 4:  A boxplot of single variable (e.g. NPP)... note that I've cut off some outlying observations, note to self:  Be sure and note this in caption.  
points$study <- factor(points$study,levels=c("wstates_grid", "nafsn", "wsf"))
one_variable_box <- ggplot(points, aes(study, modis_npp, fill=study)) +
  geom_boxplot(position="dodge") + 
  scale_fill_manual(name = "", values = c('wstates_grid' = 'goldenrod4', 'nafsn' = 'firebrick4', 'wsf' = 'darkolivegreen'), labels=c('Western states forestland', 'NAFSN', 'This study')) +
  #scale_x_discrete("") + 
  scale_y_continuous(bquote('Annual NPP ('*'g'~C~m^-2*')'), breaks=c(0, 5000, 10000), limits=c(0,12000)) +
  theme_bw(base_size = 22, base_family = "Arial") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
ggsave(path="/YOURPATH", "one_variable_box.png", width = 9.5, height = 7.5, units = "in", dpi=500)

#A t-test is probably not the right statistical test in this case (because of potential spatial autocorrelation of the western states grid).  But we'll try it anyway.  I also tried a spaMM model (see below).  But that model wouldn't run on 900,000+ observations.  I tried a random subset of the data and got the same basic result as a t-test (the random subset didn't exhibit autocorrelation).  
#Subset data
grid <- filter(points, study=="wstates_grid")
wfi  <- points %>% filter(study == "wsf")
nafsn  <- filter(points, study == "nafsn")
#Run t-tests
t.test(grid$modis_npp, wfi$modis_npp)
t.test(grid$prism_tmean_c, wfi$prism_tmean_c)
t.test(grid$prism_tmax_c, wfi$prism_tmax_c)
#According to the t-test, there is strong evidence of a difference in mean NPP, temp, and precip between WFI and the systematic sample of westside environmental space.

#Games-Howell test: pairwise comparisons without assuming equal variances.  This is a better approach than a t-test, but still doesn't resolve the spatial dependence issue.  
#Some cleaning to make the Games Howell function work
points_clean <- points %>%
  filter(!is.na(study), study %in% c("wsf", "nafsn", "wstates_grid")) %>%
  filter(!is.na(modis_npp), is.finite(modis_npp)) %>%
  mutate(study = droplevels(factor(study)))
#Games Howell test for npp.... All pairwise differences in mean NPP are statistically significant, even when accounting for unequal group sizes (wsf = 36, nafsn = 1666, wstates_grid = 900k+) and unequal variances... our study sites have much higher mean NPP than both nafsn and the general grid.
gamesHowellTest(modis_npp ~ study, data = points_clean)
#Games Howell test for tmax:  Present study has a significantly lower maximum monthly temperature than NAFSN sites....  Present study vs western states grid:  No statistically significant difference in maximum temperature (p = 0.095), although this suggests suggest slightly cooler conditions at WFI sites.
gamesHowellTest(prism_tmax_c ~ study, data = points_clean)
#Games Howell test for precipitation: Significant difference between all groups.
gamesHowellTest(prism_ppt_mm ~ study, data = points_clean)

#But the above approaches don't really account for spatial dependence.  So here's what we're going to do:
#Calculate the observed mean for each of three environmental variables — annual precipitation, maximum monthly temperature, and net primary productivity (NPP) — across the 36 WFI sites.
#Randomly select 36 sites from the western states grid dataset (without replacement) and computed a mean value for each environmental variable across that sample.
#Repeat this procedure 10,000 times to generate a null distribution of permuted means for each variable (I used 10,000 permutations in the paper).
#Calculate a two-tailed empirical p-value as the proportion of permuted means that were as or more extreme than the observed mean for sites reported in the present study relative to the permutation distribution.
#Function for a spatial permutation test
permutation_test <- function(var_name, reference_df, label = "Permutation Test", n_reps = 10000, seed = 42) {
  set.seed(seed)
  #Extract observed and reference values
  obs_vals <- wfi[[var_name]]
  ref_vals <- reference_df[[var_name]]
  #Drop missing values
  obs_vals <- obs_vals[is.finite(obs_vals)]
  ref_vals <- ref_vals[is.finite(ref_vals)]
  obs_mean <- mean(obs_vals)
  n_obs <- length(obs_vals)
  #Permutation loop
  perm_means <- replicate(n_reps, {
    mean(sample(ref_vals, n_obs), na.rm = TRUE)
  })
  #Two-tailed empirical p-value
  p_val <- mean(abs(perm_means - mean(perm_means)) >= abs(obs_mean - mean(perm_means)))
  #Optional:  Plot
  #hist(perm_means, breaks = 50, col = "lightgray",
  #     main = paste(label, ":", var_name),
  #     xlab = paste("Permuted Mean", var_name))
  #abline(v = obs_mean, col = "red", lwd = 2)
  #legend("topright", legend = paste("Observed =", round(obs_mean, 2)), col = "red", lwd = 2)
  #Return result
  tibble(
    variable = var_name,
    observed_mean = obs_mean,
    permuted_mean = mean(perm_means),
    p_value = p_val
  )
}

#Run function on western states grid
grid_results <- bind_rows(
  permutation_test("prism_ppt_mm", grid, "WFI vs Grid"),
  permutation_test("modis_npp", grid, "WFI vs Grid"),
  permutation_test("prism_tmax_c", grid, "WFI vs Grid")
)
grid_results

#Run function on NAFSN sites
nafsn_results <- bind_rows(
  permutation_test("prism_ppt_mm", nafsn, "WFI vs NAFSN"),
  permutation_test("modis_npp", nafsn, "WFI vs NAFSN"),
  permutation_test("prism_tmax_c", nafsn, "WFI vs NAFSN")
)
nafsn_results

#Environmental space graphic pt. 2

#Read environmental variables file, select HUC8 names, reduce to most common HUC8 name per site, and join to site coordinates
env_vars <- read.csv("env_data.csv")
env_vars <- env_vars %>% select(site, site_name, huc8_name)
env_vars <- env_vars %>% group_by(site, site_name) %>% count(huc8_name) %>% slice(which.max(n))
env_vars$site <- as.character(env_vars$site)
wfi <- wfi %>% left_join(env_vars) 

#Plot it (Panel B of Figure 4)... may take a minute to plot (lots of western states points)
env_dense <- ggplot(grid, aes(prism_tmax_c, prism_ppt_mm)) +
  #geom_point(color="darkolivegreen", shape=1, size=.025, alpha=.15) +
  stat_density2d(aes(fill = after_stat(level)), geom="polygon", bins=15, alpha=1) +
  stat_density2d(bins=15, color="black", linewidth=.25, alpha=.75) +
  geom_point(data=nafsn, aes(prism_tmax_c, prism_ppt_mm), color="red", shape=24, size=.95, alpha=.75) +
  geom_point(data=wfi, aes(prism_tmax_c, prism_ppt_mm, color=huc8_name), shape=22, fill="black", size=3.25, stroke=1.75) +
  geom_text(
    data = wfi,
    aes(prism_tmax_c, prism_ppt_mm, label = site),
    color = "black",
    size = 4,                         
    hjust = 0.0,
    nudge_x = 0.175,                   
    vjust = 0.0,
    nudge_y = 0.0375,                 
    family = "Arial",
    fontface = "bold") +
  scale_fill_gradient("", low = "grey90", high = "darkolivegreen") +
  scale_x_continuous("TMax (C)") +
  scale_y_continuous("PPT (mm)") +
  theme_bw(base_size = 22, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none"  
  )
ggsave(path="/YOURPATH", "env_dense.png", width = 9.5, height = 8.5, units = "in", dpi=500)


#########################################
######### Create fire histories #########
#########################################

#Clear workspace and set working directory
rm(list=ls())
setwd("/YOURPATH")

#Libraries for this section of code
library(tidyverse)
library(tidytext)
library(leaflet)
library(KernSmooth)
library(ggplot2)
library(ggnewscale)

#Read cross-sections file, select pertinent columns, and drop series with no first or last year
xsections <- read.csv("xsections.csv")
xsections <- xsections[,c("site", "sample_id", "species", "sample_ht_cm", "d10_mm", "pith","outer_ring", "inner_ring", "latitude", "longitude")]
xsections <- xsections %>% drop_na(outer_ring) 
xsections <- xsections %>% drop_na(inner_ring) 
nrow(xsections)

#Make all samples 6 digits (to remove, A, B, C, D, etc. because we are combining histories from the same tree), and convert numeric columns to numeric
xsections$sample_id <- substr(xsections$sample_id, 1, 6)
xsections$pith <- as.numeric(xsections$pith)
xsections$outer_ring <- as.numeric(xsections$outer_ring)
xsections$inner_ring <- as.numeric(xsections$inner_ring)
xsections$sample_ht_cm <- as.numeric(xsections$sample_ht_cm)
xsections$d10_mm <- as.numeric(xsections$d10_mm)
nrow(xsections)
length(unique(xsections$sample_id))

#Collapse data from A, B, C, D, etc. of same sample and make result a dataframe (for some reason the tibble output seems to choke the establishment data model).
xsections_collapse <- xsections %>% 
  group_by(site, sample_id, species) %>% 
  summarize(
    latitude=collapse::fmin(latitude, na.rm = TRUE),
    longitude=collapse::fmin(longitude, na.rm = TRUE),
    sample_ht_cm=collapse::fmin(sample_ht_cm, na.rm = TRUE),
    d10_mm=collapse::fmin(d10_mm, na.rm = TRUE),
    outer_ring=collapse::fmax(outer_ring, na.rm = TRUE),
    inner_ring=collapse::fmin(inner_ring, na.rm = TRUE),
    pith=collapse::fmin(pith, na.rm = TRUE))
xsections_collapse <- data.frame(xsections_collapse)
str(xsections_collapse)
target_xsections <- xsections_collapse 
unique(target_xsections$site)
nrow(target_xsections)

#Read fires file, remove non-fire injuries, and select pertinent columns. 
events <- read.csv("events.csv")
events <- events %>% filter(type=="FS") 
events <- events[,c("site", "sample_id", "event_year")]
target_events <- events 

#Change continuous variables to numeric
target_xsections$pith <- as.numeric(target_xsections$pith)
target_xsections$outer_ring <- as.numeric(target_xsections$outer_ring)
target_xsections$inner_ring <- as.numeric(target_xsections$inner_ring)

#Read height calibration data
hts <- read.csv("ht_calibration.csv")

#Add dummy species variable to xsections file to match the column name in the height calibration file and change these dummy species names that aren't in the ht_calibration file (so as not to alter the species as it will be coded below)
target_xsections$spp <- target_xsections$species
target_xsections$spp <- ifelse(target_xsections$spp=="" | target_xsections$spp=="THPL" | target_xsections$spp=="PIMO" | target_xsections$spp=="TSME" | target_xsections$spp=="UNKN" | target_xsections$spp=="LAOC" | target_xsections$spp=="ABPR?" | target_xsections$spp=="TSHE?", "PSME", target_xsections$spp)
target_xsections$d10_mm <- as.numeric(target_xsections$d10_mm)
target_xsections$sample_ht_cm <- as.numeric(target_xsections$sample_ht_cm)

#Create a linear model for years to mineral soil
lm1 <- lm(rings_to_soil ~ sample_ht_cm + spp + d10_mm, data=hts)

#Predict years to mineral soil for xsections, bind that to xsections file, and make any predictions that are negative numbers 0
predict_years <- predict(lm1, newdata = target_xsections)
target_xsections <- cbind(target_xsections, predict_years)
target_xsections$predict_years <- ifelse(target_xsections$predict_years < 0, 0, target_xsections$predict_years)
target_xsections$pith <- as.numeric(target_xsections$pith)
target_xsections$predict_years <- round(target_xsections$predict_years)
target_xsections$estab_year <- target_xsections$pith-target_xsections$predict_years

#Cohort detection function (note new arrangement for padding)
cohort_function <- function(x, yearsvec, nosims, siglevel){
  estab <- x[,yearsvec]
  minest <- min(estab)
  maxest <- max(estab)
  minest_pad <- min(estab)-100
  maxest_pad <- max(estab)+100
  notrees <- length(estab)
  nosims <- nosims
  unif_resampfunct <- function(){
    simestab <- sample(minest:maxest, size = notrees, replace = TRUE)
    return(simestab)
  }
  unif_resampfunct_pad <- function(){
    simestab <- sample(minest_pad:maxest_pad, size = notrees, replace = TRUE)
    return(simestab)
  }
  sim_data <- data.frame(replicate(n = nosims, 
                                   expr = if(max(estab)-min(estab) >= 250)
                                     unif_resampfunct() else unif_resampfunct_pad()))
  names(sim_data) <- gsub(x = names(sim_data), pattern = "\\X", replacement = "sim") 
  all_data <- cbind(sim_data, estab)
  bw_wide <- dpik(estab)/1.25
  bw_narrow <- dpik(estab)/1.5
  all_density <- apply(all_data, 2, bkde, bandwidth=if(max(estab)-min(estab) >= 150) bw_narrow else bw_wide)  
  long_density <- do.call(rbind.data.frame, all_density)
  long_density$run <- rownames(long_density)
  long_density$group <- substr(long_density$run, 1, 3)
  sim_dens <- long_density[(long_density$group=="sim"),]
  est_dens <- long_density[(long_density$group=="est"),]
  sim_dens_sig <- quantile(sim_dens$y, siglevel)
  if(any(est_dens$y >= sim_dens_sig)){    
    estab_sig <- est_dens[(est_dens$y >= sim_dens_sig),]
    estab_sig$lower_bound <- rep(sim_dens_sig, nrow(estab_sig))
    estab_sig <- estab_sig[,c("x", "y", "lower_bound")]
    estab_sig$run_no <- as.numeric(sub(".*b.", "", 
                                       rownames(estab_sig)))
    rownames(estab_sig) <- 1:nrow(estab_sig)
    estab_sig$rows <- as.numeric(rownames(estab_sig))
    ints <- c(0, which(diff(estab_sig$run_no) != 1), length(estab_sig$run_no))
    estab_sig$cohort_no <- cut(estab_sig$rows, breaks=ints, labels=FALSE)
    estab_sig$rows <- NULL
  } else {
    estab_sig <- data.frame(x = 0, y = 0, lower_bound = 0, run_no = 0, cohort_no="NA")
  }
  cohort_table <- estab_sig %>%
    group_by(cohort_no) %>%
    summarize(min_cohort_year = round(min(x),0),
              max_cohort_year = round(max(x),0))
  cohort_table
}

#Create a cross-sections file for the cohort detection procedure
xsections_cohorts <- target_xsections %>% drop_na(estab_year)
xsections_cohorts <- xsections_cohorts %>% select(site, sample_id, estab_year) 

#Split the xsections dataframe into a list and run function on every element of the list.  Important note:  The "map" function is a purrr function which seems to work without loading purrr explicitly.  This library has given me problems in the past (an old version is loaded as of this writing).  If the code breaks down here, it's probably a problem with purrr. The 1,000 simulations specified here and at other points in the code (the simulation procedure below is repeated in other sections) should take less than minute.  The 10,000 simulations used in the paper will take longer.  
set.seed(1234)
site_list <- split(xsections_cohorts, xsections_cohorts$site)
results_list <- map(site_list, cohort_function, yearsvec="estab_year", nosims=10000, siglevel=.99)

#Convert list to dataframe with site column and convert 0 to NA and add unique site/cohort identifier
wfi_cohorts <- do.call(rbind, unname(Map(cbind, site = names(results_list), results_list)))
wfi_cohorts[wfi_cohorts == 0] <- NA
wfi_cohorts[order(-wfi_cohorts$min_cohort_year),] 
wfi_cohorts$site_cohort <- paste(wfi_cohorts$site, wfi_cohorts$cohort_no, sep=" no")

#Optional:  Fix issue with site
wfi_cohorts <- wfi_cohorts[!wfi_cohorts$site_cohort == "24 no1",]

#Get all the years between the first and last cohort year for each site, pad the years by two, and rename "year" to "match_year" to match with target_xsections file
cohort_seq <- wfi_cohorts %>% select(site, min_cohort_year, max_cohort_year)
cohort_seq$min_cohort_year <- cohort_seq$min_cohort_year-2
cohort_seq$max_cohort_year <- cohort_seq$max_cohort_year+2
cohort_seq <- cohort_seq %>% drop_na()
cohort_years <- cohort_seq %>%
  transmute(site, year = map2(min_cohort_year, max_cohort_year, `:`)) %>%
  unnest(cols = c(year))
cohort_years$site <- as.integer(cohort_years$site)
cohort_years$cohort <- "TRUE"
colnames(cohort_years)[colnames(cohort_years) == 'year'] <- 'match_year'

#Join with target_xsections file (make site into factor), and change NA to "False".
target_xsections$match_year <- target_xsections$estab_year
cohort_years$site <- as.factor(cohort_years$site)
target_xsections$site <- as.factor(target_xsections$site)
cohort_years$match_year <- as.integer(cohort_years$match_year)
target_xsections$match_year <- as.integer(target_xsections$match_year)
target_xsections <- target_xsections %>% left_join(cohort_years, by=c("site", "match_year"))
target_xsections$cohort <- ifelse(target_xsections$cohort=="TRUE", "True", "False")
target_xsections$cohort <- target_xsections$cohort %>% replace_na("False")
head(target_xsections)

#Make a new column with the oldest date present on either sample, which is either the establishment year or the inner ring year.  This new column will be used to order samples within sites below.
target_xsections <- target_xsections %>% select(site, sample_id, species, outer_ring, inner_ring, estab_year, cohort)
target_xsections$last_graph <- with(target_xsections, pmin(estab_year, inner_ring, na.rm = TRUE))

#Read environmental variables file (for ordering samples by elevation)
env_vars <- read.csv("env_data.csv")
unique(env_vars$site)
#Get average snow duration day by site
env_vars <- env_vars %>% select(site, huc8_name, sdd_average)
forsort <- env_vars %>% group_by(site, huc8_name)  %>% 
  summarize(mean_sdd=mean(sdd_average))
forsort <- forsort %>% arrange(huc8_name, -mean_sdd)
print(forsort, n=50)
data.frame(forsort)

#Order sites by watershed (north to south) and then snow day disappearance (highest to lowest)
#unique(target_xsections$site)
target_xsections$site <- factor(target_xsections$site, levels=c("55", "58", "54", "56", "57", "62", "53", "59", "64", "52", "61", "60", "68", "71", "72", "85", "90", "89", "92", "93", "88", "87", "91", "83", "86", "84", "81", "82", "96", "98", "26", "24",  "22", "23", "25", "28"))
target_events$site <- factor(target_events$site, levels=c("55", "58", "54", "56", "57", "62", "53", "59", "64", "52", "61", "60", "68", "71", "72", "85", "90", "89", "92", "93", "88", "87", "91", "83", "86", "84", "81", "82", "96", "98", "26", "24",  "22", "23", "25", "28"))
wfi_cohorts$site <- factor(wfi_cohorts$site, levels=c("55", "58", "54", "56", "57", "62", "53", "59", "64", "52", "61", "60", "68", "71", "72", "85", "90", "89", "92", "93", "88", "87", "91", "83", "86", "84", "81", "82", "96", "98", "26", "24",  "22", "23", "25", "28"))

#Optional:  Order sites by stand age instead
#unique(target_xsections$site)
#target_xsections$site <- factor(target_xsections$site,levels=c("84", "86", "22", "85", "28", "61", "60", "62", "64", "98", "96", "57", "55", "24", "71", "59", "81", "72", "82", "88", "23", "89", "25", "87", "90", "26", "54", "83", "58", "56", "68", "52", "53", "92", "93", "91"))
#target_events$site <- factor(target_events$site, levels=c("84", "86", "22", "85", "28", "61", "60", "62", "64", "98", "96", "57", "55", "24", "71", "59", "81", "72", "82", "88", "23", "89", "25", "87", "90", "26", "54", "83", "58", "56", "68", "52", "53", "92", "93", "91"))
#wfi_cohorts$site <- factor(wfi_cohorts$site, levels=c("84", "86", "22", "85", "28", "61", "60", "62", "64", "98", "96", "57", "55", "24", "71", "59", "81", "72", "82", "88", "23", "89", "25", "87", "90", "26", "54", "83", "58", "56", "68", "52", "53", "92", "93", "91"))

#Order sample IDs within site facets by first year
target_xsections <- target_xsections %>%
  group_by(site) %>%
  ungroup() %>%
  mutate(site=as.factor(site),
         sample_id=reorder_within(site, last_graph, sample_id)) 
#Give the reordered sample IDs a new name, create the old column name ("sample_id") from the new column name, and use it join to the fires column (so the fires sample_id column can also be reordered)
target_xsections <- rename(target_xsections, "sample_id_match"="sample_id")
target_xsections$sample_id <- stringi::stri_sub(target_xsections$sample_id_match, -6)
target_xsections_to_match <- target_xsections %>%
  select(sample_id, sample_id_match)
target_events <- left_join(target_events, target_xsections_to_match, by = "sample_id")

#Determine if any fires preceded overlapped minimum cohort year... first create list of all years 7 years before and 3 years after first cohort year by site and cohort, then join to fires, then match... this is really clunky code...  might be simpler to do it manually
wfi_cohorts$fire_search_min <- wfi_cohorts$min_cohort_year-10
wfi_cohorts$fire_search_max <- wfi_cohorts$min_cohort_year+5
wfi_cohorts <- wfi_cohorts[complete.cases(wfi_cohorts), ]
fire_search <- wfi_cohorts %>%
  transmute(site, cohort_no, cohort_initiate = map2(fire_search_min, fire_search_max, `:`)) %>%
  unnest(cols = c(cohort_initiate))
fire_search$cohort_no <- NULL
cohort_fires <- target_events %>%
  left_join(fire_search) %>%
  group_by(site, sample_id_match) %>%
  filter(event_year %in% cohort_initiate) %>%
  distinct(site, event_year)
cohort_fires

#Get unique cohort initiating fires
cohort_fires_unique <- cohort_fires %>% 
  group_by(site) %>% 
  distinct(event_year) 
cohort_fires_unique

#Make segment lengths match the breaks
first_year_graph <- 1150
last_year_graph <- 2022
target_xsections$outer_ring <- ifelse(target_xsections$outer_ring < first_year_graph, first_year_graph, target_xsections$outer_ring)
target_xsections$inner_ring <- ifelse(target_xsections$inner_ring > last_year_graph, last_year_graph, target_xsections$inner_ring)

#Rename species
target_xsections$species <- ifelse(target_xsections$species=="ABAM", "Other", ifelse(target_xsections$species=="ABGR", "Other", ifelse(target_xsections$species=="PIMO", "Other", ifelse(target_xsections$species=="PIPO", "Other", ifelse(target_xsections$species=="THPL", "Other", ifelse(target_xsections$species=="UNKN", "Other", target_xsections$species))))))
#Order species
target_xsections$species <- factor(target_xsections$species, levels=c("PSME", "TSHE", "TSME", "ABPR", "Other"))

#Create dataframe for boxes indicating watersheds
boxes <- forsort
unique(boxes$huc8_name)
boxes$minyear <- first_year_graph
boxes$maxyear <- last_year_graph
#Lower Columbia-Sandy box
sandbox <- boxes %>% filter(huc8_name=="Lower Columbia-Sandy")
sandbox <- sandbox %>% select(site, minyear, maxyear)
#Clackamas box
clackbox <- boxes %>% filter(huc8_name=="Clackamas")
clackbox <- clackbox %>% select(site, minyear, maxyear)
#Northern Santiam box
nsantbox <- boxes %>% filter(huc8_name=="North Santiam")
nsantbox <- nsantbox %>% select(site, minyear, maxyear)
#Southern Santiam box
ssantbox <- boxes %>% filter(huc8_name=="South Santiam")
ssantbox <- ssantbox %>% select(site, minyear, maxyear)
#Mckenzie box
macbox <- boxes %>% filter(huc8_name=="Mckenzie")
macbox <- macbox %>% select(site, minyear, maxyear)
#Middle Fork box
midbox <- boxes %>% filter(huc8_name=="Middle Fork Willamette")
midbox <- midbox %>% select(site, minyear, maxyear)

#Graph all samples from all sites
fire_chart <- ggplot() +
  geom_rect(data=sandbox, aes(xmin=minyear, xmax=maxyear, ymin=-Inf, ymax=Inf), color="green2", linewidth=0.025, fill="green2", alpha=.1) +
  geom_rect(data=clackbox, aes(xmin=minyear, xmax=maxyear, ymin=-Inf, ymax=Inf), color="orangered", linewidth=0.025, fill="orangered", alpha=.1) +
  geom_rect(data=nsantbox, aes(xmin=minyear, xmax=maxyear, ymin=-Inf, ymax=Inf), color="cyan2", linewidth=0.025, fill="cyan2", alpha=.1) +
  geom_rect(data=ssantbox, aes(xmin=minyear, xmax=maxyear, ymin=-Inf, ymax=Inf), color="darkorchid2", linewidth=0.025, fill="darkorchid2", alpha=.1) +
  geom_rect(data=macbox, aes(xmin=minyear, xmax=maxyear, ymin=-Inf, ymax=Inf), color="goldenrod1", linewidth=0.025, fill="goldenrod1", alpha=.1) +
  geom_rect(data=midbox, aes(xmin=minyear, xmax=maxyear, ymin=-Inf, ymax=Inf), color="violetred1", linewidth=0.025, fill="violetred1", alpha=.1) +
  geom_segment(data=target_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=.4, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=target_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=.7) +
  geom_point(data=target_events, aes(x = event_year, y = sample_id_match), size=.4, shape = 23, color="darkred", fill="darkred") +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey45')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000), limits=c(first_year_graph, last_year_graph)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.25, .25)) +
  facet_grid(site ~ ., scales = "free", space = "free") +
  theme_bw(base_size=12, base_family = "Arial") +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(strip.background = element_blank(), panel.border = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(strip.text.y.right = element_text(angle = 0), panel.spacing.y = unit(0.0, "lines")) 
ggsave(path="/YOURPATH", "fire_chart.png", width = 5.25, height = 8, units = "in", dpi=800)

#Graph individual site (site 64) for Figure 1
site_xsections <- target_xsections %>% filter(site==64)
site_events <- target_events %>% filter(site==64)
site_cohort <- wfi_cohorts %>% filter(site==64)
site_box <- clackbox %>% filter(site==64)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1950, 1950, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1450, 1450, site_xsections$inner_ring)
site64_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1450, xmax=1950, ymin=-Inf, ymax=Inf), color="orangered", linewidth=0.025, fill="orangered", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.4) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=.9) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=.9, shape = 23, color="darkred", fill="darkred") +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey45')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1600, 1700, 1800, 1900), limits=c(1450, 1950), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.05, 0.05)) +
  theme_bw(base_size=15, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site64_chart.png", width = 1.7, height = 1.9, units = "in", dpi=800)

#Graph individual site (site 28) for Figure 1
site_xsections <- target_xsections %>% filter(site==28)
site_events <- target_events %>% filter(site==28)
site_cohort <- wfi_cohorts %>% filter(site==28)
site_box <- midbox %>% filter(site==28)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1950, 1950, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1450, 1450, site_xsections$inner_ring)
site28_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1450, xmax=1950, ymin=-Inf, ymax=Inf), color="violetred1", linewidth=0.025, fill="violetred1", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.4) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=.9) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=.9, shape = 23, color="darkred", fill="darkred") +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey45')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1600, 1700, 1800, 1900), limits=c(1450, 1950), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.05, 0.05)) +
  theme_bw(base_size=15, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site28_chart.png", width = 1.7, height = 1.9, units = "in", dpi=800)

#Graph individual site (site 93) for Figure 1
site_xsections <- target_xsections %>% filter(site==93)
site_events <- target_events %>% filter(site==93)
site_cohort <- wfi_cohorts %>% filter(site==93)
site_box <- ssantbox %>% filter(site==93)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1950, 1950, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1450, 1450, site_xsections$inner_ring)
site93_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1450, xmax=1950, ymin=-Inf, ymax=Inf), color="darkorchid2", linewidth=0.025, fill="darkorchid2", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year+20, xmax=max_cohort_year-30, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.4) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=.9) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=.9, shape = 23, color="darkred", fill="darkred") +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey45')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1600, 1700, 1800, 1900), limits=c(1450, 1950), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.05, 0.05)) +
  theme_bw(base_size=15, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_blank()) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site93_chart.png", width = 1.7, height = 1.9, units = "in", dpi=800)

#Graph individual site for large fires graphic (site 54)
site_xsections <- target_xsections %>% filter(site==54)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1560, 1560, site_xsections$outer_ring)
site_events <- target_events %>% filter(site==54)
site_box <- sandbox %>% filter(site==54)
site_cohort <- wfi_cohorts %>% filter(site==54)
site54_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1480, xmax=1560, ymin=-Inf, ymax=Inf), color="green2", linewidth=0.025, fill="green2", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.8, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=1.5, stroke=1.5) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("54"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1520, 1540), limits=c(1480, 1560), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site54_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for large fires graphic (site 56)
site_xsections <- target_xsections %>% filter(site==56)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1560, 1560, site_xsections$outer_ring)
site_events <- target_events %>% filter(site==56)
site_cohort <- wfi_cohorts %>% filter(site==56)
site_box <- sandbox %>% filter(site==56)
site56_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1480, xmax=1560, ymin=-Inf, ymax=Inf), color="green2", linewidth=0.025, fill="green2", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.8, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=1.5, stroke=1.5) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("56"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1520, 1540), limits=c(1480, 1560), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site56_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for large fires graphic (site 58)
site_xsections <- target_xsections %>% filter(site==58)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1560, 1560, site_xsections$outer_ring)
site_events <- target_events %>% filter(site==58)
site_cohort <- wfi_cohorts %>% filter(site==58)
site_box <- sandbox %>% filter(site==58)
site58_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1480, xmax=1560, ymin=-Inf, ymax=Inf), color="green2", linewidth=0.025, fill="green2", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.7, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=1.5, stroke=1.5) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("58"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1520, 1540), limits=c(1480, 1560), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site58_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for large fires graphic (site 55)
site_xsections <- target_xsections %>% filter(site==55)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1560, 1560, site_xsections$outer_ring)
site_events <- target_events %>% filter(site==55)
site_cohort <- wfi_cohorts %>% filter(site==55)
site_box <- sandbox %>% filter(site==55)
site55_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1480, xmax=1560, ymin=-Inf, ymax=Inf), color="green2", linewidth=0.025, fill="green2", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.8, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=1.5, stroke=1.5) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("55"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1520, 1540), limits=c(1480, 1560), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site55_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for large fires graphic (site 64)
site_xsections <- target_xsections %>% filter(site==64)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1750, 1750, site_xsections$outer_ring)
site_events <- target_events %>% filter(site==64)
site_box <- clackbox %>% filter(site==64)
site_cohort <- wfi_cohorts %>% filter(site==64)
site64large_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1680, xmax=1760, ymin=-Inf, ymax=Inf), color="orangered", linewidth=0.025, fill="orangered", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.7, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=1.5, stroke=1) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("64"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1680, 1700, 1720, 1740, 1760), limits=c(1680, 1760), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site64large_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for large fires graphic (site 62)
site_xsections <- target_xsections %>% filter(site==62)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1750, 1750, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1680, 1680, site_xsections$inner_ring)
site_events <- target_events %>% filter(site==62)
site_box <- clackbox %>% filter(site==62)
site_cohort <- wfi_cohorts %>% filter(site==62)
site62_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1680, xmax=1760, ymin=-Inf, ymax=Inf), color="orangered", linewidth=0.025, fill="orangered", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.7, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=2.25, stroke=1.25) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("62"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1680, 1700, 1720, 1740, 1760), limits=c(1680, 1760), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site62_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for large fires graphic (site 61)
site_xsections <- target_xsections %>% filter(site==61)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1750, 1750, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1680, 1680, site_xsections$inner_ring)
site_events <- target_events %>% filter(site==61)
site_box <- clackbox %>% filter(site==61)
site_cohort <- wfi_cohorts %>% filter(site==61)
site61_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1680, xmax=1760, ymin=-Inf, ymax=Inf), color="orangered", linewidth=0.025, fill="orangered", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.7, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=2.25, stroke=1.25) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("61"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1680, 1700, 1720, 1740, 1760), limits=c(1680, 1760), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site61_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for large fires graphic (site 60)
site_xsections <- target_xsections %>% filter(site==60)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1750, 1750, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1680, 1680, site_xsections$inner_ring)
site_events <- target_events %>% filter(site==60)
site_box <- clackbox %>% filter(site==60)
site_cohort <- wfi_cohorts %>% filter(site==60)
site60_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1680, xmax=1760, ymin=-Inf, ymax=Inf), color="orangered", linewidth=0.025, fill="orangered", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.7, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=2.25, stroke=1.25) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("60"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1680, 1700, 1720, 1740, 1760), limits=c(1680, 1760), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site60_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for large fires graphic (site 52)
site_xsections <- target_xsections %>% filter(site==52)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1750, 1750, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1680, 1680, site_xsections$inner_ring)
site_events <- target_events %>% filter(site==52)
site_box <- clackbox %>% filter(site==52)
site_cohort <- wfi_cohorts %>% filter(site==52)
site52_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1680, xmax=1760, ymin=-Inf, ymax=Inf), color="orangered", linewidth=0.025, fill="orangered", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.7, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=2.25, stroke=1.25) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.75, shape = 23, color="darkred", fill="darkred") +
  annotate(geom="text", x = Inf, y = Inf, label = paste("52"), vjust = 1.55, hjust = 7.1, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1680, 1700, 1720, 1740, 1760), limits=c(1680, 1760), expand = c(0.05, 0.05)) +
  scale_y_discrete("", breaks = NULL, expand = c(0.1, 0.025)) +
  theme_bw(base_size=19, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site52_chart.png", width = 2.28, height = 2.36, units = "in", dpi=500)

#Graph individual site for cohort demonstration graphic (site 96)
site_xsections <- target_xsections %>% filter(site==96)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1855, 1855, site_xsections$outer_ring)
#site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1495, 1495, site_xsections$inner_ring)
site_events <- target_events %>% filter(site==96)
site_box <- clackbox %>% filter(site==96)
site_cohort <- wfi_cohorts %>% filter(site==96)
site96_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1485, xmax=1855, ymin=-Inf, ymax=Inf), color="goldenrod1", linewidth=0.025, fill="goldenrod1", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=1.7, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=2.15, stroke=1.25) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.25, shape = 23, color="darkred", fill="darkred") +
  geom_text(data=site_xsections, aes(x = estab_year, y= sample_id_match, label=estab_year), size=2.5, color="black", hjust = 1.32) +
  #annotate(geom="text", x = Inf, y = Inf, label = paste("96"), vjust = 1.55, hjust = 10, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850), limits=c(1485, 1855)) +
  scale_y_discrete("", breaks = NULL) +
  theme_bw(base_size=18, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site96_chart.png", width = 3, height = 3.5, units = "in", dpi=800)

#Graph individual site for cohort demonstration graphic (site 82)
site_xsections <- target_xsections %>% filter(site==82)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1830, 1830, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1605, 1605, site_xsections$inner_ring)
site_events <- target_events %>% filter(site==82)
site_box <- clackbox %>% filter(site==82)
site_cohort <- wfi_cohorts %>% filter(site==82)
site82_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1605, xmax=1830, ymin=-Inf, ymax=Inf), color="goldenrod1", linewidth=0.025, fill="goldenrod1", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=2.45, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=2.15, stroke=1.25) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.25, shape = 23, color="darkred", fill="darkred") +
  geom_text(data=site_xsections, aes(x = estab_year, y= sample_id_match, label=estab_year), size=2.5, color="black", hjust = 1.32) +
  #annotate(geom="text", x = Inf, y = Inf, label = paste("82"), vjust = 1.55, hjust = 10, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1650, 1700, 1750, 1800), limits=c(1605, 1830)) +
  scale_y_discrete("", breaks = NULL) +
  theme_bw(base_size=18, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site82_chart.png", width = 3, height = 3.5, units = "in", dpi=800)

#Graph individual site for cohort demonstration graphic (site 26)
site_xsections <- target_xsections %>% filter(site==26)
site_xsections$outer_ring <- ifelse(site_xsections$outer_ring > 1600, 1600, site_xsections$outer_ring)
site_xsections$inner_ring <- ifelse(site_xsections$inner_ring < 1475, 1475, site_xsections$inner_ring)
site_events <- target_events %>% filter(site==26)
site_box <- midbox %>% filter(site==26)
site_cohort <- wfi_cohorts %>% filter(site==26)
site26_chart <- ggplot() +
  geom_rect(data=site_box, aes(xmin=1475, xmax=1600, ymin=-Inf, ymax=Inf), color="violetred1", linewidth=0.025, fill="violetred1", alpha=.1) +
  geom_rect(data=site_cohort, aes(xmin=min_cohort_year, xmax=max_cohort_year, ymin=-Inf, ymax=Inf), fill="grey20", alpha=.5) +
  geom_segment(data=site_xsections, aes(x = outer_ring, y = sample_id_match, xend = inner_ring, yend = sample_id_match, color=species), linewidth=2.45, alpha=1) +
  scale_color_manual(name = "", values=c('PSME' = 'darkolivegreen4', 'TSHE' = 'dodgerblue4', 'TSME' = 'slateblue4', 'ABPR' = 'lightsalmon4', 'Other'='aquamarine4'), labels=c("PSME (n=571)", "TSHE (n=30)", "TSME (n=19)", "ABPR (n=13)", "Other (n=17)")) +
  geom_point(data=site_xsections, aes(x = estab_year, y= sample_id_match, shape=cohort, fill=cohort, alpha=cohort), color="black", size=2.15, stroke=1.25) +
  geom_point(data=site_events, aes(x = event_year, y = sample_id_match), size=2.25, shape = 23, color="darkred", fill="darkred") +
  geom_text(data=site_xsections, aes(x = estab_year, y= sample_id_match, label=estab_year), size=2.5, color="black", hjust = 1.32) +
  #annotate(geom="text", x = Inf, y = Inf, label = paste("26"), vjust = 1.55, hjust = 10, size=5.25) +
  scale_fill_manual(values = c('True' = 'black', 'False' = 'grey35')) +
  scale_shape_manual(values = c('True' = 22, 'False' = 0)) +
  scale_alpha_manual(values = c('True' = 1, 'False' = .55)) +
  scale_x_continuous("", breaks = c(1500, 1550, 1600), limits=c(1475, 1600)) +
  scale_y_discrete("", breaks = NULL) +
  theme_bw(base_size=18, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, margin=margin(0,0,0,0))) +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", linewidth=.35), panel.grid.minor.x = element_blank(), panel.grid.major = element_line(colour="grey50", linewidth=0.05)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "site26_chart.png", width = 3, height = 3.5, units = "in", dpi=800)


##########################################
############ Summarize data, #############
##### calculate fire history metrics #####
########## and bootstrap MFRIs ###########
##########################################

#Clear workspace and set directory
rm(list=ls())
setwd("/YOURPATH")

#Libraries for this section
library(tidyverse)
library(officer) 

#Read cross-sections file, select pertinent columns, and drop series with no first or last year
xsections <- read.csv("xsections.csv")
xsections <- xsections[,c("site", "sample_id", "species", "sample_ht_cm", "d10_mm", "pith","outer_ring", "inner_ring", "latitude", "longitude")]
xsections <- xsections %>% drop_na(outer_ring) 
xsections <- xsections %>% drop_na(inner_ring) 
nrow(xsections)

#Make all samples 6 digits (to remove, A, B, C, D, etc. because we are combining histories from the same tree), and convert numeric columns to numeric
xsections$sample_id <- substr(xsections$sample_id, 1, 6)
xsections$pith <- as.numeric(xsections$pith)
xsections$outer_ring <- as.numeric(xsections$outer_ring)
xsections$inner_ring <- as.numeric(xsections$inner_ring)
xsections$sample_ht_cm <- as.numeric(xsections$sample_ht_cm)
xsections$d10_mm <- as.numeric(xsections$d10_mm)
nrow(xsections)
length(unique(xsections$sample_id))

#Collapse data from A, B, C, D, etc. of same sample and make result a dataframe (for some reason the tibble output seems to choke the establishment data model).
xsections_collapse <- xsections %>% 
  group_by(site, sample_id, species) %>% 
  summarize(
    latitude=collapse::fmin(latitude, na.rm = TRUE),
    longitude=collapse::fmin(longitude, na.rm = TRUE),
    sample_ht_cm=collapse::fmin(sample_ht_cm, na.rm = TRUE),
    d10_mm=collapse::fmin(d10_mm, na.rm = TRUE),
    outer_ring=collapse::fmax(outer_ring, na.rm = TRUE),
    inner_ring=collapse::fmin(inner_ring, na.rm = TRUE),
    pith=collapse::fmin(pith, na.rm = TRUE))
xsections_collapse <- data.frame(xsections_collapse)
str(xsections_collapse)
target_xsections <- xsections_collapse 
unique(target_xsections$site)
nrow(target_xsections)

#Read fires file, remove non-fire injuries, and select pertinent columns. 
events <- read.csv("events.csv")
events <- events %>% filter(type=="FS") 
events <- events[,c("site", "sample_id", "event_year")]
target_events <- events 

#Summary of total number of fires, min fire year, max fire year, total range of years between fires by site
fire_count <- target_events %>%
  group_by(site) %>%
  summarise(no_fires = n_distinct(event_year),
            first_fire=min(event_year),
            last_fire=max(event_year),
            recon_range=last_fire-first_fire) 
print(fire_count[order(fire_count$no_fires),], n=40)
mean(fire_count$no_fires)
summary(fire_count$recon_range)

#Get all the years between the first and last year range for each sample
xsections_range <- xsections_collapse %>%
  transmute(site, sample_id, year = map2(inner_ring, outer_ring, `:`)) %>%
  unnest(cols = c(year))
data.frame(xsections_range)

#Min and max year by site
xsections_years <- xsections_range %>% 
  group_by(site) %>% 
  summarize(min_year=min(year),
            max_year=max(year))
xsections_years
summary(xsections_years$min_year)

#Combine with the fire counts
sites_fires <- xsections_years %>% left_join(fire_count)
sites_fires

#Fires that burned multiple sites
#Unique fire years
fire_years_site <- target_events %>%
  group_by(site) %>%
  distinct(event_year) %>%
  arrange(site)
fire_years_site
multi_site_fires <- fire_years_site %>%
  group_by(event_year) %>%
  tally() %>% 
  arrange(n)
print(multi_site_fires %>% filter(n>1), n=50)
print(multi_site_fires %>% filter(n>2 & n<=3), n=50)

#Summarize forest types by site
ilap <- read.csv("ilap.csv")
head(ilap)
ilap_summary <- ilap %>%
  group_by(forest_type_new) %>%
  tally()
ilap_summary

#Create a file indicating every year that wood samples are potentially recording fire (all cross dated tree rings of all sample), plus how many samples were recording fire in that year ("sample_depth").
sample_depth_by_year <- xsections_range %>%
  distinct(site, sample_id, year) %>%  
  count(site, year, name = "sample_depth")  
sample_depth_by_year

#Now add fire years as a binary outcome to that file
#Drop the 'century' column from fire_years_site and add fire year flag
fire_years_clean <- fire_years_site %>%
  select(site, year = event_year) %>%
  mutate(fire_occur = 1)
#Join with sample_depth_by_year and fill fire_occur as 0 if NA
fire_history <- sample_depth_by_year %>%
  left_join(fire_years_clean, by = c("site", "year")) %>%
  mutate(fire_occur = ifelse(is.na(fire_occur), 0, fire_occur))
fire_history

#Build a table of key fire history metrics for each site... the cv_interval column is the coefficient of variation (CV) of fire return intervals for each site.  It measures how variable the time between fires is at a site, relative to the average time between fires.  It's standard deviation of intervals divided by mean of intervals.  A high CV (> 1) indicates intervals that are highly variable (some very short, some very long).  A low CV (< 1) indicates intervals that are relatively consistent.  NA indicates a site that had only one fire (so no intervals).
#Calculate recording length for each site and truncate our calculations to the first cross dated year of our wood sample records (1171) to the year 1900.
#Limit fire history data to years 1171–1900
fire_history_trimmed <- fire_history %>%
  filter(year >= 1171 & year <= 1900)
#Calculate recording length for all sites
record_length <- fire_history_trimmed %>%
  group_by(site) %>%
  summarise(record_length = n_distinct(year), .groups = "drop")
#Create fire history summary table
fire_summary_with_fires <- fire_history_trimmed %>%
  filter(fire_occur == 1) %>%
  group_by(site) %>%
  summarise(
    n_fires = n(),
    fire_years = list(sort(year)),
    .groups = "drop"
  ) %>%
  mutate(
    intervals = lapply(fire_years, function(x) diff(x)),
    min_interval = sapply(intervals, function(x) if (length(x) == 0) NA else min(x)),
    max_interval = sapply(intervals, function(x) if (length(x) == 0) NA else max(x)),
    mfri         = sapply(intervals, function(x) if (length(x) == 0) NA else mean(x)),
    sd_interval  = sapply(intervals, function(x) if (length(x) == 0) NA else sd(x)),
    cv_interval  = round(sd_interval / mfri, 2),
    first_fire   = sapply(fire_years, min),
    last_fire    = sapply(fire_years, max)
  ) %>%
  select(site, n_fires, first_fire, last_fire,
         min_interval, max_interval, mfri, cv_interval)
fire_history_summary <- record_length %>%
  left_join(fire_summary_with_fires, by = "site")

#Read the environmental data and join to ilap (to get forest series in our file)
env_data <- read.csv("env_data.csv")
env_data <- env_data %>%
  left_join(ilap %>% select(site, forest_type_new), by = "site")

#Join site-level elevation, aspect, and series to the fire history summary
env_summary <- env_data %>%
  group_by(site) %>%
  summarise(
    elevation_m = mean(elevation_m, na.rm = TRUE),
    aspect = mean(aspect, na.rm = TRUE),
    slope_percent = mean(slope_percent, na.rm = TRUE),
    series = names(which.max(table(forest_type_new))),
    .groups = "drop"
  )
fire_history_summary <- fire_history_summary %>%
  left_join(env_summary, by = "site")

#Add site size
#Read cross-sections file, select pertinent columns, and drop series with no first or last year
xsections <- read.csv("xsections.csv")
xsections <- xsections[,c("site", "sample_id", "species", "pith","outer_ring", "inner_ring", "latitude", "longitude")]
xsections <- xsections %>% drop_na(outer_ring) 
xsections <- xsections %>% drop_na(inner_ring) 

#Make all samples 6 digits (to remove, A, B, C, D, etc.) and convert to numeric
xsections$sample_id <- substr(xsections$sample_id, 1, 6)
xsections$pith <- as.numeric(xsections$pith)
xsections$outer_ring <- as.numeric(xsections$outer_ring)
xsections$inner_ring <- as.numeric(xsections$inner_ring)

#Collapse data from A, B, C, D, etc. of same sample
xsections_collapse <- xsections %>% 
  group_by(sample_id) %>% 
  summarize(outer_ring=collapse::fmax(
    outer_ring, na.rm = TRUE),
    inner_ring=collapse::fmin(inner_ring, na.rm = TRUE),
    pith=collapse::fmin(pith, na.rm = TRUE))
print(xsections_collapse, n=10)
head(xsections_collapse)

#Add back the other xsections data
addback <- xsections %>% select("site", "sample_id", "species", "latitude", "longitude")
xsections_final <- left_join(xsections_collapse, addback)
head(xsections_final)

#Create sf object from point data
xsections_sf <- xsections_final %>%
  distinct(site, latitude, longitude) %>%  
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 32610)  

# Group by site and compute minimum bounding geometry (convex hull)
site_polygons <- xsections_sf %>%
  group_by(site) %>%
  summarise(geometry = st_combine(geometry)) %>%
  mutate(geometry = st_convex_hull(geometry)) %>%
  st_as_sf()

#Compute area in hectares
site_polygons <- site_polygons %>%
  mutate(area_ha = st_area(geometry) / 10^4) %>%  # Convert from m² to ha
  st_drop_geometry()
data.frame(site_polygons)

#Join to main summary
fire_history_summary <- fire_history_summary %>%
  left_join(site_polygons, by = "site")

#Function to convert aspect in degrees to cardinal direction
aspect_to_cardinal <- function(degrees) {
  cuts <- c(0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 360)
  labels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")  # 'N' appears twice for 0–22.5 and 337.5–360
  return(cut(degrees %% 360, breaks = cuts, labels = labels, include.lowest = TRUE, right = FALSE))
}

#Apply convert aspect function to the summary table
fire_history_summary <- fire_history_summary %>%
  mutate(aspect_cardinal = aspect_to_cardinal(aspect))

#We want to add the site-level estimates of annual probability of fire and CIs to this table.  We get that from the spaMM model described below.  We haven't gotten there yet, so I'm going to provide that data here:
site_predicts_trunc <- structure(list(
  site = c(22L, 23L, 24L, 25L, 26L, 28L, 52L, 53L, 54L, 55L, 56L, 57L, 58L, 59L, 60L, 61L, 62L, 64L, 68L, 71L, 72L, 81L, 82L, 83L, 84L, 85L, 86L, 87L, 88L, 89L, 90L, 91L, 92L, 93L, 96L, 98L),
  prob_fire = c(0.010590165, 0.012802926, 0.007344811, 0.012278204, 0.007946499, 0.021164820, 0.013765197, 0.013658038, 0.016855144, 0.010872177, 0.012486599, 0.009658423, 0.012595092, 0.006695735, 0.013332486, 0.011475327, 0.006213982, 0.005700798, 0.005393982, 0.023100976, 0.009054484, 0.012181792, 0.038519034, 0.020176518, 0.009263889, 0.002569201, 0.015134791, 0.010965921, 0.006676390, 0.007478782, 0.004115546, 0.013150317, 0.004898690, 0.009357931, 0.017780992, 0.018988983),
  fixefVar_0.025 = c(0.0065690676, 0.0093878854, 0.0050838122, 0.0079427925, 0.0044725704, 0.0128714157, 0.0105566231, 0.0077042173, 0.0116581049, 0.0070727282, 0.0085025828, 0.0060983695, 0.0090582262, 0.0040782785, 0.0067539023, 0.0080580310, 0.0041533300, 0.0037356729, 0.0032779783, 0.0157422071, 0.0054032217, 0.0089179306, 0.0275772992, 0.0134908032, 0.0050183543, 0.0004669952, 0.0109528646, 0.0071295448, 0.0045980657, 0.0047032821, 0.0024666714, 0.0082449493, 0.0031212413, 0.0056991920, 0.0114114246, 0.0112607533),
  fixefVar_0.975 = c(0.017030485, 0.017438391, 0.010600662, 0.018934846, 0.014080530, 0.034614509, 0.017931309, 0.024101220, 0.024311970, 0.016678394, 0.018302927, 0.015264814, 0.017488586, 0.010974574, 0.026150182, 0.016318009, 0.009287479, 0.008690646, 0.008863763, 0.033781574, 0.015135563, 0.016620157, 0.053562951, 0.030074503, 0.017039642, 0.014002038, 0.020879716, 0.016831637, 0.009684977, 0.011872620, 0.006859051, 0.020912613, 0.007680542, 0.015329276, 0.027606607, 0.031850225)
), class = "data.frame", row.names = c(NA, -36L))
site_predicts_trunc

#Join annualized fire probabilities to the predictions by site
fire_history_summary <- fire_history_summary %>%
  left_join(
    site_predicts_trunc %>%
      select(site, prob_fire, fixefVar_0.025, fixefVar_0.975),
    by = "site"
  )

#Calculate model-based mean fire return interval as the reciprocal of the annual probability of fire
fire_history_summary <- fire_history_summary %>%
  mutate(
    mfri_est = 1 / prob_fire,
    mfri_lower = 1 / fixefVar_0.975,  # Lower bound --> upper prob
    mfri_upper = 1 / fixefVar_0.025   # Upper bound --> lower prob
  )

#View the table
print(fire_history_summary, n=50, width=Inf)
summary(fire_history_summary$n_fires)
summary(fire_history_summary$record_length)
summary(fire_history_summary$min_interval)
summary(fire_history_summary$max_interval)
summary(fire_history_summary$first_fire)
summary(fire_history_summary$cv_interval)
summary(fire_history_summary$mfri)
summary(fire_history_summary$bootstrap_mfri)
summary(fire_history_summary$mfri_est)
summary(fire_history_summary$mfri_lower)
summary(fire_history_summary$mfri_upper)

#Look at some subset of sites
#print(filter(fire_history_summary, cv_interval <= .3), width=Inf)

#Create a final table and export as Word document

#Order sites to match other graphics (reverse it below)
site_levels <- rev(c("28", "25", "23", "22", "24", "26", "98", "96", "82", "81", "84", "86", "83", "91", "87", "88", "93", "92", "89", "90", "85", "72", "71", "68", "60", "61", "52", "64", "59", "53", "62", "57", "56", "54", "58", "55"))

#Clean and format summary table
fire_table_export <- fire_history_summary %>%
  mutate(site = factor(site, levels = site_levels)) %>%
  arrange(site) %>%
  transmute(
    `Site` = as.character(site),
    `Series` = recode(as.character(series),
                      "ABGR" = "A. grandis",
                      "PSME" = "P. menziesii",
                      "TSHE" = "T. heterophylla",
                      "ABAM" = "A. amabilis",
                      "TSME" = "T. mertensiana",
                      .default = as.character(series)),
    `Site area (ha)` = round(as.numeric(area_ha), 2),
    `Elevation (m)` = round(elevation_m),
    `Slope (%)` = round(slope_percent),
    `Aspect` = as.character(aspect_cardinal),
    `Years recording` = record_length,
    `No. of fires` = n_fires,
    `First fire` = as.integer(first_fire),
    `Last fire` = as.integer(last_fire),
    `Min. fire interval` = min_interval,
    `Max fire interval` = max_interval,
    `MFRI` = ifelse(is.na(mfri), NA, round(mfri, 1)),
    `MFRI Coeff. Var.` = ifelse(is.na(cv_interval), NA, round(cv_interval, 2)),
    `Annual prob. fire` = prob_fire,
    `Model-based MFRI` = mfri_est,
    `Model-based MFRI CIs` = ifelse(is.na(mfri_lower) | is.na(mfri_upper),
         NA,
        paste0(round(mfri_lower, 1), "–", round(mfri_upper, 1)))
  )

#Identify numeric columns
num_cols <- names(fire_table_export)[sapply(fire_table_export, is.numeric)]

#Create summary rows before rounding or coercion
mean_row <- fire_table_export %>%
  summarise(across(all_of(num_cols), ~ round(mean(.x, na.rm = TRUE), 3))) %>%
  mutate(`Site` = "Mean")

range_row <- fire_table_export %>%
  summarise(across(all_of(num_cols), ~ {
    rng <- range(.x, na.rm = TRUE)
    paste0(round(rng[1], 3), " – ", round(rng[2], 3))
  })) %>%
  mutate(`Site` = "Range")

#Fill non-numeric columns in summary rows with NA
non_num_cols <- setdiff(names(fire_table_export), num_cols)
for (col in non_num_cols) {
  mean_row[[col]] <- NA
  range_row[[col]] <- NA
}

#Reorder and coerce all to character
mean_row <- mean_row[, names(fire_table_export)] %>%
  mutate(across(everything(), as.character))

range_row <- range_row[, names(fire_table_export)] %>%
  mutate(across(everything(), as.character))

fire_table_export_chr <- fire_table_export %>%
  mutate(
    `Annual prob. fire` = round(`Annual prob. fire`, 3),
    `Model-based MFRI` = round(`Model-based MFRI`, 1),
    across(everything(), as.character)
  )

#Combine and replace missing with em-dash
fire_table_export_final <- bind_rows(
  fire_table_export_chr,
  mean_row,
  range_row
) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x) | .x == "", "–", .x)))

#Create formatted flextable... note that we are requiring flextable redundantly throughout to avoid conflicts with some packages loaded earlier.
ft <- flextable::flextable(fire_table_export_final) %>%
  flextable::autofit() %>%
  flextable::fontsize(size = 6, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::padding(padding = 3, part = "all") %>%
  flextable::italic(j = "Series", italic = TRUE, part = "body") %>%  # Italicize Series
  flextable::width(j = "Series", width = 1.75) %>%                    # Expand Series column slightly
  flextable::set_table_properties(layout = "autofit", width = .95) %>%
  flextable::hline(i = nrow(fire_table_export), border = fp_border(color = "black", width = 1))  # Line above Mean row

#Define landscape page layout
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape", width = 11, height = 8.5),
  type = "continuous"
)

#Export as Word to Desktop
desktop_path <- file.path(Sys.getenv("HOME"), "Desktop", "fire_summary_table.docx")
doc <- read_docx() %>%
  body_set_default_section(sect_properties) %>%
  flextable::body_add_flextable(ft)
print(doc, target = desktop_path)


#####################################
#### Create fire occurrence file ####
#####################################

#The goal of this section of code is to take data about cross-sections and data about cross-dated fire years and combine with environmental data to create a long-format data frame of all years from our fire history records with fire occurrence represented as a binomial response (0 or 1) suitable for modeling.  Running this code will result in creation of a file called fire_env which will be written to your home directory and will be loaded again in future sections.  

#Clear workspace and set directory
rm(list = ls())
setwd("/YOURPATH")

#Load libraries required for this section
library(tidyverse)

#Read cross-sections file, select pertinent columns, and drop series with no first or last year
xsections <- read.csv("xsections.csv")
xsections <- xsections[,c("site", "sample_id", "species", "sample_ht_cm", "d10_mm", "pith","outer_ring", "inner_ring", "latitude", "longitude")]
xsections <- xsections %>% drop_na(outer_ring) 
xsections <- xsections %>% drop_na(inner_ring) 

#Read events file, select pertinent columns, remove non-fire injuries, and remove unused site levels
events <- read.csv("events.csv")
events <- events %>% filter(type=="FS") 
events <- events[,c("site", "sample_id", "event_year")]
events <- droplevels(events)

#Read environmental data and select pertinent columns
env_data <- read.csv("env_data.csv")
env_data <- env_data %>% select(site, latitude:dist_accum)
env_data <- droplevels(env_data)

#Which sites are not present in the other files? It looks like we have everything.
sort(unique(xsections$site))
sort(unique(events$site))
sort(unique(env_data$site))
setdiff(env_data$site, xsections$site)
setdiff(events$site, xsections$site)
setdiff(xsections$site, events$site)
setdiff(xsections$site, env_data$site)
setdiff(xsections$site, events$site)
setdiff(events$site, env_data$site)

#Get means of environmental variables by site (mean of all samples per site) 
env_site <- env_data %>% 
  group_by(site) %>%
  summarise(across(everything(), mean))
env_site

#Get second earliest pith date (surrogate for stand age)
xsections_second_pith <- xsections %>% 
  group_by(site) %>%
  summarise(pith2 = nth(pith, 2, order_by = pith))
xsections_second_pith

#Create dataframe that consists of all years of all samples
xsections_years <- xsections %>%
  transmute(site, sample_id, year = map2(inner_ring, outer_ring, `:`)) %>%
  unnest(cols = c(year))
xsections_years
nrow(xsections_years)

#Add samples per year (i.e., sample depth)
xsections_samp_depth <- xsections_years %>%
  group_by(site, year) %>%
  tally() %>%
  rename("sample_depth"="n") 
xsections_samp_depth 
nrow(xsections_samp_depth)

#Create stand age column as the time in years between the first ring and the last ring
xsections_samp_depth$stand_age1 <- ave(xsections_samp_depth$year, xsections_samp_depth$site, FUN = seq_along)

#Join with second pith file
xsections_samp_depth <- xsections_samp_depth %>%
  left_join(xsections_second_pith)

#Create another stand age column as the time since second establishment date
xsections_samp_depth$pith2 <- as.numeric(xsections_samp_depth$pith2)
xsections_samp_depth$stand_age2 <- xsections_samp_depth$year - xsections_samp_depth$pith2
xsections_samp_depth$stand_age2[xsections_samp_depth$stand_age2 <= 0] <- NA
xsections_samp_depth$pith2 <- NULL

#Look at individual site... note that it's correctly parsed gaps in sample depth (for instance, site 84 has a hanging chronology, i.e., the chronology starts in 1202 and ends in 1587 and starts again in 1832)
testsite <- xsections_samp_depth %>% filter(site==84)
testsite %>% tibble::as_tibble() %>% print(n=400)

#Create dataframe of fire years by site
fire_years_site <- events %>% 
  group_by(site) %>% 
  distinct(event_year) %>% 
  rename("year"="event_year") %>% 
  arrange(site, year) 
fire_years_site

#Create fire occurrence indicator in fire sites file ("1) and merge unique occurrences (all of the fire sites file with no matches indicated by NA) with sample depth file.
fire_years_site$fire_occur = 1
fire_samp_depth <- merge(xsections_samp_depth, unique(fire_years_site), all.x = TRUE)
fire_samp_depth["fire_occur"][is.na(fire_samp_depth["fire_occur"])] <- 0
head(fire_samp_depth)
nrow(fire_samp_depth)

#Create time since fire variable.... this is super clunky code that I should probably improve.  
#First start a running year count by site that is reset every time fire_occur=1 (indicating a fire occurred).  Call this tsf_del.  
fire_samp_depth <- fire_samp_depth %>%
  group_by(site) %>%
  mutate(tsf_del = accumulate(fire_occur >= 1, ~ifelse(.y==TRUE, 0, .x+1)))
#Then convert 0 to NA
fire_samp_depth$tsf_del[fire_samp_depth$tsf_del==0] <- NA
#Then make the first NA enountered per site a 0
fire_samp_depth <- fire_samp_depth %>%
  group_by(site) %>%
  mutate_at(c("tsf_del"), ~replace(., row_number() == 1 & is.na(.), 0))
#Then create a column that makes any other NAs encountered (the actual year of fire) in the tsf_del column a number equivalent to the year immediately prior +1.  Call this column tsf.
fire_samp_depth <- fire_samp_depth %>%
  group_by(site)  %>%
  mutate(across(tsf_del, ~ zoo::na.locf(.x)))
fire_samp_depth$tsf <- ifelse(fire_samp_depth$fire_occur==1, fire_samp_depth$tsf_del+1, fire_samp_depth$tsf_del)
#And delete the old column
fire_samp_depth$tsf_del <- NULL

#Summarize number of years potentially recording fire per site
site_years <- fire_samp_depth %>%
  group_by(site) %>%
  count() %>%
  rename(site_rec_years = n)
print(site_years, n=50)

#Merge total site recording years 
fire_samp_depth <- fire_samp_depth %>% left_join(site_years)

#Maximum number of samples per site
no_samples <- fire_samp_depth %>%
  group_by(site) %>%
  summarize(samples_site = max(sample_depth)) %>%
  arrange(samples_site)
no_samples

#Merge maximum number of samples 
fire_samp_depth <- fire_samp_depth %>% left_join(no_samples)

#Read PDSI file, merge with PDSI file, and create one year lag PDSI column
pdsi <- read.csv("pdsi.csv")
pdsi <- pdsi[,c("YEAR", "RECON", "X20LP")]
pdsi <- pdsi %>%  
  rename("year"="YEAR",
         "gp33_pdsi"="RECON",
         "gp33_pdsi_lp"="X20LP")
fire_samp_depth <- fire_samp_depth %>% left_join(pdsi)
fire_samp_depth <- fire_samp_depth %>%                
  group_by(site) %>%
  mutate(lag_pdsi = lag(gp33_pdsi, n = 1, default = NA))

#Merge with environmental variables file
env_site$site <- as.factor(as.character(env_site$site))
fire_samp_depth$site <- as.factor(as.character(fire_samp_depth$site))
fire_env <- fire_samp_depth %>% left_join(env_site)
head(fire_env)
nrow(fire_env)
names(fire_env)

#Write .csv
write.csv(fire_env, "fire_env.csv", row.names = FALSE)


################################################
#### Create lightning ignition density file ####
################################################

#This block of code is how I created the ignition density file and so I'm including it here in case it's helpful.  But I am including the ignitions data that's created by this code, so you don't have to run this.  To run it you need the latest version of the FOD database, which is a large file that I'm not including (this block of code won't run without the FOD, which I'm not including).  You might get slightly different ignition densities using a more recent version of the FOD than I used.  

#Libraries neeeded for this section
library(sf)
library(tidyverse)

#Read site centroids data (has locations of samples for creating polygons)
site_centroids <- read.csv("site_centroids.csv")
head(site_centroids)

#Convert site centroids file into sf object
sites_sf <- site_centroids %>%
  select(site, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
head(sites_sf)

#Read FOD data... I'm not including this data here.  You'll have to download it yourself (see https://www.fs.usda.gov/rds/archive/catalog/RDS-2013-0009.6).  It is a large file.  
fod <- st_read("/YOURPATHFOD/FPA_FOD_20210617.gpkg", quiet = TRUE)
fod_lightning <- fod %>% filter(NWCG_CAUSE_CLASSIFICATION=="Natural")
unique(fod_lightning$NWCG_CAUSE_CLASSIFICATION)

#Buffer sites by 10km
sites_buff_10km <- st_buffer(sites_sf, 10000)

#Get counts of lightning ignitions within site buffers and calculate ignition density per km2
st_crs(fod_lightning) = st_crs(sites_buff_10km)
sites_buff_10km$ignition_count <- lengths(st_intersects(sites_buff_10km, fod_lightning))
sites_buff_10km$ignition_density_km2 <- sites_buff_10km$ignition_count / 314.16
ignitions <- sites_buff_10km %>% select(site, ignition_density_km2)

#Again, you might get different ignition densities using a newer or older version of the FOD than I used.  What I got was:
structure(list(
  site = c(22L, 23L, 24L, 25L, 26L, 28L, 52L, 53L, 54L, 55L, 56L, 57L, 58L, 59L, 60L, 61L, 62L, 64L, 68L, 71L, 72L, 81L, 82L, 83L, 84L, 85L, 86L, 87L, 88L, 89L, 90L, 91L, 92L, 93L, 96L, 98L), ignition_density_km2 = c(0.178253119, 0.155971480, 0.658899924, 0.079577285, 0.235548765, 0.216450216, 0.031830914, 0.149605297, 0.009549274, 0.019098549, 0.009549274, 0.120957474, 0.022281640, 0.280112045, 0.015915457, 0.133689840, 0.076394194, 0.098675834, 0.222816399, 0.254647313, 0.280112045, 
0.124140565, 0.410618793, 0.582505730, 0.184619302, 0.178253119, 0.222816399, 0.057295646, 0.149605297, 0.171886937, 0.200534759, 0.101858925, 0.187802394, 0.127323657, 0.143239114, 0.136872931)), class = "data.frame", row.names = c(NA, -36L))
#Drop geometry and (optiponal) write file
ignitions <- sf::st_drop_geometry(ignitions)
#write.csv(ignitions, "ignitions.csv", row.names = FALSE)

#Some quick and dirty analysis of this data:

#Count of unique fire years by site
events <- read.csv("events.csv")
events <- events %>% filter(type=="FS") 
site_counts <- events %>% 
  group_by(site) %>% 
  arrange(site, event_year, sample_id) %>% 
  summarise(count = n_distinct(event_year)) %>% 
  rename("no_fires"="count")
site_counts

#Join files
final_ignition_data <- sites_buff_10km %>% left_join(site_counts)
head(final_ignition_data)

#Read environmental data, summarize elevation, and join to ignition data
env_data <- read.csv("env_data.csv")
elev <- env_data %>% group_by(site) %>%
  summarize(elev=mean(elevation_m))
final_ignition_data <- final_ignition_data %>% left_join(elev)

#Plot relationship between total fire and lightning density
formula <- y ~ x
quartz(width=10, height=8)
fire_ignitions <- ggplot(final_ignition_data, aes(ignition_density_km2, no_fires, label=site)) +
  geom_smooth(method = "lm", formula = formula, linewidth=.75, alpha=.2) +
  ggpmisc::stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                        formula = formula, parse = TRUE, size = 4) +  geom_point(fill="darkred", color="black", shape=24, size=3, alpha=.75) +
  ggrepel::geom_text_repel(size=3, nudge_x = 0.01, nudge_y = 0.01) +
  scale_x_continuous(bquote(Lightning~ignitions~km^2)) +
  scale_y_continuous("Total reconstructed fires") +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position="none") +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) 
fire_ignitions
#ggsave(path="/YOURPATH", "fire_ignitions.png", width = 10, height = 8, units = "in", dpi=500)

#Plot relationship between lightning density and elevation
formula <- y ~ x
quartz(width=10, height=8)
elev_ignitions <- ggplot(final_ignition_data, aes(elev, ignition_density_km2, label=site)) +
  geom_smooth(method = "lm", formula = formula, linewidth=.75, alpha=.2) +
  ggpmisc::stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                        formula = formula, parse = TRUE, size = 4) +  geom_point(fill="darkred", color="black", shape=24, size=3, alpha=.75) +
  ggrepel::geom_text_repel(size=3, nudge_x = 0.01, nudge_y = 0.01) +
  scale_x_continuous("Elevation (m)") +
  scale_y_continuous(bquote(Lightning~ignitions~km^2)) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(legend.position="none") +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) 
elev_ignitions
#ggsave(path="/YOURPATH", "elev_ignitions.png", width = 10, height = 8, units = "in", dpi=500)

#This quick and dirty analysis suggests that there is a relationship between total fires and lightning ignition density, but it's a pretty weak relationship.  There appears to be stronger relationship between lightning ignition density and elevation (which itself is strongly related to snow disappearance day, see below).  So I'm not really excited about including lightning ignition density in the model, although I do that below as a post hoc exercise.  


#######################################
### spaMM model for fire occurrence ###
#######################################

#The goal of this section is to build a statistical model of the key influences on the binary fire occurrence response.  Some of this code will take a while to run, and I try and flag those parts below.  This section requires the fire_env file which is created above.  

rm(list=ls())
setwd("/YOURPATH")

#Load libraries required for this section
library(tidyverse)
library(spaMM)
library(lme4)
library(DHARMa)
library(sp)
library(sf)

### Organize the data for modeling ###

#Read site-level fire occurence data created in an earlier section
fire_env <- read.csv("fire_env.csv", header = TRUE) 
head(fire_env)

#Change latitude and longitude to UTMs by first converting to sf object, then transforming to a UTM CRS, which should provide more intuitive interpretation of variogram distances and perhaps a more logical spatial random effect term
fire_env_sf <- st_as_sf(x = fire_env, coords = c("longitude", "latitude"), crs = 4269)
fire_env_utm <- fire_env_sf %>% sf::st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")
fire_env <- fire_env_utm %>%
  dplyr::mutate(easting = sf::st_coordinates(.)[,1],
                northing = sf::st_coordinates(.)[,2])
#Confirm we've got the coordinates right
ggplot() + geom_sf(data = fire_env) + geom_text(data = fire_env, aes(x=easting, y=northing, label = site))
#Drop geometry
fire_env <- st_drop_geometry(fire_env)

#Summarize fire occurrence by site
unique(fire_env$site)
table(fire_env$site, fire_env$fire_occur)

#Filter data to years prior to the modern period and after the period for which we have decent sample depth (see text)
min_year <- 1400
max_year <- 1900
fire_env <- fire_env %>% filter(year>=min_year & year<=max_year)
head(fire_env)
nrow(fire_env)

#Correlation of different variables... we can't have highly correlated terms in the same model.  You could try and fit a model with many other terms that are included with this data (bonus environmental variables).  But sdd_average will prove to be a slightly better fit of other variables (such as the below) that are basically speaking to fuel drying potential of landscape setting.  
corr_data <- fire_env %>% select(elevation_m, sdd_average, vpd_max, temp_max_c, precip_mm)
PerformanceAnalytics::chart.Correlation(corr_data, histogram=TRUE, pch=19)

### Initial model and diagnostics ###

#The first step is to try and figure out which of strongly correlated variables that influence fuel moisture should be used in a final model.  We're going to fit these models without the spatial correlation terms, which is the same fit we could achieve with glmer.  Just as a reminder of how to interpret binomial regression outputs:  The estimate is the coefficient associated with the variable. It is the estimated amount by which the log odds of fire occurrence would increase if that variable (i.e., sdd_average) were one unit higher.  When a coefficient is negative, the higher the value of the variable, the lower the probability of fire.  

#First model: Elevation
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm1.1 <- fitme(fire_occur ~ 
                  gp33_pdsi + 
                    elevation_m + 
                  slope_percent +
                  curvature +
                  stand_age1 + 
                  tsf +
                  (1|site),
                data=fire_env, 
                family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)

#Second model: Snow disappearance day.  
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm1.2 <- fitme(fire_occur ~ 
                    gp33_pdsi + 
                    sdd_average + 
                    slope_percent +
                    curvature +
                    stand_age1 + 
                    tsf +
                    (1|site),
                  data=fire_env, 
                  family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)

#Third model: VPD
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm1.3 <- fitme(fire_occur ~ 
                    gp33_pdsi + 
                    vpd_max + 
                    slope_percent +
                    curvature +
                    stand_age1 + 
                    tsf +
                    (1|site),
                  data=fire_env, 
                  family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)

#Fourth model: Max temp
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm1.4 <- fitme(fire_occur ~ 
                    gp33_pdsi + 
                    temp_max_c + 
                    slope_percent +
                    curvature +
                    stand_age1 + 
                    tsf +
                    (1|site),
                  data=fire_env, 
                  family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)

#Fifth model: Precipitation
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm1.5 <- fitme(fire_occur ~ 
                    gp33_pdsi + 
                    precip_mm + 
                    slope_percent +
                    curvature +
                    stand_age1 + 
                    tsf +
                    (1|site),
                  data=fire_env, 
                  family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)

#Comparisons
summary(spamm1.1)
summary(spamm1.2)
summary(spamm1.3)
summary(spamm1.4)
summary(spamm1.5)
AIC(spamm1.1)
AIC(spamm1.2)
AIC(spamm1.3)
AIC(spamm1.4)
AIC(spamm1.5)

#The model with snow disappearance day is slightly better.  There's not actually strong correlation with precipitation, so we could try a model with snow disappearance day and precipitation... but this isn't better, so we are going to use spamm1.2, which is the model with just snow disappearance day.  
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm1.6 <- fitme(fire_occur ~ 
                    gp33_pdsi + 
                    sdd_average + 
                    precip_mm +
                    slope_percent +
                    curvature +
                    stand_age1 + 
                    tsf +
                    (1|site),
                  data=fire_env, 
                  family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)
summary(spamm1.6)
AIC(spamm1.6)
AIC(spamm1.2)

#Some diagnostics using the very helpful DHARMa package... There does not seem to be overdispersion or zero inflation... the residuals to predicted fits seem good.  There does seem to be some outliers, but we can safely ignore them for now (see note about outliers below)
testDispersion(spamm1.2)
testZeroInflation(spamm1.2)
simres <- simulateResiduals(spamm1.2, plot = TRUE)
simres
outl <- testOutliers(simres, type="bootstrap")
outl

#FYI:  Model residuals for spaMM models are stored here:
head(simres$scaledResiduals)

### Tests for temporal and spatial dependence ###

#First, we're going to test for temporal dependence 
#First, bind simulated residuals to data for testing
fortest <- fire_env %>% dplyr::select(site, year)
fortest$res <- simres$scaledResiduals
head(fortest)

#PACF for individual sites... first select a site and test for PACF and ACF... you would need to repeat this for every site
unique(fortest$site)
a_site <- fortest %>% filter(site==98)
pacf(a_site$res)
acf(a_site$res)

#Graph residuals over time... this doesn't seem to indicate significant temporal dependence.
resid_time <- ggplot(filter(fortest), aes(x = year, y = res)) +
  facet_wrap(vars(site)) +
  geom_line() +
  theme_bw()
resid_time

#Test for temporal autocorrelation by looping over sites  If many more than 0.05 * n_sites of the p values are below 0.05, this would indicate general autocorrelation. You could also plot some PACFs using pacf(res_i) if you wanted another look. Basically, we want to inspect the time series of residuals for each tree.  Inside the loop, you are trying to assign a p-value to the ith position in tests. This procedure does not indicate significant temporal dependence, so I think we can ignore temporal correlation moving forward
sites <- unique(fortest$site)
tests <- double(length(unique(fortest$site)))
for(i in 1:length(sites)){
  res_i <- subset(fortest, sites == sites[i])$res
  sum_lm <- summary(lm(res_i[-1] ~ res_i[-length(res_i)]))
  tests[i] <- sum_lm$coefficients[8]
}
sort(tests)
length(tests)
length(tests[tests < 0.05])
mean(tests < 0.05)
hist(tests, breaks=26)
abline(v = .05, col="red")

#Test for spatial dependence 
#Bind simulated residuals to data that includes northing and easting for testing
fortest_spat <- fire_env %>% dplyr::select(site, year, northing, easting)
fortest_spat$res <- simres$scaledResiduals
head(fortest_spat)

#Subset a year and convert that year file to a spatial points dataframe and plot variogram.  I'm not running this here because the variogram call relies on the gstat package which conflicts with other packages that I use throughout this code.  But you could start a new session, load library(gstat), and run this and it should work just fine.  
#Years that burned multiple sites:  1499, 1648, 1652, 1670, 1689, 1693, 1700, 1706, 1714, 1730, 1738, 1740, 1765, 1768, 1773, 1794, 1798, 1807, 1835, 1883, 1896, 1690, 1783, 1812, 1831, 1836, 1844, 1846, 1522 
#fortest_year <- fortest_spat %>% filter(year==1522)
#coordinates(fortest_year) <- ~northing+easting
#vario <- variogram(res~1, fortest_year)
#vario.fit <- autofitVariogram(res~1,
#                    fortest_year,
#                    model = c("Sph"),
#                    kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
#                    fix.values = c(NA, NA, NA),
#                   start_vals = c(NA,NA,NA),
#                    verbose = T)

#Another approach... this uses the automap library which also has conflicts with previously loaded packages.  But you could start a new session and load library(automap) and give this a try.  
#fortest_year <- fortest_spat %>% filter(year==1522)
#coordinates(fortest_year) <- ~northing+easting
#vario.fit <- autofitVariogram(res~1, fortest_year)
#plot(vario, vario.fit$var_model, main = "Fitted variogram")

#Also, graph residuals in space... hard to see any clear patterns
fortest_year <- fortest_spat %>% filter(year==1522)
ggplot(fortest_year, aes(x = easting, y = northing, size = res)) +
  geom_point() +
  scale_size_continuous(range = c(1,10))

#A formal test with Moran's I
fortest_year <- fortest_spat %>% filter(year==1522)
testSpatialAutocorrelation(fortest_year$res, x = fortest_year$easting, y = fortest_year$northing, plot = TRUE)

#Test for spatial autocorrelation by looping over sites with a Moran's I test.  
fortest_list <- fortest_spat %>%
  group_split(year)
moransi_out <- lapply(fortest_list, function(x) testSpatialAutocorrelation(x$res, x$easting, x$northing, plot = FALSE)) 
thepvalues <- lapply(moransi_out, function (x) x['p.value'])
sig_spatial_auto <- unlist(lapply(thepvalues, function(x) length(which(x<=0.05))))
sum(sig_spatial_auto > 0) / length(sig_spatial_auto)
#So about 5% of the years exhibit spatial autocorrelation, which probably justifies a spatial model, which we will fit below

### Final models and diagnostics ###

#Now we'll fit the model with snow disapperance day as above with a spatially correlated random effect term (easting and northing).  This approach is justified by the spatial correlation observed above, plus, this is a better model than our best model to date as measured by AIC.  
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm2.1 <- fitme(fire_occur ~ 
                  gp33_pdsi + 
                  sdd_average + 
                  slope_percent +
                  curvature +
                  stand_age1 + 
                  tsf +
                  (1|site) +
                  Matern(1|easting+northing), 
                data=fire_env, 
                family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)
#The coefficient estimate in the output indicate the average change in the log odds of the response variable associated with a one unit increase in each predictor variable.  Put another way:  the estimates are the estimated amount by which the log odds of fire occurrence would increase if the fixed effect, i.e., sdd_average were one unit higher. The log odds of fire occurrence when fixed effects are set at O is the intercept. 
summary(spamm2.1)
AIC(spamm1.2)
AIC(spamm2.1)

#We have an interest in how fire occurence may vary by forest type, so we'll add forest type to the data
ilap <- read.csv("ilap.csv")
ilap <- ilap %>% select(site, forest_type_new)
fire_env <- fire_env %>% left_join(ilap) 
head(fire_env)

#Add a forest type term to the model above
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm2.2 <- fitme(fire_occur ~ 
                    gp33_pdsi + 
                    sdd_average + 
                    slope_percent +
                    curvature +
                    stand_age1 + 
                    tsf +
                    forest_type_new +
                    (1|site) +
                    Matern(1|easting+northing), 
                  data=fire_env, 
                  family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)
summary(spamm2.2)

#Now a post hoc model with lightning ignition data

#Read lightning data and join with fire-environmental variables data
ignitions <- read.csv("ignitions.csv")
head(ignitions)
fire_env <- fire_env %>% left_join(ignitions) 

#Model with lightning ignitions... I also tried some interactions , but none of them was a better model.  
avail_thr <- parallel::detectCores(logical=FALSE) - 1L 
spamm2.3 <- fitme(fire_occur ~ 
                    gp33_pdsi +
                    sdd_average +                    
                    slope_percent +
                    curvature +
                    stand_age1 + 
                    tsf +
                    forest_type_new +
                    ignition_density_km2 +
                    (1|site) +
                    Matern(1|easting+northing), 
                  data=fire_env, 
                  family=binomial(), control.HLfit=list(NbThreads=max(avail_thr, 1L), max.iter.mean=1000, spaMM.options(separation_max=500)))
options(scipen=999)
summary(spamm2.3)

#Compare AIC
AIC(spamm1.1)
AIC(spamm1.2)
AIC(spamm1.3)
AIC(spamm1.4)
AIC(spamm1.5)
AIC(spamm2.1)
AIC(spamm2.3)

#Evaluation of AIC provides little reason to include forest type or lightning ignition density in the model.  We do have an interest in whether fire occurrence varies as a function of forest type, so we will estimate effect of forest type in steps below.   

#To summarize where we're at:  We're using a spatial generalized linear mixed model with a binomial (logit) link to model fire_occur as a function of:
#Fixed effects: gp33_pdsi, sdd_average, slope_percent, curvature, stand_age1, tsf, (and forest_type_new in one iteration of the model)
#Random effects: (1|site) and Matern(1|easting + northing) for spatial autocorrelation
#This setup is the most appropriate way I can think of to model site-level fire occurrence with spatial dependence, and the model converges successfully.

#Some diagnostics... There doesn't seem to be overdisperson or zero inflation.  KS test p-value greater than .05 provides no evidence that residuals are not normally distributed.  There were some outlier issues with n=10.  From Dharma documentation:  "Because no data was simulated in the range of the observed value, we don't know "how strongly" these values deviate from the model expectation, so the term "outlier" should be used with a grain of salt. It is not a judgment about the magnitude of the residual deviation, but simply a dichotomous sign that we are outside the simulated range. Moreover, the number of outliers will decrease as we increase the number of simulations."  And wouldn't you know it, increasing the number of simulations resolves the outlier issue.  So apparently there are no meaningful outliers.  But running the testOutliers function with n=1,000 will take almost several minutes at least.
testDispersion(spamm2.1)
testZeroInflation(spamm2.1)
simres <- simulateResiduals(spamm2.1, plot = TRUE, n=1000)
outl <- testOutliers(simres, type="bootstrap", nBoot=1000)
outl

#There's interesting stuff in the spaMM summary not found in other R mixed model outputs.  There's the fixed effects (betas) which are the estimated regression parameters (slopes). Then the correlation parameter nu and rho which represent the strength and the speed of decay in the spatial effect respectively, which we can interpret as the actual spatial correlation effect by plotting the estimated correlation between locations against their distance (this takes a while to run with the full dataset... we could use just use one year).  With this data plotted, we can say something like:  The correlation of fire occurence is below .1 at distance greater than 20km, or something along those lines. 
#Spatial correlation of fire occurrence at sites
summary(spamm2.1)
sliced_obs <- fire_env %>% slice_sample(n=1000, replace = FALSE)
unique(sliced_obs$site)
dd <- dist(sliced_obs[,c("easting","northing")])
mm <- MaternCorr(dd, nu = 16.6666666667, rho = 0.0008748874)

#Coerce these "dist" objects into dataframes
dd_mat <- as.matrix(dd)
dd_df <- reshape2::melt(dd_mat, c("x", "y"), value.name = "dd")
mm_mat <- as.matrix(mm)
mm_df <- reshape2::melt(mm_mat, c("x", "y"), value.name = "mm")
mm_dd_df <- left_join(dd_df, mm_df)
mm_dd_df <- mm_dd_df[,3:4]
nrow(mm_dd_df)
mm_dd_df <- mm_dd_df %>% filter(dd!=0 & mm!=0)
nrow(mm_dd_df)
summary(mm_dd_df)
head(mm_dd_df)

#Plot these objects with ggplot... this will take several minutes.
dist_corr_plot <- ggplot(mm_dd_df, aes(x=dd/1000, y=mm)) +
  geom_point(shape=21, alpha=.75, color="black", fill="darkred", size=2) +
  scale_x_continuous("Distance between pairs of sites (km)", breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), limits=c(0,max(mm_dd_df$dd/1000))) +
  scale_y_continuous("Estimated correlation") +
  theme_bw(base_size = 10, base_family = "Arial")
dist_corr_plot
ggsave(path="/YOURPATH", "dist_corr.png", width = 6, height = 4, units = "in", dpi=500)

### Effects of random and fixed variables on fire probability ###

#We can compute confidence intervals for the estimates directly, but be warned, this procedure takes a really long time (an hour or more).  We're going to calculate CIs for the estimates below in a different step and I find this most helpful for evaluating significance of these variables.  Here we see that there's just a tiny bit of overlap with zero for PDSI CIs (these are log odds, so this overlap is pretty insignificant), no overlap with zero for snow disappearance day CIs, and no overlap with zero for time since fire CIs.  CIs for other fixed effects all overlap with zero quite a bit and it's hard to imagine there's a strongly predictive relationship between fire occurence and these variables.  As noted above, I tried a number of interactions and those didn't improve model performance.    
#Note that "forest_type_newABAM" is the intercept so this procedure doesn't calculate a CI for that forest type directly.
#summary(spamm2.1)
#pdsi_ci <- confint(spamm2.1,"gp33_pdsi") 
#sdd_ci <- confint(spamm2.1,"sdd_average")
#slope_ci <- confint(spamm2.1,"slope_percent") 
#curv_ci <- confint(spamm2.1,"curvature") 
#age_ci <- confint(spamm2.1,"stand_age1") 
#tsf_ci <- confint(spamm2.1,"tsf") 
#summary(spamm2.3)
#pdsi_ci <- confint(spamm2.3,"gp33_pdsi") 
#sdd_ci <- confint(spamm2.3,"sdd_average")
#slope_ci <- confint(spamm2.3,"slope_percent") 
#curv_ci <- confint(spamm2.3,"curvature") 
#age_ci <- confint(spamm2.3,"stand_age1") 
#tsf_ci <- confint(spamm2.3,"tsf") 
#lightning_ci <- confint(spamm2.3,"ignition_density_km2")
#abgr_ci <- confint(spamm2.3, "forest_type_newABGR") 
#tshe_ci <- confint(spamm2.3, "forest_type_newTSHE") 
#tsme_ci <- confint(spamm2.3, "forest_type_newTSME") 

#Prediction for sites
forprediction <- fire_env %>% select(northing, easting, gp33_pdsi, sdd_average, slope_percent, curvature, stand_age1, tsf, site)
forprediction_means <- forprediction %>%
  group_by(site) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  )
forprediction_means$prob_fire <- as.numeric(predict(spamm2.1, forprediction_means)) 

#Get 95% confidence intervals for site predictions
site_predicts <- cbind(forprediction_means, get_intervals(spamm2.1, newdata = forprediction_means, intervals = "fixefVar"))
site_predicts$site <- as.character(site_predicts$site)
summary(site_predicts$prob_fire)
summary(site_predicts$fixefVar_0.025)
summary(site_predicts$fixefVar_0.975)

#Reorder factors:
site_predicts <- site_predicts %>%
  mutate(site = forcats::fct_reorder(site, prob_fire)) 
site_predicts$site <- factor(site_predicts$site, levels=c("28", "25", "23", "22", "24", "26", "98", "96", "82", "81", "84", "86", "83", "91", "87", "88", "93", "92", "89", "90", "85", "72", "71", "68", "60", "61", "52", "64", "59", "53", "62", "57", "56", "54", "58", "55"))

#Read environmental variables file, select HUC8 names, reduce to most common HUC8 name per site, and join to site coordinates
env_vars <- read.csv("env_data.csv")
env_vars <- env_vars %>% select(site, huc8_name)
env_vars <- env_vars %>% group_by(site) %>% count(huc8_name) %>% slice(which.max(n))
env_vars$site <- as.factor(env_vars$site)
site_predicts <- site_predicts %>% left_join(env_vars) 
head(site_predicts)

#Graph by site
site_prob_cis <- ggplot(site_predicts, aes(prob_fire, site, color=huc8_name)) +
  geom_pointrange(aes(xmin = fixefVar_0.975, xmax = fixefVar_0.025), show.legend = FALSE) +
  geom_point(shape=22, fill="black", size=2.5, stroke=1.75) +
  #geom_point(aes(prob_fire_future, site), color="midnightblue", shape=8, size=.75) +
  scale_color_manual(name = "", values=c("Lower Columbia-Sandy" = "green2", "Clackamas"="orangered", "North Santiam"="cyan2", "South Santiam"="darkorchid2", "Mckenzie"="goldenrod1", "Middle Fork Willamette"="violetred1"), guide = "none") +  
  scale_x_continuous("Annual probability of fire", expand = c(0.004, .004)) +
  scale_y_discrete("") +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(path="/YOURPATH", "site_prob_cis.png", width = 2.4, height = 7, units = "in", dpi=500)

#Now, the same thing as above, but we're going to predict for the model that includes forest type (series).
forprediction_type <- fire_env %>% select(northing, easting, gp33_pdsi, sdd_average, slope_percent, curvature, stand_age1, tsf, forest_type_new, site)
forprediction_type_means <- forprediction_type %>%
  group_by(site) %>%
  summarise(
    across(where(is.numeric), mean),
    forest_type_new = first(forest_type_new),  
    .groups = "drop"
  )
forprediction_type_means$prob_fire <- as.numeric(predict(spamm2.1, forprediction_type_means)) 

#Get 95% confidence intervals for forest type predictions
type_predicts <- cbind(forprediction_type_means, get_intervals(spamm2.2, newdata = forprediction_type_means, intervals = "fixefVar"))

#Recode series name and reorder levels
type_predicts <- type_predicts %>% mutate(forest_type_new=recode(forest_type_new, 'ABGR'='Grand fir series (n=3)', 'TSHE'='W. hemlock series (n=20)', 'ABAM'='Silver fir series (n=8)', 'TSME'='Mtn. hemlock series (n=5)'))
type_predicts$forest_type_new <- factor(type_predicts$forest_type_new, levels=c('Grand fir series (n=3)', 'W. hemlock series (n=20)', 'Silver fir series (n=8)', 'Mtn. hemlock series (n=5)'))
#Graph forest type (series)
forest_type_prob_graph <- ggplot(data=type_predicts, aes(x=forest_type_new, y=prob_fire, fill=forest_type_new)) +
  geom_boxplot(color="black", width=0.5) +
  scale_fill_manual(name = "", values=c('Grand fir series (n=3)' = 'darkgreen', 'W. hemlock series (n=20)'='dodgerblue4', 'Silver fir series (n=8)' = 'deeppink3', 'Mtn. hemlock series (n=5)' = 'slateblue4')) +
  scale_x_discrete("", labels = scales::label_wrap(15)) +
  scale_y_continuous("Annual prob. fire") +
  theme_bw(base_size = 17, base_family = "Arial") +
  theme(legend.position = "none")
ggsave(path="/YOURPATH", "forest_type_prob_graph.png", width = 7, height = 3.15, units = "in", dpi=600)

#Compare difference between forest types
#Test assuming equal variances 
type_predicts %>%
  rstatix::anova_test(prob_fire ~ forest_type_new,
                      detailed = TRUE)

#A test that does not assume equal variances
(oneway <- stats::oneway.test(prob_fire ~ forest_type_new,
                              data = type_predicts,
                              var.equal = FALSE))
#Effect size
effectsize::eta_squared(oneway)
#Effect size rank epsilon squared
effectsize::rank_epsilon_squared(type_predicts$prob_fire ~ type_predicts$forest_type_new)

#Post hoc comparisons of two groups at a time, aka "pairwise comparisons" using Kruskall-Wallis test
type_predicts %>% rstatix::kruskal_test(prob_fire ~ forest_type_new)

#Post-hoc test for non-parametric data
phdt <- type_predicts %>% rstatix::dunn_test(prob_fire ~ forest_type_new)
data.frame(phdt)

#Now we're going to predict fixed effects while controlling for the spatial and site random effects
#Generate prediction grid for PDSI
summary(fire_env$gp33_pdsi)
pdsi_newdat <- expand.grid(
  easting = mean(fire_env$easting), 
  northing = mean(fire_env$northing),
  gp33_pdsi = seq(-3.95, 5, length.out = 90),
  sdd_average = mean(fire_env$sdd_average), 
  slope_percent = mean(fire_env$slope_percent),
  curvature = mean(fire_env$curvature),
  stand_age1 = mean(fire_env$stand_age1), 
  tsf = mean(fire_env$tsf), 
  site = unique(fire_env$site)
)

#Predict fire probability (marginal: fixed effects only)
pdsi_newdat$prob_fire <- as.numeric(predict(spamm2.1, pdsi_newdat, re.form = NA))

#Add 95% confidence intervals for fixed effects
pdsi_newdat <- cbind(
  pdsi_newdat,
  get_intervals(spamm2.1, newdata = pdsi_newdat, intervals = "fixefVar", re.form = NA)
)

#Collapse predictions to one line per gp33_pdsi value (i.e., average over sites)
pdsi_plot_data <- pdsi_newdat %>%
  group_by(gp33_pdsi) %>%
  summarise(
    prob_fire = mean(prob_fire, na.rm = TRUE),
    lower = mean(fixefVar_0.025, na.rm = TRUE),
    upper = mean(fixefVar_0.975, na.rm = TRUE),
    .groups = "drop"
  )

#Plot the smoothed prediction curve
pdsi_graph <- ggplot(pdsi_plot_data, aes(x = gp33_pdsi, y = prob_fire)) +
  geom_line(color = "darkred", linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey50") +
  scale_y_continuous("Annual prob. fire") +
  scale_x_continuous("Reconstructed PDSI (SDs)") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(plot.margin = unit(c(.5, 1, .5, 1), "lines"))
pdsi_graph

#Generate prediction grid for snow disappearance day
summary(fire_env$sdd_average)
sdd_newdat <- expand.grid(
  easting = mean(fire_env$easting), 
  northing = mean(fire_env$northing),
  gp33_pdsi = mean(fire_env$gp33_pdsi),
  sdd_average = seq(72.82, 232.04, length.out = 90), 
  slope_percent = mean(fire_env$slope_percent),
  curvature = mean(fire_env$curvature),
  stand_age1 = mean(fire_env$stand_age1), 
  tsf = mean(fire_env$tsf), 
  site = unique(fire_env$site)
)

#Predict fire probability (marginal: fixed effects only)
sdd_newdat$prob_fire <- as.numeric(predict(spamm2.1, sdd_newdat, re.form = NA))

#Add 95% confidence intervals for fixed effects
sdd_newdat <- cbind(
  sdd_newdat,
  get_intervals(spamm2.1, newdata = sdd_newdat, intervals = "fixefVar", re.form = NA)
)

#Collapse predictions to one line per sdd_average value (i.e., average over sites)
sdd_plot_data <- sdd_newdat %>%
  group_by(sdd_average) %>%
  summarise(
    prob_fire = mean(prob_fire, na.rm = TRUE),
    lower = mean(fixefVar_0.025, na.rm = TRUE),
    upper = mean(fixefVar_0.975, na.rm = TRUE),
    .groups = "drop"
  )

#Plot the smoothed prediction curve
sdd_graph <- ggplot(sdd_plot_data, aes(x = sdd_average, y = prob_fire)) +
  geom_line(color = "darkred", linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey50") +
  scale_y_continuous("Annual prob. fire") +
  scale_x_continuous("Snow disappearance day") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(plot.margin = unit(c(.5, 1, .5, 1), "lines"))
sdd_graph

#Generate prediction grid for slope
summary(fire_env$slope_percent)
slope_newdat <- expand.grid(
  easting = mean(fire_env$easting), 
  northing = mean(fire_env$northing),
  gp33_pdsi = mean(fire_env$gp33_pdsi),
  sdd_average = mean(fire_env$sdd_average),
  slope_percent = seq(0, 60.73, length.out = 90), 
  curvature = mean(fire_env$curvature),
  stand_age1 = mean(fire_env$stand_age1), 
  tsf = mean(fire_env$tsf), 
  site = unique(fire_env$site)
)

#Predict fire probability (marginal: fixed effects only)
slope_newdat$prob_fire <- as.numeric(predict(spamm2.1, slope_newdat, re.form = NA))

#Add 95% confidence intervals for fixed effects
slope_newdat <- cbind(
  slope_newdat,
  get_intervals(spamm2.1, newdata = slope_newdat, intervals = "fixefVar", re.form = NA)
)

#Collapse predictions to one line per slope_percent value (i.e., average over sites)
slope_plot_data <- slope_newdat %>%
  group_by(slope_percent) %>%
  summarise(
    prob_fire = mean(prob_fire, na.rm = TRUE),
    lower = mean(fixefVar_0.025, na.rm = TRUE),
    upper = mean(fixefVar_0.975, na.rm = TRUE),
    .groups = "drop"
  )

#Plot the smoothed prediction curve
slope_graph <- ggplot(slope_plot_data, aes(x = slope_percent, y = prob_fire)) +
  geom_line(color = "darkred", linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey50") +
  scale_y_continuous("Annual prob. fire") +
  scale_x_continuous("Slope percent") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(plot.margin = unit(c(.5, 1, .5, 1), "lines"))
slope_graph

#Generate prediction grid for curvature
summary(fire_env$curvature)
curve_newdat <- expand.grid(
  easting = mean(fire_env$easting), 
  northing = mean(fire_env$northing),
  gp33_pdsi = mean(fire_env$gp33_pdsi),
  sdd_average = mean(fire_env$sdd_average),
  slope_percent = mean(fire_env$slope_percent), 
  curvature = seq(-0.74651, 1.69604, length.out = 90), 
  stand_age1 = mean(fire_env$stand_age1), 
  tsf = mean(fire_env$tsf), 
  site = unique(fire_env$site)
)

#Predict fire probability (marginal: fixed effects only)
curve_newdat$prob_fire <- as.numeric(predict(spamm2.1, curve_newdat, re.form = NA))

#Add 95% confidence intervals for fixed effects
curve_newdat <- cbind(
  curve_newdat,
  get_intervals(spamm2.1, newdata = curve_newdat, intervals = "fixefVar", re.form = NA)
)

#Collapse predictions to one line per curvature value (i.e., average over sites)
curve_plot_data <- curve_newdat %>%
  group_by(curvature) %>%
  summarise(
    prob_fire = mean(prob_fire, na.rm = TRUE),
    lower = mean(fixefVar_0.025, na.rm = TRUE),
    upper = mean(fixefVar_0.975, na.rm = TRUE),
    .groups = "drop"
  )

#Plot the smoothed prediction curve
curve_graph <- ggplot(curve_plot_data, aes(x = curvature, y = prob_fire)) +
  geom_line(color = "darkred", linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey50") +
  scale_y_continuous("Annual prob. fire") +
  scale_x_continuous("Curvature") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(plot.margin = unit(c(.5, 1, .5, 1), "lines"))
curve_graph

#Generate prediction grid for stand age
summary(fire_env$stand_age1)
standage_newdat <- expand.grid(
  easting = mean(fire_env$easting), 
  northing = mean(fire_env$northing),
  gp33_pdsi = mean(fire_env$gp33_pdsi),
  sdd_average = mean(fire_env$sdd_average),
  slope_percent = mean(fire_env$slope_percent), 
  curvature = mean(fire_env$curvature),  
  stand_age1 = seq(1, 730, length.out = 90), 
  tsf = mean(fire_env$tsf), 
  site = unique(fire_env$site)
)

#Predict fire probability (marginal: fixed effects only)
standage_newdat$prob_fire <- as.numeric(predict(spamm2.1, standage_newdat, re.form = NA))

#Add 95% confidence intervals for fixed effects
standage_newdat <- cbind(
  standage_newdat,
  get_intervals(spamm2.1, newdata = standage_newdat, intervals = "fixefVar", re.form = NA)
)

#Collapse predictions to one line per stand age value (i.e., average over sites)
standage_plot_data <- standage_newdat %>%
  group_by(stand_age1) %>%
  summarise(
    prob_fire = mean(prob_fire, na.rm = TRUE),
    lower = mean(fixefVar_0.025, na.rm = TRUE),
    upper = mean(fixefVar_0.975, na.rm = TRUE),
    .groups = "drop"
  )

#Plot the smoothed prediction curve
standage_graph <- ggplot(standage_plot_data, aes(x = stand_age1, y = prob_fire)) +
  geom_line(color = "darkred", linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey50") +
  scale_y_continuous("Annual prob. fire") +
  scale_x_continuous("Stand age") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(plot.margin = unit(c(.5, 1, .5, 1), "lines"))
standage_graph

#Generate prediction grid for time since fire
summary(fire_env$tsf)
tsf_newdat <- expand.grid(
  easting = mean(fire_env$easting), 
  northing = mean(fire_env$northing),
  gp33_pdsi = mean(fire_env$gp33_pdsi),
  sdd_average = mean(fire_env$sdd_average),
  slope_percent = mean(fire_env$slope_percent), 
  curvature = mean(fire_env$curvature),  
  stand_age1 = mean(fire_env$stand_age1), 
  tsf = seq(0, 420, length.out = 90), 
  site = unique(fire_env$site)
)

#Predict fire probability (marginal: fixed effects only)
tsf_newdat$prob_fire <- as.numeric(predict(spamm2.1, tsf_newdat, re.form = NA))

#Add 95% confidence intervals for fixed effects
tsf_newdat <- cbind(
  tsf_newdat,
  get_intervals(spamm2.1, newdata = tsf_newdat, intervals = "fixefVar", re.form = NA)
)

#Collapse predictions to one line per tsf value (i.e., average over sites)
tsf_plot_data <- tsf_newdat %>%
  group_by(tsf) %>%
  summarise(
    prob_fire = mean(prob_fire, na.rm = TRUE),
    lower = mean(fixefVar_0.025, na.rm = TRUE),
    upper = mean(fixefVar_0.975, na.rm = TRUE),
    .groups = "drop"
  )

#Plot the smoothed prediction curve
tsf_graph <- ggplot(tsf_plot_data, aes(x = tsf, y = prob_fire)) +
  geom_line(color = "darkred", linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey50") +
  scale_y_continuous("Annual prob. fire") +
  scale_x_continuous("Time since fire") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(plot.margin = unit(c(.5, 1, .5, 1), "lines"))
tsf_graph

#Graph what our visual and quantitative assessment (see confidence intervals above) demonstrate to be significant relationships.  
significant_predicts <- ggpubr::ggarrange(sdd_graph, pdsi_graph, tsf_graph, ncol = 1, nrow = 3)
ggsave(path="/YOURPATH", "significant_predicts.png", width = 3, height = 7, units = "in", dpi=800)


####################################
#### Evaluate watershed cohorts ####
####################################

#Here we're going to run the same simulation procedure that we used to identify tree establishment cohorts within stands at the scale of all tree establishment dates within the six large river drainages in which study sites are located.  The simulation procedure part of this code is redundant with a previous section, just implemented on a coarser aggregation of the tree establishment data set.  

#Clear workspace and set directory
rm(list=ls())
setwd("/YOURPATH")

#Load libraries for this section
library(tidyverse)
library(tidytext)
library(ggplot2)
library(leaflet)
library(KernSmooth)
library(ggpubr)

#Read cross-sections file, select pertinent columns, drop series with no first or last year, and convert to numeric
xsections <- read.csv("xsections.csv")
xsections <- xsections[,c("site", "sample_id", "species", "sample_ht_cm", "d10_mm", "pith","outer_ring", "inner_ring", "latitude", "longitude")]
xsections <- xsections %>% drop_na(outer_ring) 
xsections <- xsections %>% drop_na(inner_ring) 
xsections$pith <- as.numeric(xsections$pith)
xsections$outer_ring <- as.numeric(xsections$outer_ring)
xsections$inner_ring <- as.numeric(xsections$inner_ring)
str(xsections)

#Make all samples 6 digits (to remove, A, B, C, D, etc.) and convert to numeric
xsections$sample_id <- substr(xsections$sample_id, 1, 6)

#Collapse data from A, B, C, D, etc. of same sample
xsections_collapse <- xsections %>% 
  group_by(sample_id) %>% 
  summarize(outer_ring=collapse::fmax(outer_ring, na.rm = TRUE),
            inner_ring=collapse::fmin(inner_ring, na.rm = TRUE),
            pith=collapse::fmin(pith, na.rm = TRUE))
print(xsections_collapse, n=57)

#Add back the other xsections data
addback <- xsections %>% select("site", "sample_id", "species", "sample_ht_cm", "d10_mm", "latitude", "longitude")
xsections_final <- left_join(xsections_collapse, addback)
target_xsections <- xsections_final 
head(target_xsections)

#Get counts of samples per site
sample_count <- target_xsections %>% 
  group_by(site) %>% 
  tally()
data.frame(sample_count[order(-sample_count$n),])
min(sample_count$n)
max(sample_count$n)
mean(sample_count$n)
sum(sample_count$n)

#Read height calibration data
hts <- read.csv("ht_calibration.csv")

#Add dummy species variable to xsections file to match the column name in the height calibration file and change these dummy species names that aren't in the ht_calibration file (so as not to alter the species as it will be coded below)
target_xsections$spp <- target_xsections$species
target_xsections$spp <- ifelse(target_xsections$spp=="" | target_xsections$spp=="THPL" | target_xsections$spp=="PIMO" | target_xsections$spp=="TSME" | target_xsections$spp=="UNKN" | target_xsections$spp=="LAOC" | target_xsections$spp=="ABPR?" | target_xsections$spp=="TSHE?", "PSME", target_xsections$spp)
target_xsections$d10_mm <- as.numeric(target_xsections$d10_mm)
target_xsections$sample_ht_cm <- as.numeric(target_xsections$sample_ht_cm)

#Read environmental variables file, select HUC8 names, reduce to most common HUC8 name per site, and join to cross-sections
env_vars <- read.csv("env_data.csv")
env_vars <- env_vars %>% select(site, huc8_name)
env_vars$site <- as.character(env_vars$site)
env_vars <- env_vars %>% group_by(site) %>% count(huc8_name) %>% slice(which.max(n))
target_xsections$site <- as.character(target_xsections$site)
target_xsections <- target_xsections %>% left_join(env_vars) 
head(target_xsections)

#Create a linear model for years to mineral soil
lm1 <- lm(rings_to_soil ~ sample_ht_cm + spp + d10_mm, data=hts)

#Predict years to mineral soil for xsections, bind that to xsections file, and make any predictions that are negative numbers 0
predict_years <- predict(lm1, newdata = target_xsections)
target_xsections <- cbind(target_xsections, predict_years)
target_xsections$predict_years <- ifelse(target_xsections$predict_years < 0, 0, target_xsections$predict_years)
target_xsections$pith <- as.numeric(target_xsections$pith)
target_xsections$predict_years <- round(target_xsections$predict_years)
target_xsections$estab_year <- target_xsections$pith-target_xsections$predict_years

#Cohort detection function
cohort_function <- function(x, yearsvec, nosims, siglevel){
  estab_vec <- x[,yearsvec]
  estab_vec_pad <- c(min(estab_vec)-50, estab_vec, max(estab_vec)+50)
  estab <- if(max(estab_vec)-min(estab_vec) >= 150) 
    estab_vec else estab_vec_pad
  minest <- min(estab)
  maxest <- max(estab)
  notrees <- length(estab)
  nosims <- nosims
  unif_resampfunct <- function(){
    simestab <- sample(minest:maxest, size = notrees, replace = TRUE)
    return(simestab)
  }
  sim_data <- data.frame(replicate(n = nosims, 
                                   expr = unif_resampfunct()))
  names(sim_data) <- gsub(x = names(sim_data), pattern = "\\X",
                          replacement = "sim") 
  all_data <- cbind(sim_data, estab)
  bw <- dpik(estab)
  all_density <- apply(all_data, 2, bkde, bandwidth=bw/2.75)  
  long_density <- do.call(rbind.data.frame, all_density)
  long_density$run <- rownames(long_density)
  long_density$group <- substr(long_density$run, 1, 3)
  sim_dens <- long_density[(long_density$group=="sim"),]
  est_dens <- long_density[(long_density$group=="est"),]
  sim_dens_sig <- quantile(sim_dens$y, siglevel)
  if(any(est_dens$y >= sim_dens_sig)){    
    estab_sig <- est_dens[(est_dens$y >= sim_dens_sig),]
    estab_sig$lower_bound <- rep(sim_dens_sig, nrow(estab_sig))
    estab_sig <- estab_sig[,c("x", "y", "lower_bound")]
    estab_sig$run_no <- as.numeric(sub(".*b.", "", 
                                       rownames(estab_sig)))
    rownames(estab_sig) <- 1:nrow(estab_sig)
    estab_sig$rows <- as.numeric(rownames(estab_sig))
    ints <- c(0, which(diff(estab_sig$run_no) != 1), length(estab_sig$run_no))
    estab_sig$cohort_no <- cut(estab_sig$rows, breaks=ints, labels=FALSE)
    estab_sig$rows <- NULL
  } else {
    estab_sig <- data.frame(x = 0, y = 0, lower_bound = 0, run_no = 0, cohort_no="NA")
  }
  return(list(crit_value=sim_dens_sig, sim_dens=sim_dens, est_dens=est_dens, estab_sig=estab_sig))
}

#Create a file to run the cohort detection procedure on,  create strata column (for splitting... could just be just one if you don't want to split it), and add site HUCs
xsections_cohorts <- target_xsections %>% drop_na(estab_year) %>% select(sample_id, estab_year, huc8_name)
head(xsections_cohorts)

#Split the establishment data into a list split by watershed 
strata_list <- split(xsections_cohorts, xsections_cohorts$huc8_name)

#Run the cohort detection function on every element of the list (run the function for each strata).  "no-sims" sets the number of simulated establishment years (for each site) and siglevel sets the critical threshold for evaluating significant of actual establishment data.  1,000 simulations will yield a reasonable estimate of significant cohorts and will run on most computers in less than 30 seconds.  I used 10,000 simulations in the paper.  
set.seed(123)
results_list <- map(strata_list, cohort_function, yearsvec="estab_year", nosims=1000, siglevel=.99)

#View the resulting list structure to get an idea how the function reports out results
data.tree::FromListSimple(results_list)

#Create a dataframe with critical values for cohort significance
crit_values <- unname(unlist(map(results_list, 1)))
crit_values_df <- data.frame(strata=substr(names(results_list), 1, 2), crit_value=crit_values) 

#Create data frame of real establishment determined to represent tree cohorts (significant given critical threshold). 
estab_sig <- map(results_list, 4)
estab_sig <- do.call(rbind.data.frame, estab_sig)
estab_sig$strata <- substr(rownames(estab_sig), 1, 2)
estab_sig$strata <- gsub("[[:punct:]]", "", estab_sig$strata)

#Some necessary housekeeping:  Remove the critical values and significant cohorts data from the results list
results_list <- lapply(results_list, function(x) x[-1])
results_list <- lapply(results_list, function(x) x[-3])

#Create data frame with establishment and simulated establishment kernel density estimate curves
est_sim <- do.call(rbind, unlist(results_list, recursive = FALSE))
est_sim$strata <- substr(rownames(est_sim), 1, 2)
est_sim$strata <- gsub("[[:punct:]]", "", est_sim$strata)
rownames(est_sim) <- NULL

#Split this dataframe into actual and simulated establishment
est_dens <- est_sim %>% filter(group=="est")
sim_dens <- est_sim %>% filter(group=="sim")

#Min and max of each cohort (for the estab graphic)
cohort_min_max <- estab_sig %>% 
  group_by(strata, cohort_no) %>% 
  summarize(min_cohort=min(x),
            max_cohort=max(x))
colnames(cohort_min_max)[which(colnames(cohort_min_max) == 'strata')] <- 'huc8_name'

#Rename and recode strata to match estab graphic files
cohort_min_max$huc8_name <- dplyr::recode(cohort_min_max$huc8_name,
        "Lower Columbia-Sandy" = "Lo",
        "Clackamas" = "Cl",
        "North Santiam" = "No",
        "South Santiam" = "So",
        "Mckenzie" = "Mc",
        "Middle Fork Willamette" = "Mi")
cohort_min_max

#Create graphical file for the cohort detection output 
estab_sig$strata <- factor(estab_sig$strata, levels = c("Lo", "Cl", "No", "So", "Mc", "Mi"))
sim_dens$strata <- factor(sim_dens$strata, levels = c("Lo", "Cl", "No", "So", "Mc", "Mi"))
est_dens$strata <- factor(est_dens$strata, levels = c("Lo", "Cl", "No", "So", "Mc", "Mi"))
crit_values_df$strata <- factor(crit_values_df$strata, levels = c("Lo", "Cl", "No", "So", "Mc", "Mi"))
establishment_cohorts <- ggplot() +
  geom_ribbon(data=estab_sig, aes(x=x, ymax=y, ymin=lower_bound, group=cohort_no, fill=cohort_no), fill="darkolivegreen", alpha=0.6) +
  geom_path(data=sim_dens, aes(x=x, y=y), color="snow4", linewidth=.025, alpha=.1) +
  geom_path(data=est_dens, aes(x=x, y=y), color="darkolivegreen", linewidth=.65, alpha=1) +
  geom_hline(data=crit_values_df, aes(yintercept=crit_value), linetype='longdash', color = "black", linewidth=.25) +
  scale_x_continuous("", breaks = c(1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900), limits=c(1150, 1950)) + 
  scale_y_continuous("Kernel density estimate", breaks = c(.0, .02), position = "right") + 
  facet_wrap(~ strata, ncol = 1, strip.position="left") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank(), panel.grid.minor = element_line(linewidth = 0.1), panel.grid.major = element_line(linewidth = .25), panel.spacing.y = unit(.35, "lines"), plot.margin=unit(c(.25,.25,.25,.25), 'cm')) 
ggsave(path="/YOURPATH", "establishment_cohorts.png", width = 4.25, height = 7.25, units = "in", dpi=800)

#Tally establishment in the demonstration sites and rename the resulting file
xsections_cohorts$decade_bin <- plyr::round_any(xsections_cohorts$estab_year, 10)
estab_tally <- xsections_cohorts %>%
  group_by(huc8_name, decade_bin) %>%
  count()
estab_tally

#Rename and recode strata to match estab graphic files
estab_tally$huc8_name <- dplyr::recode(estab_tally$huc8_name,
                  "Lower Columbia-Sandy" = "Lo",
                  "Clackamas" = "Cl",
                  "North Santiam" = "No",
                  "South Santiam" = "So",
                  "Mckenzie" = "Mc",
                  "Middle Fork Willamette" = "Mi")
estab_tally

#scale_color_manual(name = "", values=c("Lower Columbia-Sandy" = "green2", "Clackamas"="orangered", "North Santiam"="cyan2", "South Santiam"="darkorchid2", "Mckenzie"="goldenrod1", "Middle Fork Willamette"="violetred1"), guide = "none")

#Histogram of the demonstration sites
estab_tally$huc8_name <- factor(estab_tally$huc8_name, levels=c("Lo", "Cl", "No", "So", "Mc", "Mi"))
cohort_min_max$huc8_name <- factor(cohort_min_max$huc8_name, levels=c("Lo", "Cl", "No", "So", "Mc", "Mi"))
#quartz(width=7, height=7)
estab_histogram <- ggplot() +
  geom_rect(data=cohort_min_max, aes(xmin=min_cohort, xmax=max_cohort, ymin=-Inf, ymax=Inf), color=NA, fill="darkolivegreen", alpha=.2) +
  geom_bar(data=estab_tally, aes(x=decade_bin, y=n, fill=huc8_name), color="black", stat = "identity", position = "stack", width=8, linewidth=.175) +
  scale_fill_manual(name = "", values=c("Lo" = "green2", "Cl"="orangered", "No"="cyan2", "So"="darkorchid2", "Mc"="goldenrod1", "Mi"="violetred1"), guide = "none") +
  scale_x_continuous("", breaks = c(1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900), limits=c(1150, 1950)) + 
  scale_y_continuous("Count of tree establishment by decade", breaks = c(0, 20), limits=c(0, 32)) +
  facet_wrap(~ huc8_name, ncol = 1, strip.position="right") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), strip.background = element_blank(), panel.grid.minor = element_line(linewidth = 0.1), panel.grid.major = element_line(linewidth = .25), panel.spacing.y = unit(.35, "lines"), plot.margin=unit(c(.25,.25,.25,.25), 'cm')) +
  theme(legend.position="none")
estab_histogram
ggsave(path="/YOURPATH", "estab_histogram.png", width = 4.25, height = 7.25, units = "in", dpi=800)

#Optional:  Combine both graphs and print 
#combine_graphs <- ggarrange(agestructurehistogram, #cohorts_graph, heights = c(3, 2), labels = c("", ""), ncol = 1, nrow = 2, align = "v")
#ggsave("combine_graphs.png", dpi=500)


################################
##### Evaluate large fires #####
################################

#Here we're going to undertake extensive exploratory analysis using maps, a permutation procedure, and evaluation of overlap of tree cohorts to attempt to determine if our data can be used to reconstruct all or part of the footprint of large fires that occurred in the distant past.  This is a pretty involved procedure that will take quite a bit of time to run.  

#Clear workspace and set directory
rm(list=ls())
setwd("/YOURPATH")

#Load libraries needed for this section
library(tidyverse)
library(RcppAlgos) 
library(geodist)
library(elevatr)
library(terra)
library(sf)
library(ggnewscale)
library(KernSmooth)
library(igraph)
library(geosphere)

## Evaluate synchrony in fire years ##

#Read fires file, remove non-fire injuries, and select pertinent columns, 
events <- read.csv("events.csv")
events <- events %>% filter(type=="FS") 
events <- events[,c("site", "sample_id", "event_year")]
unique(events$site)

#Read site centroids file and change site to character
site_coords <- read.csv("site_centroids.csv")
site_coords$site <- as.character(site_coords$site)
head(site_coords)

#Calculate distances between sites (run geodist function, turn 0s to NA, and find min, max and mean by column)
site_dist <- geodist(site_coords, measure = 'geodesic')
site_dist[site_dist == 0] <- NA
apply(site_dist,2,min, na.rm=TRUE)
min(apply(site_dist,2,min, na.rm=TRUE))
mean(apply(site_dist,2,min, na.rm=TRUE))
min(apply(site_dist,2,max, na.rm=TRUE))
max(apply(site_dist,2,max, na.rm=TRUE))

#Summarize unique fire years
fire_recons <- events %>% distinct(event_year) %>% arrange(event_year)
fire_recons %>% filter(event_year>=1900)

#Summarize unique fire years by site and convert site to character
fire_years_site <- events %>% 
  group_by(site) %>% 
  distinct(event_year) %>% 
  rename("fire_year"="event_year")  %>% 
  arrange(site, fire_year)  
fire_years_site$site <- as.character(fire_years_site$site)
fire_years_site

#Join unique fire years by site to coordinates
fire_years_site <- fire_years_site %>%
  left_join(site_coords)
fire_years_site

#Permute fire years and sites function... this one is a bear...  see all the way down towards the bottom of this section for the really simple way to do this...  I couldn't use the really simple way because (among othe reasons) I felt like I needed to make sure there was a maximum and minimum number of fire years to stuff into each site when permuting to match the real data... In the more complex permutation procedure coded below, we're creating a file for joining back the site coordinates, and reducing the fire years file to fires that burned more than one site in the same year.  Then, the heart of the function: I use the partitionsSample function to assign random group sizes.  Then I use the sample_degseq function to randomly assign values to groups.  Then I run that custom random split function on the data and map it to a dataframe.  Then I replicate the whole thing n times.  Then do some cleanup of the resulting files.  This is a somewhat hacky piece of code that can and should be cleaned up but it does what it's supposed to do and I haven't had time to return to it.  Note that we're going to split it again when we run it, which is a little silly, but it works.  I'm sure I'll return to this procedure and streamline it at some point, so write to me before you get to frustrated with the hacky parts. 
permute_function <- function(x, noperms){
  for_join <- x  %>%
    select(site, longitude, latitude) %>% 
    distinct() %>% 
    arrange(site)
  for_join$site <- as.character(for_join$site)
  shuffle <- function(x){
    real_data <- x %>% 
      group_by(fire_year) %>% 
      filter(n() > 1) %>%
      arrange(fire_year)
    random_split <- function(values, groups, size.min, size.max) {
      R <- table(values)
      S <- partitionsSample( # assign random group sizes
        size.min:size.max, length(groups), TRUE, target = length(values), n = 1
      )[1,]
      d <- length(R) - length(S)
      with(
        as_data_frame(
          sample_degseq( # randomly assign values to groups
            c(R, integer(max(0, -d))),
            c(sample(S), integer(max(0, d))),
            "simple.no.multiple.uniform"
          )
        ),
        split(as(names(R), class(values))[from], groups[to])
      )
    }
    split_results <- random_split(real_data$fire_year, sample(unique(x$site), length(unique(real_data$site))), min(table(real_data$site)), max(table(real_data$site)))
    perm_df <- map_df(split_results, ~as.data.frame(.x), .id="site") %>%
      rename("fire_year"=".x") 
  }
  perm_list <- replicate(n=noperms, expr=shuffle(x), simplify = FALSE)
  all_perms <- dplyr::bind_rows(perm_list, .id = 'run')
  all_perms <- all_perms %>% left_join(for_join, by="site")
}

#Run permutation procedure... this takes about 15 seconds with 1,000 permutations.
perm_results <- permute_function(fire_years_site, 1000)
head(perm_results)
tail(perm_results)
nrow(perm_results)

#Create file with real sites and all real fire years that burned more than one site (this is what the first part of the function does)
fire_years_dup <- fire_years_site %>% 
  group_by(fire_year) %>% 
  filter(n() > 1) %>%
  arrange(fire_year)
fire_years_dup
unique(fire_years_dup$fire_year)
length(unique(fire_years_dup$fire_year))

#Split the real results into a list by fire year
real_split <- split(fire_years_dup, fire_years_dup$fire_year)

#Map geodist function to elements of the real fire list, convert to dataframe and do some manipulation of that dataframe
real_fire_dist <- map(real_split, geodist, measure = 'geodesic')
real_fire_dist_df <- data.frame(unlist(real_fire_dist))
real_fire_dist_df <- real_fire_dist_df %>% 
  rename("dist_m" = "unlist.real_fire_dist.") %>% 
  filter(dist_m > 0)
real_fire_dist_df$fire_year <- substr(rownames(real_fire_dist_df), 1, 4)
head(real_fire_dist_df)

#Summarize distances for each real fire year (divide by one thousand to put in terms of kilometers)
real_fire_dist_summary <- real_fire_dist_df %>%
  group_by(fire_year) %>%
  summarize(mean_dist_km = mean(dist_m)/1000)
real_fire_dist_summary

#Split the permutation results into a list by fire year
perm_split <- split(perm_results, perm_results$fire_year)

#Map geodist function to elements of the permutation list, convert to dataframe and do some manipulation of that dataframe... this may take up to 3-5 minutes.  
perm_fire_dist <- map(perm_split, geodist, measure = 'geodesic')
perm_fire_dist_df <- data.frame(unlist(perm_fire_dist))
perm_fire_dist_df <- perm_fire_dist_df %>% 
  rename("dist_m" = "unlist.perm_fire_dist.") %>% 
  filter(dist_m > 0)
perm_fire_dist_df$fire_year <- substr(rownames(perm_fire_dist_df), 1, 4)
head(perm_fire_dist_df)

#Summarize distances for each permuted fire year (divide by one thousand to get in terms of kilometers)
perm_fire_dist_summary <- perm_fire_dist_df %>%
  group_by(fire_year) %>%
  summarize(mean_dist_km = mean(dist_m)/1000)
perm_fire_dist_summary

#Plot distances... we're looking for fire years that burned more than one site that have inter-site distances well below simulated inter-site fire distances.... this may take quite a bit of time to run.
dist_box <- ggplot(perm_fire_dist_df, aes(x=fire_year, y=dist_m/1000)) +
  geom_boxplot(position="dodge", color="black", fill="steelblue4") + 
  geom_point(data=real_fire_dist_summary, aes(x=fire_year, y=mean_dist_km), size=3, shape = 17, color="darkred") +
  scale_x_discrete("Fire year") +
  scale_y_continuous("Inter-site distance (km)") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) 
dist_box
ggsave(path="/YOURPATH", "dist_box.png", width = 20, height = 12, units = "cm", dpi=800)

#So, we have the following fire years to evaluate that appear to have occurred closer together in space than would have been predicted by chance:  1499, 1768, 1783, 1798, 1831, 1846, 1883, and 1958... probably should also check out 1807, 1812, and 1844.
flist <- c(1499, 1768, 1783, 1798, 1807, 1812, 1831, 1844, 1846, 1883, 1958)
#Upon reflection, there's little chance that most of these fires were part of the same fire perimeter.  Let's look at instead just:
flist <- c(1499, 1768, 1783, 1798, 1831, 1846, 1958)
examine <- fire_years_dup %>% filter(fire_year %in% flist)
examine$group <- ifelse(examine$fire_year <1800, 1, ifelse(examine$fire_year >=1883, 3, 2))

#Map the oldest close in space fire years
filterexam <- examine %>% filter(group==1)
early_fires <- ggplot(data=filterexam, aes(x=longitude, y=latitude)) +
  geom_point(aes(color=as.factor(fire_year))) +
  geom_text(aes(label=site), color="black", size=4, hjust = 0.0, nudge_x = 0.005, vjust = 0.0, nudge_y = 0.005, family = "Arial", check_overlap = TRUE) +
  geom_polygon(aes(fill = as.factor(fire_year),
                   colour = as.factor(fire_year)),
               alpha = 0.3) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude)) +
  theme_bw()
early_fires
ggsave(path="/YOURPATH", "early_fires.png", width = 16, height = 16, units = "cm", dpi=800)

#Map the next oldest close-in-space fire years
filterexam <- examine %>% filter(group==2)
#quartz(width=8, height=8)
mid_fires <- ggplot(data=filterexam, aes(x=longitude, y=latitude)) +
  geom_point(aes(color=as.factor(fire_year))) +
  geom_text(aes(label=site), color="black", size=4, hjust = 0.0, nudge_x = 0.005, vjust = 0.0, nudge_y = 0.005, family = "Arial", check_overlap = TRUE) +
  geom_polygon(aes(fill = as.factor(fire_year),
                   colour = as.factor(fire_year)),
               alpha = 0.3) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude)) +
  theme_bw()
mid_fires
ggsave(path="/YOURPATH", "mid_fires.png", width = 16, height = 16, units = "cm", dpi=800)

#Map the youngest close in space fire years
filterexam <- examine %>% filter(group==3)
#quartz(width=8, height=8)
late_fires <- ggplot(data=filterexam, aes(x=longitude, y=latitude)) +
  geom_point(aes(color=as.factor(fire_year))) +
  geom_text(aes(label=site), color="black", size=4, hjust = 0.0, nudge_x = 0.005, vjust = 0.0, nudge_y = 0.005, family = "Arial", check_overlap = TRUE) +
  geom_polygon(aes(fill = as.factor(fire_year),
                   colour = as.factor(fire_year)),
               alpha = 0.3) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude)) +
  theme_bw()
late_fires
ggsave(path="/YOURPATH", "late_fires.png", width = 16, height = 16, units = "cm", dpi=800)

#Map all fire years
#quartz(width=8, height=8)
all_fires <- ggplot(data=examine, aes(x=longitude, y=latitude)) +
  geom_point(aes(color=as.factor(fire_year))) +
  geom_text(aes(label=site), color="black", size=4, hjust = 0.0, nudge_x = 0.005, vjust = 0.0, nudge_y = 0.005, family = "Arial", check_overlap = TRUE) +
  geom_polygon(aes(fill = as.factor(fire_year),
                   colour = as.factor(fire_year)),
               alpha = 0.3) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude)) +
  theme_bw()
all_fires
ggsave(path="/YOURPATH", "all_fires.png", width = 16, height = 16, units = "cm", dpi=800)

### Evaluate cohorts ###

#Here we're going to run the simulation procedure we've already run.  Below I'm going to run 1,000 simulations, but I used 1,000 in the paper. 

#Read cross-sections file, select pertinent columns, and drop series with no first or last year
xsections <- read.csv("xsections.csv")
xsections <- xsections[,c("site", "sample_id", "species", "sample_ht_cm", "d10_mm", "pith","outer_ring", "inner_ring", "latitude", "longitude")]
xsections <- xsections %>% drop_na(outer_ring) 
xsections <- xsections %>% drop_na(inner_ring) 

#Make all samples 6 digits (to remove, A, B, C, D, etc. because we are combining histories from the same tree), and convert numeric columns to numeric
xsections$sample_id <- substr(xsections$sample_id, 1, 6)
xsections$pith <- as.numeric(xsections$pith)
xsections$outer_ring <- as.numeric(xsections$outer_ring)
xsections$inner_ring <- as.numeric(xsections$inner_ring)
xsections$sample_ht_cm <- as.numeric(xsections$sample_ht_cm)
xsections$d10_mm <- as.numeric(xsections$d10_mm)
nrow(xsections)
length(unique(xsections$sample_id))

#Collapse data from A, B, C, D, etc. of same sample and make result a dataframe (for some reason the tibble output seems to choke the establishment data model).
xsections_collapse <- xsections %>% 
  group_by(site, sample_id, species) %>% 
  summarize(
    latitude=collapse::fmin(latitude, na.rm = TRUE),
    longitude=collapse::fmin(longitude, na.rm = TRUE),
    sample_ht_cm=collapse::fmin(sample_ht_cm, na.rm = TRUE),
    d10_mm=collapse::fmin(d10_mm, na.rm = TRUE),
    outer_ring=collapse::fmax(outer_ring, na.rm = TRUE),
    inner_ring=collapse::fmin(inner_ring, na.rm = TRUE),
    pith=collapse::fmin(pith, na.rm = TRUE))
xsections_collapse <- data.frame(xsections_collapse)
str(xsections_collapse)

#Add back the xsections data
addback <- xsections %>% select("site", "sample_id", "species", "sample_ht_cm", "d10_mm", "latitude", "longitude")
xsections_final <- left_join(xsections_collapse, addback)
target_xsections <- xsections_final 
head(target_xsections)

#Change continuous variables to numeric
target_xsections$pith <- as.numeric(target_xsections$pith)
target_xsections$outer_ring <- as.numeric(target_xsections$outer_ring)
target_xsections$inner_ring <- as.numeric(target_xsections$inner_ring)

#Get counts of samples per site
sample_count <- target_xsections %>% 
  group_by(site) %>% 
  tally()
data.frame(sample_count[order(-sample_count$n),])
min(sample_count$n)
max(sample_count$n)
mean(sample_count$n)
sum(sample_count$n)

#Read site centroids file and filter to Oregon national forest sites, change site to character
site_coords <- read.csv("site_centroids.csv")
site_coords$site <- as.character(site_coords$site)
head(site_coords)

#Read height calibration data
hts <- read.csv("ht_calibration.csv")

#Add dummy species variable to xsections file to match the column name in the height calibration file and change these dummy species names that aren't in the ht_calibration file (so as not to alter the species as it will be coded below)
target_xsections$spp <- target_xsections$species
target_xsections$spp <- ifelse(target_xsections$spp=="" | target_xsections$spp=="THPL" | target_xsections$spp=="PIMO" | target_xsections$spp=="TSME" | target_xsections$spp=="UNKN" | target_xsections$spp=="LAOC" | target_xsections$spp=="ABPR?" | target_xsections$spp=="TSHE?", "PSME", target_xsections$spp)
target_xsections$d10_mm <- as.numeric(target_xsections$d10_mm)
target_xsections$sample_ht_cm <- as.numeric(target_xsections$sample_ht_cm)

#Create a linear model for years to mineral soil
lm1 <- lm(rings_to_soil ~ sample_ht_cm + spp + d10_mm, data=hts)

#Predict years to mineral soil for xsections, bind that to xsections file, and make any predictions that are negative numbers 0
predict_years <- predict(lm1, newdata = target_xsections)
target_xsections <- cbind(target_xsections, predict_years)
target_xsections$predict_years <- ifelse(target_xsections$predict_years < 0, 0, target_xsections$predict_years)
target_xsections$pith <- as.numeric(target_xsections$pith)
target_xsections$predict_years <- round(target_xsections$predict_years)
target_xsections$estab_year <- target_xsections$pith-target_xsections$predict_years

#Cohort detection function
cohort_function <- function(x, yearsvec, nosims, siglevel){
  estab_vec <- x[,yearsvec]
  estab_vec_pad <- c(min(estab_vec)-50, estab_vec, max(estab_vec)+50)
  estab <- if(max(estab_vec)-min(estab_vec) >= 150) 
    estab_vec else estab_vec_pad
  minest <- min(estab)
  maxest <- max(estab)
  notrees <- length(estab)
  nosims <- nosims
  unif_resampfunct <- function(){
    simestab <- sample(minest:maxest, size = notrees, replace = TRUE)
    return(simestab)
  }
  sim_data <- data.frame(replicate(n = nosims, 
                                   expr = unif_resampfunct()))
  names(sim_data) <- gsub(x = names(sim_data), pattern = "\\X",
                          replacement = "sim") 
  all_data <- cbind(sim_data, estab)
  bw <- dpik(estab)
  all_density <- apply(all_data, 2, bkde, bandwidth=bw)  
  long_density <- do.call(rbind.data.frame, all_density)
  long_density$run <- rownames(long_density)
  long_density$group <- substr(long_density$run, 1, 3)
  sim_dens <- long_density[(long_density$group=="sim"),]
  est_dens <- long_density[(long_density$group=="est"),]
  sim_dens_sig <- quantile(sim_dens$y, siglevel)
  if(any(est_dens$y >= sim_dens_sig)){    
    estab_sig <- est_dens[(est_dens$y >= sim_dens_sig),]
    estab_sig$lower_bound <- rep(sim_dens_sig, nrow(estab_sig))
    estab_sig <- estab_sig[,c("x", "y", "lower_bound")]
    estab_sig$run_no <- as.numeric(sub(".*b.", "", 
                                       rownames(estab_sig)))
    rownames(estab_sig) <- 1:nrow(estab_sig)
    estab_sig$rows <- as.numeric(rownames(estab_sig))
    ints <- c(0, which(diff(estab_sig$run_no) != 1), length(estab_sig$run_no))
    estab_sig$cohort_no <- cut(estab_sig$rows, breaks=ints, labels=FALSE)
    estab_sig$rows <- NULL
  } else {
    estab_sig <- data.frame(x = 0, y = 0, lower_bound = 0, run_no = 0, cohort_no="NA")
  }
  cohort_table <- estab_sig %>%
    group_by(cohort_no) %>%
    summarize(min_cohort_year = round(min(x),0),
              max_cohort_year = round(max(x),0))
  cohort_table
}

#Create a cross-sections file for the cohort detection procedure
xsections_cohorts <- target_xsections %>% drop_na(estab_year) 
#xsections_touse <- droplevels(xsections_touse)

#Split the xsections dataframe into a list and run function on every element of the list.  IMPORTANT NOTE:  The "map" function is a purrr function which seems to work without loading purrr explicitly.  This library has given me problems in the past (and old version is loaded as of this writing).  If the code breaks down here, it's probably a problem with purrr. 
set.seed(1234)
site_list <- split(xsections_cohorts, xsections_cohorts$site)
results_list <- map(site_list, cohort_function, yearsvec="estab_year", nosims=1000, siglevel=.99)

#Convert list to dataframe with site column, remove 0s, add unique site/cohort identifier, and pad all cohorts by a couple of years
wfi_cohorts <- do.call(rbind, unname(Map(cbind, site = names(results_list), results_list)))
wfi_cohorts <- wfi_cohorts %>% filter(min_cohort_year>0)
wfi_cohorts[order(-wfi_cohorts$min_cohort_year),] 
wfi_cohorts$site_cohort <- paste(wfi_cohorts$site, wfi_cohorts$cohort_no, sep=" no")
wfi_cohorts$min_cohort_year <- wfi_cohorts$min_cohort_year-2
wfi_cohorts$max_cohort_year <- wfi_cohorts$max_cohort_year+2

#Add column indicating cohorts that overlap with another cohort
wfi_cohorts_range <- wfi_cohorts %>%
  mutate(Range = ivs::iv(min_cohort_year, max_cohort_year), .keep = "unused")
wfi_cohorts_range <- wfi_cohorts_range %>%
  mutate(Overlap = ivs::iv_count_overlaps(Range, Range) > 1L)
wfi_cohorts <- wfi_cohorts %>% left_join(wfi_cohorts_range)
head(wfi_cohorts)

#Join site coordinates to cohort file
wfi_cohorts <- wfi_cohorts %>% left_join(site_coords)

#Graph
cohorts_graph <- ggplot() +
  geom_segment(data=wfi_cohorts, aes(x = min_cohort_year, y = reorder(site_cohort, min_cohort_year), xend = max_cohort_year, yend = factor(site_cohort), color=Overlap), linewidth=6, alpha=1) +
  scale_x_continuous("", breaks=seq(1250,2000,50), limits=c(1250,2000)) +
  scale_y_discrete("Site") +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(path="/YOURPATH", "cohorts_graph.png", width = 20, height = 16, units = "cm", dpi=800)

#The following cohorts that overlap are close enough together in space (based on the analysis of fire years that likely constitue the same fire perimeter, see below) that they may have originated as a consequence of the same fire perimeter:

#To map:  
#Sites 60, 61, 62, 64, and 71 had cohorts initiated in 1700, 1708, 1703, 1698, and 1700 ***
#Sites 56 and 58 had cohorts initiated in 1500 and 1498 respectively and sites 54 and 58 had a fire year recorded in 1499.  Map all of these Bull Run sites. **
#Sites 87 and 90 had cohorts initiated in 1522 and 1518 respectively
#Sites 88 and 89 had cohorts initiated in 1569 and 1571, and sites 89 and 93 had cohorts initiated in 1720 and 1713 respectively
#Sites 26 and 83 had cohorts initiated in 1531 and 1526
#Sites 25 and 23 had cohorts initiated in 1558 and 1564
#Sites 22 and 28 had cohorts initiated in 1783 and 1785
#Sites 28 and 82, fire year=1768
#Sites 28, 26, and 96, fire year=1783
#Sites 96 and 83, fire year=1798
clack <- c(60, 61, 62, 64, 71)
bullrun <- c(54, 56, 58)
mac <- c(84, 86)
santiam1 <- c(87, 90)
santiam2 <- c(88, 89)
mf1 <- c(26, 83)
mf2 <- c(25, 23)
mf3 <- c(28, 22)

#Extract the different cohorts and bind
forbind <- wfi_cohorts
forbind$Range <- NULL
bull_cohorts <- forbind %>% filter(site %in% bullrun)
bull_cohorts$group <- "bullrun"
clack_cohorts <- forbind %>% filter(site %in% clack)
clack_cohorts$group <- "clack"
mac_cohorts <- forbind %>% filter(site %in% mac)
mac_cohorts$group <- "mckenzie"
santiam1_cohorts <- forbind %>% filter(site %in% santiam1)
santiam1_cohorts$group <- "santiam1"
santiam2_cohorts <- forbind %>% filter(site %in% santiam2)
santiam2_cohorts$group <- "santiam2"
mf1_cohorts <- forbind %>% filter(site %in% mf1)
mf1_cohorts$group <- "mf1"
mf2_cohorts <- forbind %>% filter(site %in% mf2)
mf2_cohorts$group <- "mf2"
mf3_cohorts <- forbind %>% filter(site %in% mf3)
mf3_cohorts$group <- "mf3"
allcohs <- do.call("rbind", list(bull_cohorts, clack_cohorts, santiam1_cohorts, santiam2_cohorts, mf1_cohorts, mf2_cohorts, mf3_cohorts))
allcohs1 <- allcohs %>% filter(cohort_no==1)
allcohs2 <- allcohs %>% filter(cohort_no==2)

#Quick and dirty map of one potential fire... these next two maps are just for exploratory purposes.
thefire <- mf3_cohorts
onefirequickmap <- ggplot() +
  geom_point(data=thefire, aes(x=longitude, y=latitude, color=as.factor(site_cohort))) +
  geom_label(data=thefire, aes(x=longitude, y=latitude, label=site_cohort), color="black", size=3, hjust = 0.0, nudge_x=0.005, vjust=0.0, nudge_y=0.006, family = "Arial") +
  geom_text(data=thefire, aes(x=longitude, y=latitude, label=min_cohort_year), color="black", size=3, hjust = 0.0, nudge_x=0.08, vjust=-1.1, nudge_y=0.005, family = "Arial", check_overlap = TRUE) +
  geom_polygon(data=thefire, aes(x=longitude, y=latitude, fill=group, colour=group), linewidth=1, alpha = 0.3) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude)) +
  theme_bw()
ggsave(path="/YOURPATH", "onefirequickmap.png", width = 20, height = 20, units = "cm", dpi=300)
thefire

#Quick and dirty map of all
quickmap <- ggplot() +
  geom_point(data=allcohs1, aes(x=longitude, y=latitude, color=as.factor(site_cohort))) +
  geom_label(data=allcohs1, aes(x=longitude, y=latitude, label=site_cohort), color="black", size=3, hjust = 0.0, nudge_x=0.005, vjust=0.0, nudge_y=0.006, family = "Arial") +
  geom_text(data=allcohs1, aes(x=longitude, y=latitude, label=min_cohort_year), color="black", size=3, hjust = 0.0, nudge_x=0.08, vjust=-1.1, nudge_y=0.005, family = "Arial", check_overlap = TRUE) +
  geom_label(data=allcohs2, aes(x=longitude, y=latitude, label=site_cohort), color="black", size=3, hjust = 0.0, nudge_x=0.005, vjust=0, nudge_y=-0.04, family = "Arial") +
  geom_text(data=allcohs2, aes(x=longitude, y=latitude, label=min_cohort_year), color="black", size=3, hjust = 0.0, nudge_x=0.08, vjust=.9, nudge_y=0, family = "Arial", check_overlap = TRUE) +
  geom_polygon(data=allcohs2, aes(x=longitude, y=latitude, fill=group, colour=group), linewidth=3, linetype = "dashed", alpha = 0.3) +
  geom_polygon(data=allcohs1, aes(x=longitude, y=latitude, fill=group, colour=group), linewidth=1, alpha = 0.3) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude)) +
  theme_bw()
quickmap
ggsave(path="/YOURPATH", "quickmap.png", width = 20, height = 20, units = "cm", dpi=300)

### Make nicer maps for publication ###

#First, read environmental variables file, select HUC8 names, reduce to most common HUC8 name per site, and join to fire_years_dup file
head(fire_years_dup)
env_vars <- read.csv("env_data.csv")
env_vars <- env_vars %>% select(site, huc8_name)
env_vars <- env_vars %>% group_by(site) %>% count(huc8_name) %>% slice(which.max(n))
env_vars$site <- as.character(env_vars$site)
fire_years_dup <- fire_years_dup %>% 
  left_join(env_vars) %>%  
  select(site, fire_year, longitude, latitude, huc8_name)
head(fire_years_dup)
#Join huc8 name to site coordinates
site_coords <- site_coords %>% left_join(env_vars)

#Load data for all maps
#Define path to GeoPackages
gpkg_study <- "wfi_layers.gpkg"
gpkg_region <- "wfi_layers_region.gpkg"

#Read spatial layers for the study area map
study_area_nfs <- st_read(gpkg_study, layer = "national_forests", quiet = TRUE)
ugbs <- st_read(gpkg_study, layer = "urban_growth_boundaries", quiet = TRUE)
rivers <- st_read(gpkg_study, layer = "large_rivers", quiet = TRUE)
lakes_res <- st_read(gpkg_study, layer = "lakes_reservoirs", quiet = TRUE)
large_fires <- st_read(gpkg_study, layer = "large_fires", quiet = TRUE)

#Map of large fire in Bull Run watershed

#Download DEM for Bull Run extent, convert to raster and then a dataframe, estimate slope, aspect, orientation, hillshade, and combine those effects.
fireyear_extent <- data.frame(x = c(-122.17, -121.765), y = c(45.41, 45.53))
fireyear_extent
prj_dd <- "EPSG:4326"
fireyear_dem <- get_elev_raster(fireyear_extent, prj = prj_dd, z = 10)
fireyear_rast <- rast(fireyear_dem) 
fireyear_rast_df <- as.data.frame(fireyear_rast, xy = TRUE)
names(fireyear_rast_df)[3] <- "elev_m"
slope <- terrain(fireyear_rast, "slope", unit="radians")
aspect <- terrain(fireyear_rast, "aspect", unit="radians")
fireyear_hill <- shade(slope, aspect, 
                       angle = 45, 
                       direction = 300,
                       normalize= TRUE)
fireyear_hill_df <- as.data.frame(fireyear_hill, xy=TRUE)
head(fireyear_hill_df)

#Bull run map
bullrun_map <- ggplot() +
  geom_raster(data = fireyear_hill_df,
              aes(x, y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_raster(data = fireyear_rast_df,
              aes(x, y, fill = elev_m),
              alpha = .25) +
  scale_fill_gradientn(colours = c("seashell", "cornsilk2", "darkkhaki", "darkolivegreen4", "darkolivegreen", "darkolivegreen", "khaki4", "navajowhite4", "navajowhite3", "whitesmoke", "grey100"), breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200)) +
  #geom_sf(data=study_area_nfs, color="black", linewidth=.55, fill=NA) +
  geom_sf(data = ugbs, color="black", fill="grey70", linewidth=.1, alpha=.5) +
  geom_sf(data=rivers, color="navy", linewidth=.65, alpha=.85) +
  geom_sf(data=lakes_res, fill="navy", color="navy", alpha=.95) +
  #geom_sf_text(data = rivers, aes(label = Name)) +
  #geom_sf_text(data = ugbs, aes(label = instName)) +
  geom_sf(data = large_fires, color="darkred", fill="darkred", linewidth=.05, alpha=.2) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude, color=huc8_name), shape=22, fill="black", size=8, stroke=5) +
  scale_color_manual(name = "", values=c("Middle Fork Willamette"="violetred1", "Clackamas"="orangered", "Lower Columbia-Sandy" = "green2", "North Santiam"="cyan2", "Mckenzie"="goldenrod1", "South Santiam"="darkorchid2"), guide = "none") +
  geom_text(data=site_coords, aes(x=longitude, y=latitude, label=site), color="black", size=9, hjust = 0.0, nudge_x = 0.006, vjust = 0.0, nudge_y = 0.006, family = "Arial", check_overlap = TRUE) +
  #geom_text(data=fireyear1499, aes(x=longitude, y=latitude, label=fire_year), color="firebrick4", size=3.25, hjust = 0.5, nudge_x = 0, vjust = 2.5, nudge_y = 0, family = "Arial", check_overlap = TRUE) +
  scale_x_continuous("", breaks = c(-122, -121.8), position="top") +
  scale_y_continuous("", breaks = c(45.5)) +
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey10", "white"),
    text_family = "Arial",
    text_cex = 2,
    line_width = 1,
    height = unit(.55, "cm")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(2, "cm"), 
    width = unit(2, "cm"),
    pad_x = unit(0.01, "cm"), pad_y = unit(0.15, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "grey80"),
      line_col = "grey10",
      text_family = "Arial"
    )) +
  guides(fill = guide_colorsteps(barwidth = 20,
                                 barheight = .65,
                                 title.position = "right")) +
  labs(fill = "m") +
  coord_sf(xlim = c(-122.17, -121.765), y = c(45.41, 45.53)) +
  theme(axis.text = element_text(size=24)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "bullrun_map.png", width = 30, height = 20, units = "cm", dpi=500)

#Map of large fire in Clackamas watershed

#Download DEM for upper Clackamas River extent, convert to raster and then a dataframe, estimate slope, aspect, orientation, hillshade, and combine those effects.
fireyear_extent <- data.frame(x = c(-122.26, -121.99), y = c(45.249, 44.835))
fireyear_extent
prj_dd <- "EPSG:4326"
fireyear_dem <- get_elev_raster(fireyear_extent, prj = prj_dd, z = 10)
fireyear_rast <- rast(fireyear_dem) 
fireyear_rast_df <- as.data.frame(fireyear_rast, xy = TRUE)
names(fireyear_rast_df)[3] <- "elev_m"
slope <- terrain(fireyear_rast, "slope", unit="radians")
aspect <- terrain(fireyear_rast, "aspect", unit="radians")
fireyear_hill <- shade(slope, aspect, 
                       angle = 45, 
                       direction = 300,
                       normalize= TRUE)
fireyear_hill_df <- as.data.frame(fireyear_hill, xy=TRUE)
head(fireyear_hill_df)

#Clackamas map
clackamas_map <- ggplot() +
  geom_raster(data = fireyear_hill_df,
              aes(x, y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_raster(data = fireyear_rast_df,
              aes(x, y, fill = elev_m),
              alpha = .25) +
  scale_fill_gradientn(colours = c("seashell", "cornsilk2", "darkkhaki", "darkolivegreen4", "darkolivegreen", "darkolivegreen", "khaki4", "navajowhite4", "navajowhite3", "whitesmoke", "grey100"), breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200)) +
  #geom_sf(data=study_area_nfs, color="black", linewidth=.55, fill=NA) +
  geom_sf(data = ugbs, color="black", fill="grey70", linewidth=.1, alpha=.5) +
  geom_sf(data=rivers, color="navy", linewidth=.65, alpha=.85) +
  geom_sf(data=lakes_res, fill="navy", color="navy", alpha=.95) +
  #geom_sf_text(data = rivers, aes(label = Name)) +
  #geom_sf_text(data = ugbs, aes(label = instName)) +
  geom_sf(data = large_fires, color="darkred", fill="darkred", linewidth=.05, alpha=.2) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude, color=huc8_name), shape=22, fill="black", size=8, stroke=5) +
  scale_color_manual(name = "", values=c("Middle Fork Willamette"="violetred1", "Clackamas"="orangered", "Lower Columbia-Sandy" = "green2", "North Santiam"="cyan2", "Mckenzie"="goldenrod1", "South Santiam"="darkorchid2"), guide = "none") +
  geom_text(data=site_coords, aes(x=longitude, y=latitude, label=site), color="black", size=9, hjust = 0.0, nudge_x = 0.009, vjust = 0.0, nudge_y = 0.009, family = "Arial", check_overlap = TRUE) +
  #geom_text(data=fireyear1499, aes(x=longitude, y=latitude, label=fire_year), color="firebrick4", size=3.25, hjust = 0.5, nudge_x = 0, vjust = 2.5, nudge_y = 0, family = "Arial", check_overlap = TRUE) +
  scale_x_continuous("", breaks = c(-122, -122.2)) +
  scale_y_continuous("", breaks = c(45, 45.2)) +
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey10", "white"),
    text_family = "Arial",
    text_cex = 2,
    line_width = 4,
    height = unit(0.55, "cm")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(2, "cm"), 
    width = unit(2, "cm"),
    pad_x = unit(0.01, "cm"), pad_y = unit(0.15, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "grey80"),
      line_col = "grey10",
      text_family = "Arial"
    )) +
  guides(fill = guide_colorsteps(barwidth = 20,
                                 barheight = .65,
                                 title.position = "right")) +
  labs(fill = "m") +
  coord_sf(xlim = c(-122.25, -121.99), y = c(45.249, 44.848)) +
  theme(axis.text = element_text(size=24)) +
  theme(legend.position="none")
ggsave(path="/YOURPATH", "clackamas_map.png", width = 20, height = 30, units = "cm", dpi=500)

#Make maps of potential large fires in the Middle Fork drainage (this map was not included in the final paper)

#Download DEM for Middle Fork (and south McKenzie) watershed extent, convert to raster and then a dataframe, estimate slope, aspect, orientation, hillshade, and combine those effects.
fireyear_extent <- data.frame(x = c(-122.65, -121.98), y = c(43.683, 43.968))
fireyear_extent
prj_dd <- "EPSG:4326"
fireyear_dem <- get_elev_raster(fireyear_extent, prj = prj_dd, z = 10)
fireyear_rast <- rast(fireyear_dem) 
fireyear_rast_df <- as.data.frame(fireyear_rast, xy = TRUE)
names(fireyear_rast_df)[3] <- "elev_m"
slope <- terrain(fireyear_rast, "slope", unit="radians")
aspect <- terrain(fireyear_rast, "aspect", unit="radians")
fireyear_hill <- shade(slope, aspect, 
                       angle = 45, 
                       direction = 300,
                       normalize= TRUE)
fireyear_hill_df <- as.data.frame(fireyear_hill, xy=TRUE)
head(fireyear_hill_df)

#Middle Fork map
middleforkmap <- ggplot() +
  geom_raster(data = fireyear_hill_df,
              aes(x, y, fill = hillshade),
              show.legend = FALSE) +
  scale_fill_distiller(palette = "Greys") +
  new_scale_fill() +
  geom_raster(data = fireyear_rast_df,
              aes(x, y, fill = elev_m),
              alpha = .25) +
  scale_fill_gradientn(colours = c("seashell", "cornsilk2", "darkkhaki", "darkolivegreen4", "darkolivegreen", "darkolivegreen", "khaki4", "navajowhite4", "navajowhite3", "whitesmoke", "grey100"), breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200)) +
  #geom_sf(data=study_area_nfs, color="black", linewidth=.55, fill=NA) +
  geom_sf(data = ugbs, color="black", fill="grey70", linewidth=.1, alpha=.5) +
  geom_sf(data=rivers, color="navy", linewidth=.35, alpha=.85) +
  geom_sf(data=lakes_res, fill="navy", color="navy", alpha=.95) +
  #geom_sf_text(data = rivers, aes(label = Name)) +
  #geom_sf_text(data = ugbs, aes(label = instName)) +
  geom_sf(data = large_fires, color="darkred", fill="darkred", linewidth=.05, alpha=.2) +
  geom_point(data=site_coords, aes(x=longitude, y=latitude, color=huc8_name), shape=22, fill="black", size=3.25, stroke=2.5) +
  scale_color_manual(name = "", values=c("Middle Fork Willamette"="violetred1", "Clackamas"="orangered", "Lower Columbia-Sandy" = "green2", "North Santiam"="cyan2", "Mckenzie"="goldenrod1", "South Santiam"="darkorchid2"), guide = "none") +
  #geom_text(data=site_coords, aes(x=longitude, y=latitude, label=site), color="black", size=5.5, hjust = 0.0, nudge_x = 0.005, vjust = 0.0, nudge_y = 0.005, family = "Arial", check_overlap = TRUE) +
  #geom_text(data=fireyear1783, aes(x=longitude, y=latitude, label=fire_year), color="firebrick4", size=3.25, hjust = 0.5, nudge_x = 0, vjust = 2.5, nudge_y = 0, family = "Arial", check_overlap = TRUE) +
  #geom_text(data=fireyear1768, aes(x=longitude, y=latitude, label=fire_year), color="firebrick4", size=3.25, hjust = 0.5, nudge_x = 0, vjust = 4.5, nudge_y = 0, family = "Arial", check_overlap = TRUE) +
  scale_x_continuous("", breaks = c(-122.5)) +
  scale_y_continuous("", breaks = c(43.75)) +
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey10", "white"),
    text_family = "Arial",
    text_cex = 1,
    line_width = .8,
    height = unit(0.35, "cm")
  ) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    height = unit(1.15, "cm"), 
    width = unit(1.15, "cm"),
    pad_x = unit(0.0, "cm"), pad_y = unit(0.15, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "grey80"),
      line_col = "grey10",
      text_family = "Arial"
    )) +
  guides(fill = guide_colorsteps(barwidth = 20,
                                 barheight = .5,
                                 title.position = "right")) +
  labs(fill = "m") +
  coord_sf(xlim = c(-122.65, -121.98), y = c(43.683, 43.968)) +
  theme(legend.position = "bottom", axis.text = element_text(size=12))
middleforkmap
ggsave(path="/YOURPATH", "middleforkmap.png", width = 25, height = 25, units = "cm", dpi=500)

### Other optional stuff which I didn't end up using but that may be helpful ###

#Create a centroid functon and calculate centroides by site
cntrd <- function(x) {
  data.frame(centroid(as.matrix(x[,c("longitude", "latitude")])))
}
xsections <- xsections %>% tidyr::drop_na(latitude)
xsections <- xsections %>% filter(!longitude=="gps")
xsections <- xsections %>% filter(!longitude=="")
xsections$latitude <- as.numeric(xsections$latitude)
xsections$longitude <- as.numeric(xsections$longitude)
site_coords <- xsections %>%
  dplyr::select(site, latitude, longitude) %>%
  group_by(site) %>%
  do(cntrd(.))
site_coords

#Alternative version of permutation procedure:
set.seed(123)
values <- c(2499,2499,2522,2522,2522,2522,2648,2648,2652,2652,2670,2670,2689,2689,2690,2690,2693,2693,2700,2700,2706,2706,2714,2714,2730,2730,2738,2738,2740,2740,2765,2765,2768,2768,2773,2773,2783,2783,2794,2794,2798,2798,2807,2807,2812,2812,2831,2831,2831,2835,2835,2836,2836,2836,2844,2844,2844,2846,2846,2846,2883,2883,2964,2964)
str(values)
groups <- 1:26
str(groups)

#The simpler permutation function
random_split <- function(values, groups, size.min, size.max) {
  R <- table(values)
  S <- tabulate(
    sample(
      rep(groups, size.max - size.min),
      length(values) - size.min*length(groups)
    ),
    length(groups)
  ) + size.min
  d <- length(R) - length(S)
  with(
    as_data_frame(
      sample_degseq( # randomly assign values to groups
        c(R, integer(max(0, -d))),
        c(sample(S), integer(max(0, d))),
        "simple.no.multiple.uniform"
      )
    ),
    split(as(names(R), class(values))[from], groups[to])
  )
}

#The simpler version of the permutation procedure
fire_years_dup <- fire_years_site %>% 
  group_by(fire_year) %>% 
  filter(n() > 1) %>%
  arrange(fire_year)
fire_years_dup

#Permutation function (reshuffle fire years among sites n times)
permute_function <- function(x, noperms){
  perm_fire <- function(){
    transform(x, random_fire=sample(fire_year)) %>%
      select(site, random_fire, longitude, latitude)
  }
  perm_list <- replicate(noperms, perm_fire(), simplify = FALSE )
  #all_perms <- dplyr::bind_rows(perm_list, .id = 'run')
  perm_list
}

#Run function on fire years duplicates file
permute_results <- permute_function(fire_years_dup, 5)
permute_results

#Split the permutation results into a list by fire year
perm_split <- lapply(permute_results, function(df) split(df, df$random_fire))

#Calculate pairwise geodesic distance matrices for each fire year within each permutation
perm_fire_dist <- lapply(perm_split, function(fireyear_list) {
  lapply(fireyear_list, function(df) {
    if (nrow(df) > 1) {
      geodist(df[, c("longitude", "latitude")], measure = "geodesic")
    } else {
      NA  
    }
  })
})

perm_fire_dist_df <- data.frame(unlist(perm_fire_dist))
perm_fire_dist_df <- perm_fire_dist_df %>% 
  rename("dist_m" = "unlist.perm_fire_dist.") %>% 
  filter(dist_m > 0)
perm_fire_dist_df$random_fire_year <- substr(rownames(perm_fire_dist_df), 1, 4)
head(perm_fire_dist_df)

#Summarize distances for each fire year
perm_fire_dist_summary <- perm_fire_dist_df %>%
  group_by(random_fire_year) %>%
  summarize(mean_dist_km = mean(dist_m)/1000)
perm_fire_dist_summary

#Plot distances
dist_box <- ggplot(perm_fire_dist_df, aes(x=random_fire_year, y=dist_m/1000)) +
  geom_boxplot(position="dodge", color="black", fill="steelblue4") + 
  #scale_fill_manual(name = "", values = c('Western States' = 'darkolivegreen4', 'NAFSN' = 'firebrick3', 'This study' = 'firebrick3')) +
  scale_x_discrete("") +
  scale_y_continuous("Distance (km)") +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) 
dist_box

#Calculate mean distance between sites
site_dist_matrix <- geodist(site_coords, measure = 'geodesic' )/1000 
colnames(site_dist_matrix) <- site_coords$site
rownames(site_dist_matrix) <- site_coords$site
site_dist_dist_df <- data.frame(unlist(site_dist_matrix))
mean(as.matrix(site_dist_dist_df))
site_dist_dist_df[site_dist_dist_df<1] <- NA
mean(as.matrix(site_dist_dist_df), na.rm=TRUE)
