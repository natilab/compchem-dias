## DIAS Processing
# Processing of Distortion-Interaction Activation Strain Model data.

# STEP 3: 
# Generate DIAS plots (distortion-interaction diagrams)

# author: Natalia Labadie
#--------------------------------------------------------------------

# This scripts generates all plots for the reactions studied, for use in PhD thesis

#--------------------------------------------------------------------

# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(scales)
library(transformr)
library(cowplot)

options(OutDec= ",") # for THESIS, comma as decimal separator



# --------------------------------------
## Import Data

dias_data <- readRDS("data/dias-2022_alldata_clean.rds")



# --------------------------------------
## AUX FUNCTIONS
# --------------------------------------

# Functions to filter and tidy data

filter_data_total <- function(data, selec_rxns, selec_subs, selec_lev, selec_sv) {
  "Filters whole dataset and returns selected reactions, substrates, 
  lev_theory and solvents. Selects total energies with variables of interest and orders in long format."
  
  filtered_data <- filter(data,
                          reaction %in% selec_rxns &
                            substrate %in% selec_subs &
                            level_theory %in% selec_lev &
                            solvent %in% selec_sv) %>% 
    select(DistanceFinal, substrate, reaction, Etotal, EInt, EDistT, 
           Ets) %>% 
    mutate(ts = !is.na(Ets)) %>% 
    pivot_longer(cols = c(Etotal, EInt, EDistT), names_to = "type", 
                 values_to = "energy") %>% 
    mutate(ts_energy = ifelse(ts, energy, NA)) %>% 
    select(-Ets) %>% 
    filter(!is.na(energy))
  
  return(filtered_data)
}

filter_data_dist <- function(data, selec_rxns, selec_subs, selec_lev, selec_sv) {
  "Filters whole dataset and returns selected reactions, substrates, 
  lev_theory and solvents. Selects fragment distortion energies and variables 
  of interest and orders in long format."
  
  filtered_data <- filter(data,
                          reaction %in% selec_rxns &
                            substrate %in% selec_subs &
                            level_theory %in% selec_lev &
                            solvent %in% selec_sv) %>% 
    select(DistanceFinal, substrate, reaction, Etotal, EDist1, EDist2, 
           Ets) %>% 
    mutate(ts = !is.na(Ets)) %>% 
    pivot_longer(cols = c(Etotal, EDist1, EDist2), names_to = "type", 
                 values_to = "energy") %>% 
    mutate(ts_energy = ifelse(ts, energy, NA)) %>% 
    select(-Ets) %>% 
    filter(!is.na(energy))
  
  return(filtered_data)
}


# --------------------------------------
## PLOTTING FUNCTIONS
# --------------------------------------

# Function for plotting total energies.

plot_total_energies_rxn <- function(data, selec_rxns, selec_subs, 
                                    selec_lev, selec_sv, 
                                    colorpalette, energy_labels,
                                    textsizes, legkeysize,
                                    xtitle, ytitle, 
                                    rxn_labels,
                                    linesize = 0.5, xlims = c(3.4, 1.9),
                                    xbreaks = seq(3.4, 2.0, -0.2),
                                    ylims = c(-75, 75),
                                    ybreaks = seq(-60, 60, 20)) {
  "Filters data and returns plot with total energies.
  colorpalette: vector with colors for total energy, dist energy and int energy.
  energy_labels: vector with labels for total energy, dist energy and int energy.
  textsizes is named vector wih names legend, xtitle, ytitle, axis.
  legkeysize is unit() element."
  
  filtered_data <- filter_data_total(data, selec_rxns, selec_subs, selec_lev, selec_sv)
  
  plot <- filtered_data %>% 
    filter(between(DistanceFinal, xlims[2], xlims[1])) %>% 
    ggplot(aes(x = DistanceFinal, y = energy, 
               group = interaction(reaction, type))) +
    geom_line(aes(y = 0), size = 0.5, linetype = 4, colour = "gray") +
    geom_line(aes(color = type, linetype = reaction)) +
    geom_point(aes(y = ts_energy, color = type)) +
    scale_x_reverse(xtitle, expand=c(0,0), limits = xlims, 
                    breaks = xbreaks) + 
    scale_y_continuous(ytitle, expand=c(0,0), limits = ylims, 
                       breaks = ybreaks) +
    theme_classic() +
    scale_color_manual(values = colorpalette, breaks = c("Etotal", "EDistT", "EInt"), 
                       labels = energy_labels) +
    scale_linetype_manual(values = c(1, 2),
                          breaks = selec_rxns,
                          labels = rxn_labels) +
    
    guides(color = "none") +
    theme(legend.title = element_blank(), 
          legend.key.size = legkeysize,
          legend.text = element_markdown(size = textsizes["legend"], hjust = 0),
          legend.position = c(0.3, 0.8),
          axis.title.x = element_text(size = textsizes["xtitle"], 
                                      hjust = 0.5, face = "bold.italic"),
          axis.title.y = element_text(size = textsizes["ytitle"], 
                                      hjust = 0.5),
          axis.text = element_text(size = textsizes["axis"])) +
    annotate("text", x = 2.35, y = -30, label = energy_labels[3], color = colorpalette[3],
             size = 6) +
    annotate("text", x = 2.1, y = 25, label = energy_labels[1], color = colorpalette[1],
             size = 6) +
    annotate("text", x = 2.35, y = 50, label = energy_labels[2], color = colorpalette[2],
             size = 6) 
  
  return(plot)
}



# Function for plotting distortion energies

plot_dist_energies_rxn <- function(data, selec_rxns, selec_subs, 
                                   selec_lev, selec_sv, 
                                   colorpalette, energy_labels,
                                   textsizes, legkeysize,
                                   xtitle, ytitle, 
                                   rxn_labels,
                                   linesize = 0.5, xlims = c(3.4, 1.9),
                                   xbreaks = seq(3.4, 2.0, -0.2),
                                   ylims = c(-20, 60),
                                   ybreaks = seq(0, 60, 20)) {
  "Filters data and returns plot with total energies.
  colorpalette: vector with colors for total energy, dist energy of frag1 and dist energy of frag2.
  energy_labels: vector with labels for total energy, dist energy of frag1 and dist energy of frag2.
  textsizes is named vector wih names legend, xtitle, ytitle, axis.
  legkeysize is unit() element."
  
  filtered_data <- filter_data_dist(data, selec_rxns, selec_subs, selec_lev, selec_sv)
  
  plot <- filtered_data %>% 
    filter(between(DistanceFinal, xlims[2], xlims[1])) %>% 
    ggplot(aes(x = DistanceFinal, y = energy, 
               group = interaction(reaction, type))) +
    geom_line(aes(y = 0), size = 0.5, linetype = 4, colour = "gray") +
    geom_line(aes(color = type, linetype = reaction)) +
    geom_point(aes(y = ts_energy, color = type)) +
    scale_x_reverse(xtitle, expand=c(0,0), limits = xlims, 
                    breaks = xbreaks) + 
    scale_y_continuous(ytitle, expand=c(0,0), limits = ylims, 
                       breaks = ybreaks) +
    theme_classic() +
    scale_color_manual(values = colorpalette, breaks = c("Etotal", "EDist1", "EDist2"), 
                       labels = energy_labels) +
    scale_linetype_manual(values = c(1, 2),
                          breaks = selec_rxns,
                          labels = rxn_labels) +
    
    guides(color = "none") +
    theme(legend.title = element_blank(), 
          legend.key.size = legkeysize,
          legend.text = element_markdown(size = textsizes["legend"], hjust = 0),
          legend.position = c(0.3, 0.8),
          axis.title.x = element_text(size = textsizes["xtitle"], 
                                      hjust = 0.5, face = "bold.italic"),
          axis.title.y = element_text(size = textsizes["ytitle"], 
                                      hjust = 0.5),
          axis.text = element_text(size = textsizes["axis"])) +
    annotate("text", x = 2.3, y = 30, label = energy_labels[2], color = colorpalette[2],
             size = 6) +
    annotate("text", x = 2.6, y = 18, label = energy_labels[1], color = colorpalette[1],
             size = 6) +
    annotate("text", x = 2.23, y = 2, label = energy_labels[3], color = colorpalette[3],
             size = 6) 
  
  return(plot)
  
}


plot_total_energies_sub <- function(data, selec_rxns, selec_subs, 
                                    selec_lev, selec_sv, 
                                    colorpalette, energy_labels,
                                    textsizes, legkeysize,
                                    xtitle, ytitle, 
                                    rxn_labels,
                                    linesize = 0.5, xlims = c(3.4, 1.9),
                                    xbreaks = seq(3.4, 2.0, -0.2),
                                    ylims = c(-75, 75),
                                    ybreaks = seq(-60, 60, 20)) {
  "Filters data and returns plot with total energies.
  colorpalette: vector with colors for total energy, dist energy and int energy.
  energy_labels: vector with labels for total energy, dist energy and int energy.
  textsizes is named vector wih names legend, xtitle, ytitle, axis.
  legkeysize is unit() element."
  
  filtered_data <- filter_data_total(data, selec_rxns, selec_subs, selec_lev, selec_sv)
  
  plot <- filtered_data %>% 
    filter(between(DistanceFinal, xlims[2], xlims[1])) %>% 
    ggplot(aes(x = DistanceFinal, y = energy, 
               group = interaction(substrate, type))) +
    geom_line(aes(y = 0), size = 0.5, linetype = 4, colour = "gray") +
    geom_line(aes(color = type, linetype = substrate)) +
    geom_point(aes(y = ts_energy, color = type)) +
    scale_x_reverse(xtitle, expand=c(0,0), limits = xlims, 
                    breaks = xbreaks) + 
    scale_y_continuous(ytitle, expand=c(0,0), limits = ylims, 
                       breaks = ybreaks) +
    theme_classic() +
    scale_color_manual(values = colorpalette, breaks = c("Etotal", "EDistT", "EInt"), 
                       labels = energy_labels) +
    scale_linetype_manual(values = c(1, 2),
                          breaks = selec_subs,
                          labels = rxn_labels) +
    
    guides(color = "none") +
    theme(legend.title = element_blank(), 
          legend.key.size = legkeysize,
          legend.text = element_markdown(size = textsizes["legend"], hjust = 0),
          legend.position = c(0.3, 0.8),
          axis.title.x = element_text(size = textsizes["xtitle"], 
                                      hjust = 0.5, face = "bold.italic"),
          axis.title.y = element_text(size = textsizes["ytitle"], 
                                      hjust = 0.5),
          axis.text = element_text(size = textsizes["axis"])) +
    annotate("text", x = 2.35, y = -30, label = energy_labels[3], color = colorpalette[3],
             size = 6) +
    annotate("text", x = 2.1, y = 25, label = energy_labels[1], color = colorpalette[1],
             size = 6) +
    annotate("text", x = 2.35, y = 50, label = energy_labels[2], color = colorpalette[2],
             size = 6) 
  
  return(plot)
}



# Function for plotting distortion energies

plot_dist_energies_sub <- function(data, selec_rxns, selec_subs, 
                                   selec_lev, selec_sv, 
                                   colorpalette, energy_labels,
                                   textsizes, legkeysize,
                                   xtitle, ytitle, 
                                   rxn_labels,
                                   linesize = 0.5, xlims = c(3.4, 1.9),
                                   xbreaks = seq(3.4, 2.0, -0.2),
                                   ylims = c(-20, 60),
                                   ybreaks = seq(0, 60, 20)) {
  "Filters data and returns plot with total energies.
  colorpalette: vector with colors for total energy, dist energy of frag1 and dist energy of frag2.
  energy_labels: vector with labels for total energy, dist energy of frag1 and dist energy of frag2.
  textsizes is named vector wih names legend, xtitle, ytitle, axis.
  legkeysize is unit() element."
  
  filtered_data <- filter_data_dist(data, selec_rxns, selec_subs, selec_lev, selec_sv)
  
  plot <- filtered_data %>% 
    filter(between(DistanceFinal, xlims[2], xlims[1])) %>% 
    ggplot(aes(x = DistanceFinal, y = energy, 
               group = interaction(substrate, type))) +
    geom_line(aes(y = 0), size = 0.5, linetype = 4, colour = "gray") +
    geom_line(aes(color = type, linetype = substrate)) +
    geom_point(aes(y = ts_energy, color = type)) +
    scale_x_reverse(xtitle, expand=c(0,0), limits = xlims, 
                    breaks = xbreaks) + 
    scale_y_continuous(ytitle, expand=c(0,0), limits = ylims, 
                       breaks = ybreaks) +
    theme_classic() +
    scale_color_manual(values = colorpalette, breaks = c("Etotal", "EDist1", "EDist2"), 
                       labels = energy_labels) +
    scale_linetype_manual(values = c(1, 2),
                          breaks = selec_subs,
                          labels = rxn_labels) +
    
    guides(color = "none") +
    theme(legend.title = element_blank(), 
          legend.key.size = legkeysize,
          legend.text = element_markdown(size = textsizes["legend"], hjust = 0),
          legend.position = c(0.3, 0.8),
          axis.title.x = element_text(size = textsizes["xtitle"], 
                                      hjust = 0.5, face = "bold.italic"),
          axis.title.y = element_text(size = textsizes["ytitle"], 
                                      hjust = 0.5),
          axis.text = element_text(size = textsizes["axis"])) +
    annotate("text", x = 2.3, y = 30, label = energy_labels[2], color = colorpalette[2],
             size = 6) +
    annotate("text", x = 2.6, y = 18, label = energy_labels[1], color = colorpalette[1],
             size = 6) +
    annotate("text", x = 2.23, y = 2, label = energy_labels[3], color = colorpalette[3],
             size = 6) 
  
  return(plot)
  
}

plot_both_rxn <- function(data, selec_rxns, selec_subs, 
                          selec_lev, selec_sv, 
                          textsizes, legkeysize,
                          xtitle, ytitle, 
                          rxn_labels,
                          linesize = 0.5, xlims = c(3.4, 1.9),
                          xbreaks = seq(3.4, 2.0, -0.2),
                          colorpalette_tot, energy_labels_tot,
                          ylims_tot = c(-75, 75),
                          ybreaks_tot = seq(-60, 60, 20),
                          colorpalette_dist, energy_labels_dist,
                          ylims_dist = c(-20, 60),
                          ybreaks_dist = seq(0, 60, 20)) {
  
  plot1 <- plot_total_energies_rxn(data, selec_rxns = selec_rxns, selec_subs = selec_subs,
                                   selec_lev = selec_lev, selec_sv = selec_sv, 
                                   textsizes = textsizes, legkeysize = legkeysize,
                                   xtitle = xtitle, ytitle = ytitle, 
                                   rxn_labels = rxn_labels, linesize = linesize,
                                   xlims = xlims, xbreaks = xbreaks, 
                                   colorpalette = colorpalette_tot, 
                                   energy_labels = energy_labels_tot,
                                   ylims = ylims_tot, ybreaks = ybreaks_tot)
  
  plot2 <- plot_dist_energies_rxn(data, selec_rxns = selec_rxns, selec_subs = selec_subs,
                                  selec_lev = selec_lev, selec_sv = selec_sv, 
                                  textsizes = textsizes, legkeysize = legkeysize,
                                  xtitle = xtitle, ytitle = ytitle, 
                                  rxn_labels = rxn_labels, linesize = linesize,
                                  xlims = xlims, xbreaks = xbreaks, 
                                  colorpalette = colorpalette_dist, 
                                  energy_labels = energy_labels_dist,
                                  ylims = ylims_dist, ybreaks = ybreaks_dist)
  
  return(plot_grid(plot1, plot2, labels = c('A', 'B'), label_size = 20))
}


plot_both_sub <- function(data, selec_rxns, selec_subs, 
                          selec_lev, selec_sv, 
                          textsizes, legkeysize,
                          xtitle, ytitle, 
                          rxn_labels,
                          linesize = 0.5, xlims = c(3.4, 1.9),
                          xbreaks = seq(3.4, 2.0, -0.2),
                          colorpalette_tot, energy_labels_tot,
                          ylims_tot = c(-75, 75),
                          ybreaks_tot = seq(-60, 60, 20),
                          colorpalette_dist, energy_labels_dist,
                          ylims_dist = c(-20, 60),
                          ybreaks_dist = seq(0, 60, 20)) {
  
  plot1 <- plot_total_energies_sub(data, selec_rxns = selec_rxns, selec_subs = selec_subs,
                                   selec_lev = selec_lev, selec_sv = selec_sv, 
                                   textsizes = textsizes, legkeysize = legkeysize,
                                   xtitle = xtitle, ytitle = ytitle, 
                                   rxn_labels = rxn_labels, linesize = linesize,
                                   xlims = xlims, xbreaks = xbreaks, 
                                   colorpalette = colorpalette_tot, 
                                   energy_labels = energy_labels_tot,
                                   ylims = ylims_tot, ybreaks = ybreaks_tot)
  
  plot2 <- plot_dist_energies_sub(data, selec_rxns = selec_rxns, selec_subs = selec_subs,
                                  selec_lev = selec_lev, selec_sv = selec_sv, 
                                  textsizes = textsizes, legkeysize = legkeysize,
                                  xtitle = xtitle, ytitle = ytitle, 
                                  rxn_labels = rxn_labels, linesize = linesize,
                                  xlims = xlims, xbreaks = xbreaks, 
                                  colorpalette = colorpalette_dist, 
                                  energy_labels = energy_labels_dist,
                                  ylims = ylims_dist, ybreaks = ybreaks_dist)
  
  return(plot_grid(plot1, plot2, labels = c('A', 'B'), label_size = 20))
}


#------------------------------------------------------------------------------------
# GENERATE ALL PLOTS
#-----------------------------------------------------------------------------------


# formatos generales
legkeysize <- unit(0.5, "inch")
textsizes <- c(legend = 13, xtitle = 14, ytitle = 14, axis = 12)
linesize <- 0.8

xtitle <- expression(bolditalic(paste("r (C\U00B7\U00B7\U00B7",  "C) / ", ring(A), sep = "")))
ytitle <- expression(bolditalic(paste(Delta, "E / kcal ", " mol "^-1)))

totE <- expression(italic(paste(Delta,"E")))
distE <- expression(italic(paste(Delta, "E"["dist"]))) 
intE <- expression(italic(paste(Delta, "E"["int"]))) 
energy_labels_tot <- c(totE, distE, intE)

dist1 <- expression(italic(paste(Delta, "E"["dist, Cp"]))) 
dist2 <- expression(italic(paste(Delta, "E"["dist, Df"]))) 
energy_labels_dist <- c(totE, dist1, dist2)

colorpalette <- c("#000000", "#0000ff", "#ff0000")
colorpalette_dist <- c("#000000", "#880044", "#157F1F")

arrow <- "\u2192"


# args de cada grafico

args1 <- list(rxn_labels = c(paste("**3.4a** + **3.2** ", arrow, " **3.5aN**"),
                             paste("**3.4a** + **3.2** ", arrow, " **3.6aE**")),
              selec_rxns = c("Product 5N", "Product 6E"),
              selec_subs = c("AllenylBPin"))


args2 <- list(rxn_labels = c(paste("**3.4a** + **3.2** ", arrow, " **3.5aX**"),
                             paste("**3.4a** + **3.2** ", arrow, " **3.6aZ**")),
              selec_rxns = c("Product 5X", "Product 6Z"),
              selec_subs = c("AllenylBPin"))

args3 <- list(rxn_labels = c(paste("**3.4b** + **3.2** ", arrow, " **3.5bN**"),
                             paste("**3.4b** + **3.2** ", arrow, " **3.6bE**")),
              selec_rxns = c("Product 5N", "Product 6E"),
              selec_subs = c("AllenylCO2Me"))


args4 <- list(rxn_labels = c(paste("**3.4b** + **3.2** ", arrow, " **3.5bX**"),
                             paste("**3.4b** + **3.2** ", arrow, " **3.6bZ**")),
              selec_rxns = c("Product 5X", "Product 6Z"),
              selec_subs = c("AllenylCO2Me"))


args5 <- list(rxn_labels = c(paste("**3.4a** + **3.2** ", arrow, " **3.5aN**"),
                             paste("**3.4b** + **3.2** ", arrow, " **3.5bN**")),
              selec_rxns = c("Product 5N"),
              selec_subs = c("AllenylBPin", "AllenylCO2Me"))

args6 <- list(rxn_labels = c(paste("**3.4a** + **3.2** ", arrow, " **3.5aX**"),
                             paste("**3.4b** + **3.2** ", arrow, " **3.5bX**")),
              selec_rxns = c("Product 5X"),
              selec_subs = c("AllenylBPin", "AllenylCO2Me"))

args7 <- list(rxn_labels = c(paste("**3.4a** + **3.2** ", arrow, " **3.6aE**"),
                             paste("**3.4b** + **3.2** ", arrow, " **3.6bE**")),
              selec_rxns = c("Product 6E"),
              selec_subs = c("AllenylBPin", "AllenylCO2Me"))

args8 <- list(rxn_labels = c(paste("**3.4a** + **3.2** ", arrow, " **3.6aZ**"),
                             paste("**3.4b** + **3.2** ", arrow, " **3.6bZ**")),
              selec_rxns = c("Product 6Z"),
              selec_subs = c("AllenylBPin", "AllenylCO2Me"))

args9 <- list(rxn_labels = c(paste("**3.1a** + **3.2** ", arrow, " **3.3aN**"),
                             paste("**3.1b** + **3.2** ", arrow, " **3.3bN**")),
              selec_rxns = c("Product 3N"),
              selec_subs = c("VinylBPin", "Methylacrylate"))

args10 <- list(rxn_labels = c(paste("**3.1a** + **3.2** ", arrow, " **3.3aX**"),
                              paste("**3.1b** + **3.2** ", arrow, " **3.3bX**")),
               selec_rxns = c("Product 3X"),
               selec_subs = c("VinylBPin", "Methylacrylate"))

args11 <- list(rxn_labels = c(paste("**3.4a** + **3.2** ", arrow, " **3.5aN**"),
                              paste("**3.1a** + **3.2** ", arrow, " **3.3aN**")),
               selec_rxns = c("Product 5N", "Product 3N"),
               selec_subs = c("AllenylBPin", "VinylBPin"))

args12 <- list(rxn_labels = c(paste("**3.4a** + **3.2** ", arrow, " **3.5aX**"),
                              paste("**3.1a** + **3.2** ", arrow, " **3.3aX**")),
               selec_rxns = c("Product 5X", "Product 3X"),
               selec_subs = c("AllenylBPin", "VinylBPin"))

args13 <- list(rxn_labels = c(paste("**3.4b** + **3.2** ", arrow, " **3.5bN**"),
                              paste("**3.1b** + **3.2** ", arrow, " **3.3bN**")),
               selec_rxns = c("Product 5N", "Product 3N"),
               selec_subs = c("AllenylCO2Me", "Methylacrylate"))

args14 <- list(rxn_labels = c(paste("**3.4b** + **3.2** ", arrow, " **3.5bX**"),
                              paste("**3.1b** + **3.2** ", arrow, " **3.3bX**")),
               selec_rxns = c("Product 5X", "Product 3X"),
               selec_subs = c("AllenylCO2Me", "Methylacrylate"))

args_rxn <- purrr::map(1:4, \(x) get(paste0("args", x)))

args_sub <- purrr::map(5:14, \(x) get(paste0("args", x)))

# for (i in 1:2) {
#   print(args[[i]]$rxn_labels)
# }

for (i in 1:length(args_rxn)) {
  plots <- plot_both_rxn(data = dias_data, selec_rxns = args_rxn[[i]]$selec_rxns, 
                         selec_subs = args_rxn[[i]]$selec_subs, 
                         selec_lev = c("M062X/6-31+G*"), selec_sv = c("toluene"), 
                         textsizes, legkeysize,
                         xtitle, ytitle, 
                         rxn_labels = args_rxn[[i]]$rxn_labels,
                         colorpalette_tot = colorpalette, energy_labels_tot = energy_labels_tot,
                         colorpalette_dist = colorpalette_dist, energy_labels_dist = energy_labels_dist)
  
  ggsave(paste0("out_plots/dias_plot_rxn_", i, ".png"), plots, height = 6,
         width = 12,
         units = "in",
         dpi = 600)
}

for (i in 1:length(args_sub)) {
  plots <- plot_both_sub(data = dias_data, selec_rxns = args_sub[[i]]$selec_rxns, 
                         selec_subs = args_sub[[i]]$selec_subs, 
                         selec_lev = c("M062X/6-31+G*"), selec_sv = c("toluene"), 
                         textsizes, legkeysize,
                         xtitle, ytitle, 
                         rxn_labels = args_sub[[i]]$rxn_labels,
                         colorpalette_tot = colorpalette, energy_labels_tot = energy_labels_tot,
                         colorpalette_dist = colorpalette_dist, energy_labels_dist = energy_labels_dist)
  
  ggsave(paste0("out_plots/dias_plot_sub_", i, ".png"), plots, height = 6,
         width = 12,
         units = "in",
         dpi = 600)
}

# EN WINDOWS QUEDAN BIEN LOS GRAFICOS PQ USA ARIAL PARA LAS EXPRESSIONS
