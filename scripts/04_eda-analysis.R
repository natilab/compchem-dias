## EDA Plots and analysis
# Analysing Energy Decomposition Analysis results

# Generate plots from processed computational EDA results and analyze results

# author: Natalia Labadie
#--------------------------------------------------------------------

# This scripts generates all plots for the reactions studied, for use in PhD thesis

#--------------------------------------------------------------------

# Libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggsci)
library(ggrepel)
library(cowplot)

options(OutDec= ",", scipen = 999)


# --------------------------------------
## Import and Clean Data

datos_eda_gp <- read_excel("data/allenyl-eda-results.xlsx", 
                           sheet = "EDA-geom23A-gasphase-multiwfn",
                           skip = 3) %>% 
  select(!...10) %>% 
  rename(EDA_Int = `Delta E Total Int`,
         EOrb = `Delta E Orbital`,ESteric = `Delta E Steric`) %>% 
  mutate(sub_letter = if_else(substrate %in% c("AllenylBPin", "VinylBPin"), "a", "b"),
         prod_lab = paste0("3.", str_extract(reaction, "\\d"), sub_letter,
                           str_extract(reaction, "(?<=\\d)\\w")),
         substrate = case_when(substrate == "AllenylBPin" ~ "3.4a",
                               substrate == "AllenylCO2Me" ~ "3.4b",
                               substrate == "VinylBPin" ~ "3.1a",
                               substrate == "Methylacrylate" ~ "3.1b"))

datos_eda_gp_long <- datos_eda_gp %>% 
  pivot_longer(cols = c(EOrb, ESteric, EDistT, EDist1, EDist2, EInt), 
               names_to = "energy_type",
               values_to = "energy_value")


# --------------------------------------
## Analyse EDA and make plot for all systems

lm_orb_int_all <- lm(EDA_Int ~ energy_value,
                     filter(datos_eda_gp_long,
                            energy_type %in% c("EOrb")))

cor_orb_all_data <- filter(datos_eda_gp_long,
                           energy_type %in% c("EOrb"))

qqnorm(cor_orb_all_data$EDA_Int); qqline(cor_orb_all_data$EDA_Int, col = 2) # feito en valores pequeños
shapiro.test(cor_orb_all_data$EDA_Int) # pvalue 0.5418

qqnorm(cor_orb_all_data$energy_value); qqline(cor_orb_all_data$energy_value, col = 2) # en extremos se aleja bastante de la recta
shapiro.test(cor_orb_all_data$energy_value) # pvalue 0.5215

cor_orb_all <- cor.test(cor_orb_all_data$EDA_Int, cor_orb_all_data$energy_value)
# r = 0.9369614 
# p-value = 7.049e-06
pvalue_orb <- cor_orb_all$p.value
r_orb <- cor_orb_all$estimate[[1]]

plot1 <- ggplot(data = filter(datos_eda_gp_long,
                              energy_type %in% c("EOrb")),
                aes(x = energy_value, y = EDA_Int)) +
  geom_smooth(method = "lm", se = FALSE, 
              show.legend = FALSE, col = "black") +
  geom_point(aes(color = substrate), size = 3) + 
  geom_text_repel(aes(label=prod_lab), size = 5) +
  ylab(expression(italic(paste(Delta, "E"["int"]," / kcal ", "mol"^-1)))) +
  xlab(expression(italic(paste(Delta, "E"["orb"]," / kcal ", "mol"^-1)))) +
  scale_y_continuous(limits = c(-13, 0)) +
  scale_color_d3() +
  labs(color="sustrato") +
  guides(colour = guide_legend(nrow = 1)) +
  geom_label(data = data.frame(y = -13, x = -41),
             aes(x, y, 
                 label = paste0("r = ", round(r_orb, 4),
                                " (p-value = ", 
                                round(pvalue_orb, 5), ")")),
             size = 6, color = "black", fill = "white", 
             fontface = "bold") + 
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.35, 0.9),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",
                                             linetype = "dashed"),
        legend.text = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16))



# --------------------------------------
## Analyse EDA and make plot for only allenes

lm_orb_int_allene <- lm(EDA_Int ~ energy_value,
                        filter(datos_eda_gp_long,
                               energy_type %in% c("EOrb"),
                               substrate %in% c("3.4a", "3.4b")))

cor_orb_allene_data <- filter(datos_eda_gp_long,
                              energy_type %in% c("EOrb"),
                              substrate %in% c("3.4a", "3.4b"))

qqnorm(cor_orb_allene_data$EDA_Int); qqline(cor_orb_allene_data$EDA_Int, col = 2)
shapiro.test(cor_orb_allene_data$EDA_Int) # pvalue 0.7212

qqnorm(cor_orb_allene_data$energy_value); qqline(cor_orb_allene_data$energy_value, col = 2)
shapiro.test(cor_orb_allene_data$energy_value) # pvalue 0.732

cor_orb_allene <- cor.test(cor_orb_allene_data$EDA_Int, cor_orb_allene_data$energy_value)
pvalue_orb_allene <- cor_orb_allene$p.value
r_orb_allene <- cor_orb_allene$estimate[[1]]
# r = 0.9784529
# p-value = 2.461e-05


plot3 <- ggplot(data = filter(datos_eda_gp_long,
                              energy_type %in% c("EOrb"),
                              substrate %in% c("3.4a", "3.4b")),
                aes(x = energy_value, y = EDA_Int)) +
  geom_smooth(method = "lm", se = FALSE, 
              show.legend = FALSE, col = "black") +
  geom_point(aes(color = substrate), size = 3) + 
  geom_text_repel(aes(label=prod_lab), size = 5) +
  ylab(expression(italic(paste(Delta, "E"["int"]," / kcal ", "mol"^-1)))) +
  xlab(expression(italic(paste(Delta, "E"["orb"]," / kcal ", "mol"^-1)))) +
  scale_y_continuous(limits = c(-13, 0)) +
  scale_color_d3() +
  labs(color="sustrato") +
  guides(colour = guide_legend(nrow = 1)) +
  geom_label(data = data.frame(y = -13, x = -41),
             aes(x, y, 
                 label = paste0("r = ", round(r_orb_allene, 4),
                                " (p-value = ",
                                round(pvalue_orb_allene, 5), ")")),
             size = 6, color = "black", 
             fill = "white", fontface = "bold") + 
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.9),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", linetype = "dashed"),
        legend.text = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16))



# --------------------------------------
## Make combined plot

top_row <- plot_grid(plot1, plot3, labels = c('A', 'B'), 
                     label_size = 14)

bottom_row <- plot_grid(NULL, plot2, NULL, labels = c('', 'C', ''),
                        label_size = 14, 
                        rel_widths = c(0.5, 1, 0.5),
                        nrow = 1)
plot_final <-  plot_grid(top_row, bottom_row, #labels = NULL, 
                         ncol = 1)

ggsave("out_plots/scatter_eda_tesis.png", plot = plot_final, width = 4500, height = 4000, units = c("px"))
# no me funciona bien, uso opciones de RStudio

