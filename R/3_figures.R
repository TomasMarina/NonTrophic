## Ecological networks of an Antarctic ecosystem: a full description of non-trophic interactions
## Date: September 2022
## Authors: Vanesa Salinas, Tomás Ignacio Marina, Georgina Cordone, Fernando Momo

# 3. FIGURES


# Load pkgs ----

packages <- c("dplyr", "tidyr", "ggplot2", "scales", "ggthemes")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load data ----

load("data/igraph_multiplex_data.rda")
load("results/complexity_&_sppattr.rda")


# Figure 5 ----
# Caption: Distribution of trophic and non-trophic interactions among basal (without prey), intermediate (with prey and predator) 
# and top (without predator) species. Interaction: amensalistic (-/0), commensalistic (+/0), competitive (-/-), mutualistic (+/+) 
# and trophic (+/-).

sp_attr_all

troph_list <- sp_attr_all %>% 
  dplyr::select(TrophicSpecies, Degree.troph, TrophicLevel, Type) %>% 
  slice(rep(seq_len(n()), Degree.troph)) %>% 
  dplyr::select(-Degree.troph) %>% 
  mutate(Interaction = "trophic")
mut_list <- sp_attr_all %>% 
  dplyr::select(TrophicSpecies, Degree.mut, TrophicLevel, Type) %>% 
  slice(rep(seq_len(n()), Degree.mut)) %>% 
  dplyr::select(-Degree.mut) %>% 
  mutate(Interaction = "mutualistic")
comp_list <- sp_attr_all %>% 
  dplyr::select(TrophicSpecies, Degree.comp, TrophicLevel, Type) %>% 
  slice(rep(seq_len(n()), Degree.comp)) %>% 
  dplyr::select(-Degree.comp) %>% 
  mutate(Interaction = "competitive")
com_list <- sp_attr_all %>% 
  dplyr::select(TrophicSpecies, Degree.com, TrophicLevel, Type) %>% 
  slice(rep(seq_len(n()), Degree.com)) %>% 
  dplyr::select(-Degree.com) %>% 
  mutate(Interaction = "commensalistic")
am_list <- sp_attr_all %>% 
  dplyr::select(TrophicSpecies, Degree.am, TrophicLevel, Type) %>% 
  slice(rep(seq_len(n()), Degree.am)) %>% 
  dplyr::select(-Degree.am) %>% 
  mutate(Interaction = "amensalistic")

all_int <- bind_rows(troph_list, mut_list, comp_list, com_list, am_list) %>% 
  group_by(TrophicSpecies, Type, Interaction) %>% 
  summarise(n = n())

# RGB code
# Trophic [255, 0, 255] (+/-)
# Mutualistic [255,0, 51] (+/+)
# Commensalistic [255, 204, 102] (+/0)
# Amensalistic [102, 255, 0] (-/0)
# Competitive [0, 51, 204] (-/-)

troph_col <- rgb(255,0,255, maxColorValue = 255)
mut_col <- rgb(255,0,51, maxColorValue = 255)
com_col <- rgb(255,204,102, maxColorValue = 255)
am_col <- rgb(102,255,0, maxColorValue = 255)
comp_col <- rgb(0,51,204, maxColorValue = 255)

fig5 <- all_int %>% 
  mutate(Interaction = factor(Interaction, levels = c("trophic", "mutualistic", "commensalistic",
                                   "amensalistic", "competitive"))) %>% 
  ggplot(aes(x = Type, y = n, colour = Interaction, fill = Interaction)) +
  geom_boxplot(alpha = 0.2) +
  scale_colour_manual(labels = c("(+/-)", "(+/+)", "(+/0)", "(-/0)", "(-/-)"), 
                      values = c(troph_col, mut_col, com_col, am_col, comp_col)) +
  scale_fill_manual(labels = c("(+/-)", "(+/+)", "(+/0)", "(-/0)", "(-/-)"), 
                      values = c(troph_col, mut_col, com_col, am_col, comp_col)) +
  scale_x_discrete(labels = c("Basal", "Intermediate", "Top")) +
  labs(x = "Species", y = "Number of interactions", colour = "Interaction") +
  theme_classic()
fig5


# Save results ----
save(all_int, fig5,
     file = "results/figure5.rda")
