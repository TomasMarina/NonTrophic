## Ecological networks of an Antarctic ecosystem: a full description of non-trophic interactions
## Date: September 2022
## Authors: Vanesa Salinas, Tomás Ignacio Marina, Georgina Cordone, Fernando Momo

# 1. DATA LOADING


# Load pkgs ----

packages <- c("dplyr", "reshape2", "igraph")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load data ----

trophic <- read.csv("data/Red_Trofica.csv")
mutualistic <- read.csv("data/Red_Mutualista.csv")
competitive_pred <- read.csv("data/Red_Competencia.csv")
competitive_res <- read.csv("data/edgelist_competitive_resources.csv")
commensalistic <- read.csv("data/edgelist_comensalism.csv")
amensalistic <- read.csv("data/edgelist_amensalism.csv")


# Convert to igraph ----

## Trophic ----
troph <- trophic[-c(1)]
adj_troph <- as.matrix(troph[, -1])
row.names(adj_troph) <- troph[, 1]
edge_troph <- melt(adj_troph) %>% 
  filter(value > 0)
g_troph <- graph_from_edgelist(as.matrix(edge_troph[,1:2]), directed = TRUE)

## Mutualistic ----
adj_mut <- as.matrix(mutualistic[, -1])
row.names(adj_mut) <- mutualistic[, 1]
edge_mut <- melt(adj_mut) %>% 
  filter(value > 0)
g_mut <- graph_from_edgelist(as.matrix(edge_mut[,1:2]), directed = FALSE)

## Competitive ----
# Only btw predators sharing prey
adj_comp_pred <- as.matrix(competitive_pred[, -1])
row.names(adj_comp_pred) <- competitive_pred[, 1]
edge_comp_pred <- melt(adj_comp_pred) %>% 
  filter(value > 0)
g_comp_pred <- graph_from_edgelist(as.matrix(edge_comp_pred[,1:2]), directed = FALSE)
# Convert and save as edge list
edgelist_comp_pred <- as_edgelist(g_comp_pred, names = TRUE)
edgelist_comp_pred_df <- data.frame(edgelist_comp_pred)
colnames(edgelist_comp_pred_df) <- c("Predator_1","Predator_2")
#write.csv(edgelist_comp_pred_df, file = "data/edgelist_competitive_predators.csv")

# Bind competitive lists (predators and resources)
edgelist_comp_pred_ok <- edgelist_comp_pred_df %>% 
  mutate(TypeCompetition = "Prey") %>% 
  rename(Competitor_1 = Predator_1, Competitor_2 = Predator_2)
edgelist_comp_res_ok <- competitive_res %>% 
  mutate(TypeCompetition = "OtherResources") %>% 
  dplyr::select(Competitor_1, Competitor_2, TypeCompetition)
edgelist_comp_final <- bind_rows(edgelist_comp_pred_ok, edgelist_comp_res_ok)

g_comp <- graph_from_edgelist(as.matrix(edgelist_comp_final[,1:2]), directed = FALSE)

## Commensalistic ----
g_com <- graph_from_edgelist(as.matrix(commensalistic[,1:2]), directed = FALSE)

## Amensalistic ----
g_am <- graph_from_edgelist(as.matrix(amensalistic[,1:2]), directed = FALSE)

## Multiplex ---
multi_troph <- edge_troph %>% 
  rename(from = Var1, to = Var2, layer = value) %>% 
  mutate(layer = "trophic")
multi_mut <- edge_mut %>% 
  rename(from = Var1, to = Var2, layer = value) %>% 
  mutate(layer = "mutualistic")
multi_comp <- edgelist_comp_final %>% 
  dplyr::select(Competitor_1, Competitor_2) %>% 
  rename(from = Competitor_1, to = Competitor_2) %>% 
  mutate(layer = "competitive")
multi_com <- commensalistic %>% 
  rename(from = from, to = to) %>% 
  mutate(layer = "commensalism")
multi_am <- amensalistic %>% 
  rename(from = from, to = to) %>% 
  mutate(layer = "amensalism")
edge_multi <- bind_rows(multi_troph, multi_mut, multi_comp, multi_com, multi_am)

g_multi <- graph_from_edgelist(as.matrix(edge_multi[,1:2]), directed = FALSE)
# Add interaction type as edge attribute
g_multi <- g_multi %>%
  set_edge_attr(., name = 'type', index = E(g_multi), value = edge_multi[,3])
edge_attr(g_multi)


# Save data ----

save(g_troph, g_mut, g_comp, g_com, g_am, edge_multi, g_multi,
     file = "data/igraph_multiplex_data.rda")
