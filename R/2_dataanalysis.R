## Ecological networks of an Antarctic ecosystem: a full description of non-trophic interactions
## Date: September 2022
## Authors: Vanesa Salinas, Tomás Ignacio Marina, Georgina Cordone, Fernando Momo

# 1. DATA ANALYSIS


# Load pkgs ----

packages <- c("igraph", "tidyr", "multiweb", "NetIndices", "intergraph")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load data ----

load("data/igraph_multiplex_data.rda")


# Network analysis ----

## Complexity ----
g_list <- list(g_troph, g_mut, g_comp, g_com, g_am)
names(g_list) <- c("trophic", "mutualistic", "competitive", "commensalistic", "amensalistic")

complexity <- lapply(g_list, calc_topological_indices)
df_complex <- bind_rows(complexity) %>% 
  mutate(Network = c("Trophic", "Mutualistic", "Competitive", "Commensalistic", "Amensalistic")) %>% 
  dplyr::select(Network, Size, Links, LD, Connectance, Components)


## Degree ----
# Trophic
V(g_troph)$degree <- degree(g_troph, v = V(g_troph), mode = "total")
# Mutualistic
V(g_mut)$degree <- degree(g_mut, v = V(g_mut), mode = "total")
# Competitive
V(g_comp)$degree <- degree(g_comp, v = V(g_comp), mode = "total")
# Commensalistic
V(g_com)$degree <- degree(g_com, v = V(g_com), mode = "total")
# Amensalistic
V(g_am)$degree <- degree(g_am, v = V(g_am), mode = "total")
# Multiplex
V(g_multi)$degree <- degree(g_multi, v = V(g_multi), mode = "total")

## Trophic level ----
# Trophic
adj_troph <- get.adjacency(g_troph, sparse = FALSE)
tl_attr <- tibble(round(TrophInd(adj_troph, Dead=c("Necromass", "FreshDetritus", "AgedDetritus")), digits = 3)) %>% 
  mutate(id = V(g_troph)$name) %>% transform(TL = as.numeric(TL), OI = as.numeric(OI))
V(g_troph)$TL <- tl_attr$TL
# Mutualistic
df_mut <- igraph::as_data_frame(g_mut, 'both')
df_mut$vertices <- df_mut$vertices %>% 
  left_join(tl_attr, c('name' = 'id'))
g_mut <- graph_from_data_frame(df_mut$edges,
                               directed = F,
                               vertices = df_mut$vertices)
g_mut <- delete_vertex_attr(g_mut, "OI")  # discard omnivory
# Competitive
df_comp <- igraph::as_data_frame(g_comp, 'both')
df_comp$vertices <- df_comp$vertices %>% 
  left_join(tl_attr, c('name' = 'id'))
g_comp <- graph_from_data_frame(df_comp$edges,
                                directed = F,
                                vertices = df_comp$vertices)
g_comp <- delete_vertex_attr(g_comp, "OI")
# Commensalistic
df_com <- igraph::as_data_frame(g_com, 'both')
df_com$vertices <- df_com$vertices %>% 
  left_join(tl_attr, c('name' = 'id'))
g_com <- graph_from_data_frame(df_com$edges,
                               directed = F,
                               vertices = df_com$vertices)
g_com <- delete_vertex_attr(g_com, "OI")
# Amensalistic
df_am <- igraph::as_data_frame(g_am, 'both')
df_am$vertices <- df_am$vertices %>% 
  left_join(tl_attr, c('name' = 'id'))
g_am <- graph_from_data_frame(df_am$edges,
                              directed = F,
                              vertices = df_am$vertices)
g_am <- delete_vertex_attr(g_am, "OI")
# Multiplex
df_multi <- igraph::as_data_frame(g_multi, 'both')
df_multi$vertices <- df_multi$vertices %>% 
  left_join(tl_attr, c('name' = 'id'))
g_multi <- graph_from_data_frame(df_multi$edges,
                                 directed = F,
                                 vertices = df_multi$vertices)
g_multi <- delete_vertex_attr(g_multi, "OI")


# Species attributes ----

## By interaction ----
#Trophic
sp_attr_troph <- bind_cols(TrophicSpecies = V(g_troph)$name, Degree.troph = V(g_troph)$degree, 
                           TrophicLevel = V(g_troph)$TL, IntType = "trophic")
# Mutualistic
sp_attr_mut <- bind_cols(TrophicSpecies = V(g_mut)$name, Degree.mut = V(g_mut)$degree, 
                           TrophicLevel = V(g_mut)$TL, IntType = "mutualistic")
# Competitive
sp_attr_comp <- bind_cols(TrophicSpecies = V(g_comp)$name, Degree.comp = V(g_comp)$degree, 
                         TrophicLevel = V(g_comp)$TL, IntType = "competitive")
# Commensalistic
sp_attr_com <- bind_cols(TrophicSpecies = V(g_com)$name, Degree.com = V(g_com)$degree, 
                          TrophicLevel = V(g_com)$TL, IntType = "commensalistic")
# Amensalistic
sp_attr_am <- bind_cols(TrophicSpecies = V(g_am)$name, Degree.am = V(g_am)$degree, 
                         TrophicLevel = V(g_am)$TL, IntType = "amensalistic")
# Multiplex
sp_attr_multi <- bind_cols(TrophicSpecies = V(g_multi)$name, Degree = V(g_multi)$degree, 
                        TrophicLevel = V(g_multi)$TL, IntType = "multiplex")

## All interactions ----
sp_attr_deg <- sp_attr_troph[1:3] %>% 
  left_join(sp_attr_mut[1:2]) %>% 
  left_join(sp_attr_comp[1:2]) %>% 
  left_join(sp_attr_com[1:2]) %>% 
  left_join(sp_attr_am[1:2]) %>% 
  rowwise() %>% 
  mutate(Degree.total = sum(c_across(starts_with("Degree")), na.rm = T)) %>% 
  dplyr::select(TrophicSpecies, TrophicLevel, everything())

## Species type ----
# basal, intermediate or top
V(g_troph)$in.degree <- degree(g_troph, mode = "in")
V(g_troph)$out.degree <- degree(g_troph, mode = "out")

sp_type <- asDF(g_troph)[["vertexes"]] %>% 
  mutate(Type = case_when(out.degree == 0 ~ "top",
                          in.degree == 0 ~ "basal",
                          TRUE ~ "intermediate")) %>% 
  rename(TrophicSpecies = name) %>%
  dplyr::select(TrophicSpecies, Type)

## All attributes ----
sp_attr_all <- sp_attr_deg %>%
  left_join(sp_type) %>%
  replace(is.na(.), 0)


# Save results ----
save(df_complex, g_troph, g_mut, g_comp, g_com, g_com, sp_attr_all,
     file = "results/complexity_&_sppattr.rda")
