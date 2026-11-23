##***************************
## Assignment 2
##
## Group's Partner: Anthony 
##
## Name: Minh-Phuc Pham
##
## 2026-11-23
##
##***************************

## load packages -----
rm(list = ls())

library(tidyverse)

## import data -----

df_s_raw <- read_csv("data/trichinella_spiralis_te.csv")
df_zim_raw <- read_csv("data/trich_z.csv")

## filtering data for DNA transposons for T. zimbabwensis -----
dna_zim <- df_zim_raw |>
  pull(matching_class) |>
  unique() |>
  str_subset("^(DNA/|RC/)")

df_zim <- df_zim_raw |>
  filter(matching_class %in% dna_zim) |>
  mutate(species = "Trichinella zimbabwensis")

## filtering data for DNA transposons for T. spiralis -----
dna_s <- df_s_raw |>
  pull(matching_class) |>
  unique() |>
  str_subset("^(DNA/|RC/)")

df_s <- df_s_raw |>
  filter(matching_class %in% dna_s) |>
  mutate(species = "Trichinella spiralis")

## identifying and filtering the common DNA transposons between the two species -----

#using intersect() to automate the commonly shared dna transposons b/w the two species rather than manually typing it
dna_shared <- intersect(dna_s, dna_zim) 

## Filter shared DNA transposon classes ----
df_s_shared <- df_s |>
  filter(matching_class %in% dna_shared) 

df_zim_shared <- df_zim |>
  filter(matching_class %in% dna_shared)

# Retain all observations rather than selecting only the longest sequences to examine the full range of sequence lengths.

#Combine both species for comparison
df_both <- bind_rows(df_s_shared, df_zim_shared)


## Edit 2: Validate changes in data selection ----
# Compare Ali's original filtering with our revised approach.

## Reproduce Anthony group's original distinct() filtering, with my updated naming variable ----
df_old_s <- df_s_shared |>
  group_by(matching_class) |>
  select(matching_class, qry_width, perc_insert, species) |>
  distinct(perc_insert, .keep_all = TRUE) |>
  distinct(qry_width, .keep_all = TRUE)

df_old_zim <- df_zim_shared |>
  group_by(matching_class) |>
  select(matching_class, qry_width, perc_insert, species) |>
  distinct(perc_insert, .keep_all = TRUE) |>
  distinct(qry_width, .keep_all = TRUE)


## Compare the number of observations ----

# All observations retained in our revised analysis
n_full <- nrow(df_both)

# Observations retained after Ali's distinct() filtering
n_distinct <- nrow(df_old_s) + nrow(df_old_zim)

# Observations retained after Ali's top-five selection
n_top5 <- nrow(df_old_s |> top_n(5, qry_width)) +
  nrow(df_old_zim |> top_n(5, qry_width))

# Observations retained after selecting the largest per class
n_max <- nrow(df_old_s |> slice(which.max(qry_width))) +
  nrow(df_old_zim |> slice(which.max(qry_width)))

# Calculate the number and percentage removed by distinct()
n_removed <- n_full - n_distinct
pct_removed <- round(n_removed / n_full * 100, 1)


## Print results ----

print(c(
  n_full = n_full,
  n_distinct = n_distinct,
  n_top5 = n_top5,
  n_max = n_max,
  n_removed = n_removed,
  pct_removed = pct_removed
))

## Compare sequence length and insertion percentage ----
p_both <- ggplot(
  df_both,
  aes(x = qry_width, y = perc_insert, colour = species)
) +
  geom_point(alpha = 0.5, size = 2) +
  facet_wrap(~ matching_class, scales = "free_x") +
  labs(
    title = "Sequence Length and Insertion Percentage",
    subtitle = "Shared DNA transposon classes in two Trichinella species",
    x = "Query sequence length (bp)",
    y = "Insertion percentage",
    colour = "Species"
  ) +
  theme_light()

print(p_both)

ggsave(
  "shared_dna_comparison.png",
  path = "output",
  plot = p_both
)






