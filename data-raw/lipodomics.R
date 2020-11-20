## code to prepare `lipodomics` dataset goes here



count_matrix <- count_file_loader("lipid_NR.xlsx")
setup <- load_metadata("setup.xlsx")

setup <- setup %>%
  unite(Group, Genotype, NR, sep = "_", remove = F)

count_matrix <- select_sufficient_counts(count_matrix,setup,2)


#as data is normalized pr info from Mesut, they should be directly comparable
#dir.create(here::here("data/figures"), showWarnings = F)
#Quality_control_plots(count_matrix,setup)

#based on QC, 330 is excluded (was steatotic)

setup <- setup %>%
  filter(!ID=="330")
colnames(count_matrix) <- as.integer(colnames(count_matrix))
count_matrix <- count_matrix %>%
  dplyr::select(-"330")

all(setup$ID==colnames(count_matrix))

lipodomics_results <- DEG_analysis(count_matrix, setup)

#heatmap<-heatmap(count_matrix,setup)


