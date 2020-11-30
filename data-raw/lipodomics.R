## code to prepare `lipodomics` dataset goes here



count_matrix <- count_file_loader("lipid_NR.xlsx")
setup <- load_metadata("setup.xlsx")

setup <- setup %>%
  unite(Group, Genotype, NR, sep = "_", remove = F)


#Might not be beneficial to remove low counts as data are percentile
#count_matrix <- select_sufficient_counts(count_matrix,setup,2)

#split the file into the two subgroups
count_matrix_groups <- count_matrix[217:238,]
count_matrix_sub <- count_matrix[1:216,]

#sum(count_matrix_groups/26)
#sum(count_matrix_sub/26)
#now data are organized as 100% pr column


#as data is normalized pr info from Mesut, they should be directly comparable
#dir.create(here::here("data/figures"), showWarnings = F)
#Quality_control_plots(count_matrix_group,setup)
#rename folder so it doesn't overwrite
#Quality_control_plots(count_matrix_sub,setup)

#based on QC, 330 is excluded (was steatotic)

setup <- setup %>%
  filter(!ID=="330")
colnames(count_matrix) <- as.integer(colnames(count_matrix))

count_matrix <- count_matrix %>%
  dplyr::select(-"330")

count_matrix_groups <- count_matrix[217:238,]
count_matrix_sub <- count_matrix[1:216,]


all(setup$ID==colnames(count_matrix))

#LPC O- causes problems for analysis so remove it
count_matrix_groups <- count_matrix_groups %>%
  dplyr::filter(rownames(count_matrix_groups)!="LPC O-")

lipodomics_results <- DEG_analysis(count_matrix_groups, setup)


count_matrix_sub <- select_sufficient_counts(count_matrix_sub,setup,2)

lipodomics_results_sub <- DEG_analysis(count_matrix_sub, setup)

#heatmap<-heatmap(count_matrix,setup)


