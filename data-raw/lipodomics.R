## code to prepare `lipodomics` dataset goes here

usethis::use_data(lipodomics, overwrite = TRUE)


count_matrix <- count_file_loader("lipid_NR.xlsx")
setup <- load_metadata("setup.xlsx")

