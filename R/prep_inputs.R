crp_file <- "data/CRP/CRP_all_SNPs.txt"

crp <- crp_file |>
  vroom::vroom()
crp_small <- crp |>
  dplyr::filter(P < 1e-8) 
crp_small |>
  vroom::vroom_write(file = "CRP_1e-8.csv")

# ARMD

armd_file <- "data/ARMD/GCST90086112.h.tsv.gz"

armd <- armd_file |>
  vroom::vroom()

### 
armd_small <- armd |>
  dplyr::filter(rs_id_all %in% crp_small$SNP)

armd_small |>
  vroom::vroom_write("ARMD_small.csv")
