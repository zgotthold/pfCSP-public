#setwd("C:/Users/zgottho2/OneDrive - Johns Hopkins/Collaboration with Dr. Wesolowski/CSP project/")
setwd("/Users/gotthold/OneDrive - Johns Hopkins/Collaboration with Dr. Wesolowski/CSP PAPER/")         

### prevalence filter within each pop of >1 %
### wsaf of > 1%
## alleles < 5

library(readr)
library(tidyverse)
library(dplyr)
library(Biostrings)
#library(malaria.em)

set.seed(1)

### DO THIS ON THE AA LEVEL
V31B_amp = "QGHNMPNDPNRNVDENANANSAVKNNNNEEPSDKHIKEYLNKIQNSLSTEWSPCSVTCGNGIQ"
V32_amp = "TCGNGIQVRIKPGSANKPKDELDYANDIEKKICKMEKCSSVFNVVNSSIGLIMVLSFLFL"
DV10 = substring(V31B_amp, 293-281+1, 302 - 281+1)
Th2R = substring(V31B_amp, 311-281+1, 327 - 281+1)
Th3R = substring(V32_amp, 352-337+1, 363-337+1)


check_match <- function(string, correct_seq){
  ans <- ""
  for(char_ix in 1:nchar(correct_seq)){
    if(substring(correct_seq, char_ix, char_ix)==substring(string, char_ix, char_ix)){
      char = "_"
    }
    else {char = substring(string, char_ix, char_ix)}
    ans <- paste(ans, char, sep="")
  }
  return(ans)
}

mhap <- read.csv("DATA/ExportAllelesCSP_ASV_FILTERED.csv")

mhap_AA <- mhap %>%
  dplyr::select(specimen_id, target_id, seq, read_count, WSAF, My_Regions, COUNTRY, form.section2.age) %>%
  mutate(
    seq = case_when(target_id == "Pf3D7_03_v3-221463-221712-1B" ~ 
                     as.character(
                       Biostrings::translate(reverseComplement(DNAStringSet(substr(seq,3,193-2))))),
                   target_id == "Pf3D7_03_v3-221295-221532-2" ~
                     as.character(
                       Biostrings::translate(reverseComplement(DNAStringSet(substr(seq,1,182-2)))))))

#### MAIN PHASING

# Import data -------------------------------------------------------------

province_list <- mhap_AA %>% group_by(My_Regions) %>% 
  dplyr::summarize(sample_count = n_distinct(specimen_id)) %>%
  arrange(sample_count) %>% 
  mutate(percent.removed = NA,
         number.removed = NA,
         final_size = NA)



full_results <- data.frame(gt_id = c(), 
                           target_id = c(),
                           seq = c(),
                           freq = c(),
                           freq_se = c(),
                           province = c())

all_removed <- c()
all_kept <- c()
all_seqs <- c()


filtered_mhap <- data.frame()
for(ix in rep(1:nrow(province_list))){
  province = unlist(province_list[ix, "My_Regions"])
  cat("PROVINCE = ", province, "\n")

  mhap_province <- mhap_AA %>% 
    filter(My_Regions==province) %>% 
    dplyr::select(c(specimen_id, target_id, read_count, seq, WSAF))
  n_samples <- n_distinct(mhap_province$specimen_id)
  cat("N SAMPLES", n_samples)

  
  ## removing samples with > 5 unique alleles
  sample_diversity_high <- mhap_province %>% group_by(specimen_id, target_id) %>% 
    mutate(n_samples = n_distinct(seq)) %>% 
    dplyr::summarize(n_samples = min(n_samples)) %>%
    group_by(specimen_id) %>% filter(n_samples > 5)

  problematic_specimens <- mhap_province %>% filter(specimen_id %in% sample_diversity_high$specimen_id)
  
  province_list[ix, "percent.removed"] <- length(unique(problematic_specimens$specimen_id))/ n_samples*100
  province_list[ix, "number.removed"] <- length(unique(problematic_specimens$specimen_id))
  
  cat("REMOVING NUMBER OF SAMPLES:", length(unique(problematic_specimens$specimen_id)),"this is", length(unique(problematic_specimens$specimen_id))/ n_samples*100, "%\n")
  all_removed <- append(all_removed, problematic_specimens$specimen_id)
  
  mhap_province <- mhap_province %>% filter(!(specimen_id %in% problematic_specimens$specimen_id))
  province_list[ix, "final_size"] <- length(unique(mhap_province$specimen_id))
  all_kept <- append(all_kept, unique(mhap_province$specimen_id))
  
  filtered_mhap <- rbind(filtered_mhap, mhap_province %>% mutate(province = province))
  
  mhap_wide <- mhap_province |>
    dplyr::select(-c(read_count, WSAF)) |>
    group_by(specimen_id, target_id) |>
    pivot_wider(names_from = target_id,
                values_from = seq,
                values_fn = ~ paste(unique(.), collapse = " ")) |>
    column_to_rownames(var = "specimen_id") 
  
  locus_name <- mhap_province |> distinct(target_id) |> rowid_to_column("id") |> mutate(hap_id = paste0("loc-", id)) |> dplyr::select(-id)
  sample_name <- mhap_province |> distinct(specimen_id) |> rowid_to_column("ids") 
  
  
  matrix_EM <- mhap_wide |> 
    drop_na() |>
    as.matrix()

  max_needed <- num.unique.allele.locus(matrix_EM)
  print(paste("Will run until", max_needed))

   if(max_needed > 0){
     next
   }
  output_EM <- malaria.em(matrix_EM, sizes=1:max_needed, locus.label=colnames(matrix_EM))
  
  # Population-level multilocus genotype frequency estimate + standard error 
  gt_freq_summary <- output_EM$haplo.prob.tab |> 
    as.data.frame() |>
    rowid_to_column("gt_id") |>
    pivot_longer(cols = -c("gt_id", "hap.prob", "hap.prob.std"), 
                 names_to = "hap_id", 
                 values_to = "seq") |>
    mutate(province = province) |>
    dplyr::select (gt_id, 
            target_id = hap_id, 
            seq, 
            freq = hap.prob, 
            freq_se = hap.prob.std,
            province)

  full_results <- rbind(full_results, gt_freq_summary)

}


full_results <- full_results %>% group_by(seq) %>% 
  mutate(id_col = cur_group_id(), 
         freq_se = as.numeric(freq_se),
         freq = as.numeric(freq))

full_results <- full_results %>% arrange(desc(freq))

write.csv(full_results, "DATA/phasing_results/ALL_phasing_AA.csv")
write.csv(province_list, "DATA/phasing_results/malaria.em.size_AA.csv")

removed_samples <- data.frame(all_removed)
write.csv(removed_samples, "DATA/phasing_results/malaria.em.removed2.csv", row.names = FALSE)


#############################
##### Studying Children #####
#############################
to_remove <- read.csv("DATA/phasing_results/malaria.em.removed.csv")
mhap_AA_children <- mhap_AA %>% filter(form.section2.age == 1)

full_results_children <- data.frame(gt_id = c(), 
                           target_id = c(),
                           seq = c(),
                           freq = c(),
                           freq_se = c(),
                           province = c())

for(ix in rep(1:nrow(province_list))){
  province = unlist(province_list[ix, "My_Regions"])
  cat("PROVINCE = ", province, "CHILDREN ONLY", "\n")
  
  mhap_child <- mhap_AA_children %>% 
    filter(!(specimen_id %in% to_remove$all_removed)) %>% 
    filter(My_Regions==province) %>% 
    dplyr::select(c(specimen_id, target_id, read_count, seq, WSAF))
  n_samples <- n_distinct(mhap_child$specimen_id)
  cat("N SAMPLES", n_samples)

  
  mhap_wide <- mhap_child |>
    dplyr::select(-c(read_count, WSAF)) |>
    group_by(specimen_id, target_id) |>
    pivot_wider(names_from = target_id,
                values_from = seq,
                values_fn = ~ paste(unique(.), collapse = " ")) |>
    column_to_rownames(var = "specimen_id") 
  
  locus_name <- mhap_child |> distinct(target_id) |> rowid_to_column("id") |> mutate(hap_id = paste0("loc-", id)) |> dplyr::select(-id)
  sample_name <- mhap_child |> distinct(specimen_id) |> rowid_to_column("ids") 
  
  
  matrix_EM <- mhap_wide |> 
    drop_na() |>
    as.matrix()
  
  max_needed <- num.unique.allele.locus(matrix_EM)
  print(paste("Will run until", max_needed))
  
  if(max_needed > 5){
    next
  }
  output_EM <- malaria.em(matrix_EM, sizes=1:max_needed, locus.label=colnames(matrix_EM))
  
  # Population-level multilocus genotype frequency estimate + standard error 
  gt_freq_summary <- output_EM$haplo.prob.tab |> 
    as.data.frame() |>
    rowid_to_column("gt_id") |>
    pivot_longer(cols = -c("gt_id", "hap.prob", "hap.prob.std"), 
                 names_to = "hap_id", 
                 values_to = "seq") |>
    mutate(province = province) |>
    dplyr::select (gt_id, 
                   target_id = hap_id, 
                   seq, 
                   freq = hap.prob, 
                   freq_se = hap.prob.std,
                   province)
  
  full_results_children <- rbind(full_results_children, gt_freq_summary)
  
}

full_results_children <- full_results_children %>% group_by(seq) %>% 
  mutate(id_col = cur_group_id(), 
         freq_se = as.numeric(freq_se),
         freq = as.numeric(freq))

full_results_children <- full_results_children %>% arrange(desc(freq))

write.csv(full_results_children, "DATA/phasing_results/ALL_phasing_children_under_5_AA.csv")



##############################
##### COUNTRY-LEVEL DATA #####
##############################


full_results_country <- data.frame(gt_id = c(), 
                                    target_id = c(),
                                    seq = c(),
                                    freq = c(),
                                    freq_se = c(),
                                    province = c())

for(c in c("Eswatini", "Namibia", "Angola", "South Africa", "Zambia")){
  print(c)
  
  mhap_country <- mhap_AA %>% 
    filter(!(specimen_id %in% to_remove$all_removed)) %>%
    filter(COUNTRY==c) %>% 
    dplyr::select(c(specimen_id, target_id, read_count, seq, WSAF))
  n_samples <- n_distinct(mhap_country$specimen_id)
  cat("N SAMPLES", n_samples)
  
  mhap_wide <- mhap_country |>
    dplyr::select(-c(read_count, WSAF)) |>
    group_by(specimen_id, target_id) |>
    pivot_wider(names_from = target_id,
                values_from = seq,
                values_fn = ~ paste(unique(.), collapse = " ")) |>
    column_to_rownames(var = "specimen_id") 
  
  locus_name <- mhap_country |> distinct(target_id) |> rowid_to_column("id") |> mutate(hap_id = paste0("loc-", id)) |> dplyr::select(-id)
  sample_name <- mhap_country |> distinct(specimen_id) |> rowid_to_column("ids") 
  
  
  matrix_EM <- mhap_wide |> 
    drop_na() |>
    as.matrix()
  
  max_needed <- num.unique.allele.locus(matrix_EM)
  print(paste("Will run until", max_needed))
  
  if(max_needed > 5){
    next
  }
  output_EM <- malaria.em(matrix_EM, sizes=1:max_needed, locus.label=colnames(matrix_EM))
  
  # Population-level multilocus genotype frequency estimate + standard error 
  gt_freq_summary <- output_EM$haplo.prob.tab |> 
    as.data.frame() |>
    rowid_to_column("gt_id") |>
    pivot_longer(cols = -c("gt_id", "hap.prob", "hap.prob.std"), 
                 names_to = "hap_id", 
                 values_to = "seq") |>
    mutate(province = c) |>
    dplyr::select (gt_id, 
                   target_id = hap_id, 
                   seq, 
                   freq = hap.prob, 
                   freq_se = hap.prob.std,
                   province)
  
  full_results_country <- rbind(full_results_country, gt_freq_summary)
  
}

 full_results_country <- full_results_country %>% group_by(seq) %>% 
  mutate(id_col = cur_group_id(), 
         freq_se = as.numeric(freq_se),
         freq = as.numeric(freq))

full_results_country <- full_results_country %>% arrange(desc(freq))

#write.csv(full_results_country, "DATA/phasing_results/ALL_phasing_country_level_AA.csv")



#########################
##### EVERYONE DATA #####
#########################
full_results_all <- data.frame(gt_id = c(), 
                                   target_id = c(),
                                   seq = c(),
                                   freq = c(),
                                   freq_se = c(),
                                   province = c())

for(c in c("All Regions")){
  print(c)
  
  mhap_all <- mhap_AA %>% 
    filter(!(specimen_id %in% to_remove$all_removed)) %>%
    dplyr::select(c(specimen_id, target_id, read_count, seq, WSAF))
  n_samples <- n_distinct(mhap_all$specimen_id)
  cat("N SAMPLES", n_samples)
  
  mhap_wide <- mhap_all |>
    dplyr::select(-c(read_count, WSAF)) |>
    group_by(specimen_id, target_id) |>
    pivot_wider(names_from = target_id,
                values_from = seq,
                values_fn = ~ paste(unique(.), collapse = " ")) |>
    column_to_rownames(var = "specimen_id") 
  
  locus_name <- mhap_all |> distinct(target_id) |> rowid_to_column("id") |> mutate(hap_id = paste0("loc-", id)) |> dplyr::select(-id)
  sample_name <- mhap_all |> distinct(specimen_id) |> rowid_to_column("ids") 
  
  
  matrix_EM <- mhap_wide |> 
    drop_na() |>
    as.matrix()
  
  max_needed <- num.unique.allele.locus(matrix_EM)
  print(paste("Will run until", max_needed))
  
  if(max_needed > 5){
    next
  }
  output_EM <- malaria.em(matrix_EM, sizes=1:max_needed, locus.label=colnames(matrix_EM))
  
  # Population-level multilocus genotype frequency estimate + standard error 
  gt_freq_summary <- output_EM$haplo.prob.tab |> 
    as.data.frame() |>
    rowid_to_column("gt_id") |>
    pivot_longer(cols = -c("gt_id", "hap.prob", "hap.prob.std"), 
                 names_to = "hap_id", 
                 values_to = "seq") |>
    mutate(province = c) |>
    dplyr::select (gt_id, 
                   target_id = hap_id, 
                   seq, 
                   freq = hap.prob, 
                   freq_se = hap.prob.std,
                   province)
  
  full_results_all <- rbind(full_results_all, gt_freq_summary)
  
}

full_results_all <- full_results_all %>% group_by(seq) %>% 
  mutate(id_col = cur_group_id(), 
         freq_se = as.numeric(freq_se),
         freq = as.numeric(freq))

full_results_all <- full_results_all%>% arrange(desc(freq))
#write.csv(full_results_all, "DATA/phasing_results/ALL_phasing_single_level_AA.csv")

###############################
##### METADATA PROCESSING #####
###############################
library(gtable)
library(gtsummary)
##all removed
all_removed_df <- data.frame(removed = all_removed)
#write.csv(all_removed_df, "DATA/all_removed.csv", row.names = F)

##all kept
all_kept_df <- data.frame(kept = all_kept)
#write.csv(all_kept_df, "DATA/all_kept.csv", row.names = F)

## now we compare what is missing with "all removed"


### THIS RELIES ON OBJECTS FROM STEP 1 ####

metadata_pulled_em <- metadata_summarized_df %>% 
  mutate(status = ifelse(form.participant_id %in% to_remove$all_removed, 
                         "Removed for malaria.em",
                         "Kept for malaria.em"))
metadata_pulled_kept = metadata_pulled_em %>% filter(status == "Kept for malaria.em") %>% dplyr::select(-status)
metadata_compare <- rbind(metadata_pulled_kept %>% mutate(status = "Kept for malaria.em"),
                          metadata_summarized_df %>% mutate(status = "All samples")
  
)

metadata_compare_summary <- metadata_compare %>% 
  tbl_summary(include = c(`Age Group`, Profession, Gender), by =status)  %>% as_gt() %>%
  gt::gtsave(filename = "FIGURES/TABLES/metadata_malaria.em.docx")

# 
# ### code for running a simulation
#   
# p1_time <- c() #moxico
# p2_time <- c() #benguela
# p3_time <- c() #KwaZulu-Natal
# 
# 
# 
# test_timing <- function(province_test, max_size = 6){
#   p_time <- c()
#   for(test_coi in rep(1:max_size)){
#     
#     my_vals <- system.time({
#     #for(ix in c(8)){
#     province = province_test
#     
#     mhap <- mhap_meta %>% 
#       filter(My_Regions==province) %>% 
#       dplyr::select(c(specimen_id, target_id, read_count, seq, WSAF))
#     n_samples <- n_distinct(mhap$specimen_id)
#     sample_diversity <- mhap %>% group_by(specimen_id, target_id) %>% 
#       mutate(n_samples = n_distinct(seq)) %>% 
#       filter(WSAF > 0.02) %>%
#       dplyr::summarize(n_samples = min(n_samples),
#                        n_samples2 = n_distinct(seq)) %>%
#       group_by(specimen_id) %>% filter(n_samples > test_coi)
#     #ggplot(sample_diversity) 
# 
#     
#     
#     #problematic_specimens <<- unique((mhap %>% filter(seq %in% not_top_prevalence$seq))$specimen_id)
#     problematic_specimens <- mhap %>% filter(specimen_id %in% sample_diversity$specimen_id)
#     #problematic_specimens2 <- mhap %>% filter(specimen_id %in% sample_diversity2$specimen_id)
#     #problematic_alleles <- unique(low_prevalence$seq)
#     province_list[ix, "percent.removed"] <- length(unique(problematic_specimens$specimen_id))/ n_samples*100
#    # province_list[ix, "percent.removed2"] <- length(unique(problematic_specimens2$specimen_id))/ n_samples*100
#     province_list[ix, "number.removed"] <- length(unique(problematic_specimens$specimen_id))
#    # province_list[ix, "number.removed2"] <- length(unique(problematic_specimens2$specimen_id))
#     
#     cat("REMOVING NUMBER OF SAMPLES:", length(unique(problematic_specimens$specimen_id)),"this is", length(unique(problematic_specimens$specimen_id))/ n_samples*100, "%\n")
#     all_removed <- append(all_removed, problematic_specimens$specimen_id)
#     old_mhap <- mhap
#     mhap <- mhap %>% filter(!(specimen_id %in% problematic_specimens$specimen_id))
#     province_list[ix, "final_size"] <- length(unique(mhap$specimen_id))
#     
#     mhap_wide <<- mhap |>
#       dplyr::select(-c(read_count, WSAF)) |>
#       group_by(specimen_id, target_id) |>
#       pivot_wider(names_from = target_id,
#                   values_from = seq,
#                   values_fn = ~ paste(unique(.), collapse = " ")) |>
#       column_to_rownames(var = "specimen_id") 
#     
#     
#     old_mhap_wide <- old_mhap |>
#       dplyr::select(-c(read_count, WSAF)) |>
#       group_by(specimen_id, target_id) |>
#       pivot_wider(names_from = target_id,
#                   values_from = seq,
#                   values_fn = ~ paste(unique(.), collapse = " ")) |>
#       column_to_rownames(var = "specimen_id") 
#     
#     
#     locus_name <- mhap |> distinct(target_id) |> rowid_to_column("id") |> mutate(hap_id = paste0("loc-", id)) |> dplyr::select(-id)
#     sample_name <- mhap |> distinct(specimen_id) |> rowid_to_column("ids") 
#     
#     
#     laos_mat_2 <- mhap_wide |> 
#       drop_na() |>
#       as.matrix()
#     
#     old_laos_mat_2 <- old_mhap_wide |> 
#       drop_na() |>
#       as.matrix()
#     
#     # laos_mat <- as.matrix(mhap_wide[, 2:4])
#     
#     max_needed <- num.unique.allele.locus(laos_mat_2)
#     province_list[ix, "m.em.size.old"] <- num.unique.allele.locus(old_laos_mat_2)
#     province_list[ix, "m.em.size"] <- max_needed
#     
#     print(paste("Will run until", max_needed))
#     
#     #exit
#     
#     # 
#     if(max_needed > 6){
#       next
#     }
#     # province_list[ix, "m.em.size.highprev"] <- max_needed
#     output_laos <- malaria.em(laos_mat_2, sizes=1:max_needed, locus.label=colnames(laos_mat_2))
#     #print.malaria.em(output_laos)
#   })
#   print(my_vals[3])
#   p_time <- c(p_time, my_vals[3])
#   
#   }
#   return(p_time)
# }
# 
# test_frame <- data.frame(province = c(), sizes = c(), time = c())
# for(provinces in c("Moxico", "Benguela", "Zaire", 
# "Luapula Province", "Eswatini", "KwaZulu-Natal")){
#   test_values <- test_timing(provinces)
#   test_frame <- rbind(test_frame, data.frame(province = provinces, sizes = rep(1:length(test_values)),
#                                              time = test_values))
# }
# 
# 
# ggplot(test_frame, aes(x = sizes/60, y = time, color= province)) + geom_point() +
#   ylab("Time in Minutes") +
#   xlab("Max Unique Alleles")
# scale_y_log10()  
