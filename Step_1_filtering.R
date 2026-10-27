#setwd("C:/Users/zgottho2/OneDrive - Johns Hopkins/Collaboration with Dr. Wesolowski/CSP project/")
#setwd("C:/Users/zoeor/OneDrive - Johns Hopkins/Collaboration with Dr. Wesolowski/CSP PAPER/")
#setwd("./JHU - OneDrive/OneDrive - Johns Hopkins/Collaboration with Dr. Wesolowski/CSP PAPER")
setwd("/Users/gotthold/OneDrive - Johns Hopkins/Collaboration with Dr. Wesolowski/CSP PAPER/")     

### prevalence filter within each pop of >1 %
### wsaf of > 1%
## alleles < 5

library(readr)
library(tidyverse)
library(dplyr)
library(Biostrings)
#library(malaria.em)


# ### DO THIS ON THE AA LEVEL
# V31B_amp = "QGHNMPNDPNRNVDENANANSAVKNNNNEEPSDKHIKEYLNKIQNSLSTEWSPCSVTCGNGIQ"
# V32_amp = "TCGNGIQVRIKPGSANKPKDELDYANDIEKKICKMEKCSSVFNVVNSSIGLIMVLSFLFL"
# DV10 = substring(V31B_amp, 293-281+1, 302 - 281+1)
# Th2R = substring(V31B_amp, 311-281+1, 327 - 281+1)
# Th3R = substring(V32_amp, 352-337+1, 363-337+1)
# 
# 
# check_match <- function(string, correct_seq){
#   ans <- ""
#   for(char_ix in 1:nchar(correct_seq)){
#     if(substring(correct_seq, char_ix, char_ix)==substring(string, char_ix, char_ix)){
#       char = "_"
#     }
#     else {char = substring(string, char_ix, char_ix)}
#     ans <- paste(ans, char, sep="")
#   }
#   return(ans)
# }

alleles_export <- read_tsv("DATA/ExportAllelesCSP_ASV.txt", show_col_types = FALSE)



# Import data -------------------------------------------------------------

country_key <- read.csv("DATA/country_key_map.csv")

contaminated_samples <- read.csv("DATA/exclusion_16Feb2024.csv")

exclusion_clusters <- read.csv("DATA/exclusion_clusters.tsv", sep = "\t") %>%
  filter(group!="nogroup")

facility_provinces = read.csv("DATA/facility_metadata.tsv", sep = "\t")
meta <- read.csv("DATA/All_GenE8_Countries_Epi_Harmonized.csv") %>%
  left_join(facility_provinces %>% dplyr::select(facility_id ,province_name), by = join_by(form.facility == facility_id)) %>%
  dplyr::select(-facility_province) %>%
  dplyr::rename(facility_province = province_name) %>%
  left_join(country_key, by = join_by(facility_province==Provinces_Contained_SNPs))

#write.csv(meta, "DATA/All_GenE8_Countries_Epi_Harmonized_SA.csv")



alleles_export <- alleles_export %>% dplyr::rename(specimen_id = FieldID,
                        target_id = Locus,
                        read_count = Reads,
                        WSAF = Reads.Norm,
                        seq = ASV) %>%
  filter(target_id != "Pf3D7_03_v3-221185-221419-1B") %>% # only phasing the first two loci
  dplyr::select(c(specimen_id, target_id, read_count, seq, WSAF)) 


alleles_export_meta <- alleles_export %>% left_join(meta, by = join_by(specimen_id == form.participant_id)) %>%
  filter(is.na(My_Regions) | My_Regions != "Lunda Sul") # only removing lunda sul

alleles_export_meta <- alleles_export_meta %>% filter(!is.na(My_Regions))

sample_filtering <- data.frame(category = "initial samples from our study regions",
                               n = n_distinct(alleles_export_meta$specimen_id))

alleles_export_meta <- alleles_export_meta %>%
  filter(!(specimen_id %in% contaminated_samples$sample)) # removing contaminated samples


sample_filtering <- sample_filtering %>% rbind(list(category = "after removal of contaminated samples",
                                                    n = n_distinct(alleles_export_meta$specimen_id)))

### removing HF and duplicated samples:
alleles_export_meta <- alleles_export_meta %>%
  filter(!tolower(facility_name) %in% c(
      "zemba health centre",
      "chibolya health post",
      "chibansa urban health centre",
      "maheba d rural health centre",
      "kazungula health post"
    )) 

sample_filtering <- sample_filtering %>% rbind(list(category = "after removal of excluded health facilities",
                                                    n = n_distinct(alleles_export_meta$specimen_id)))


alleles_export_meta <- alleles_export_meta %>% filter(!(specimen_id %in% exclusion_clusters$sample)) # removing exclusion clusters

sample_filtering <- sample_filtering %>% rbind(list(category = "after removal of exclusion clusters",
                                                    n = n_distinct(alleles_export_meta$specimen_id)))


table((alleles_export_meta %>% distinct(specimen_id, facility_province))$facility_province, useNA = "ifany")



alleles_export_meta <- alleles_export_meta %>% group_by(My_Regions) %>%
  mutate(n_samples = n_distinct(specimen_id)) %>% group_by(My_Regions, seq) %>%
  mutate(n_seqs = n_distinct(specimen_id)) %>% ungroup() %>%
  mutate(population_prevalence = n_seqs/n_samples) %>% filter(population_prevalence > 0.01)


old_size <- n_distinct(alleles_export_meta$specimen_id)
print("removing samples that did not amplify on both targets")

alleles_export_meta <- alleles_export_meta %>% group_by(specimen_id) %>%
  filter("Pf3D7_03_v3-221463-221712-1B" %in% target_id,
         "Pf3D7_03_v3-221295-221532-2" %in% target_id)
new_size <- n_distinct(alleles_export_meta$specimen_id)
cat("We lost", old_size - new_size, "samples from", old_size)

sample_filtering <- sample_filtering %>% rbind(list(category = "after removing samples that did not amplify on both targets",
                                                    n = n_distinct(alleles_export_meta$specimen_id)))


COI_csp_cals <- alleles_export_meta %>% group_by(specimen_id, target_id) %>%
  dplyr::summarize(n = n_distinct(seq)) %>%
  pivot_wider(names_from = target_id, values_from = n) %>%
  mutate(estimated_coi = pmax(`Pf3D7_03_v3-221295-221532-2`, `Pf3D7_03_v3-221463-221712-1B`)) 

#write.csv(COI_csp_cals, "DATA/coi_table.csv")

# 
# ggplot(COI_csp_cals, aes(x = estimated_coi, y = coi)) + geom_count() + 
#   xlab("COI (by CSP estimation)") +
#   ylab("COI (by Alfred)") + xlim(c(0, 12)) + ylim(c(0, 12)) + geom_abline(slope = 1, intercept = 0)
metadata_pulled <-meta %>% filter(form.participant_id %in% alleles_export_meta$specimen_id)
#write.csv(alleles_export_meta, "DATA/ExportAllelesCSP_ASV_FILTERED.csv")


original_SNPs_file <- read.csv("DATA/csp_SNPsconsolidated_250418.txt.gz", sep="\t")
filtered_SNPs_file <- original_SNPs_file %>% filter(FieldID %in% alleles_export_meta$specimen_id) %>%
  left_join(meta, by = join_by(FieldID == form.participant_id)) %>%
  filter(is.na(My_Regions) | My_Regions != "Lunda Sul")


sample_filtering <- sample_filtering %>% rbind(list(category = "final",
                                                    n = n_distinct(alleles_export_meta$specimen_id)))

EXIT


#write.csv(filtered_SNPs_file, "DATA/csp_SNPsconsolidated_250418_FILTERED.csv")

###############################
##### METADATA PROCESSING #####
###############################
library(gtable)
library(gtsummary)

metadata_summarized_df <- metadata_pulled %>%
  mutate(`Age Group`= case_match(form.section2.age,
                                 1 ~ "<5",
                                 2 ~ "5-14",
                                 3 ~ "15-24",
                                 4 ~ "25-39",
                                 5 ~ "40-59",
                                 6 ~ "60+"),
         `Age Group` = factor(`Age Group`, levels=c("<5", "5-14", "15-24", "25-39", "40-59", "60+")),
         `Profession` = case_match(form.section2.occupation, 
                                   1 ~ "Farming/Agriculture",
                                   2 ~ "Business (other than farming)",
                                   3 ~ "Student",
                                   4 ~ "Child or minor, but not student",
                                   5 ~ "Mining",
                                   6 ~ "Unemployed",
                                   7 ~ "Driver", 
                                   8 ~ "Self-employed", 
                                   88 ~ "Other", ),
         `Gender` = case_match(form.section2.gender, 
                               1 ~ "Male",
                               2 ~ "Female",
                               NA ~ NA,
                               99 ~ NA,
                               .default = "Other"))

metadata_summarized <- metadata_summarized_df %>%
  tbl_summary(include = c(`Age Group`, Profession, Gender), by =COUNTRY) %>%
  add_p(test = everything() ~ "chisq.test",
    test.args = all_tests("chisq.test") ~ list(simulate.p.value=TRUE)) %>% as_gt() # %>%
  # gt::gtsave(filename = "FIGURES/TABLES/metadata_country.docx")

metadata_summarized

write.csv(sample_filtering, "FIGURES/TABLES/sample_filtering.csv")

library(ggstats)

age_summary <- metadata_summarized_df %>%
  filter(!is.na(`Age Group`)) %>%
  group_by(COUNTRY, `Age Group`) %>%
  dplyr::summarise(count = n()) %>%
  group_by(COUNTRY) %>%
  mutate(proportion = count / sum(count)) %>%
  mutate(label = paste0(count, " (", round(proportion*100, 1), "%)")) 


age_by_country = ggplot(age_summary, aes(x = COUNTRY, fill = `Age Group`, y = proportion)) +
  geom_bar(stat= "identity") +
  geom_text(aes(label = label), position = position_fill(0.5)) +
  ylab("Proportion of samples") +
  xlab(NULL) + theme_minimal()

age_summary_total <- metadata_summarized_df %>%
  filter(!is.na(`Age Group`)) %>%
  group_by(`Age Group`) %>%
  dplyr::summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  mutate(label = paste0(count, " (", round(proportion*100, 1), "%)")) %>%
  mutate(COUNTRY = "OVERALL")

age_overall = 
  ggplot(age_summary_total, aes(x = COUNTRY, fill = `Age Group`, y = proportion)) +
  geom_bar(stat= "identity") +
  geom_text(aes(label = label), position = position_fill(0.5)) +
  ylab("Proportion of samples") +
  xlab(NULL) + theme_minimal()

(age_overall + age_by_country + plot_layout(widths = c(1, 5), guides = "collect", axes = "collect")) +
  ggtitle("Age Distribution by Country")


##chi-squared tests for each country
age_summary %>% filter(COUNTRY == "Angola") %>% pull(count)
age_summary_total %>% pull(count)

x = age_summary %>% filter(COUNTRY == "Angola") %>% pull(count)
y = age_summary_total %>% pull(count)
test_matrix = as.matrix(rbind(x, y))


###### SEQUENCING QUALITY #####
get_iqr <- function(vec){
  p25 = quantile(vec, 0.25)
  p75 = quantile(vec, 0.75)
  cat(p25, "-", p75)
  cat("\nthe median is", median(vec))
}

get_iqr(alleles_export_meta %>%
          filter(target_id == "Pf3D7_03_v3-221463-221712-1B") %>%
          pull(read_count))
get_iqr(alleles_export_meta %>%
          filter(target_id == "Pf3D7_03_v3-221295-221532-2") %>%
          pull(read_count))


##alleles per sample table

alleles_export_meta %>% group_by(target_id) %>%
  summarize(total = n_distinct(seq))

alleles_per_sample <- alleles_export_meta %>% group_by(specimen_id, target_id, My_Regions) %>%
  dplyr::summarize(n_alleles = n_distinct(seq),
                   monoclonal = n_alleles == 1) %>% group_by(target_id, My_Regions) %>%
  dplyr::summarize(mean = mean(n_alleles),
            median = median(n_alleles),
            per_mono = sum(monoclonal) / n_distinct(specimen_id),
            IQR = paste0(quantile(n_alleles, 0.25), ",", quantile(n_alleles, 0.75)))

alleles_per_sample_overall <- alleles_export_meta %>% group_by(specimen_id, target_id) %>%
  dplyr::summarize(n_alleles = n_distinct(seq),
                   monoclonal = n_alleles == 1) %>% group_by(target_id) %>%
  dplyr::summarize(mean = mean(n_alleles),
            median = median(n_alleles),
            per_mono = sum(monoclonal) / n_distinct(specimen_id),
            IQR = paste0(quantile(n_alleles, 0.25), ",", quantile(n_alleles, 0.75))) %>%
  mutate(My_Regions = "OVERALL")

#calculating monoclonal at both
alleles_monoclonal <- alleles_export_meta %>% group_by(specimen_id, target_id) %>%
  dplyr::summarize(n_alleles = n_distinct(seq),
                   monoclonal = n_alleles == 1) %>% ungroup() %>% group_by(specimen_id) %>%
  dplyr::summarise(max_alleles = max(n_alleles),
                   min_alleles = min(n_alleles)) 
cat("Out of", nrow(alleles_monoclonal), "samples, only", sum(alleles_monoclonal$max_alleles == 1), "are monoclonal at both sites.",
    sum(alleles_monoclonal$min_alleles == 1), "samples are monoclonal at either site (naively phasable),",
    sum(alleles_monoclonal$max_alleles == 1)/nrow(alleles_monoclonal),
    sum(alleles_monoclonal$min_alleles == 1)/nrow(alleles_monoclonal))

num_alleles_naive <- alleles_export_meta %>% filter(specimen_id %in% 
                                                      (alleles_monoclonal %>% filter(min_alleles == 1) %>%
                                                         pull(specimen_id)))



out <- num_alleles_naive %>%
  dplyr::select(target_id, seq, specimen_id) %>%
  group_by(specimen_id) %>%   # use your ID columns here if you have them
  reframe({
    t1 <- seq[target_id == "Pf3D7_03_v3-221463-221712-1B"]
    t2 <- seq[target_id == "Pf3D7_03_v3-221295-221532-2"]
    
    if (length(t1) == 1 && length(t2) >= 1) {
      tibble(seq = paste0(t1[1], "_",t2))
    } else if (length(t2) == 1 && length(t1) >= 1) {
      tibble(seq = paste0(t1, "_", t2[1]))
    } else {
      stop("Expected exactly one of T1/T2 to have a single seq value.")
    }
  })
cat("We have", n_distinct(out$seq), "unique haplotypes among naively phasable samples")
write.csv(out, "DATA/naive_phasing.csv")

alleles_per_sample <- rbind(alleles_per_sample, alleles_per_sample_overall) %>%
  mutate(target_id_clean = case_when(target_id == "Pf3D7_03_v3-221463-221712-1B" ~ "T1",
                                     target_id == "Pf3D7_03_v3-221295-221532-2" ~ "T2")) %>%
  ungroup() %>%
  mutate(mean = round(mean, 2),
         median = paste0(median, " [", IQR, "]"),
         per_mono = round(per_mono *100, 1)) %>%
  dplyr::select(-c(target_id, IQR)) %>% pivot_wider(names_from = target_id_clean, values_from = c(mean, median, per_mono),
                                               names_glue = "{.value}_{target_id_clean}")

cleaned_gtable_alleles_export <- alleles_per_sample %>%
  gt(rowname_col = "My_Regions") %>%
  tab_spanner(label = "DV10-Th2R", columns = c(mean_T1, median_T1, per_mono_T1)) %>%
  tab_spanner(label = "DV10-Th3R", columns = c(mean_T2, median_T2, per_mono_T2)) %>%
  cols_label(
    mean_T1 = "mean",
    median_T1 = "median [IQR]",
    mean_T2 = "mean",
    median_T2 = "median [IQR]",
    per_mono_T1 = "monoclonal (%)",
    per_mono_T2 = "monoclonal (%)"
  )

#cleaned_gtable_alleles_export %>% gtsave("FIGURES/TABLES/alleles_per_sample_by_province.docx")

#write.csv(alleles_per_sample,"FIGURES/TABLES/alleles_per_sample_by_province.csv", row.names = FALSE)
