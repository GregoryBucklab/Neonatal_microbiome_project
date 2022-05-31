library(outliers)
library('dplyr')
library("ggpubr")
#rm(list = ls.str(mode = 'numeric'))
#rm(list = ls.str(mode = 'character'))
#rm(list = ls.str(mode = 'logical'))

##### get metadata and 16s, cytokine reads table #####
# run 'get_pop2_ptb47_metadata.R' 'get_lipidomics.R'

setwd('/Users/binzhu/Desktop/Mom-baby/')
##### data preparation #####
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Lachnospiraceae_G-9_bacterium_HMT_924')] = 'BVAB1'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Lachnospiraceae_G-9_bacterium_HMT_924_BT')] = 'BVAB1_BT'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Sneathia_amnii_Not_Validly_Published')] = 'Sneathia_amnii'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Sneathia_amnii_Not_Validly_Published_BT')] = 'Sneathia_amnii_BT'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == '1-68_OTU_12')] = 'Tissierellaceae_1_68_OTU_12'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'WAL_1855D_OTU_17')] = 'Tissierellaceae_WAL_1855D_OTU_17'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'WAL_1855D_OTU_17_BT')] = 'Tissierellaceae_WAL_1855D_OTU_17_BT'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Rs-045_OTU_3')] = 'TM7_H1'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'ph2_OTU_4')] = 'Tissierellaceae_ph2_OTU_4'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == '1-68_OTU_17')] = 'Tissierellaceae_1_68_OTU_17'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'WAL_1855D_OTU_29')] = 'Tissierellaceae_WAL_1855D_OTU_29'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == '1-68_OTU_6')] = 'Tissierellaceae_1_68_OTU_6'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == '258ds10_OTU_4_BT')] = 'Fibrobacteria_258ds10_OTU_4_BT'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == '1-68_OTU_16')] = 'Tissierellaceae_1_68_OTU_16'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == '1-68_OTU_18')] = 'Tissierellaceae_1_68_OTU_18'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'WAL_1855D_OTU_19')] = 'Tissierellaceae_WAL_1855D_OTU_19'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'WAL_1855D_OTU_27')] = 'Tissierellaceae_WAL_1855D_OTU_27'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'WAL_1855D_OTU_8')] = 'Tissierellaceae_WAL_1855D_OTU_8'
row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'ph2_OTU_13')] = 'Tissierellaceae_ph2_OTU_13'

taxonomy$V3[which(taxonomy$V3 == 'Lachnospiraceae_G-9_bacterium_HMT_924')] = 'BVAB1'
taxonomy$V3[which(taxonomy$V3 == 'Lachnospiraceae_G-9_bacterium_HMT_924_BT')] = 'BVAB1_BT'
taxonomy$V3[which(taxonomy$V3 == 'Sneathia_amnii_Not_Validly_Published')] = 'Sneathia_amnii'
taxonomy$V2[which(taxonomy$V3 == 'Sneathia_amnii')] = 'k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;g__Sneathia;s__amnii'
taxonomy$V3[which(taxonomy$V3 == 'Sneathia_amnii_Not_Validly_Published_BT')] = 'Sneathia_amnii_BT'
taxonomy$V2[which(taxonomy$V3 == 'Sneathia_amnii_BT')] = 'k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;g__Sneathia;s__amnii'
taxonomy$V3[which(taxonomy$V3 == '1-68_OTU_12')] = 'Tissierellaceae_1_68_OTU_12'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_1_68_OTU_12')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'WAL_1855D_OTU_17')] = 'Tissierellaceae_WAL_1855D_OTU_17'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_WAL_1855D_OTU_17')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'WAL_1855D_OTU_17_BT')] = 'Tissierellaceae_WAL_1855D_OTU_17_BT'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_WAL_1855D_OTU_17_BT')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'Rs-045_OTU_3')] = 'TM7_H1'
taxonomy$V2[which(taxonomy$V3 == 'TM7_H1')] = 'k__Bacteria;p__TM7;c__TM7-3;o__;f__;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'ph2_OTU_4')] = 'Tissierellaceae_ph2_OTU_4'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_ph2_OTU_4')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == '1-68_OTU_17')] = 'Tissierellaceae_1_68_OTU_17'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_1_68_OTU_17')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'WAL_1855D_OTU_29')] = 'Tissierellaceae_WAL_1855D_OTU_29'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_WAL_1855D_OTU_29')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == '1-68_OTU_6')] = 'Tissierellaceae_1_68_OTU_6'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_1_68_OTU_6')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == '258ds10_OTU_4_BT')] = 'Fibrobacteria_258ds10_OTU_4_BT'
taxonomy$V2[which(taxonomy$V3 == 'Fibrobacteria_258ds10_OTU_4_BT')] = 'k__Bacteria;p__Fibrobacteres;c__Fibrobacteria;o__;f__;g__;s__'
taxonomy$V3[which(taxonomy$V3 == '1-68_OTU_16')] = 'Tissierellaceae_1_68_OTU_16'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_1_68_OTU_16')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'WAL_1855D_OTU_19')] = 'Tissierellaceae_WAL_1855D_OTU_19'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_WAL_1855D_OTU_19')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'WAL_1855D_OTU_27')] = 'Tissierellaceae_WAL_1855D_OTU_27'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_WAL_1855D_OTU_27')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'WAL_1855D_OTU_8')] = 'Tissierellaceae_WAL_1855D_OTU_8'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_WAL_1855D_OTU_8')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == 'ph2_OTU_13')] = 'Tissierellaceae_ph2_OTU_13'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_ph2_OTU_13')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
taxonomy$V3[which(taxonomy$V3 == '1-68_OTU_18')] = 'Tissierellaceae_1_68_OTU_18'
taxonomy$V2[which(taxonomy$V3 == 'Tissierellaceae_1_68_OTU_18')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'

# give flags
setwd('/Users/binzhu/Desktop/Mom-baby/')

# remove baby more than 2 day
metadata_16s$days_rel2birth <- as.numeric(as.character(metadata_16s$days_rel2birth))
keep = metadata_16s$days_rel2birth > 2 & !is.na(metadata_16s$days_rel2birth)
sum(keep)
keep = !keep
reads_table_16s = reads_table_16s[,keep]
metadata_16s = metadata_16s[keep,]
metadata_16s$days_rel2birth[metadata_16s$KitType == 'Baby_Birth' & is.na(metadata_16s$days_rel2birth)] = 0

# flag
color1=c('#E43B2D','#0070FF','#4DAF4A','#FEF840','#F27E33','#F180BF')

metadata_16s$Flag = NA
metadata_16s$Flag[metadata_16s$SampleType == "BCKD" & metadata_16s$days_rel2birth == 2 & !is.na(metadata_16s$days_rel2birth)] = "BCKD_birth_2_days"
metadata_16s$Flag[metadata_16s$SampleType == "BRCD" & metadata_16s$days_rel2birth == 2 & !is.na(metadata_16s$days_rel2birth)] = "BRCD_birth_2_days"
metadata_16s$Flag[metadata_16s$SampleType == "BS1D" & metadata_16s$days_rel2birth == 2 & !is.na(metadata_16s$days_rel2birth)] = "BS1D_birth_2_days"

metadata_16s$Flag[metadata_16s$SampleType == "BCKD" & metadata_16s$days_rel2birth == 1 & !is.na(metadata_16s$days_rel2birth)] = "BCKD_birth_1_days"
metadata_16s$Flag[metadata_16s$SampleType == "BRCD" & metadata_16s$days_rel2birth == 1 & !is.na(metadata_16s$days_rel2birth)] = "BRCD_birth_1_days"
metadata_16s$Flag[metadata_16s$SampleType == "BS1D" & metadata_16s$days_rel2birth == 1 & !is.na(metadata_16s$days_rel2birth)] = "BS1D_birth_1_days"

metadata_16s$Flag[metadata_16s$SampleType == "BCKD" & metadata_16s$KitType == 'Baby_Birth' & metadata_16s$days_rel2birth == 0 & !is.na(metadata_16s$days_rel2birth)] = 'BCKD_birth_0_day'
metadata_16s$Flag[metadata_16s$SampleType == "BRCD" & metadata_16s$KitType == 'Baby_Birth' & metadata_16s$days_rel2birth == 0 & !is.na(metadata_16s$days_rel2birth)] = 'BRCD_birth_0_day'
metadata_16s$Flag[metadata_16s$SampleType == "BS1D" & metadata_16s$KitType == 'Baby_Birth' & metadata_16s$days_rel2birth == 0 & !is.na(metadata_16s$days_rel2birth)] = 'BS1D_birth_0_day'

for (a in 1:nrow(metadata_16s)) {
  if (!is.na(metadata_16s$Flag[a])) {
    next
  }
  
  if (metadata_16s$weeks_pregnant[a] <= 12 & !is.na(metadata_16s$weeks_pregnant[a])) {
    metadata_16s$Flag[a] = paste(metadata_16s$SampleType[a],'1st_trimester',sep = '_')
  } 
  if (metadata_16s$weeks_pregnant[a] > 12 & metadata_16s$weeks_pregnant[a] <= 26 & !is.na(metadata_16s$weeks_pregnant[a])) {
    metadata_16s$Flag[a] = paste(metadata_16s$SampleType[a],'2nd_trimester',sep = '_')
  } 
  if (metadata_16s$weeks_pregnant[a] > 26 & !is.na(metadata_16s$weeks_pregnant[a])) {
    metadata_16s$Flag[a] = paste(metadata_16s$SampleType[a],'3rd_trimester',sep = '_')
  }

}

keep = !is.na(metadata_16s$Flag)
reads_table_16s = reads_table_16s[,keep]
metadata_16s = metadata_16s[keep,]

 metadata_16s$SampleID_ori = metadata_16s$SampleID
 metadata_16s$SampleID = paste('Psudo_sample_ID',c(1:nrow(metadata_16s)), sep = '_')  ### run when output contains sample ID or participant ID
 colnames(reads_table_16s) = paste('Psudo_sample_ID',c(1:nrow(metadata_16s)), sep = '_')  ### run when output contains sample ID or participant ID

 
 for (a in 1: ncol(reads_table_cytokines)) {
   n = which(metadata_16s$SampleID_ori == colnames(reads_table_cytokines)[a])
   if (length(n) ==1) {
     colnames(reads_table_cytokines)[a] = metadata_16s$SampleID[n]
   }
 }
 keep = colnames(reads_table_cytokines) %in% metadata_16s$SampleID
 reads_table_cytokines = reads_table_cytokines[,keep]
 
 for (a in 1: ncol(reads_table_lipid)) {
   n = which(metadata_16s$SampleID_ori == colnames(reads_table_lipid)[a])
   if (length(n) ==1) {
     colnames(reads_table_lipid)[a] = metadata_16s$SampleID[n]
   }
 }
 keep = colnames(reads_table_lipid) %in% metadata_16s$SampleID
 reads_table_lipid = reads_table_lipid[,keep]
 
 for (a in 1: ncol(reads_table_human_MTG)) {
   n = which(metadata_16s$SampleID_ori == colnames(reads_table_human_MTG)[a])
   if (length(n) ==1) {
     colnames(reads_table_human_MTG)[a] = metadata_16s$SampleID[n]
   }
 }
 keep = colnames(reads_table_human_MTG) %in% metadata_16s$SampleID
 reads_table_human_MTG = reads_table_human_MTG[,keep]
 
 metadata_all$SampleID_2 = metadata_all$SampleID
 metadata_all$SampleID_2 = paste('Psudo_sample_ID',c(1:nrow(metadata_all)), sep = '_')  ### run when output contains sample ID or participant ID
 metadata_all$ParticipantID_2 = metadata_all$ParticipantID  ### run when output contains sample ID or participant ID
 
participant_id_unique = unique(metadata_all$ParticipantID)
for (a in 1: length(participant_id_unique)) {
  n = which(metadata_all$ParticipantID ==  participant_id_unique[a])
  metadata_all$ParticipantID[n] = paste('Psudo_participant_ID',a, sep = '_')
}

participant_id_unique = unique(metadata_16s$ParticipantID)
for (a in 1: length(participant_id_unique)) {
  n = which(metadata_16s$ParticipantID ==  participant_id_unique[a])
  metadata_16s$ParticipantID[n] = paste('Psudo_participant_ID',a, sep = '_')
}

# get main taxa list in 6 microbiomes
type_all = unique(metadata_16s$SampleType)

for (a in 1:length(type_all)) {
  keep = metadata_16s$SampleType == type_all[a]
  reads_table = reads_table_16s[,keep]
  
  # species and sample total reads thresholds and get taxonomy
  keep <- colSums(reads_table) > 5000   # sample total reads >= 5000
  reads_table = reads_table[,keep]
  
  keep <- rowSums(reads_table) >= ncol(reads_table)   # species total reads >= sample number
  sum(keep)
  reads_table <- reads_table[keep,]
  reads_table = prepare_reads_table_2(reads_table, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
  
  if (a == 1) {
    reads_table_16s_BRCD = reads_table
  } else if (a == 2) {
    reads_table_16s_BS1D = reads_table
  } else if (a == 3) {
    reads_table_16s_BCKD = reads_table
  } else if (a == 4) {
    reads_table_16s_MCKD = reads_table
  } else if (a == 5) {
    reads_table_16s_MV1D = reads_table
  } else {
    reads_table_16s_MRCD = reads_table
  } 
}

taxa_list = c(row.names(reads_table_16s_BRCD),row.names(reads_table_16s_BS1D),
              row.names(reads_table_16s_BCKD),row.names(reads_table_16s_MCKD),
              row.names(reads_table_16s_MV1D),row.names(reads_table_16s_MRCD))
taxa_list = unique(taxa_list)

# remove samples not passed threshold
passed_sample_list = c(colnames(reads_table_16s_BRCD),colnames(reads_table_16s_BS1D),
                       colnames(reads_table_16s_BCKD),colnames(reads_table_16s_MCKD),
                       colnames(reads_table_16s_MV1D),colnames(reads_table_16s_MRCD))

keep = metadata_16s$SampleID %in% passed_sample_list
sum(keep)
metadata_16s = metadata_16s[keep,]
reads_table_16s = reads_table_16s[,keep]

# remove mom samples that do not have paired baby
baby_participant_list = unique(metadata_16s$ParticipantID[metadata_16s$Mombaby == 'Baby'])
keep = metadata_16s$ParticipantID %in% baby_participant_list
sum(keep)
metadata_16s = metadata_16s[keep,]
reads_table_16s = reads_table_16s[,keep]

keep = colnames(reads_table_16s_MCKD) %in% colnames(reads_table_16s)
reads_table_16s_MCKD = reads_table_16s_MCKD[,keep]

keep = colnames(reads_table_16s_MRCD) %in% colnames(reads_table_16s)
reads_table_16s_MRCD = reads_table_16s_MRCD[,keep]

keep = colnames(reads_table_16s_MV1D) %in% colnames(reads_table_16s)
reads_table_16s_MV1D = reads_table_16s_MV1D[,keep]

# taxa threshold
keep = row.names(reads_table_16s) %in% taxa_list
sum(keep)
reads_table_16s = reads_table_16s[keep,]

# total reads threshold
keep = colSums(reads_table_16s) >= 5000
sum(keep)
reads_table_16s = reads_table_16s[,keep]
metadata_16s = metadata_16s[keep,]

# get taxonomy
taxonomy_16s = as.data.frame(row.names(reads_table_16s))
colnames(taxonomy_16s) = 'Taxa'
#taxonomy_16s$Taxa_2 = str_replace_all(taxonomy_16s$Taxa,'_BT','')
taxonomy_16s$Taxonomy = NA
for (a in 1: nrow(taxonomy_16s)) {
  n= which(taxonomy$V3 == taxonomy_16s$Taxa[a])
  taxonomy_16s$Taxonomy[a] = taxonomy$V2[n[1]]
}
rm(metadata_all,metadata_all_cha,metadata_all_num,metadata_16s_cha,metadata_16s_num, reads_table)

metadata_16s$BMI = metadata_16s$weight * 0.0283495 / ((metadata_16s$height * 0.0254)^2)

sum(metadata_16s$Flag == "BCKD_birth_0_day" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BCKD_birth_1_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BCKD_birth_2_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BRCD_birth_0_day" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BRCD_birth_1_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BRCD_birth_2_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BS1D_birth_0_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BS1D_birth_1_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BS1D_birth_2_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MCKD_1st_trimester" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MCKD_2nd_trimester" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MCKD_3rd_trimester" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MRCD_1st_trimester" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MRCD_2nd_trimester" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MRCD_3rd_trimester" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MV1D_1st_trimester" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MV1D_2nd_trimester" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MV1D_3rd_trimester" & !is.na(metadata_16s$Flag))

colnames(metadata_16s)[colnames(metadata_16s) == 'no_vdischarge'] = 'vdischarge_pregnancy'
metadata_16s$vdischarge_pregnancy[metadata_16s$vdischarge_pregnancy == 'No'] = 'xx'
metadata_16s$vdischarge_pregnancy[metadata_16s$vdischarge_pregnancy == 'Yes'] = 'No'
metadata_16s$vdischarge_pregnancy[metadata_16s$vdischarge_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_vodor'] = 'vodor_pregnancy'
metadata_16s$vodor_pregnancy[metadata_16s$vodor_pregnancy == 'No'] = 'xx'
metadata_16s$vodor_pregnancy[metadata_16s$vodor_pregnancy == 'Yes'] = 'No'
metadata_16s$vodor_pregnancy[metadata_16s$vodor_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_diarrhea '] = 'diarrhea_pregnancy'
metadata_16s$diarrhea_pregnancy[metadata_16s$diarrhea_pregnancy == 'No'] = 'xx'
metadata_16s$diarrhea_pregnancy[metadata_16s$diarrhea_pregnancy == 'Yes'] = 'No'
metadata_16s$diarrhea_pregnancy[metadata_16s$diarrhea_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_respiratory'] = 'respiratory_pregnancy'
metadata_16s$respiratory_pregnancy[metadata_16s$respiratory_pregnancy == 'No'] = 'xx'
metadata_16s$respiratory_pregnancy[metadata_16s$respiratory_pregnancy == 'Yes'] = 'No'
metadata_16s$respiratory_pregnancy[metadata_16s$respiratory_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_contractions'] = 'contractions_pregnancy'
metadata_16s$contractions_pregnancy[metadata_16s$contractions_pregnancy == 'No'] = 'xx'
metadata_16s$contractions_pregnancy[metadata_16s$contractions_pregnancy == 'Yes'] = 'No'
metadata_16s$contractions_pregnancy[metadata_16s$contractions_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_vitching'] = 'vitching_pregnancy'
metadata_16s$vitching_pregnancy[metadata_16s$vitching_pregnancy == 'No'] = 'xx'
metadata_16s$vitching_pregnancy[metadata_16s$vitching_pregnancy == 'Yes'] = 'No'
metadata_16s$vitching_pregnancy[metadata_16s$vitching_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_false_labor'] = 'false_labor_pregnancy'
metadata_16s$false_labor_pregnancy[metadata_16s$false_labor_pregnancy == 'No'] = 'xx'
metadata_16s$false_labor_pregnancy[metadata_16s$false_labor_pregnancy == 'Yes'] = 'No'
metadata_16s$false_labor_pregnancy[metadata_16s$false_labor_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_high_bp'] = 'high_bp_pregnancy'
metadata_16s$high_bp_pregnancy[metadata_16s$high_bp_pregnancy == 'No'] = 'xx'
metadata_16s$high_bp_pregnancy[metadata_16s$high_bp_pregnancy == 'Yes'] = 'No'
metadata_16s$high_bp_pregnancy[metadata_16s$high_bp_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_pets_owner'] = 'pets_owner_pregnancy'
metadata_16s$pets_owner_pregnancy[metadata_16s$pets_owner_pregnancy == 'No'] = 'xx'
metadata_16s$pets_owner_pregnancy[metadata_16s$pets_owner_pregnancy == 'Yes'] = 'No'
metadata_16s$pets_owner_pregnancy[metadata_16s$pets_owner_pregnancy == 'xx'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'no_vbleeding'] = 'vbleeding_pregnancy'
metadata_16s$vbleeding_pregnancy[metadata_16s$vbleeding_pregnancy == 'No'] = 'xx'
metadata_16s$vbleeding_pregnancy[metadata_16s$vbleeding_pregnancy == 'Yes'] = 'No'
metadata_16s$vbleeding_pregnancy[metadata_16s$vbleeding_pregnancy == 'xx'] = 'Yes'

colnames(metadata_16s)[colnames(metadata_16s) == 'mom_delivery_method'] = 'mom_delivery_method (C-section yes / vaginal no)'
colnames(metadata_16s)[colnames(metadata_16s) == 'Delivery'] = 'Preterm delivery'
colnames(metadata_16s)[colnames(metadata_16s) == 'baby_icu'] = 'baby_nicu'
colnames(metadata_16s)[colnames(metadata_16s) == 'baby_problems_none'] = 'baby_problems'

metadata_16s[metadata_16s == 'PTB' | metadata_16s == 'Yes_I_plan_to_breast_feed' | 
               metadata_16s == 'Yes_my_baby_was_put_in_the_ICU' |metadata_16s == 'Yes_I_have_been_diagnosed_with_hyperemesis' |
               metadata_16s == 'C_Section' |metadata_16s == 'Cesarean_section' |  
               metadata_16s == 'Female'| metadata_16s == 'Yes_98_100'| metadata_16s == 'Yes_101_103'| metadata_16s == 'Yes_membranes_ruptured_1_18_hours_before_labor'| metadata_16s == 'Yes_membranes_ruptured_before_37_weeks_of_pregnancy'] = 'Yes'

metadata_16s[metadata_16s == 'TB' | metadata_16s == 'No_I_have_not_been_exposed_to_any_chemicals' | metadata_16s == 'No_I_have_not_been_to_the_hospital' | metadata_16s == 'No_I_have_not_been_diagnosed_with_hyperemesis_and_I_do_not_have_symptoms' |
               metadata_16s == 'No_my_baby_was_not_put_in_the_ICU' | metadata_16s == 'No_I_do_not_plan_to_breast_feed' | metadata_16s == 'No_I_am_not_taking_any_new_medications' | metadata_16s == 'No_I_have_not_been_diagnosed_with_hyperemesis_but_I_have_symptoms' |
               metadata_16s == 'Vaginal' | metadata_16s == 'Vaginal_delivery' |  metadata_16s == 'No_I_have_not_had_a_fever_since_my_last_visit' |
               metadata_16s == 'Male' | metadata_16s == 'none'| metadata_16s == 'No_there_has_been_no_premature_rupture_of_membranes'] = 'No'

metadata_16s[metadata_16s == 'Not_Sure' | 
               metadata_16s == 'Not_sure' | 
               metadata_16s == 'Uncertain' |
               metadata_16s == 'Do_not_know'|
               metadata_16s == 'Nelson 6th Floor'|
               metadata_16s == 'Baby_NICU' | metadata_16s == 'NaN' | metadata_16s == 'Don_t_know'] = NA

write.csv(metadata_16s,'metadata_16s_ori.csv') ### run when output contains sample ID or participant ID

# smoker_second_hand
metadata_16s[metadata_16s == 'Never'] = 0
metadata_16s[metadata_16s == 'Rarely'] = 1
metadata_16s[metadata_16s == 'Almost_every_day'] = 2
metadata_16s[metadata_16s == 'Every_day'] = 3
  
# various worries
metadata_16s[metadata_16s == 'Much_Less'] = 1
metadata_16s[metadata_16s == 'A_Bit_Less'] = 2
metadata_16s[metadata_16s == 'About_Average'] = 3
metadata_16s[metadata_16s == 'A_Bit_More'] = 4
metadata_16s[metadata_16s == 'Much_More'] = 5

metadata_16s[metadata_16s == 'Never'] = 0
metadata_16s[metadata_16s == 'Almost_Never'] = 1
metadata_16s[metadata_16s == 'Sometimes'] = 2
metadata_16s[metadata_16s == 'Fairly_Often'] = 3
metadata_16s[metadata_16s == 'Very_Often'] = 4

# received_oral_frequency, douche_frequency
metadata_16s$received_oral_frequency[metadata_16s$received_oral_frequency == 'No'] = 0
metadata_16s$douche_frequency[metadata_16s$douche_frequency == 'No'] = 0

# sexual_partners_month
metadata_16s$sexual_partners_month[metadata_16s$sexual_partners_month == '6_to_10'] = 8

# yogurt, milk, cheese, ice_cream
metadata_16s[metadata_16s == '1_2_servings_per_day'] = 4
metadata_16s[metadata_16s == '3_4_servings_per_day'] = 5
metadata_16s[metadata_16s == '5_servings_per_day'] = 6
metadata_16s[metadata_16s == 'Less_than_1_serving_per_week'] = 1
metadata_16s[metadata_16s == '1_2_servings_per_week'] = 2
metadata_16s[metadata_16s == '3_4_servings_per_week'] = 3

# weight_change
metadata_16s[metadata_16s == 'My_weight_has_NOT_changed_since_my_last_visit'] = 0
metadata_16s[metadata_16s == 'I_have_LOST_weight_since_my_last_visit'] = -1
metadata_16s[metadata_16s == 'I_have_GAINED_weight_since_my_last_visit'] = 1

# weight_lost, weight_gained
metadata_16s[metadata_16s == '0_5_lbs'] = 1
metadata_16s[metadata_16s == '10_15_lbs'] = 3
metadata_16s[metadata_16s == '5_10_lbs'] = 2
metadata_16s[metadata_16s == 'More_than_15_lbs'] = 4

# physical_activity_vigorous, physical_activity_moderate, physical_activity_light
metadata_16s[metadata_16s == '0_times'] = 0
metadata_16s[metadata_16s == '1_2_times'] = 1
metadata_16s[metadata_16s == '3_4_times'] = 2
metadata_16s[metadata_16s == '5_6_times'] = 3
metadata_16s[metadata_16s == '7_times'] = 4

# uti_lifetime_number, 
metadata_16s[metadata_16s == '2_to_4'] = 3
metadata_16s[metadata_16s == '2_to_4'] = 3

# douche_frequency, vaginal_sex_frequency, vaginal_penetration_frequency, received_oral_frequency, gave_oral_frequency, anal_sex_frequency, vaginal_lubrication_frequency
metadata_16s[metadata_16s == '1_3_times_per_month'] = 2
metadata_16s[metadata_16s == '2_6_times_per_week'] = 4
metadata_16s[metadata_16s == 'Less_than_once_a_month'] = 1
metadata_16s[metadata_16s == 'Never'] = 0
metadata_16s[metadata_16s == 'Once_a_day'] = 5
metadata_16s[metadata_16s == 'Once_a_week'] = 3

# income
metadata_16s[metadata_16s == '15_000_19_999'] = 2
metadata_16s[metadata_16s == '20_000_39_999'] = 3
metadata_16s[metadata_16s == '40_000_59_999'] = 4
metadata_16s[metadata_16s == '60_000_79_999'] = 5
metadata_16s[metadata_16s == '80_000_or_more'] = 6
metadata_16s[metadata_16s == 'Under_15_000'] = 1

# education
metadata_16s[metadata_16s == '2_year_College_Degree'] = 3
metadata_16s[metadata_16s == '4_year_College_Degree'] = 4
metadata_16s[metadata_16s == 'Doctoral_or_Professional_Degree'] = 6
metadata_16s[metadata_16s == 'High_School_GED'] = 2
metadata_16s[metadata_16s == 'Less_than_High_School'] = 1
metadata_16s[metadata_16s == 'Masters_Degree'] = 5
metadata_16s[metadata_16s == 'Some_College'] = NA

# prenatal_care_start
metadata_16s[metadata_16s == '0_4_weeks_after_conception'] = 3
metadata_16s[metadata_16s == '4_weeks_to_the_end_of_the_1st_trimester'] = 2
metadata_16s[metadata_16s == '4_weeks_to_the_end_of_the_first_trimester'] = 2
metadata_16s[metadata_16s == 'Before_conception'] = 4
metadata_16s[metadata_16s == 'During_the_2nd_trimester'] = 1
metadata_16s[metadata_16s == 'During_the_second_trimester'] = 1
metadata_16s$SampleID_ori = NULL

 write.csv(reads_table_16s_BCKD,'reads_table_16s_BCKD.csv')  ### run when output contains sample ID or participant ID
 write.csv(reads_table_16s_BRCD,'reads_table_16s_BRCD.csv')  ### run when output contains sample ID or participant ID
 write.csv(reads_table_16s_BS1D,'reads_table_16s_BS1D.csv')  ### run when output contains sample ID or participant ID
 write.csv(reads_table_16s_MCKD,'reads_table_16s_MCKD.csv')  ### run when output contains sample ID or participant ID
 write.csv(reads_table_16s_MRCD,'reads_table_16s_MRCD.csv')  ### run when output contains sample ID or participant ID
 write.csv(reads_table_16s_MV1D,'reads_table_16s_MV1D.csv')  ### run when output contains sample ID or participant ID
 write.csv(reads_table_16s,'reads_table_16s_6_microbiome_combination.csv')  ### run when output contains sample ID or participant ID
 write.csv(metadata_16s,'metadata_16s.csv') ### run when output contains sample ID or participant ID
 write.csv(reads_table_cytokines,'reads_table_cytokines.csv') ### run when output contains sample ID or participant ID
 write.csv(reads_table_lipid,'reads_table_lipid.csv') ### run when output contains sample ID or participant ID

##### get paired sample ID and output to the server #####
type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby)
  keep = metadata_16s$SampleType == type1_all[a] & 
    !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
  sum(keep)
  
  metadata_baby = metadata_16s[keep,]
  metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
  keep = duplicated(metadata_baby$ParticipantID)
  metadata_baby =metadata_baby[!keep,]
  
  # MV1D
  Participant_list = metadata_baby$ParticipantID
  metadata_mom_MV = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
  metadata_mom_MV = metadata_mom_MV[metadata_mom_MV$Mombaby == 'Mom',]
  metadata_mom_MV = metadata_mom_MV[order(metadata_mom_MV$VisitNum, decreasing = T),]
  metadata_mom_MV = metadata_mom_MV[order(metadata_mom_MV$ParticipantID),]
  keep = !duplicated(metadata_mom_MV$ParticipantID)
  metadata_mom_MV = metadata_mom_MV[keep,]      # mom's metadata
  
  keep = metadata_baby$ParticipantID %in% metadata_mom_MV$ParticipantID
  metadata_baby = metadata_baby[keep,]
  
  metadata_mom_MV = metadata_mom_MV[order(metadata_mom_MV$ParticipantID),]
  metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]
  

  # MRCD
  Participant_list = metadata_baby$ParticipantID
  metadata_mom_MR = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MRCD',]
  metadata_mom_MR = metadata_mom_MR[metadata_mom_MR$Mombaby == 'Mom',]
  metadata_mom_MR = metadata_mom_MR[order(metadata_mom_MR$VisitNum, decreasing = T),]
  metadata_mom_MR = metadata_mom_MR[order(metadata_mom_MR$ParticipantID),]
  keep = !duplicated(metadata_mom_MR$ParticipantID)
  metadata_mom_MR = metadata_mom_MR[keep,]      # mom's metadata
  
  keep = metadata_baby$ParticipantID %in% metadata_mom_MR$ParticipantID
  metadata_baby = metadata_baby[keep,]
  metadata_mom_MV = metadata_mom_MV[keep,]
  
  metadata_mom_MR = metadata_mom_MR[order(metadata_mom_MR$ParticipantID),]

  
  # MCKD
  Participant_list = metadata_baby$ParticipantID
  metadata_mom_MC = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MCKD',]
  metadata_mom_MC = metadata_mom_MC[metadata_mom_MC$Mombaby == 'Mom',]
  metadata_mom_MC = metadata_mom_MC[order(metadata_mom_MC$VisitNum, decreasing = T),]
  metadata_mom_MC = metadata_mom_MC[order(metadata_mom_MC$ParticipantID),]
  keep = !duplicated(metadata_mom_MC$ParticipantID)
  metadata_mom_MC = metadata_mom_MC[keep,]      # mom's metadata
  
  keep = metadata_baby$ParticipantID %in% metadata_mom_MC$ParticipantID
  metadata_baby = metadata_baby[keep,]
  metadata_mom_MV = metadata_mom_MV[keep,]
  metadata_mom_MR = metadata_mom_MR[keep,]
  
  metadata_mom_MC = metadata_mom_MC[order(metadata_mom_MC$ParticipantID),]
  
  metadata_baby = data.frame(metadata_baby$SampleID,type1_all[a])
  colnames(metadata_baby) = c('SampleID','Type')
  metadata_mom_MC = data.frame(metadata_mom_MC$SampleID,'MCDK')
  colnames(metadata_mom_MC) = c('SampleID','Type')
  metadata_mom_MR = data.frame(metadata_mom_MR$SampleID,'MRCD')
  colnames(metadata_mom_MR) = c('SampleID','Type')
  metadata_mom_MV = data.frame(metadata_mom_MV$SampleID,'MV1D')
  colnames(metadata_mom_MV) = c('SampleID','Type')
    
  output = rbind(metadata_baby,metadata_mom_MC,metadata_mom_MR,metadata_mom_MV)
  
  setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/DADA2/mom_baby')
  write.table(output, paste0('SampleID_',type1_all[a],'.txt'), row.names = F, quote = F, col.names = F)
  setwd('/Users/binzhu/Desktop/Mom-baby/')
}






