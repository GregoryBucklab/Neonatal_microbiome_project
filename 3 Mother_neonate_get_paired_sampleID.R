library(outliers)

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

metadata_16s$Flag[metadata_16s$SampleType == "MV1D"] = 'MV'
metadata_16s$Flag[metadata_16s$SampleType == "MRCD"] = 'MR'
metadata_16s$Flag[metadata_16s$SampleType == "MCKD"] = 'MB'

keep = !is.na(metadata_16s$Flag)
reads_table_16s = reads_table_16s[,keep]
metadata_16s = metadata_16s[keep,]

metadata_16s$SampleID_ori = metadata_16s$SampleID

metadata_16s$SampleID_2 = metadata_16s$SampleID

participant_id_unique = unique(metadata_16s$ParticipantID)

keep = colSums(reads_table_16s) >= 5000
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

rm(metadata_all,metadata_all_cha,metadata_all_num,metadata_16s_cha,metadata_16s_num, reads_table)

# select maternal microbiomes at the last visit
keep = !is.na(metadata_16s$Flag)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]
metadata$Flag[metadata$Flag == 'MV1D_3rd_trimester'] = 'MV'
metadata$Flag[metadata$Flag == 'MV1D_2nd_trimester'] = 'MV'
metadata$Flag[metadata$Flag == 'MV1D_1st_trimester'] = 'MV'
metadata$Flag[metadata$Flag == 'MRCD_3rd_trimester'] = 'MR'
metadata$Flag[metadata$Flag == 'MRCD_2nd_trimester'] = 'MR'
metadata$Flag[metadata$Flag == 'MRCD_1st_trimester'] = 'MR'
metadata$Flag[metadata$Flag == 'MCKD_3rd_trimester'] = 'MB'
metadata$Flag[metadata$Flag == 'MCKD_2nd_trimester'] = 'MB'
metadata$Flag[metadata$Flag == 'MCKD_1st_trimester'] = 'MB'

metadata_2 = metadata
metadata = metadata[metadata$Mombaby == "Baby",]
metadata_2 = metadata_2[metadata_2$Mombaby == "Mom",]

metadata_2$keep = paste(metadata_2$ParticipantID, metadata_2$SampleType, sep = '_')
metadata_2 = metadata_2[order(metadata_2$VisitNum, decreasing = T),]
keep = duplicated(metadata_2$keep)
sum(keep)
metadata_2 = metadata_2[!keep,]

output_table = data.frame(participant = unique(metadata_2$ParticipantID), MCKD= NA, MRCD=NA, MV1D = NA,
                          BCKD_birth_0_day = NA, BCKD_birth_1_days = NA,BCKD_birth_2_days = NA,
                          BRCD_birth_0_day = NA, BRCD_birth_1_days = NA,BRCD_birth_2_days = NA,
                          BS1D_birth_0_day = NA, BS1D_birth_1_days = NA,BS1D_birth_2_days = NA)
for (a in 1:nrow(metadata_2)) {
  x = which(output_table$participant == metadata_2$ParticipantID[a])
  y = which(colnames(output_table) == metadata_2$SampleType[a])
  output_table[x,y] = metadata_2$SampleID[a]
}
for (a in 1:nrow(metadata)) {
  x = which(output_table$participant == metadata$ParticipantID[a])
  y = which(colnames(output_table) == metadata$Flag[a])
  output_table[x,y] = metadata$SampleID[a]
}

n = which(output_table$participant == 'P38290')  
output_table = output_table[-n,]

sample_list = output_table
sample_list$participant = NULL
sample_list = gather(output_table)
sample_list = sample_list[sample_list$key != 'participant',]
sample_list = sample_list[!is.na(sample_list$value),]

reads_table_16s = reads_table_16s[,colnames(reads_table_16s) %in% sample_list$value]
metadata_16s = metadata_16s[metadata_16s$SampleID %in% sample_list$value, ]

reads_table_cytokines = reads_table_cytokines[,colnames(reads_table_cytokines) %in% sample_list$value]
reads_table_lipid = reads_table_lipid[,colnames(reads_table_lipid) %in% sample_list$value]

setwd('/Users/binzhu/secure/godel/gpfs_fs/home/bzhu/dbGAP_mom_baby')
sample_list = sample_list$value
write.table(sample_list,'Mom_baby_16s_sample_list.txt', row.names = F, quote = F, col.names = F)
write.table(colnames(reads_table_lipid),'Mom_baby_lipid_sample_list.txt', row.names = F, quote = F, col.names = F)
write.table(colnames(reads_table_cytokines),'Mom_baby_cytokine_sample_list.txt', row.names = F, quote = F, col.names = F)
x = colnames(metadata_2_imputation)[c(16:ncol(metadata_2_imputation))]
x = x[-c(4,127)]
write.table(x,'Mom_baby_metadata_list.txt', row.names = F, quote = F, col.names = F)
setwd('/Users/binzhu/Desktop/Mom-baby/')
