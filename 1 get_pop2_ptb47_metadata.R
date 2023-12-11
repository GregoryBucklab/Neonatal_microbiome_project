####### library ##########
library(stringr)
library(parallel)
#library(FEAST)
######################################## PTB47 metadata ################################################
####### get PTB47 mapping metadata #######################
# get mapping data
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/mappings/')
metadata <- read.table('ptb47.Participant-Kit-Sample-Mapping.internal.20170724.txt', header =T)
setwd('/Users/binzhu/Desktop/Mom-baby')

metadata$SampleID <- as.character(metadata$SampleID )

metadata$SampleType <- as.character(metadata$SampleType )

metadata$ParticipantID_Mombaby = paste(metadata$ParticipantID,substr(metadata$SampleID, start = 1, stop = 1), sep= "_")
metadata$Mombaby = substr(metadata$SampleID, start = 1, stop = 1)
metadata$Mombaby[metadata$Mombaby == 'M'] = 'Mom'
metadata$Mombaby[metadata$Mombaby == 'B'] = 'Baby'



####### get PTB47 tags information ######################
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/tags/')
data <- read.delim('ptb47.subjectsPlusKits.20170724.txt')
setwd('/Users/binzhu/Desktop/Mom-baby')
data = data[!duplicated(data$ParticipantID),]

for (a in 1: nrow(metadata)) {
  n =which(data$ParticipantID == metadata$ParticipantID[a])
  metadata$Group[a] = data$Group[n]
}

col_number <- ncol(metadata)

####### get PTB47 clinical information ####### 
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/clinical/')
data <- read.delim('ptb47.clinicalData.20170825.txt')
setwd('/Users/binzhu/Desktop/Mom-baby')

sample_name <- unique(data$kitId)

clinical <- unique(data$VariableName)

reads_table <- matrix(0, nrow = length(sample_name), ncol = length(clinical))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = sample_name 
colnames(reads_table) =clinical 

data$kitId <- as.character(data$kitId)
data$VariableName <- as.character(data$VariableName)
data$Value <- as.character(data$Value)

for (a in 1: dim(data)[1]) {
  
  kitID_num <- which(sample_name == data[a,2])
  clinical_num <- which(clinical == data[a,5])
  
  reads_table[kitID_num,clinical_num] = data$Value[a]
  
}

reads_table <- as.data.frame(reads_table,stringsAsFactors=FALSE)

reads_table$KitID <- row.names(reads_table)
#reads_table %>% mutate_all(as.character)
#metadata <- merge(metadata,reads_table , 'KitID')
data <- matrix(data = NA, ncol= ncol(reads_table), nrow= nrow(metadata))
colnames(data) <- colnames(reads_table)

metadata <- cbind(metadata,data)

for (a in 1: nrow(reads_table)) {
  n=which(metadata$KitID == reads_table$KitID[a])
  
  if (length(n) == 0) {
    next
  }
  
  metadata[n,(col_number+1):ncol(metadata)] = reads_table[a,]
}

metadata[,ncol(metadata)] <- NULL
col_number <- ncol(metadata)






metadata$clinic<- as.character(metadata$clinic)

# manager clinic data 
metadata[metadata == 0] = NA

metadata <- metadata[order(metadata$VisitNum, decreasing = T),]
metadata <- metadata[order(metadata$ParticipantID_Mombaby),]

unique_name <- unique(metadata$ParticipantID_Mombaby)
metadata_2<- matrix(data=NA, nrow =length(unique_name), ncol = 22)
row.names(metadata_2) = unique_name
colnames(metadata_2) = c('average_height','average_weight','average_pulse','average_systolic','average_diastolic','average_vaginal_ph',
                         'average_reported_ga','average_day_real2birth','average_true_ga','average_ga_at_delivery','average_fundal_height',
                         'last_height','last_weight','last_pulse','last_systolic','last_diastolic','last_vaginal_ph',
                         'last_reported_ga','last_day_real2birth','last_true_ga','last_ga_at_delivery','last_fundal_height')
metadata_2 <- as.data.frame(metadata_2)

for (a in 1:nrow(metadata_2)) {
  n =which(metadata$ParticipantID_Mombaby  == row.names(metadata_2)[a] )
  
  for (b in 1:11){
    metadata_2[a,b] = mean(as.numeric(as.character(metadata[n,b+11])), na.rm = TRUE)
    
    for (c in 1: length(n)) {
      if (!is.na(metadata[n[c],b+11])) {
        metadata_2[a,b+11] = metadata[n[c],b+11]
      }
    }
  }
}

metadata_2$ParticipantID_Mombaby = row.names(metadata_2)
metadata_2 = metadata_2[,c(-7,-8,-9,-10,-19,-20)]

data <- matrix(data = NA, ncol= ncol(metadata_2), nrow= nrow(metadata))
colnames(data) <- colnames(metadata_2)

metadata <- cbind(metadata,data)

for (a in 1: nrow(metadata_2)) {
  n=which(metadata$ParticipantID_Mombaby == metadata_2$ParticipantID_Mombaby[a])
  
  if (length(n) == 0) {
    next
  }
  
  metadata[n,(col_number+1):ncol(metadata)] = metadata_2[a,]
}

metadata[,ncol(metadata)] <- NULL
col_number <- ncol(metadata)

#metadata[,c(7,12:38)] <- lapply(metadata[,c(7,12:38)], as.numeric)

####### get PTB47 survey information ###########
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/survey/')
data <- read.delim('ptb47.surveyData.20170927.txt')
setwd('/Users/binzhu/Desktop/Mom-baby')

sample_name <- unique(data$KitID)
survey <- unique(data$VariableName)

reads_table <- matrix(NA, nrow = length(sample_name), ncol = length(survey))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = sample_name 
colnames(reads_table) =survey 

data$Value <- as.character(data$Value)
data$KitID <- as.character(data$KitID)
data$VariableName <- as.character(data$VariableName)

for (a in 1: nrow(data)) {
  kitID_num <- which(sample_name == data[a,2])
  survey_num <- which(survey == data[a,5])
  reads_table[kitID_num,survey_num] = data[a,6]
}


survey_list_for_not_repeat = read.table('survey_list_for_not_repeat.txt')

keep = colnames(reads_table) %in% survey_list_for_not_repeat$V1

reads_table_1 = reads_table[,keep]  ############ not repeat 
reads_table_2 = reads_table[,!keep]  ##############  repeat 

###
data <- matrix(data = NA, ncol= ncol(reads_table_1), nrow= nrow(metadata))
colnames(data) <- colnames(reads_table_1)

metadata <- cbind(metadata,data)

for (a in 1: nrow(reads_table_1)) {
  n=which(metadata$KitID == row.names(reads_table_1)[a])
  
  if (length(n) == 0) {
    next
  }
  
  for (b in 1: length(n)) {
    metadata[n[b],(col_number+1):ncol(metadata)] = reads_table_1[a,]
  }
  
}

col_number <- ncol(metadata)

###
data <- matrix(data = NA, ncol= ncol(reads_table_2), nrow= nrow(metadata))
colnames(data) <- colnames(reads_table_2)

for (a in 1: ncol(reads_table_2)) {
  for (b in 1: nrow(reads_table_2)) {
    if (is.na(reads_table_2[b,a])) {
      next
    }
    
    m = which(metadata$KitID == row.names(reads_table_2)[b])
    n = which(metadata$ParticipantID == metadata$ParticipantID[m[1]])
    
    data[n,a] = reads_table_2[b,a]
  }
}

metadata <- cbind(metadata,data)

col_number <- ncol(metadata)

###
metadata$mom_delivery_method[metadata$mom_delivery_method == 'C_Section'] = 'Cesarean_section'
metadata$mom_delivery_method[metadata$mom_delivery_method == 'C_section'] = 'Cesarean_section'
metadata$mom_delivery_method[metadata$mom_delivery_method == 'Not_Sure'] = 'Do_not_know'
metadata$mom_delivery_method[is.na(metadata$mom_delivery_method)] = 'Do_not_know'

metadata_PTB <- metadata 
rm(metadata, reads_table_1, reads_table_2, reads_table, survey_list_for_not_repeat)



######################################## TB metadata ################################################
####### get TB mapping metadata #######################
# get mapping data
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/mappings/')
metadata <- read.table('pop2.Participant-Kit-Sample-Mapping.internal.20170605.txt', header =T)
setwd('/Users/binzhu/Desktop/Mom-baby')

metadata$SampleID <- as.character(metadata$SampleID )

metadata$SampleType <- as.character(metadata$SampleType )

metadata$ParticipantID_Mombaby = paste(metadata$ParticipantID,substr(metadata$SampleID, start = 1, stop = 1), sep= "_")
metadata$Mombaby = substr(metadata$SampleID, start = 1, stop = 1)
metadata$Mombaby[metadata$Mombaby == 'M'] = 'Mom'
metadata$Mombaby[metadata$Mombaby == 'B'] = 'Baby'

####### get TB tags information ######################
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/tags/')
data <- read.delim('pop2.SubjectsPlusKits.20160617.txt')
setwd('/Users/binzhu/Desktop/Mom-baby')
data = data[!duplicated(data$ParticipantID),]
data$Group = NA

for (a in 1: nrow(metadata)) {
  n =which(data$ParticipantID == metadata$ParticipantID[a])
  metadata$Group[a] = data$Group[n]
}

col_number <- ncol(metadata)

####### get TB clinical information ####### 
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/clinical/')
data <- read.delim('pop2.clinicalData.20170928.txt')
setwd('/Users/binzhu/Desktop/Mom-baby')

sample_name <- unique(data$kitId)

clinical <- unique(data$VariableName)

reads_table <- matrix(0, nrow = length(sample_name), ncol = length(clinical))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = sample_name 
colnames(reads_table) =clinical 

data$kitId <- as.character(data$kitId)
data$VariableName <- as.character(data$VariableName)
data$Value <- as.character(data$Value)

for (a in 1: dim(data)[1]) {
  
  kitID_num <- which(sample_name == data[a,2])
  clinical_num <- which(clinical == data[a,5])
  
  reads_table[kitID_num,clinical_num] = data$Value[a]
  
}

reads_table <- as.data.frame(reads_table,stringsAsFactors=FALSE)

reads_table$KitID <- row.names(reads_table)
#reads_table %>% mutate_all(as.character)
#metadata <- merge(metadata,reads_table , 'KitID')
data <- matrix(data = NA, ncol= ncol(reads_table), nrow= nrow(metadata))
colnames(data) <- colnames(reads_table)

metadata <- cbind(metadata,data)

for (a in 1: nrow(reads_table)) {
  n=which(metadata$KitID == reads_table$KitID[a])
  
  if (length(n) == 0) {
    next
  }
  
  metadata[n,(col_number+1):ncol(metadata)] = reads_table[a,]
}

metadata[,ncol(metadata)] <- NULL
col_number <- ncol(metadata)






metadata$clinic<- as.character(metadata$clinic)

# manager clinic data 
metadata[metadata == 0] = NA

metadata <- metadata[order(metadata$VisitNum, decreasing = T),]
metadata <- metadata[order(metadata$ParticipantID_Mombaby),]

unique_name <- unique(metadata$ParticipantID_Mombaby)
metadata_2<- matrix(data=NA, nrow =length(unique_name), ncol = 22)
row.names(metadata_2) = unique_name
colnames(metadata_2) = c('average_height','average_weight','average_pulse','average_systolic','average_diastolic','average_vaginal_ph',
                         'average_reported_ga','average_day_real2birth','average_true_ga','average_ga_at_delivery','average_fundal_height',
                         'last_height','last_weight','last_pulse','last_systolic','last_diastolic','last_vaginal_ph',
                         'last_reported_ga','last_day_real2birth','last_true_ga','last_ga_at_delivery','last_fundal_height')
metadata_2 <- as.data.frame(metadata_2)

for (a in 1:nrow(metadata_2)) {
  n =which(metadata$ParticipantID_Mombaby  == row.names(metadata_2)[a] )
  
  for (b in 1:11){
    metadata_2[a,b] = mean(as.numeric(as.character(metadata[n,b+11])), na.rm = TRUE)
    
    for (c in 1: length(n)) {
      if (!is.na(metadata[n[c],b+11])) {
        metadata_2[a,b+11] = metadata[n[c],b+11]
      }
    }
  }
}

metadata_2$ParticipantID_Mombaby = row.names(metadata_2)
metadata_2 = metadata_2[,c(-7,-8,-9,-10,-19,-20)]

data <- matrix(data = NA, ncol= ncol(metadata_2), nrow= nrow(metadata))
colnames(data) <- colnames(metadata_2)

metadata <- cbind(metadata,data)

for (a in 1: nrow(metadata_2)) {
  n=which(metadata$ParticipantID_Mombaby == metadata_2$ParticipantID_Mombaby[a])
  
  if (length(n) == 0) {
    next
  }
  
  metadata[n,(col_number+1):ncol(metadata)] = metadata_2[a,]
}

metadata[,ncol(metadata)] <- NULL
col_number <- ncol(metadata)

#metadata[,c(7,12:38)] <- lapply(metadata[,c(7,12:38)], as.numeric)

####### get TB survey information ###########
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/survey/')
data <- read.delim('pop2.surveyData.20170927.txt')
setwd('/Users/binzhu/Desktop/Mom-baby')

sample_name <- unique(data$KitID)
survey <- unique(data$VariableName)

reads_table <- matrix(NA, nrow = length(sample_name), ncol = length(survey))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = sample_name 
colnames(reads_table) =survey 

data$Value <- as.character(data$Value)
data$KitID <- as.character(data$KitID)
data$VariableName <- as.character(data$VariableName)

for (a in 1: nrow(data)) {
  kitID_num <- which(sample_name == data[a,2])
  survey_num <- which(survey == data[a,5])
  reads_table[kitID_num,survey_num] = data[a,6]
}


survey_list_for_not_repeat = read.table('survey_list_for_not_repeat.txt')

keep = colnames(reads_table) %in% survey_list_for_not_repeat$V1

reads_table_1 = reads_table[,keep]  ############ not repeat 
reads_table_2 = reads_table[,!keep]  ##############  repeat 

###
data <- matrix(data = NA, ncol= ncol(reads_table_1), nrow= nrow(metadata))
colnames(data) <- colnames(reads_table_1)

metadata <- cbind(metadata,data)

for (a in 1: nrow(reads_table_1)) {
  n=which(metadata$KitID == row.names(reads_table_1)[a])
  
  if (length(n) == 0) {
    next
  }
  
  for (b in 1: length(n)) {
    metadata[n[b],(col_number+1):ncol(metadata)] = reads_table_1[a,]
  }
  
}

col_number <- ncol(metadata)

###
data <- matrix(data = NA, ncol= ncol(reads_table_2), nrow= nrow(metadata))
colnames(data) <- colnames(reads_table_2)

for (a in 1: ncol(reads_table_2)) {
  for (b in 1: nrow(reads_table_2)) {
    if (is.na(reads_table_2[b,a])) {
      next
    }
    
    m = which(metadata$KitID == row.names(reads_table_2)[b])
    n = which(metadata$ParticipantID == metadata$ParticipantID[m[1]])
    
    data[n,a] = reads_table_2[b,a]
  }
}

metadata <- cbind(metadata,data)

col_number <- ncol(metadata)

###
metadata$mom_delivery_method[metadata$mom_delivery_method == 'C_Section'] = 'Cesarean_section'
metadata$mom_delivery_method[metadata$mom_delivery_method == 'C_section'] = 'Cesarean_section'
metadata$mom_delivery_method[metadata$mom_delivery_method == 'Not_Sure'] = 'Do_not_know'
metadata$mom_delivery_method[is.na(metadata$mom_delivery_method)] = 'Do_not_know'

metadata_TB <- metadata 
rm(metadata, reads_table_1, reads_table_2, reads_table, survey_list_for_not_repeat)



######################################## combind TB and PTB and make metadata columns as numeric if only contains numbers and NA ########
metadata_all = rbind(metadata_PTB,metadata_TB)

keep = duplicated(metadata_all$SampleID)
sum(keep)
metadata_all = metadata_all[!keep,]

metadata_num = as.data.frame(matrix(data = NA, nrow = nrow(metadata_all), ncol =0))
metadata_cha = as.data.frame(matrix(data = NA, nrow = nrow(metadata_all), ncol =0))
n =1
m =1
metadata_all[metadata_all == NaN] = NA
metadata_all[,1:ncol(metadata_all)] <- lapply(metadata_all[,1:ncol(metadata_all)],as.character)

c = c('0','1','2','3','4','5','6','7','8','9','.')

for (a in 1:ncol(metadata_all)) {
  b = metadata_all[,a]
  b = b[!is.na(b)]
  b = strsplit(b,'*')
  b = unlist(b)
  b = unique(b)
  keep = b %in% c
  
  if (sum(!keep) == 0) {
    metadata_all[,a] <- as.numeric(metadata_all[,a])
    metadata_num <- cbind(metadata_num,metadata_all[,a] )
    colnames(metadata_num)[n] = colnames(metadata_all)[a]
    n=n+1
  } else {
    metadata_cha <- cbind(metadata_cha,metadata_all[,a])
    colnames(metadata_cha)[m] = colnames(metadata_all)[a]
    m=m+1
  }
}

metadata_all[,12:39][metadata_all[,12:39] == 0 & !is.na(metadata_all[,12:39])] = NA
metadata_num[,2:29][metadata_num[,2:29] == 0 & !is.na(metadata_num[,2:29])] = NA

metadata_all$weeks_pregnant[metadata_all$weeks_pregnant > 42] = NA

metadata_all_num = metadata_num
metadata_all_cha = metadata_cha
metadata_all_num$weeks_pregnant[metadata_all_num$weeks_pregnant > 42] = NA

rm(alpha.shannon,alpha.shannon_diversity,cts,data,metadata_2,metadata_PTB,metadata_TB,reads_table,
   reads_table_cytokines_PTB,reads_table_cytokines_TB,shannon_rarefaction,threshold,metadata_num,metadata_cha)
####### get 16s reads table and metadata #########
setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/New_16s_rRNA_database')
taxonomy = read.table('gg_HOMD.V1_V3.NR97.txt',sep ='\t')

setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/new_reads_by_new_db')
data = read.table('16s.momspi.new_reads.txt', header = T, sep = '\t')

setwd('/Users/binzhu/Desktop/Mom-baby')
data = data[data$Threshold.Status == 'AT',]
data = data[data$SampleID %in% metadata_all$SampleID,]

taxa_list = unique(data$Taxa)
taxa_list = taxa_list[order(taxa_list)]
sample_list = unique(data$SampleID)

{
  trials = c(1: length(sample_list))
  
  func_1 = function(trial) {
    keep = data$SampleID == sample_list[trial]
    data_2 = data[keep,]
    data_2 = data_2[,c(3,5)]
    keep = duplicated(data_2$Taxa)
    
    if (sum(keep) >0) {
      data_2 = data_2[order(data_2$Taxa),]
      keep = duplicated(data_2$Taxa)
      keep = which(keep)
      data_2$No_of_Reads[keep-1] = data_2$No_of_Reads[keep-1] + data_2$No_of_Reads[keep]
      data_2 = data_2[-keep,]
    }
    
    keep = taxa_list %in% data_2$Taxa
    sum(keep)
    data_3 = taxa_list[!keep]
    data_3 = as.data.frame(data_3)
    colnames(data_3) = 'Taxa'
    data_3$No_of_Reads = 0
    
    data_2 = rbind(data_2,data_3)
    data_2 = data_2[order(data_2$Taxa),]
    colnames(data_2) = c('Taxa',sample_list[trial])
    data_2$Taxa = NULL
    return(data_2)
  }
  
  data_3 = mclapply(trials, func_1, mc.cores = 8)
  
  
  reads_table_16s = as.data.frame(matrix(data = 0, ncol = length(sample_list) , nrow = length(taxa_list)))
  colnames(reads_table_16s) = sample_list
  row.names(reads_table_16s) = taxa_list
  
  for (a in 1: length(data_3)) {
    data_4 = data_3[a]
    data_3[a] = NA
    data_4 = unlist(data_4)
    reads_table_16s[,a] = as.numeric(as.character(data_4))
  }
}

# get 16s metadata
keep = metadata_all$SampleID %in% colnames(reads_table_16s); sum(keep)
metadata_16s = metadata_all[keep,]

keep = colnames(reads_table_16s) %in% metadata_all$SampleID
reads_table_16s = reads_table_16s[,keep]
reads_table_16s = reads_table_16s[metadata_16s$SampleID]

metadata_16s$Group[metadata_16s$Group == 'control'] = 'No'
metadata_16s$Group[metadata_16s$Group == 'preterm'] = 'Yes'
colnames(metadata_16s)[colnames(metadata_16s) == 'Group'] = 'Preterm'

rm(data,data_2,data_3, data_list)


# get 16s reads table and taxonomy
keep = rowSums(reads_table_16s) != 0;sum(keep)
reads_table_16s = reads_table_16s[keep,]

data = as.data.frame(row.names(reads_table_16s))
colnames(data) = 'rowname_16s'

data$V2 = str_remove_all(data$rowname_16s,'_AT')
data$V2 = str_remove_all(data$V2,'_BT')
data$V3 = str_sub(data$rowname_16s,-2,-1)
data$Taxa_name = NA
data$Taxonomy = NA
data$reads = rowSums(reads_table_16s)

keep = data$V3 == 'BT'; sum(keep)
data = data[!keep,]
reads_table_16s = reads_table_16s[!keep,]
length(unique(data$V2))

keep = taxonomy$V1 %in% data$V2
sum(keep)
data_2 = taxonomy[keep,]
data_2 = data_2[match(data$V2,data_2$V1),]

taxa_list = unique(data_2$V3)
reads_table_16s$taxa = data_2$V3

data = as.data.frame(table(reads_table_16s$taxa))
reads_table_16s_1 = reads_table_16s[reads_table_16s$taxa %in% data$Var1[data$Freq ==1],]
reads_table_16s_2 = reads_table_16s[reads_table_16s$taxa %in% data$Var1[data$Freq >1],]

{
  data_list <- split(reads_table_16s_2, f = reads_table_16s_2$taxa)
  
  taxa_list = names(data_list) 
  length(taxa_list)
  
  trials = c(1: length(taxa_list))
  func_1 = function(trial) {
    data = data_list[[trial]]
    data$taxa = NULL
    data = colSums(data)
    data = as.numeric(data)

    return(data)
  }
  data = mclapply(trials, func_1, mc.cores = 8)
  
  reads_table_16s_2$taxa = NULL
  for (a in 1:length(taxa_list)) {
    reads_table_16s_2[a,] = data[[a]]
  }
  reads_table_16s_2 = reads_table_16s_2[c(1:length(taxa_list)),]
  row.names(reads_table_16s_2) = taxa_list
  
  row.names(reads_table_16s_1) = reads_table_16s_1$taxa
  reads_table_16s_1$taxa = NULL
  
  reads_table_16s = rbind(reads_table_16s_1,reads_table_16s_2)
}

row.names(reads_table_16s) = str_replace_all(row.names(reads_table_16s), '_clade','_OTU')
taxonomy$V3 = str_replace_all(taxonomy$V3, '_clade','_OTU')

taxonomy_16s = taxonomy[taxonomy$V3 %in% row.names(reads_table_16s),]
taxonomy_16s = taxonomy_16s[!duplicated(taxonomy_16s$V3),]
rm(x,reads_table,data,data_2,data_list,output,reads_table_16s_1,reads_table_16s_2, metadata_16s_cha,metadata_16s_num,metadata_all,metadata_all_cha,metadata_all_num, taxonomy)

####### get cytokines #########
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/cytokines')
data <- read.table('ptb47.cytokines.20171103.txt',sep='\t', header = T)
setwd('/Users/binzhu/Desktop/mom-baby')
sample_name <- unique(data$SampleID)
species_name <- unique(data$Cytokine)

reads_table <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = species_name 
colnames(reads_table) = sample_name

for (a in 1: nrow(data)) {
  
  column_num <- which(sample_name == data$SampleID[a])
  row_num <- which(species_name == data$Cytokine[a])
  
  reads_table[row_num,column_num] = as.numeric(as.character(data$NormConc.pg_cyto.mg_prot.[a]))
  
}
reads_table_cytokines_PTB = as.data.frame(reads_table)
colnames(reads_table_cytokines_PTB) = str_replace_all(colnames(reads_table_cytokines_PTB),'MVAX','MV1D')
reads_table_cytokines_PTB = reads_table_cytokines_PTB[,!is.na(colSums(reads_table_cytokines_PTB))]

keep = c("IL-1b","Eotaxin","IL-8","TNF-a","IL-17A","MIP-1b","IL-6","IP-10", "RANTES")
reads_table_cytokines_PTB = reads_table_cytokines_PTB[row.names(reads_table_cytokines_PTB) %in% keep,]

rm(metadata_2,data)


setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/cytokines')
data <- read.table('pop2.cytokines.20171019.txt',sep='\t', header = T)
setwd('/Users/binzhu/Desktop/mom-baby')
sample_name <- unique(data$SampleID)
species_name <- unique(data$Cytokine)

reads_table <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = species_name 
colnames(reads_table) = sample_name

for (a in 1: nrow(data)) {
  
  column_num <- which(sample_name == data$SampleID[a])
  row_num <- which(species_name == data$Cytokine[a])
  
  reads_table[row_num,column_num] = as.numeric(as.character(data$NormConc.pg_cyto.mg_prot.[a]))
  
}
reads_table_cytokines_TB = as.data.frame(reads_table)
colnames(reads_table_cytokines_TB) = str_replace_all(colnames(reads_table_cytokines_TB),'MVAX','MV1D')
reads_table_cytokines_TB = reads_table_cytokines_TB[,!is.na(colSums(reads_table_cytokines_TB))]

keep = c("IL-1b","Eotaxin","IL-8","TNF-a","IL-17A","MIP-1b","IL-6","IP-10", "RANTES")
reads_table_cytokines_TB = reads_table_cytokines_TB[row.names(reads_table_cytokines_TB) %in% keep,]

rm(metadata_2,data)


reads_table_cytokines_PTB = reads_table_cytokines_PTB[order(row.names(reads_table_cytokines_PTB)),]
reads_table_cytokines_TB = reads_table_cytokines_TB[order(row.names(reads_table_cytokines_TB)),]
reads_table_cytokines = cbind(reads_table_cytokines_PTB,reads_table_cytokines_TB)

keep = duplicated(colnames(reads_table_cytokines))
sum(keep)
reads_table_cytokines = reads_table_cytokines[,!keep]

rm(alpha.shannon,alpha.shannon_diversity,cts,data,metadata_2,metadata_PTB,metadata_TB,reads_table,
   reads_table_cytokines_PTB,reads_table_cytokines_TB,shannon_rarefaction,threshold,metadata_num,metadata_cha)

####### get lipidomics #########
setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/momspi/Lipidomics/all_lipids')
data <- read.delim('MOMS_PI_lipidomics.20160108.txt')

setwd('/Users/binzhu/Desktop/Mom-baby/')
data$Expression[data$Expression == 'ND'] = 0

sample_name <- unique(data$SampleID)

clinical <- unique(data$Lipid)

reads_table <- matrix(0, nrow = length(sample_name), ncol = length(clinical))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = sample_name 
colnames(reads_table) =clinical 

data$SampleID <- as.character(data$SampleID)
data$Lipid <- as.character(data$Lipid)
data$Expression <- as.numeric(as.character(data$Expression))

for (a in 1: dim(data)[1]) {
  
  SampleID_num <- which(sample_name == data[a,1])
  clinical_num <- which(clinical == data[a,3])
  
  reads_table[SampleID_num,clinical_num] = data$Expression[a]
  
}

reads_table <- as.data.frame(reads_table,stringsAsFactors=FALSE)
row.names(reads_table) = str_replace_all(row.names(reads_table), 'MVAL','MV1D')
reads_table_lipid = as.data.frame(t(reads_table))
