setwd('/Users/binzhu/Desktop/Mom-baby/')

#####  20 sets of samples, combine three neonatal sites #####
# find covariance
{
  # get mom metadata and microbiome alpha & beta diversity
  {
    keep = metadata_16s$Mombaby == 'Baby'
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    colnames(metadata_baby)[colnames(metadata_baby) == 'BMI'] = 'baby_BMI'
    colnames(metadata_baby)[colnames(metadata_baby) == 'pulse'] = 'baby_pulse'
    colnames(metadata_baby)[colnames(metadata_baby) == 'height'] = 'baby_height'
    colnames(metadata_baby)[colnames(metadata_baby) == 'weight'] = 'baby_weight'
    colnames(metadata_baby)[colnames(metadata_baby) == "days_rel2birth"] = "time_after_birth"
    
    Participant_list = metadata_baby$ParticipantID
    
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list,]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]
    metadata_baby = metadata_baby[match(metadata_mom$ParticipantID, metadata_baby$ParticipantID),]
    metadata_baby$BMI = metadata_mom[,which(colnames(metadata_mom) == 'BMI')]
    metadata_baby$pulse = metadata_mom[,which(colnames(metadata_mom) == 'pulse')]
    metadata_baby$height = metadata_mom[,which(colnames(metadata_mom) == 'height')]
    metadata_baby$weight = metadata_mom[,which(colnames(metadata_mom) == 'weight')]
    
    metadata_baby[metadata_baby == 'Yes'] =1
    metadata_baby[metadata_baby == 'No'] =0
    
    metadata_mom[metadata_mom == 'Yes'] =1
    metadata_mom[metadata_mom == 'No'] =0
    
    metadata_mom_MV1D = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom_MV1D = metadata_mom_MV1D[metadata_mom_MV1D$Mombaby == 'Mom',]
    metadata_mom_MV1D = metadata_mom_MV1D[order(metadata_mom_MV1D$VisitNum, decreasing = T),]
    metadata_mom_MV1D = metadata_mom_MV1D[order(metadata_mom_MV1D$ParticipantID),]
    keep = !duplicated(metadata_mom_MV1D$ParticipantID)
    metadata_mom_MV1D = metadata_mom_MV1D[keep,]
    
    metadata_mom_MRCD = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MRCD',]
    metadata_mom_MRCD = metadata_mom_MRCD[metadata_mom_MRCD$Mombaby == 'Mom',]
    metadata_mom_MRCD = metadata_mom_MRCD[order(metadata_mom_MRCD$VisitNum, decreasing = T),]
    metadata_mom_MRCD = metadata_mom_MRCD[order(metadata_mom_MRCD$ParticipantID),]
    keep = !duplicated(metadata_mom_MRCD$ParticipantID)
    metadata_mom_MRCD = metadata_mom_MRCD[keep,]
    
    metadata_mom_MCKD = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MCKD',]
    metadata_mom_MCKD = metadata_mom_MCKD[metadata_mom_MCKD$Mombaby == 'Mom',]
    metadata_mom_MCKD = metadata_mom_MCKD[order(metadata_mom_MCKD$VisitNum, decreasing = T),]
    metadata_mom_MCKD = metadata_mom_MCKD[order(metadata_mom_MCKD$ParticipantID),]
    keep = !duplicated(metadata_mom_MCKD$ParticipantID)
    metadata_mom_MCKD = metadata_mom_MCKD[keep,]
    
    # get alpha & beta diversity
    reads_table_mom_MV1D = reads_table_16s_MV1D[,colnames(reads_table_16s_MV1D) %in% metadata_mom_MV1D$SampleID] 
    reads_table_mom_MV1D = reads_table_mom_MV1D[metadata_mom_MV1D$SampleID]
    reads_table_baby_2 = reads_table_mom_MV1D
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MV1D <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MV1D <- as.matrix(alpha.shannon_diversity_MV1D/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MV1D <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MV1D
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MV1D = pic
    colnames(tsne_MV1D) = c('MV_t_SNE1','MV_t_SNE2')
    
    
    
    reads_table_mom_MRCD = reads_table_16s_MRCD[,colnames(reads_table_16s_MRCD) %in% metadata_mom_MRCD$SampleID] 
    reads_table_mom_MRCD = reads_table_mom_MRCD[metadata_mom_MRCD$SampleID]
    reads_table_baby_2 = reads_table_mom_MRCD
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MRCD <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MRCD <- as.matrix(alpha.shannon_diversity_MRCD/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MRCD <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MRCD
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MRCD = pic
    colnames(tsne_MRCD) = c('MR_t_SNE1','MR_t_SNE2')
    
    
    reads_table_mom_MCKD = reads_table_16s_MCKD[,colnames(reads_table_16s_MCKD) %in% metadata_mom_MCKD$SampleID] 
    reads_table_mom_MCKD = reads_table_mom_MCKD[metadata_mom_MCKD$SampleID]
    reads_table_baby_2 = reads_table_mom_MCKD
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MCKD <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MCKD <- as.matrix(alpha.shannon_diversity_MCKD/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MCKD <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MCKD
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MCKD = pic
    colnames(tsne_MCKD) = c('MB_t_SNE1','MB_t_SNE2')
    
    data = matrix(data = NA, ncol = 15, nrow = nrow(metadata_baby))
    row.names(data) = row.names(metadata_baby)
    metadata_baby = cbind(data,metadata_baby)
    colnames(metadata_baby)[1:15] = c('Shannon_MV','Evenness_MV','Number_of_taxa_MV','MV_t_SNE1','MV_t_SNE2',
                                      'Shannon_MR','Evenness_MR','Number_of_taxa_MR','MR_t_SNE1','MR_t_SNE2',
                                      'Shannon_MB','Evenness_MB','Number_of_taxa_MB','MB_t_SNE1','MB_t_SNE2')
    
    for (a in 1:nrow(metadata_baby)) {
      n = which(metadata_mom_MV1D$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,1] = alpha.shannon_diversity_MV1D[n]
        metadata_baby[a,2] = alpha.evenness_MV1D[n]
        metadata_baby[a,3] = alpha.ovserved_OTU_MV1D[n]
        metadata_baby[a,4] = tsne_MV1D[n,1]
        metadata_baby[a,5] = tsne_MV1D[n,2]
        
      }
      
      n = which(metadata_mom_MRCD$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,6] = alpha.shannon_diversity_MRCD[n]
        metadata_baby[a,7] = alpha.evenness_MRCD[n]
        metadata_baby[a,8] = alpha.ovserved_OTU_MRCD[n]
        metadata_baby[a,9] = tsne_MRCD[n,1]
        metadata_baby[a,10] = tsne_MRCD[n,2]
      }
      
      n = which(metadata_mom_MCKD$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,11] = alpha.shannon_diversity_MCKD[n]
        metadata_baby[a,12] = alpha.evenness_MCKD[n]
        metadata_baby[a,13] = alpha.ovserved_OTU_MCKD[n]
        metadata_baby[a,14] = tsne_MCKD[n,1]
        metadata_baby[a,15] = tsne_MCKD[n,2]
      }
    }
    
    metadata_baby_2 = metadata_baby[,str_detect(colnames(metadata_baby), 'Participant|t_SNE|Shannon|Evenness|Number_of_taxa|
                                                baby_height|baby_weight|baby_pulse|baby_BMI|time_after_birth')]
    metadata_baby_2 = metadata_baby_2[,!str_detect(colnames(metadata_baby_2), 'Mombaby|child|fundal|average|last|change|lost|gain')]
    
  }
  
  # get sample set with matched NB, NR, and NS, totally 20 sets
  {
    {
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
    }
    {
      keep = is.na(output_table$BCKD_birth_0_day) & is.na(output_table$BCKD_birth_1_days) &is.na(output_table$BCKD_birth_2_days) &is.na(output_table$BRCD_birth_0_day) & is.na(output_table$BRCD_birth_1_days) &is.na(output_table$BRCD_birth_2_days) &is.na(output_table$BS1D_birth_0_day) & is.na(output_table$BS1D_birth_1_days) &is.na(output_table$BS1D_birth_2_days)
      sum(keep)
      
      output_table_2 = as.data.frame(matrix(data = NA, ncol = 5, nrow = 0))
      for (a in 1: nrow(output_table)) {
        if (!is.na(output_table$MCKD[a]) & !is.na(output_table$MRCD[a]) & !is.na(output_table$MV1D[a]) & 
            (!is.na(output_table$BCKD_birth_0_day[a]) | !is.na(output_table$BCKD_birth_1_days[a]) | !is.na(output_table$BCKD_birth_2_days[a])) &
            (!is.na(output_table$BRCD_birth_0_day[a]) | !is.na(output_table$BRCD_birth_1_days[a]) | !is.na(output_table$BRCD_birth_2_days[a]))) {
          output_table_3 = as.data.frame(matrix(data = NA, ncol = 5, nrow = 1))
          output_table_3$V1 = output_table$MCKD[a]; output_table_3$V2 = output_table$MRCD[a]; output_table_3$V3 = output_table$MV1D[a]
          
          if (!is.na(output_table$BCKD_birth_2_days[a])) {
            output_table_3$V4 = output_table$BCKD_birth_2_days[a]
          } else if (!is.na(output_table$BCKD_birth_1_days[a])) {
            output_table_3$V4 = output_table$BCKD_birth_1_days[a]
          } else {
            output_table_3$V4 = output_table$BCKD_birth_0_day[a]
          }
          
          if (!is.na(output_table$BRCD_birth_2_days[a])) {
            output_table_3$V5 = output_table$BRCD_birth_2_days[a]
          } else if (!is.na(output_table$BRCD_birth_1_days[a])) {
            output_table_3$V5 = output_table$BRCD_birth_1_days[a]
          } else {
            output_table_3$V5 = output_table$BRCD_birth_0_day[a]
          }
          
          output_table_2 = rbind(output_table_2, output_table_3)
        } 
      }
      
      #setwd('/Users/binzhu/secure/godel/gpfs_fs/home/bzhu/works/MGS_samples')
      #write.csv(output_table_2 , 'Sample_table_paired_mom_baby_135.csv')
      {
        keep = metadata_16s$SampleType == "MV1D"
        metadata = metadata_16s[keep,]
        
        metadata = metadata[order(metadata$VisitNum, decreasing = F),]
        keep = !duplicated(metadata$ParticipantID)
        sum(keep)
        
        metadata = metadata[keep,]
        metadata = metadata[metadata$VisitNum == 1 | metadata$VisitNum == 2,]
        sample_list_mom_1st_sample = metadata$SampleID
        #setwd('/Users/binzhu/secure/godel/gpfs_fs/home/bzhu/works/MGS_samples')
        #write.csv(sample_list_mom_1st_sample , 'sample_list_mom_1st_sample.csv')
      }
    }
    keep = ((!is.na(output_table$BCKD_birth_0_day) | !is.na(output_table$BCKD_birth_1_days)| !is.na(output_table$BCKD_birth_2_days)) & 
              (!is.na(output_table$BRCD_birth_0_day) | !is.na(output_table$BRCD_birth_1_days)| !is.na(output_table$BRCD_birth_2_days)) & 
              (!is.na(output_table$BS1D_birth_1_days)| !is.na(output_table$BS1D_birth_2_days)))
    output_table = output_table[keep,]
    Participant_list = output_table[,1]
    output_table = output_table[,c(5:13)]
    
    for (a in 1:nrow(output_table)) {
      if (!is.na(output_table[a,1])) {
        output_table[a,2] = NA
        output_table[a,3] = NA
      } else if ((!is.na(output_table[a,2]))) {
        output_table[a,3] = NA
      }
      
      if (!is.na(output_table[a,4])) {
        output_table[a,5] = NA
        output_table[a,6] = NA
      } else if ((!is.na(output_table[a,5]))) {
        output_table[a,6] = NA
      }
      
      if (!is.na(output_table[a,7])) {
        output_table[a,8] = NA
        output_table[a,9] = NA
      } else if ((!is.na(output_table[a,8]))) {
        output_table[a,9] = NA
      }
    }
    
    output_table_output = output_table
    row.names(output_table_output) = Participant_list
    write.csv(output_table_output,'output_table_20sets.csv')
    output_table = gather(output_table)
    output_table$ParticipantID = rep(Participant_list,9)
    output_table = output_table[!is.na(output_table$value),]
    
    Sample_list = output_table
    Sample_list$NB = Sample_list$key
    Sample_list$NB[str_detect(Sample_list$NB,'BCKD')] = 1
    Sample_list$NB[str_detect(Sample_list$NB,'B')] = 0
    Sample_list$NR = Sample_list$key
    Sample_list$NR[str_detect(Sample_list$NR,'BRCD')] = 1
    Sample_list$NR[str_detect(Sample_list$NR,'B')] = 0
    Sample_list$NS = Sample_list$key
    Sample_list$NS[str_detect(Sample_list$NS,'BS1D')] = 1
    Sample_list$NS[str_detect(Sample_list$NS,'B')] = 0
    
    Sample_list$key = NULL
    colnames(Sample_list)[1] = 'SampleID'
    
    metadata = merge(Sample_list, metadata_baby_2, 'ParticipantID')
    metadata = merge(metadata, metadata_mom, 'ParticipantID')
    
    Sample_list = Sample_list[order(Sample_list$SampleID),]
    metadata = metadata[order(metadata$SampleID.x),]
  }
  
  # get numeric factors 
  {
    metadata_2 = metadata
    metadata_2[,1:ncol(metadata_2)] <- lapply(metadata_2[,1:ncol(metadata_2)],as.character)
    c2 = c('0','1','2','3','4','5','6','7','8','9','.','-')
    for (a2 in 1:ncol(metadata_2)) {
      b2 = metadata_2[,a2]
      b2 = b2[!is.na(b2)]
      b2 = strsplit(b2,'*')
      b2 = unlist(b2)
      b2 = unique(b2)
      keep = b2 %in% c2
      
      if (sum(!keep) == 0) {
        metadata_2[,a2] <- as.numeric(metadata_2[,a2])
      } else {
        metadata_2[,a2] <- NA
      }
    }
  }
  
  # filter low quality metadata
  {
    metadata_2 = apply(metadata_2, 2, as.numeric)
    
    # remove input with less than 2 levels
    keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[,j])) >= 2 ))
    metadata_2 = metadata_2[,keep]
    
    metadata_2 = metadata_2[, !str_detect(colnames(metadata_2),'average|pregnancy|last_weight|last_pulse|last_ga_at_delivery|last_weight|weight_change') ]
    
    involved_factors = read.csv('involved_factors.csv', header = F)   # in the file, metadata with less than 10% positive cases in any of the matched neonatal microbiomes are removed.
    involved_factors_list = c(involved_factors$V1, c('Shannon_MV','Evenness_MV','Number_of_taxa_MV','MV_t_SNE1','MV_t_SNE2',
                                                     'Shannon_MR','Evenness_MR','Number_of_taxa_MR','MR_t_SNE1','MR_t_SNE2',
                                                     'Shannon_MB','Evenness_MB','Number_of_taxa_MB','MB_t_SNE1','MB_t_SNE2'))
    metadata_2 = metadata_2[,(colnames(metadata_2) %in% involved_factors_list)]
    
    missing_values_metadata_2 = colSums(is.na(metadata_2))
    missing_values_metadata_2 = as.matrix(missing_values_metadata_2)
    write.csv(missing_values_metadata_2,'mom_baby_shannon_missing_values.csv')
    
    keep = colSums(is.na(metadata_2)) <= nrow(metadata_2)/4
    sum(keep)
    metadata_2 = metadata_2[,keep]
    
    keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[,j][!is.na(metadata_2[,j])])) >= 2 ))
    sum(keep)
    metadata_2 = metadata_2[,keep]
    
    
  }
  
  # Imputation by mice
  {
    library(mice)
    data = as.data.frame(metadata_2)
    md.pattern(data, rotate.names = T)
    x.impmi<- mice(data, m = 2, printFlag = FALSE)
    x.impmi_2 = x.impmi$imp
    
    x.impmi_name = names(data)  
    for (b in 1:length(x.impmi_2)) {
      x.impmi_3 = x.impmi_2[[b]]
      colnames(x.impmi_3) = c('A','B')
      x.impmi_3$result = (x.impmi_3$A + x.impmi_3$B)/2
      
      if (length(x.impmi_3$result) == 0) {
        next
      }
      
      n = which(colnames(data) == x.impmi_name[b])
      
      for (c in 1: nrow(x.impmi_3)) {
        m = as.numeric(row.names(x.impmi_3)[c])
        data[m,n]= x.impmi_3$result[c]
      }
    }
    
    metadata_2_imputation = data
    
    for (a in 1: ncol(metadata_2_imputation)) {
      metadata_2_imputation[,a][is.na(metadata_2_imputation[,a])] <-  median(metadata_2_imputation[,a][!is.na(metadata_2_imputation[,a])])
    }
    
    metadata_2_imputation = as.data.frame(scale(metadata_2_imputation))
    
    metadata_2_imputation = as.data.frame(metadata_2_imputation)
    metadata_2_imputation$body_site = NA
    
    for (a in 1: dim(metadata_2_imputation)) {
      if (metadata$NB[a] == '1') {
        metadata_2_imputation$body_site[a] = 1
      } else if (metadata$NR[a] == '1') {
        metadata_2_imputation$body_site[a] = 2
      } else {
        metadata_2_imputation$body_site[a] = 3
      }
    }
  }
}

# beta adonis
{
  # get reads_table
  keep = colnames(reads_table_16s) %in% Sample_list$SampleID
  reads_table = reads_table_16s[,keep]
  reads_table = reads_table[Sample_list$SampleID]
  reads_table = prepare_reads_table_2(reads_table, total_reads_threshold = 5000, species_threshold = 0.0001, mc.cores = 8) 
  
  reads_table = as.data.frame(t(reads_table))
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
  write.csv(reads_table,'mom_baby_beta_reads_table_20set.csv')
  write.csv(metadata_2_imputation,'mom_baby_metadata_2_imputation_20set.csv')
#  write.csv(p_cluster,'mom_baby_p_cluster.csv')
  setwd('/Users/binzhu/Desktop/Mom-baby/')
  ### below should be run on Fenn server ###
  $ conda activate ALDEx2
  #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
  $ R
  
  {
    library(vegan)
    library(parallel)
    reads_table = read.csv('mom_baby_beta_reads_table_20set.csv', row.names = 1)
    metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation_20set.csv', row.names = 1)
    
    pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(metadata_2_imputation)))
    row.names(pvalue_all) = colnames(metadata_2_imputation)
    
    trials = c(1: ncol(metadata_2_imputation))
    func_1 = function(trial) {
      metadata = metadata_2_imputation[,trial]
      metadata = as.data.frame(metadata)
      pvalue <- adonis2(reads_table ~ ., data = metadata, 
                        method = "bray", parallel = 60)
      pvalue <- c(pvalue[1,1],pvalue[1,2],pvalue[1,3],pvalue[1,4],pvalue[1,5])
      return(pvalue)
    }
    pvalue = mclapply(trials, func_1, mc.cores = 60)
    
    for (a in 1:(nrow(pvalue_all))) {
      pvalue_all[a,] = pvalue[[a]]
    }
    colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P-value')
    
    write.csv(pvalue_all,'Adonis_20sets.csv')
    
    pvalue_all = pvalue_all[order(pvalue_all$`P-value`, decreasing = F),]
    metadata_2_imputation = metadata_2_imputation[row.names(pvalue_all)]
    
    metadata = as.data.frame(metadata_2_imputation)
    pvalue <- adonis2(reads_table ~ ., data = metadata,method = "bray", parallel = 60)
    pvalue
    write.csv(pvalue,'multi_variable_analysis_adonis_pvalue_20set2.csv')

  
  }

  
  setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
  data = read.csv('Adonis_20sets.csv')
  data_2 = read.csv('multi_variable_analysis_adonis_pvalue_20set2.csv')
  setwd('/Users/binzhu/Desktop/Mom-baby/')
  write.csv(data,'Adonis_20sets.csv')
  write.csv(data_2,'multi_variable_analysis_adonis_pvalue_20set2.csv')
}

#####  all samples, combine three neonatal sites #####
# find covariance
{
  # get mom metadata and microbiome alpha & beta diversity
  {
    keep = metadata_16s$Mombaby == 'Baby'
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    colnames(metadata_baby)[colnames(metadata_baby) == 'BMI'] = 'baby_BMI'
    colnames(metadata_baby)[colnames(metadata_baby) == 'pulse'] = 'baby_pulse'
    colnames(metadata_baby)[colnames(metadata_baby) == 'height'] = 'baby_height'
    colnames(metadata_baby)[colnames(metadata_baby) == 'weight'] = 'baby_weight'
    colnames(metadata_baby)[colnames(metadata_baby) == "days_rel2birth"] = "time_after_birth"
    
    Participant_list = metadata_baby$ParticipantID
    
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list,]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]
    metadata_baby = metadata_baby[match(metadata_mom$ParticipantID, metadata_baby$ParticipantID),]
    metadata_baby$BMI = metadata_mom[,which(colnames(metadata_mom) == 'BMI')]
    metadata_baby$pulse = metadata_mom[,which(colnames(metadata_mom) == 'pulse')]
    metadata_baby$height = metadata_mom[,which(colnames(metadata_mom) == 'height')]
    metadata_baby$weight = metadata_mom[,which(colnames(metadata_mom) == 'weight')]
    
    metadata_baby[metadata_baby == 'Yes'] =1
    metadata_baby[metadata_baby == 'No'] =0
    
    metadata_mom[metadata_mom == 'Yes'] =1
    metadata_mom[metadata_mom == 'No'] =0
    
    metadata_mom_MV1D = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom_MV1D = metadata_mom_MV1D[metadata_mom_MV1D$Mombaby == 'Mom',]
    metadata_mom_MV1D = metadata_mom_MV1D[order(metadata_mom_MV1D$VisitNum, decreasing = T),]
    metadata_mom_MV1D = metadata_mom_MV1D[order(metadata_mom_MV1D$ParticipantID),]
    keep = !duplicated(metadata_mom_MV1D$ParticipantID)
    metadata_mom_MV1D = metadata_mom_MV1D[keep,]
    
    metadata_mom_MRCD = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MRCD',]
    metadata_mom_MRCD = metadata_mom_MRCD[metadata_mom_MRCD$Mombaby == 'Mom',]
    metadata_mom_MRCD = metadata_mom_MRCD[order(metadata_mom_MRCD$VisitNum, decreasing = T),]
    metadata_mom_MRCD = metadata_mom_MRCD[order(metadata_mom_MRCD$ParticipantID),]
    keep = !duplicated(metadata_mom_MRCD$ParticipantID)
    metadata_mom_MRCD = metadata_mom_MRCD[keep,]
    
    metadata_mom_MCKD = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MCKD',]
    metadata_mom_MCKD = metadata_mom_MCKD[metadata_mom_MCKD$Mombaby == 'Mom',]
    metadata_mom_MCKD = metadata_mom_MCKD[order(metadata_mom_MCKD$VisitNum, decreasing = T),]
    metadata_mom_MCKD = metadata_mom_MCKD[order(metadata_mom_MCKD$ParticipantID),]
    keep = !duplicated(metadata_mom_MCKD$ParticipantID)
    metadata_mom_MCKD = metadata_mom_MCKD[keep,]
    
    # get alpha & beta diversity
    reads_table_mom_MV1D = reads_table_16s_MV1D[,colnames(reads_table_16s_MV1D) %in% metadata_mom_MV1D$SampleID] 
    reads_table_mom_MV1D = reads_table_mom_MV1D[metadata_mom_MV1D$SampleID]
    reads_table_baby_2 = reads_table_mom_MV1D
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MV1D <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MV1D <- as.matrix(alpha.shannon_diversity_MV1D/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MV1D <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MV1D
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MV1D = pic
    colnames(tsne_MV1D) = c('MV_t_SNE1','MV_t_SNE2')
    
    
    
    reads_table_mom_MRCD = reads_table_16s_MRCD[,colnames(reads_table_16s_MRCD) %in% metadata_mom_MRCD$SampleID] 
    reads_table_mom_MRCD = reads_table_mom_MRCD[metadata_mom_MRCD$SampleID]
    reads_table_baby_2 = reads_table_mom_MRCD
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MRCD <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MRCD <- as.matrix(alpha.shannon_diversity_MRCD/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MRCD <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MRCD
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MRCD = pic
    colnames(tsne_MRCD) = c('MR_t_SNE1','MR_t_SNE2')
    
    
    reads_table_mom_MCKD = reads_table_16s_MCKD[,colnames(reads_table_16s_MCKD) %in% metadata_mom_MCKD$SampleID] 
    reads_table_mom_MCKD = reads_table_mom_MCKD[metadata_mom_MCKD$SampleID]
    reads_table_baby_2 = reads_table_mom_MCKD
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MCKD <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MCKD <- as.matrix(alpha.shannon_diversity_MCKD/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MCKD <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MCKD
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MCKD = pic
    colnames(tsne_MCKD) = c('MB_t_SNE1','MB_t_SNE2')
    
    data = matrix(data = NA, ncol = 15, nrow = nrow(metadata_baby))
    row.names(data) = row.names(metadata_baby)
    metadata_baby = cbind(data,metadata_baby)
    colnames(metadata_baby)[1:15] = c('Shannon_MV','Evenness_MV','Number_of_taxa_MV','MV_t_SNE1','MV_t_SNE2',
                                      'Shannon_MR','Evenness_MR','Number_of_taxa_MR','MR_t_SNE1','MR_t_SNE2',
                                      'Shannon_MB','Evenness_MB','Number_of_taxa_MB','MB_t_SNE1','MB_t_SNE2')
    
    for (a in 1:nrow(metadata_baby)) {
      n = which(metadata_mom_MV1D$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,1] = alpha.shannon_diversity_MV1D[n]
        metadata_baby[a,2] = alpha.evenness_MV1D[n]
        metadata_baby[a,3] = alpha.ovserved_OTU_MV1D[n]
        metadata_baby[a,4] = tsne_MV1D[n,1]
        metadata_baby[a,5] = tsne_MV1D[n,2]
        
      }
      
      n = which(metadata_mom_MRCD$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,6] = alpha.shannon_diversity_MRCD[n]
        metadata_baby[a,7] = alpha.evenness_MRCD[n]
        metadata_baby[a,8] = alpha.ovserved_OTU_MRCD[n]
        metadata_baby[a,9] = tsne_MRCD[n,1]
        metadata_baby[a,10] = tsne_MRCD[n,2]
      }
      
      n = which(metadata_mom_MCKD$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,11] = alpha.shannon_diversity_MCKD[n]
        metadata_baby[a,12] = alpha.evenness_MCKD[n]
        metadata_baby[a,13] = alpha.ovserved_OTU_MCKD[n]
        metadata_baby[a,14] = tsne_MCKD[n,1]
        metadata_baby[a,15] = tsne_MCKD[n,2]
      }
    }
    
    metadata_baby_2 = metadata_baby[,str_detect(colnames(metadata_baby), 'Participant|t_SNE|Shannon|Evenness|Number_of_taxa|
                                                baby_height|baby_weight|baby_pulse|baby_BMI|time_after_birth')]
    metadata_baby_2 = metadata_baby_2[,!str_detect(colnames(metadata_baby_2), 'Mombaby|child|fundal|average|last|change|lost|gain')]
    
  }

  metadata = merge(metadata_baby_2, metadata_mom, 'ParticipantID')
  
  # get numeric factors 
  {
    metadata_2 = metadata
    metadata_2[,1:ncol(metadata_2)] <- lapply(metadata_2[,1:ncol(metadata_2)],as.character)
    c2 = c('0','1','2','3','4','5','6','7','8','9','.','-')
    for (a2 in 1:ncol(metadata_2)) {
      b2 = metadata_2[,a2]
      b2 = b2[!is.na(b2)]
      b2 = strsplit(b2,'*')
      b2 = unlist(b2)
      b2 = unique(b2)
      keep = b2 %in% c2
      
      if (sum(!keep) == 0) {
        metadata_2[,a2] <- as.numeric(metadata_2[,a2])
      } else {
        metadata_2[,a2] <- NA
      }
    }
  }
  
  # filter low quality metadata
  {
    metadata_2 = apply(metadata_2, 2, as.numeric)
    
    # remove input with less than 2 levels
    keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[,j])) >= 2 ))
    metadata_2 = metadata_2[,keep]
    
    metadata_2 = metadata_2[, !str_detect(colnames(metadata_2),'average|pregnancy|last_weight|last_pulse|last_ga_at_delivery|last_weight|weight_change') ]
    
    involved_factors = read.csv('involved_factors.csv', header = F)   # in the file, metadata with less than 10% positive cases in any of the matched neonatal microbiomes are removed.
    involved_factors_list = c(involved_factors$V1, c('Shannon_MV','Evenness_MV','Number_of_taxa_MV','MV_t_SNE1','MV_t_SNE2',
                                                     'Shannon_MR','Evenness_MR','Number_of_taxa_MR','MR_t_SNE1','MR_t_SNE2',
                                                     'Shannon_MB','Evenness_MB','Number_of_taxa_MB','MB_t_SNE1','MB_t_SNE2'))
    metadata_2 = metadata_2[,(colnames(metadata_2) %in% involved_factors_list)]
    
    missing_values_metadata_2 = colSums(is.na(metadata_2))
    missing_values_metadata_2 = as.matrix(missing_values_metadata_2)
    write.csv(missing_values_metadata_2,'mom_baby_shannon_missing_values.csv')
    
    keep = colSums(is.na(metadata_2)) <= nrow(metadata_2)/4
    sum(keep)
    metadata_2 = metadata_2[,keep]
    
    keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[,j][!is.na(metadata_2[,j])])) >= 2 ))
    sum(keep)
    metadata_2 = metadata_2[,keep]
    
    
  }

  # Imputation by mice
  {
    library(mice)
    data = as.data.frame(metadata_2)
    md.pattern(data, rotate.names = T)
    x.impmi<- mice(data, m = 2, printFlag = FALSE)
    x.impmi_2 = x.impmi$imp
    
    x.impmi_name = names(data)  

    for (b in 1:length(x.impmi_2)) {
      x.impmi_3 = x.impmi_2[[b]]
      colnames(x.impmi_3) = c('A','B')
      x.impmi_3$result = (x.impmi_3$A + x.impmi_3$B)/2
      
      if (length(x.impmi_3$result) == 0) {
        next
      }
      
      n = which(colnames(data) == x.impmi_name[b])
      
      for (c in 1: nrow(x.impmi_3)) {
        m = as.numeric(row.names(x.impmi_3)[c])
        data[m,n]= x.impmi_3$result[c]
      }
      
    }
 
      metadata_2_imputation = data
      
      metadata_2_imputation$body_site = NA
      
      for (a in 1: dim(metadata_2_imputation)) {
        if (metadata_baby$SampleType[a] == 'BCKD') {
          metadata_2_imputation$body_site[a] = 1
        } else if (metadata_baby$SampleType[a] == 'BRCD') {
          metadata_2_imputation$body_site[a] = 2
        } else {
          metadata_2_imputation$body_site[a] = 3
        }
      }
    
  }
}

# beta adonis
{
  # get reads_table
  keep = metadata_16s$Mombaby == 'Baby'
  metadata = metadata_16s[keep,]
  reads_table = reads_table_16s[,keep]
  
  metadata_2= as.data.frame(matrix(data = NA, nrow = nrow(metadata), ncol = ncol(metadata_2_imputation)))
  colnames(metadata_2) = colnames(metadata_2_imputation)
  for (a in 1:nrow(metadata_2)) {
    n= which(metadata_baby$ParticipantID == metadata$ParticipantID[a])
    metadata_2[a,] = metadata_2_imputation[n,]
    
  }
  
  reads_table = prepare_reads_table_2(reads_table, total_reads_threshold = 0, species_threshold = 0.0001, mc.cores = 8) 
  
  reads_table = as.data.frame(t(reads_table))
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
  write.csv(reads_table,'mom_baby_beta_reads_table.csv')
  write.csv(metadata_2,'mom_baby_metadata_2_imputation.csv')
  setwd('/Users/binzhu/Desktop/Mom-baby/')
  ### below should be run on Fenn server ###
  $ conda activate ALDEx2
  #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
  $ R
  
  {
    library(vegan)
    library(parallel)
    reads_table = read.csv('mom_baby_beta_reads_table.csv', row.names = 1)
    metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation.csv', row.names = 1)
    
    pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(metadata_2_imputation)))
    row.names(pvalue_all) = colnames(metadata_2_imputation)
    
    for (a in 1:(nrow(pvalue_all))) {
      metadata = metadata_2_imputation[,a]
      metadata = as.data.frame(metadata)
      pvalue <- adonis2(reads_table ~ ., data = metadata, 
                        method = "bray", parallel = 60)
      pvalue_all[a,] = pvalue[1,]
    }
    colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
    write.csv(pvalue_all,'Adonis_all.csv')
    
#    pvalue_all$V2[row.names(pvalue_all) == 'diabetes'] =1
#    pvalue_all$V2[row.names(pvalue_all) == 'live_babies'] =1
#    pvalue_all$V2[row.names(pvalue_all) == 'Evenness_MB'] =1
#    pvalue_all$V2[row.names(pvalue_all) == 'diabetes'] =1

    pvalue_all = pvalue_all[order(pvalue_all$P_value, decreasing = F),]
    metadata_2_imputation = metadata_2_imputation[row.names(pvalue_all)]
    
    metadata = as.data.frame(metadata_2_imputation)
    pvalue <- adonis2(reads_table ~ ., data = metadata,method = "bray", parallel = 60)
    pvalue
    write.csv(pvalue,'multi_variable_analysis_adonis_pvalue_all.csv')
    
    
  }
  
  
  setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
  data = read.csv('Adonis_all.csv')
  data_2 = read.csv('multi_variable_analysis_adonis_pvalue_all.csv')
  setwd('/Users/binzhu/Desktop/Mom-baby/')
  write.csv(data,'Adonis_all.csv')
  write.csv(data_2,'multi_variable_analysis_adonis_pvalue_all.csv')
}
#####  separate three neonatal sites, alpha and beta multivariate analyses #####
# find covariance
{
  # get mom metadata and microbiome alpha & beta diversity
  {
    keep = metadata_16s$Mombaby == 'Baby'
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    colnames(metadata_baby)[colnames(metadata_baby) == 'BMI'] = 'baby_BMI'
    colnames(metadata_baby)[colnames(metadata_baby) == 'pulse'] = 'baby_pulse'
    colnames(metadata_baby)[colnames(metadata_baby) == 'height'] = 'baby_height'
    colnames(metadata_baby)[colnames(metadata_baby) == 'weight'] = 'baby_weight'
    colnames(metadata_baby)[colnames(metadata_baby) == "days_rel2birth"] = "time_after_birth"
    
    Participant_list = metadata_baby$ParticipantID
    
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list,]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]
    metadata_baby = metadata_baby[match(metadata_mom$ParticipantID, metadata_baby$ParticipantID),]
    metadata_baby$BMI = metadata_mom[,which(colnames(metadata_mom) == 'BMI')]
    metadata_baby$pulse = metadata_mom[,which(colnames(metadata_mom) == 'pulse')]
    metadata_baby$height = metadata_mom[,which(colnames(metadata_mom) == 'height')]
    metadata_baby$weight = metadata_mom[,which(colnames(metadata_mom) == 'weight')]
    
    metadata_baby[metadata_baby == 'Yes'] =1
    metadata_baby[metadata_baby == 'No'] =0
    
    metadata_mom[metadata_mom == 'Yes'] =1
    metadata_mom[metadata_mom == 'No'] =0
    
    metadata_mom_MV1D = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom_MV1D = metadata_mom_MV1D[metadata_mom_MV1D$Mombaby == 'Mom',]
    metadata_mom_MV1D = metadata_mom_MV1D[order(metadata_mom_MV1D$VisitNum, decreasing = T),]
    metadata_mom_MV1D = metadata_mom_MV1D[order(metadata_mom_MV1D$ParticipantID),]
    keep = !duplicated(metadata_mom_MV1D$ParticipantID)
    metadata_mom_MV1D = metadata_mom_MV1D[keep,]
    
    metadata_mom_MRCD = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MRCD',]
    metadata_mom_MRCD = metadata_mom_MRCD[metadata_mom_MRCD$Mombaby == 'Mom',]
    metadata_mom_MRCD = metadata_mom_MRCD[order(metadata_mom_MRCD$VisitNum, decreasing = T),]
    metadata_mom_MRCD = metadata_mom_MRCD[order(metadata_mom_MRCD$ParticipantID),]
    keep = !duplicated(metadata_mom_MRCD$ParticipantID)
    metadata_mom_MRCD = metadata_mom_MRCD[keep,]
    
    metadata_mom_MCKD = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MCKD',]
    metadata_mom_MCKD = metadata_mom_MCKD[metadata_mom_MCKD$Mombaby == 'Mom',]
    metadata_mom_MCKD = metadata_mom_MCKD[order(metadata_mom_MCKD$VisitNum, decreasing = T),]
    metadata_mom_MCKD = metadata_mom_MCKD[order(metadata_mom_MCKD$ParticipantID),]
    keep = !duplicated(metadata_mom_MCKD$ParticipantID)
    metadata_mom_MCKD = metadata_mom_MCKD[keep,]
    
    # get alpha & beta diversity
    reads_table_mom_MV1D = reads_table_16s_MV1D[,colnames(reads_table_16s_MV1D) %in% metadata_mom_MV1D$SampleID] 
    reads_table_mom_MV1D = reads_table_mom_MV1D[metadata_mom_MV1D$SampleID]
    reads_table_baby_2 = reads_table_mom_MV1D
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MV1D <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MV1D <- as.matrix(alpha.shannon_diversity_MV1D/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MV1D <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MV1D
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MV1D = pic
    colnames(tsne_MV1D) = c('MV_t_SNE1','MV_t_SNE2')
    
    
    
    reads_table_mom_MRCD = reads_table_16s_MRCD[,colnames(reads_table_16s_MRCD) %in% metadata_mom_MRCD$SampleID] 
    reads_table_mom_MRCD = reads_table_mom_MRCD[metadata_mom_MRCD$SampleID]
    reads_table_baby_2 = reads_table_mom_MRCD
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MRCD <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MRCD <- as.matrix(alpha.shannon_diversity_MRCD/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MRCD <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MRCD
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MRCD = pic
    colnames(tsne_MRCD) = c('MR_t_SNE1','MR_t_SNE2')
    
    
    reads_table_mom_MCKD = reads_table_16s_MCKD[,colnames(reads_table_16s_MCKD) %in% metadata_mom_MCKD$SampleID] 
    reads_table_mom_MCKD = reads_table_mom_MCKD[metadata_mom_MCKD$SampleID]
    reads_table_baby_2 = reads_table_mom_MCKD
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity_MCKD <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness_MCKD <- as.matrix(alpha.shannon_diversity_MCKD/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU_MCKD <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
    
    reads_table = reads_table_mom_MCKD
    reads_table = as.data.frame(t(reads_table))
    reads_table = reads_table + 0.5
    reads_table <- as.data.frame(clr(reads_table))
    tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
    pic <- tsne$Y
    tsne_MCKD = pic
    colnames(tsne_MCKD) = c('MB_t_SNE1','MB_t_SNE2')
    
    data = matrix(data = NA, ncol = 15, nrow = nrow(metadata_baby))
    row.names(data) = row.names(metadata_baby)
    metadata_baby = cbind(data,metadata_baby)
    colnames(metadata_baby)[1:15] = c('Shannon_MV','Evenness_MV','Number_of_taxa_MV','MV_t_SNE1','MV_t_SNE2',
                                      'Shannon_MR','Evenness_MR','Number_of_taxa_MR','MR_t_SNE1','MR_t_SNE2',
                                      'Shannon_MB','Evenness_MB','Number_of_taxa_MB','MB_t_SNE1','MB_t_SNE2')
    
    for (a in 1:nrow(metadata_baby)) {
      n = which(metadata_mom_MV1D$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,1] = alpha.shannon_diversity_MV1D[n]
        metadata_baby[a,2] = alpha.evenness_MV1D[n]
        metadata_baby[a,3] = alpha.ovserved_OTU_MV1D[n]
        metadata_baby[a,4] = tsne_MV1D[n,1]
        metadata_baby[a,5] = tsne_MV1D[n,2]
        
      }
      
      n = which(metadata_mom_MRCD$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,6] = alpha.shannon_diversity_MRCD[n]
        metadata_baby[a,7] = alpha.evenness_MRCD[n]
        metadata_baby[a,8] = alpha.ovserved_OTU_MRCD[n]
        metadata_baby[a,9] = tsne_MRCD[n,1]
        metadata_baby[a,10] = tsne_MRCD[n,2]
      }
      
      n = which(metadata_mom_MCKD$ParticipantID == metadata_baby$ParticipantID[a])
      if (length(n) != 0) {
        metadata_baby[a,11] = alpha.shannon_diversity_MCKD[n]
        metadata_baby[a,12] = alpha.evenness_MCKD[n]
        metadata_baby[a,13] = alpha.ovserved_OTU_MCKD[n]
        metadata_baby[a,14] = tsne_MCKD[n,1]
        metadata_baby[a,15] = tsne_MCKD[n,2]
      }
    }
    
    metadata_baby_2 = metadata_baby[,str_detect(colnames(metadata_baby), 'Participant|t_SNE|Shannon|Evenness|Number_of_taxa|baby_height|baby_weight|baby_pulse|baby_BMI|time_after_birth')]
    metadata_baby_2 = metadata_baby_2[,!str_detect(colnames(metadata_baby_2), 'Mombaby|child|fundal|average|last|change|lost|gain')]
    
  }

  # get numeric factors 
  {
    metadata_2 = merge(metadata_baby_2, metadata_mom,'ParticipantID')
    metadata_2[,1:ncol(metadata_2)] <- lapply(metadata_2[,1:ncol(metadata_2)],as.character)
    c2 = c('0','1','2','3','4','5','6','7','8','9','.','-')
    for (a2 in 1:ncol(metadata_2)) {
      b2 = metadata_2[,a2]
      b2 = b2[!is.na(b2)]
      b2 = strsplit(b2,'*')
      b2 = unlist(b2)
      b2 = unique(b2)
      keep = b2 %in% c2
      
      if (sum(!keep) == 0) {
        metadata_2[,a2] <- as.numeric(metadata_2[,a2])
      } else {
        metadata_2[,a2] <- NA
      }
    }
  }
  
  # filter low quality metadata
  {
    metadata_2 = apply(metadata_2, 2, as.numeric)
    
    # remove input with less than 2 levels
    keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[,j])) >= 2 ))
    metadata_2 = metadata_2[,keep]
    
    metadata_2 = metadata_2[, !str_detect(colnames(metadata_2),'average|pregnancy|last_weight|last_pulse|last_ga_at_delivery|last_weight|weight_change') ]
    
    involved_factors = read.csv('involved_factors.csv', header = F)
    involved_factors_list = c(involved_factors$V1, c('Shannon_MV','Evenness_MV','Number_of_taxa_MV','MV_t_SNE1','MV_t_SNE2',
                                                     'Shannon_MR','Evenness_MR','Number_of_taxa_MR','MR_t_SNE1','MR_t_SNE2',
                                                     'Shannon_MB','Evenness_MB','Number_of_taxa_MB','MB_t_SNE1','MB_t_SNE2'))
    metadata_2 = metadata_2[,(colnames(metadata_2) %in% involved_factors_list)]
    
    missing_values_metadata_2 = colSums(is.na(metadata_2))
    missing_values_metadata_2 = as.matrix(missing_values_metadata_2)
    write.csv(missing_values_metadata_2,'mom_baby_shannon_missing_values.csv')
    
    keep = colSums(is.na(metadata_2)) <= nrow(metadata_2)/4
    sum(keep)
    metadata_2 = metadata_2[,keep]
    
    keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[,j][!is.na(metadata_2[,j])])) >= 2 ))
    sum(keep)
    metadata_2 = metadata_2[,keep]
  }
  
  # spearman's correlation
  {
    p.yes.r = matrix(data = NA, ncol = ncol(metadata_2), nrow = ncol(metadata_2))
    colnames(p.yes.r) = colnames(metadata_2) 
    row.names(p.yes.r) = colnames(metadata_2) 
    
    otu.pval = matrix(data = NA, ncol = ncol(metadata_2), nrow = ncol(metadata_2))
    colnames(otu.pval) = colnames(metadata_2) 
    row.names(otu.pval) = colnames(metadata_2) 
    
    for (a in 1: (ncol(p.yes.r) -1)) {
      for (b in a : ncol(p.yes.r)) {
        
        if (a == b) {
          otu.pval[a,b] = 0
          p.yes.r[a,b] = 1
          next
        }
        data = cbind(metadata_2[,a],metadata_2[,b])
        data = data[!is.na(data[,1]) & !is.na(data[,2]),]
        pvalue = cor.test(data[,1], data[,2], method = "spearman")
        otu.pval[a,b] = pvalue$p.value
        otu.pval[b,a] = pvalue$p.value
        p.yes.r[a,b] = as.numeric(pvalue$estimate)
        p.yes.r[b,a] = as.numeric(pvalue$estimate)
        
      }
    }
    p.yes.r[is.na(p.yes.r)] = 0
    otu.pval[is.na(otu.pval)] = 1
    write.csv(p.yes.r,'metadata_spearman_adj_Rvalue.csv')

    otu.pval_adj <- adjust.p(as.matrix(otu.pval), pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
    otu.pval_adj <- otu.pval_adj$adjp
    otu.pval_adj <- otu.pval_adj$adjusted.p
    otu.pval_adj = matrix(otu.pval_adj,nrow = ncol(p.yes.r),ncol = ncol(p.yes.r))
    
    p.yes <- otu.pval_adj< 0.05  
    p.yes.r <- p.yes.r*p.yes # only select correlation values based on p-value criterion
    
    write.csv(otu.pval_adj,'metadata_spearman_adj_pvalue.csv')
  }
  
  # clustering
  {
    clustering_method = "ward.D"; clustering_distance = "canberra"; k = 20
    
    p.yes.rr = abs(p.yes.r)
#    p.yes.rr = (p.yes.r)
#    n = which(colnames(p.yes.rr) == 'new_medications')
#    m = which(colnames(p.yes.rr) == 'ovarian_cyst')
#    p.yes.rr[n,m] = 0; p.yes.rr[m,n] = 0
    
    p <- pheatmap(p.yes.rr, clustering_method = clustering_method ,clustering_distance_rows = clustering_distance,
                  clustering_distance_cols = clustering_distance)
    
    p_cluster = sort(cutree(p$tree_col, k=k))
    
    #    p_cluster[p_cluster == 7] = 2
    
    p_cluster = as.data.frame(p_cluster)
    p_cluster$p_cluster = as.factor(p_cluster$p_cluster)
    
    data = p$tree_row$labels[p$tree_row$order]
    p_cluster$x = row.names(p_cluster)
    p_cluster = p_cluster[match(data,p_cluster$x),]
    p_cluster$p_cluster = as.character(p_cluster$p_cluster)
    data = unique(p_cluster$p_cluster)
    for (a in 1: nrow(p_cluster)) {
      p_cluster$p_cluster[a] = which(data == p_cluster$p_cluster[a])
    }
    p_cluster$x = NULL
    write.csv(p_cluster,'mom_factor_cluster.csv')
    data = p$tree_row$labels[p$tree_row$order]
    data = as.data.frame(data)
    data = cbind(data,p_cluster$p_cluster)
    data$order = row.names(data)
    write.csv(data,'mom_factor_cluster_order.csv')
    
    bar_max = max(p.yes.rr,na.rm = T)
    bar_min = min(p.yes.rr,na.rm = T)
    paletteLength <- 50
    myColor <- colorRampPalette(c("white", "red"))(paletteLength)
    myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                  seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
    
    myBreaks <- unique(myBreaks)
    
    #    mycolors = list(p_cluster = c("1"="#1F78B4","2"="#33A02C","3"="#E31A1C",
    #                                  "4"="#FF7F00","5"="#6A3D9A","6"="#B15928"))
    
    p_cluster$p_cluster = as.character(p_cluster$p_cluster)
    
    pdf("Covariance_metatadat.pdf", height = 12,width = 12)
    p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = 5, 
                  treeheight_row = 50, treeheight_col = 50, annotation_col = p_cluster,annotation_row = p_cluster, 
                  clustering_method = clustering_method ,clustering_distance_rows = clustering_distance,
                  clustering_distance_cols = clustering_distance )
    dev.off()
    
  }
  
  # Imputation by mice
  {
    library(mice)
    data = as.data.frame(metadata_2)
    md.pattern(data, rotate.names = T)
    x.impmi<- mice(data, m = 2, printFlag = FALSE)
    x.impmi_2 = x.impmi$imp
    
    x.impmi_name = names(data)  
    for (b in 1:length(x.impmi_2)) {
      x.impmi_3 = x.impmi_2[[b]]
      colnames(x.impmi_3) = c('A','B')
      x.impmi_3$result = (x.impmi_3$A + x.impmi_3$B)/2
      
      if (length(x.impmi_3$result) == 0) {
        next
      }
      
      n = which(colnames(data) == x.impmi_name[b])
      
      for (c in 1: nrow(x.impmi_3)) {
        m = as.numeric(row.names(x.impmi_3)[c])
        data[m,n]= x.impmi_3$result[c]
      }
      

      
    }
    
    sum(is.na(data))
    metadata_2_imputation = data
    data = as.data.frame(data)
    
    metadata_2_imputation = as.data.frame(scale(metadata_2_imputation))
  }
}

colnames(metadata_16s)[colnames(metadata_16s) == 'mom_delivery_method (C-section yes / vaginal no)'] = 'mom_delivery_method'

type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]
library(pROC)
library(doSNOW)
library(parallel)
library(matrixStats)
library(erer) # export list
library(caret)

for (a in 1: ncol(metadata_2_imputation)) {
    keep = metadata_2_imputation[,a] > mean(metadata_2_imputation[,a])
    metadata_2_imputation[keep,a] = 'High'
    metadata_2_imputation[!keep,a] = 'Low'

}

# NB
{
  a=1
  # find paired mom factors and baby's microbiome (time closest between mom and baby)
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]

    input = cbind(metadata_mom$ParticipantID, metadata_2_imputation)
    colnames(input)[1] = 'ParticipantID'
    keep = input$ParticipantID %in% metadata_baby$ParticipantID
    sum(keep)
    input = input[keep,]
    
    input = input[order(input$ParticipantID),]
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    input = input[keep,]
    
    reads_table_baby_2 = reads_table_baby
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    
    output = alpha.shannon_diversity
    input = as.matrix(input[,-1])
  }
  # adonis
  {
    {
      # get reads_table
      reads_table = as.data.frame(t(reads_table_baby))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table <- reads_table$otu.tab.rff
      reads_table <- as.data.frame(reads_table)
      
      setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
      write.csv(reads_table,'mom_baby_beta_reads_table.csv')
      write.csv(input,'mom_baby_metadata_2_imputation.csv')

      setwd('/Users/binzhu/Desktop/Mom-baby/')
      ### below should be run on Fenn server ###
      $ conda activate ALDEx2
      #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
      $ R
      
      {
        library(vegan)
        library(parallel)
        setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
        
        reads_table = read.csv('mom_baby_beta_reads_table_NB.csv', row.names = 1)
        metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation_NB.csv', row.names = 1)
        
        pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(metadata_2_imputation)))
        row.names(pvalue_all) = colnames(metadata_2_imputation)
        
        trials = c(1: ncol(metadata_2_imputation))
        func_1 = function(trial) {
          metadata = metadata_2_imputation[,trial]
          metadata = as.data.frame(metadata)
          pvalue <- adonis2(reads_table ~ ., data = metadata, 
                            method = "bray", parallel = 60)
          pvalue <- pvalue[1,]
          return(pvalue)
        }
        pvalue = mclapply(trials, func_1, mc.cores = 60)
        
        for (a in 1:(nrow(pvalue_all))) {
          pvalue_all[a,] = pvalue[[a]]
        }
        colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
        write.csv(pvalue_all,'Adonis_NB.csv')

        pvalue_all = pvalue_all[order(pvalue_all$P_value, decreasing = F),]
        metadata_2_imputation = metadata_2_imputation[row.names(pvalue_all)]
        
        metadata = as.data.frame(metadata_2_imputation[,c(1:110)])
        pvalue <- adonis2(reads_table ~ ., data = metadata,method = "bray", parallel = 60)
        pvalue
        write.csv(pvalue,'multi_variable_analysis_adonis_pvalue_NB.csv')
        
      }
      
      
      setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
      data = read.csv('Adonis_NB.csv')
      data_2 = read.csv('multi_variable_analysis_adonis_pvalue_NB.csv')
      setwd('/Users/binzhu/Desktop/Mom-baby/')
      data4 <- adjust.p(data$P_value, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data4 <- data4$adjp
      data4 = data4$adjusted.p
      data$adj_p = data4
      write.csv(data,'Adonis_NB.csv')
      write.csv(data_2,'multi_variable_analysis_adonis_pvalue_NB.csv')
    }
  }
  
  # adonis power
  {
    # get reads_table
    reads_table = as.data.frame(t(reads_table_baby))
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
    reads_table <- reads_table$otu.tab.rff
    reads_table <- as.data.frame(reads_table)
    
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    write.csv(reads_table,'mom_baby_beta_reads_table_NB.csv')
    write.csv(input,'mom_baby_metadata_2_imputation_NB.csv')
    
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    
    ### below should be run on Fenn server ###
    $ conda activate ALDEx2
    #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
    $ R
    setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    reads_table = read.csv('mom_baby_beta_reads_table_NB.csv', row.names = 1)
    metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation_NB.csv', row.names = 1)
    
    setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa/sample_size_adonis')
    {
      for (trial in 1: ncol(metadata_2_imputation)) {
        file_name = paste0('NB_sample_size_',colnames(metadata_2_imputation)[trial],'.R')
        file.create(file_name)
        
        write("setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')", file=file_name,append=TRUE)
        write("reads_table = read.csv('mom_baby_beta_reads_table_NB.csv', row.names = 1)", file=file_name,append=TRUE)
        write("metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation_NB.csv', row.names = 1)", file=file_name,append=TRUE)
        write("adonis_power <- readRDS('adonis_power_function.RData')", file=file_name,append=TRUE)
        write(paste0('trial = ',trial), file=file_name, append=TRUE)
        write("metadata = metadata_2_imputation[,trial]", file=file_name,append=TRUE)
        write("power_size_model_lm = adonis_power(reads_table, metadata)", file=file_name,append=TRUE)
        write("setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa/sample_size_adonis')", file=file_name,append=TRUE)
        write("saveRDS(power_size_model_lm, file=paste0('NB_power_size_model_lm_',colnames(metadata_2_imputation)[trial],'.RData'))", file=file_name,append=TRUE)
      }
    }
    
  }
  
  # spearman's correlation between neonatal taxa and maternal factors
  {
    reads_table_baby_2 = reads_table_baby
    reads_table_baby_2 = as.data.frame(t(reads_table_baby_2))
    reads_table_baby_2 = reads_table_baby_2 + 0.5
    reads_table_baby_2 <- as.matrix(clr(reads_table_baby_2))     ### CLR normalization in rows
    reads_table_baby_2 = as.data.frame(t(reads_table_baby_2))
    
    {
      type = 'spearman'; pvalue = 0.05; cor_parameter= 0; style = 1; bar_max = 2; 
      bar_min = -2; pheatmap_fontsize = 5; treeheight = 50; alpha = 0.05; mc.cores =8; pvalue_adj = T
      reads_table1 = as.matrix(as.data.frame(t(input)))
      reads_table2 = as.matrix(reads_table_baby_2)
      
      correlation_data_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = 0))
      colnames(correlation_data_all) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
      
      try(for (a in 1: nrow(reads_table1)) {
        correlation_data = as.data.frame(matrix(data = NA, ncol = 5, nrow = nrow(reads_table2)))
        colnames(correlation_data) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
        correlation_data$Gene = row.names(reads_table1)[a]
        correlation_data$Taxa = row.names(reads_table2)
        
        trials = c(1: nrow(reads_table2))
        
        func_1 = function(trial) {
          data1 = as.numeric(reads_table1[a,])
          data2 = as.numeric(reads_table2[trial,])
          
          keep = (!is.na(data1)) & (!is.na(data2))
          data1 = data1[keep]
          data2 = data2[keep]
          
          if (length(data1) <8) {
            return(c(list(pvalue = NA), list(Rvalue = NA)))
          }
          
          data3 = cor.test(data1, data2, method = type)
          pvalue = data3$p.value
          Rvalue = data3$estimate
          
          return(c(list(pvalue = pvalue), list(Rvalue = Rvalue)))
        }
        
        pRvalue = mclapply(trials, func_1, mc.cores = mc.cores)
        pRvalue = (unlist(pRvalue))
        keep = seq(1, length(pRvalue), 2)
        correlation_data$Pvalue = pRvalue[keep]
        
        keep = seq(2, length(pRvalue), 2)
        correlation_data$Rvalue = pRvalue[keep]
        
        if (a ==1) {
          correlation_data_all = correlation_data
        } else {
          correlation_data_all = rbind(correlation_data_all,correlation_data)
        }
      })
      
      data4 <- adjust.p(correlation_data_all$Pvalue, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data4 <- data4$adjp
      correlation_data_all$adj_p <- data4$adjusted.p
      correlation_data_all = correlation_data_all[order(correlation_data_all$Pvalue,decreasing = F),]
      write.csv(correlation_data_all,'Correlation_neonatal_taxa_factors_NB.csv', row.names = F)
    }
  }
  
  # correlation with maternal microbiome
  {
    data = data.frame(Shannon_index_NB = output, Shannon_index_M = input[,1])
    linearMod <- lm(Shannon_index_NB ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NB)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MV microbiome', 
           y = 'Shanon index of the NR microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NB_MV.pdf',width=3, height=3)
    
    data = data.frame(Shannon_index_NB = output, Shannon_index_M = input[,4])
    linearMod <- lm(Shannon_index_NB ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NB)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MR microbiome', 
           y = 'Shanon index of the NR microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NB_MR.pdf',width=3, height=3)
    
    data = data.frame(Shannon_index_NB = output, Shannon_index_M = input[,7])
    linearMod <- lm(Shannon_index_NB ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NB)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MB microbiome', 
           y = 'Shanon index of the NR microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NB_MB.pdf',width=3, height=3)
  }
  
  # random forest model
  {
    # cross validation
    ### 
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    write.csv(input,'mom_baby_shannon_NB.filted.csv')
    write.csv(output,'mom_baby_shannon_NB.output.csv')
    write.csv(p_cluster,'mom_baby_shannon_NB.p_cluster.csv')
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    ### below should be run on Fenn server ###
    $ conda activate ALDEx2
    #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
    $ R
    {
      library(pROC)
      library(doSNOW)
      library(parallel)
      library(matrixStats)
      library(caret)
      
      input = read.csv('mom_baby_shannon_NB.filted.csv', row.names = 1)
      input = input[,!str_detect(colnames(input),'baby')]
      output = read.csv('mom_baby_shannon_NB.output.csv', row.names = 1)
      output = output$x
      p_cluster = read.csv('mom_baby_shannon_NB.p_cluster.csv', row.names = 1)
      
      pvalue_RF = vector()
      r_squared_RF = vector()
      
      p_cluster$p_cluster = (as.character(p_cluster$p_cluster))
      p_cluster_list = unique(p_cluster$p_cluster)

      # RF
      {
        trials = c(1: nrow(input))
        func_1 = function(trial) {
          input_1 = input[-trial,]
          output_1 = output[-trial]
          
          input_2 = as.data.frame(rbind(input[trial,],input[trial,]))
          
          data = cbind(output_1,input_1)
          set.seed(2022)
          mtry_number = as.integer(ncol(input_3)/3)
          rfregFit <- train(output_1 ~ ., data = data, method = "ranger",importance="permutation",
                            tuneGrid = data.frame(mtry=mtry_number,
                                                  min.node.size = ncol(input) *2,
                                                  splitrule="variance"))
          
          predicted_value_RF_x = predict(rfregFit, input_2, na.action = na.pass)[1]
          
          return(predicted_value_RF_x)
        }
        predicted_value_RF_x = mclapply(trials, func_1, mc.cores = 50)
        predicted_value_RF_x = unlist(predicted_value_RF_x)
        
        data = data.frame(Shannon_index = output, Predicted_value = predicted_value_RF_x)
        linearMod <- lm(Predicted_value ~ Shannon_index, data=data)
        linearMod = summary(linearMod)
        r_squared <- linearMod$r.squared
        pvalue <- as.numeric(linearMod$coefficients[,4][2])
        pvalue
        
        ggplot(data, aes(x = Shannon_index, y = Predicted_value)) + 
          geom_point(size = 1) +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = 'Shanon index of the NB microbiome', 
               y = 'Predicted Shannon index \n of the NB microbiome')+ 
          ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        
        ggsave('Shannon_NB_RF.pdf',width=4, height=3)
      }
      
      input = read.csv('mom_baby_shannon_NB.filted.csv', row.names = 1)
      
      for (b in 1:length(p_cluster_list)) {
        input_3 = input
        n = which(p_cluster$p_cluster == p_cluster_list[b])
        n = row.names(p_cluster)[n]
        n = which(colnames(input) %in% n)
        input_3 = input_3[,n]
        
        trials = c(1: nrow(input_3))
        func_1 = function(trial) {
          set.seed(2022+123*trial)
          input_1 = input_3[-trial,]
          output_1 = output[-trial]
          
          input_2 = as.data.frame(rbind(input_3[trial,],input_3[trial,]))
          
          data_test = cbind(output_1,input_1)
          
          mtry_number = as.integer(ncol(input_3)/3)
          if (mtry_number <1) {
            mtry_number = 1
          }
          
          rfregFit <- train(output_1 ~ ., data = data_test, method = "ranger",importance="permutation",
                            tuneGrid = data.frame(mtry= mtry_number,
                                                  min.node.size = ncol(input_3) *2,
                                                  splitrule="variance"))
          y = predict(rfregFit, input_2)[1]
          
          return(y)
        }
        predicted_value_RF_x = mclapply(trials, func_1, mc.cores = 50)
        predicted_value_RF_x = unlist(predicted_value_RF_x)
        
        data = data.frame(Shannon_index = output, Predicted_value = predicted_value_RF_x)
        linearMod <- lm(Predicted_value ~ Shannon_index, data=data)
        linearMod = summary(linearMod)
        r_squared <- linearMod$r.squared
        pvalue <- as.numeric(linearMod$coefficients[,4][2])
        
        pvalue_RF[b] = pvalue
        r_squared_RF[b] = r_squared
      }
      
      data = data.frame(Variable_name = p_cluster_list, r_squared_RF = r_squared_RF, pvalue_RF = pvalue_RF)
      write.csv(data,'mom_baby_shannon_NB.importance.csv')
      data[order(data$pvalue_RF),]
      
    }

    ### back to local ###
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    importance_RF = read.csv('mom_baby_shannon_NB.importance.csv')
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    
    importance_RF$X = NULL
    importance_RF$Importance = -log10(importance_RF$pvalue_RF)
    colnames(importance_RF)[1] = 'Factor_cluster'
#    importance_RF$Missing_values  = missing_values_metadata_2
    write.csv(importance_RF,'mom_baby_shannon_NB_importance.csv', row.names = F)
    
    data = cbind(output,input)
    set.seed(2022)
    rfregFit_final_NB <- train(output ~ ., data = data, method = "ranger",importance="permutation",
                      tuneGrid = data.frame(mtry=as.integer(ncol(input)/3),
                                            min.node.size = ncol(input) *2,
                                            splitrule="variance"))

  }
  
}

# NR
{
  a=2
  # find paired mom factors and baby's microbiome (time closest between mom and baby)
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    input = cbind(metadata_mom$ParticipantID, metadata_2_imputation)
    colnames(input)[1] = 'ParticipantID'
    keep = input$ParticipantID %in% metadata_baby$ParticipantID
    sum(keep)
    input = input[keep,]
    
    input = input[order(input$ParticipantID),]
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    input = input[keep,]
    
    reads_table_baby_2 = reads_table_baby
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    
    output = alpha.shannon_diversity
    input = as.matrix(input[,-1])
  }
  # adonis
  {
    {
      # get reads_table
      reads_table = as.data.frame(t(reads_table_baby))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table <- reads_table$otu.tab.rff
      reads_table <- as.data.frame(reads_table)
      
      setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
      write.csv(reads_table,'mom_baby_beta_reads_table.csv')
      write.csv(input,'mom_baby_metadata_2_imputation.csv')
      
      setwd('/Users/binzhu/Desktop/Mom-baby/')
      ### below should be run on Fenn server ###
      $ conda activate ALDEx2
      #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
      $ R
      
      {
        library(vegan)
        library(parallel)
        reads_table = read.csv('mom_baby_beta_reads_table.csv', row.names = 1)
        metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation.csv', row.names = 1)
        
        pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(metadata_2_imputation)))
        row.names(pvalue_all) = colnames(metadata_2_imputation)
        
        trials = c(1: ncol(metadata_2_imputation))
        func_1 = function(trial) {
          metadata = metadata_2_imputation[,trial]
          metadata = as.data.frame(metadata)
          pvalue <- adonis2(reads_table ~ ., data = metadata, 
                            method = "bray", parallel = 60)
          pvalue <- pvalue[1,]
          return(pvalue)
        }
        pvalue = mclapply(trials, func_1, mc.cores = 60)
        
        for (a in 1:(nrow(pvalue_all))) {
          pvalue_all[a,] = pvalue[[a]]
        }
        colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
        write.csv(pvalue_all,'Adonis_NR.csv')
        
        pvalue_all = pvalue_all[order(pvalue_all$P_value, decreasing = F),]
        metadata_2_imputation = metadata_2_imputation[row.names(pvalue_all)]
        
        metadata = as.data.frame(metadata_2_imputation[,c(1:85)])
        pvalue <- adonis2(reads_table ~ ., data = metadata,method = "bray", parallel = 60)
        pvalue
        write.csv(pvalue,'multi_variable_analysis_adonis_pvalue_NR.csv')
        
      }
      
      
      setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
      data = read.csv('Adonis_NR.csv')
      data_2 = read.csv('multi_variable_analysis_adonis_pvalue_NR.csv')
      setwd('/Users/binzhu/Desktop/Mom-baby/')
      data4 <- adjust.p(data$P_value, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data4 <- data4$adjp
      data4 = data4$adjusted.p
      data$adj_p = data4
      write.csv(data,'Adonis_NR.csv')
      write.csv(data_2,'multi_variable_analysis_adonis_pvalue_NR.csv')
        
      }
      
      
      setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
      data = read.csv('Adonis_NR.csv')
      data_2 = read.csv('multi_variable_analysis_adonis_pvalue_NR.csv')
      setwd('/Users/binzhu/Desktop/Mom-baby/')
      data4 <- adjust.p(data$P.value, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data4 <- data4$adjp
      data4 = data4$adjusted.p
      data$adj_p = data4
      write.csv(data,'Adonis_NR.csv')
      write.csv(data_2,'multi_variable_analysis_adonis_pvalue_NR.csv')
    }
  
  # adonis power
  {
    # get reads_table
    reads_table = as.data.frame(t(reads_table_baby))
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
    reads_table <- reads_table$otu.tab.rff
    reads_table <- as.data.frame(reads_table)
    
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    write.csv(reads_table,'mom_baby_beta_reads_table_NR.csv')
    write.csv(input,'mom_baby_metadata_2_imputation_NR.csv')
    
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    
    ### below should be run on Fenn server ###
    $ conda activate ALDEx2
    #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
    $ R
    
    setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    reads_table = read.csv('mom_baby_beta_reads_table_NR.csv', row.names = 1)
    metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation_NR.csv', row.names = 1)
    
    setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa/sample_size_adonis')
    {
      for (trial in 1: ncol(metadata_2_imputation)) {
        file_name = paste0('NR_sample_size_',colnames(metadata_2_imputation)[trial],'.R')
        file.create(file_name)
        
        write("setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')", file=file_name,append=TRUE)
        write("reads_table = read.csv('mom_baby_beta_reads_table_NR.csv', row.names = 1)", file=file_name,append=TRUE)
        write("metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation_NR.csv', row.names = 1)", file=file_name,append=TRUE)
        write("adonis_power <- readRDS('adonis_power_function.RData')", file=file_name,append=TRUE)
        write(paste0('trial = ',trial), file=file_name, append=TRUE)
        write("metadata = metadata_2_imputation[,trial]", file=file_name,append=TRUE)
        write("power_size_model_lm = adonis_power(reads_table, metadata)", file=file_name,append=TRUE)
        write("setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa/sample_size_adonis')", file=file_name,append=TRUE)
        write("saveRDS(power_size_model_lm, file=paste0('NR_power_size_model_lm_',colnames(metadata_2_imputation)[trial],'.RData'))", file=file_name,append=TRUE)
      }
    }
    # readRDS("tsne.RData")
  }
  
  # spearman's correlation between neonatal taxa and maternal factors
  {
    reads_table_baby_2 = reads_table_baby
    reads_table_baby_2 = as.data.frame(t(reads_table_baby_2))
    reads_table_baby_2 = reads_table_baby_2 + 0.5
    reads_table_baby_2 <- as.matrix(clr(reads_table_baby_2))     ### CLR normalization in rows
    reads_table_baby_2 = as.data.frame(t(reads_table_baby_2))
    
    {
      type = 'spearman'; pvalue = 0.05; cor_parameter= 0; style = 1; bar_max = 2; 
      bar_min = -2; pheatmap_fontsize = 5; treeheight = 50; alpha = 0.05; mc.cores =8; pvalue_adj = T
      reads_table1 = as.matrix(as.data.frame(t(input)))
      reads_table2 = as.matrix(reads_table_baby_2)
      
      correlation_data_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = 0))
      colnames(correlation_data_all) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
      
      try(for (a in 1: nrow(reads_table1)) {
        correlation_data = as.data.frame(matrix(data = NA, ncol = 5, nrow = nrow(reads_table2)))
        colnames(correlation_data) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
        correlation_data$Gene = row.names(reads_table1)[a]
        correlation_data$Taxa = row.names(reads_table2)
        
        trials = c(1: nrow(reads_table2))
        
        func_1 = function(trial) {
          data1 = as.numeric(reads_table1[a,])
          data2 = as.numeric(reads_table2[trial,])
          
          keep = (!is.na(data1)) & (!is.na(data2))
          data1 = data1[keep]
          data2 = data2[keep]
          
          if (length(data1) <8) {
            return(c(list(pvalue = NA), list(Rvalue = NA)))
          }
          
          data3 = cor.test(data1, data2, method = type)
          pvalue = data3$p.value
          Rvalue = data3$estimate
          
          return(c(list(pvalue = pvalue), list(Rvalue = Rvalue)))
        }
        
        pRvalue = mclapply(trials, func_1, mc.cores = mc.cores)
        pRvalue = (unlist(pRvalue))
        keep = seq(1, length(pRvalue), 2)
        correlation_data$Pvalue = pRvalue[keep]
        
        keep = seq(2, length(pRvalue), 2)
        correlation_data$Rvalue = pRvalue[keep]
        
        if (a ==1) {
          correlation_data_all = correlation_data
        } else {
          correlation_data_all = rbind(correlation_data_all,correlation_data)
        }
      })
      
      data4 <- adjust.p(correlation_data_all$Pvalue, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data4 <- data4$adjp
      correlation_data_all$adj_p <- data4$adjusted.p
      correlation_data_all = correlation_data_all[order(correlation_data_all$Pvalue,decreasing = F),]
      write.csv(correlation_data_all,'Correlation_neonatal_taxa_factors_NR.csv', row.names = F)
    }
  }
  
  # correlation with maternal microbiome
  {
    data = data.frame(Shannon_index_NR = output, Shannon_index_M = input[,1])
    linearMod <- lm(Shannon_index_NR ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NR)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MV microbiome', 
           y = 'Shanon index of the NR microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NR_MV.pdf',width=3, height=3)
    
    data = data.frame(Shannon_index_NR = output, Shannon_index_M = input[,4])
    linearMod <- lm(Shannon_index_NR ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NR)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MR microbiome', 
           y = 'Shanon index of the NR microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NR_MR.pdf',width=3, height=3)
    
    data = data.frame(Shannon_index_NR = output, Shannon_index_M = input[,7])
    linearMod <- lm(Shannon_index_NR ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NR)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MB microbiome', 
           y = 'Shanon index of the NR microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NR_MB.pdf',width=3, height=3)
  }
  
  # random forest model
  {
    # cross validation
    ### 
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    write.csv(input,'mom_baby_shannon_NR.filted.csv')
    write.csv(output,'mom_baby_shannon_NR.output.csv')
    write.csv(p_cluster,'mom_baby_shannon_NR.p_cluster.csv')
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    ### below should be run on Fenn server ###
    $ conda activate ALDEx2
    #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
    $ R
    {
      library(pROC)
      library(doSNOW)
      library(parallel)
      library(matrixStats)
      library(caret)
      
      input = read.csv('mom_baby_shannon_NR.filted.csv', row.names = 1)
      input = input[,!str_detect(colnames(input),'baby')]
      
      output = read.csv('mom_baby_shannon_NR.output.csv', row.names = 1)
      output = output$x
      p_cluster = read.csv('mom_baby_shannon_NR.p_cluster.csv', row.names = 1)
      
      pvalue_RF = vector()
      r_squared_RF = vector()
      
      p_cluster$p_cluster = (as.character(p_cluster$p_cluster))
      p_cluster_list = unique(p_cluster$p_cluster)
      
      # RF
      {
        trials = c(1: nrow(input))
        func_1 = function(trial) {
          input_1 = input[-trial,]
          output_1 = output[-trial]
          
          input_2 = as.data.frame(rbind(input[trial,],input[trial,]))
          
          data = cbind(output_1,input_1)
          set.seed(2022)
          mtry_number = as.integer(ncol(input_3)/3)
          rfregFit <- train(output_1 ~ ., data = data, method = "ranger",importance="permutation",
                            tuneGrid = data.frame(mtry=mtry_number,
                                                  min.node.size = ncol(input) *2,
                                                  splitrule="variance"))
          
          predicted_value_RF_x = predict(rfregFit, input_2, na.action = na.pass)[1]
          
          return(predicted_value_RF_x)
        }
        predicted_value_RF_x = mclapply(trials, func_1, mc.cores = 50)
        predicted_value_RF_x = unlist(predicted_value_RF_x)
        
        data = data.frame(Shannon_index = output, Predicted_value = predicted_value_RF_x)
        linearMod <- lm(Predicted_value ~ Shannon_index, data=data)
        linearMod = summary(linearMod)
        r_squared <- linearMod$r.squared
        pvalue <- as.numeric(linearMod$coefficients[,4][2])
        pvalue
        
        ggplot(data, aes(x = Shannon_index, y = Predicted_value)) + 
          geom_point(size = 1) +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = 'Shanon index of the NR microbiome', 
               y = 'Predicted Shannon index \n of the NR microbiome')+ 
          ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        
        ggsave('Shannon_NR_RF.pdf',width=4, height=3)
      }
      
      input = read.csv('mom_baby_shannon_NR.filted.csv', row.names = 1)
      for (b in 1:length(p_cluster_list)) {
        input_3 = input
        n = which(p_cluster$p_cluster == p_cluster_list[b])
        n = row.names(p_cluster)[n]
        n = which(colnames(input) %in% n)
        input_3 = input_3[,n]
        
        trials = c(1: nrow(input_3))
        func_1 = function(trial) {
          set.seed(2022+123*trial)
          input_1 = input_3[-trial,]
          output_1 = output[-trial]
          
          input_2 = as.data.frame(rbind(input_3[trial,],input_3[trial,]))
          
          data_test = cbind(output_1,input_1)
          
          mtry_number = as.integer(ncol(input_3)/3)
          if (mtry_number <1) {
            mtry_number = 1
          }
          
          rfregFit <- train(output_1 ~ ., data = data_test, method = "ranger",importance="permutation",
                            tuneGrid = data.frame(mtry= mtry_number,
                                                  min.node.size = ncol(input_3) *2,
                                                  splitrule="variance"))
          y = predict(rfregFit, input_2)[1]
          
          return(y)
        }
        predicted_value_RF_x = mclapply(trials, func_1, mc.cores = 50)
        predicted_value_RF_x = unlist(predicted_value_RF_x)
        
        data = data.frame(Shannon_index = output, Predicted_value = predicted_value_RF_x)
        linearMod <- lm(Predicted_value ~ Shannon_index, data=data)
        linearMod = summary(linearMod)
        r_squared <- linearMod$r.squared
        pvalue <- as.numeric(linearMod$coefficients[,4][2])
        
        pvalue_RF[b] = pvalue
        r_squared_RF[b] = r_squared
      }
      
      data = data.frame(Variable_name = p_cluster_list, r_squared_RF = r_squared_RF, pvalue_RF = pvalue_RF)
      write.csv(data,'mom_baby_shannon_NR.importance.csv')
      data[order(data$pvalue_RF),]
      
    }
    
    ### back to local ###
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    importance_RF = read.csv('mom_baby_shannon_NR.importance.csv')
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    
    importance_RF$X = NULL
    importance_RF$Importance = -log10(importance_RF$pvalue_RF)
    colnames(importance_RF)[1] = 'Factor_cluster'
    #    importance_RF$Missing_values  = missing_values_metadata_2
    write.csv(importance_RF,'mom_baby_shannon_NR_importance.csv', row.names = F)
    
    data = cbind(output,input)
    set.seed(2022)
    rfregFit_final_NR <- train(output ~ ., data = data, method = "ranger",importance="permutation",
                               tuneGrid = data.frame(mtry=as.integer(ncol(input)/3),
                                                     min.node.size = ncol(input) *2,
                                                     splitrule="variance"))
    
  }
  
  
  
  
}


# NS
{
  a=3
  # find paired mom factors and baby's microbiome (time closest between mom and baby)
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    input = cbind(metadata_mom$ParticipantID, metadata_2_imputation)
    colnames(input)[1] = 'ParticipantID'
    keep = input$ParticipantID %in% metadata_baby$ParticipantID
    sum(keep)
    input = input[keep,]
    
    input = input[order(input$ParticipantID),]
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    input = input[keep,]
    
    reads_table_baby_2 = reads_table_baby
    reads_table_baby_2 = t(reads_table_baby_2)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    
    output = alpha.shannon_diversity
    input = as.matrix(input[,-1])
  }
  # adonis
  {
    {
      # get reads_table
      reads_table = as.data.frame(t(reads_table_baby))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table <- reads_table$otu.tab.rff
      reads_table <- as.data.frame(reads_table)
      
      setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
      write.csv(reads_table,'mom_baby_beta_reads_table.csv')
      write.csv(input,'mom_baby_metadata_2_imputation.csv')
      
      setwd('/Users/binzhu/Desktop/Mom-baby/')
      ### below should be run on Fenn server ###
      $ conda activate ALDEx2
      #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
      $ R
      
      {
        library(vegan)
        library(parallel)
        reads_table = read.csv('mom_baby_beta_reads_table.csv', row.names = 1)
        metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation.csv', row.names = 1)
        
        pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(metadata_2_imputation)))
        row.names(pvalue_all) = colnames(metadata_2_imputation)
        
        trials = c(1: ncol(metadata_2_imputation))
        func_1 = function(trial) {
          metadata = metadata_2_imputation[,trial]
          metadata = as.data.frame(metadata)
          pvalue <- adonis2(reads_table ~ ., data = metadata, 
                            method = "bray", parallel = 60)
          pvalue <- pvalue[1,]
          return(pvalue)
        }
        pvalue = mclapply(trials, func_1, mc.cores = 60)
        
        for (a in 1:(nrow(pvalue_all))) {
          pvalue_all[a,] = pvalue[[a]]
        }
        colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
        write.csv(pvalue_all,'Adonis_NS.csv')
        
        pvalue_all = pvalue_all[order(pvalue_all$P_value, decreasing = F),]
        metadata_2_imputation = metadata_2_imputation[row.names(pvalue_all)]
        
        metadata = as.data.frame(metadata_2_imputation[,c(1:57)])
        pvalue <- adonis2(reads_table ~ ., data = metadata,method = "bray", parallel = 60)
        pvalue
        write.csv(pvalue,'multi_variable_analysis_adonis_pvalue_NS.csv')
        
      }
      
      
      setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
      data = read.csv('Adonis_NS.csv')
      data_2 = read.csv('multi_variable_analysis_adonis_pvalue_NS.csv')
      setwd('/Users/binzhu/Desktop/Mom-baby/')
      data4 <- adjust.p(data$P_value, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data4 <- data4$adjp
      data4 = data4$adjusted.p
      data$adj_p = data4
      write.csv(data,'Adonis_NS.csv')
      write.csv(data_2,'multi_variable_analysis_adonis_pvalue_NS.csv')
    }
  }
  
  # adonis power
  {
    # get reads_table
    reads_table = as.data.frame(t(reads_table_baby))
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
    reads_table <- reads_table$otu.tab.rff
    reads_table <- as.data.frame(reads_table)
    
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    write.csv(reads_table,'mom_baby_beta_reads_table_NS.csv')
    write.csv(input,'mom_baby_metadata_2_imputation_NS.csv')
    
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    
    ### below should be run on Fenn server ###
    $ conda activate ALDEx2
    #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
    $ R
    
    setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    reads_table = read.csv('mom_baby_beta_reads_table_NS.csv', row.names = 1)
    metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation_NS.csv', row.names = 1)
    
    setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa/sample_size_adonis')
    {
      for (trial in 1: ncol(metadata_2_imputation)) {
        file_name = paste0('NS_sample_size_',colnames(metadata_2_imputation)[trial],'.R')
        file.create(file_name)
        
        write("setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')", file=file_name,append=TRUE)
        write("reads_table = read.csv('mom_baby_beta_reads_table_NS.csv', row.names = 1)", file=file_name,append=TRUE)
        write("metadata_2_imputation = read.csv('mom_baby_metadata_2_imputation_NS.csv', row.names = 1)", file=file_name,append=TRUE)
        write("adonis_power <- readRDS('adonis_power_function.RData')", file=file_name,append=TRUE)
        write(paste0('trial = ',trial), file=file_name, append=TRUE)
        write("metadata = metadata_2_imputation[,trial]", file=file_name,append=TRUE)
        write("power_size_model_lm = adonis_power(reads_table, metadata)", file=file_name,append=TRUE)
        write("setwd('/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa/sample_size_adonis')", file=file_name,append=TRUE)
        write("saveRDS(power_size_model_lm, file=paste0('NS_power_size_model_lm_',colnames(metadata_2_imputation)[trial],'.RData'))", file=file_name,append=TRUE)
      }
    }
    # readRDS("tsne.RData")
  }
  
  
  # spearman's correlation between neonatal taxa and maternal factors
  {
    reads_table_baby_2 = reads_table_baby
    reads_table_baby_2 = as.data.frame(t(reads_table_baby_2))
    reads_table_baby_2 = reads_table_baby_2 + 0.5
    reads_table_baby_2 <- as.matrix(clr(reads_table_baby_2))     ### CLR normalization in rows
    reads_table_baby_2 = as.data.frame(t(reads_table_baby_2))
    
    {
      type = 'spearman'; pvalue = 0.05; cor_parameter= 0; style = 1; bar_max = 2; 
      bar_min = -2; pheatmap_fontsize = 5; treeheight = 50; alpha = 0.05; mc.cores =8; pvalue_adj = T
      reads_table1 = as.matrix(as.data.frame(t(input)))
      reads_table2 = as.matrix(reads_table_baby_2)
      
      correlation_data_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = 0))
      colnames(correlation_data_all) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
      
      try(for (a in 1: nrow(reads_table1)) {
        correlation_data = as.data.frame(matrix(data = NA, ncol = 5, nrow = nrow(reads_table2)))
        colnames(correlation_data) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
        correlation_data$Gene = row.names(reads_table1)[a]
        correlation_data$Taxa = row.names(reads_table2)
        
        trials = c(1: nrow(reads_table2))
        
        func_1 = function(trial) {
          data1 = as.numeric(reads_table1[a,])
          data2 = as.numeric(reads_table2[trial,])
          
          keep = (!is.na(data1)) & (!is.na(data2))
          data1 = data1[keep]
          data2 = data2[keep]
          
          if (length(data1) <8) {
            return(c(list(pvalue = NA), list(Rvalue = NA)))
          }
          
          data3 = cor.test(data1, data2, method = type)
          pvalue = data3$p.value
          Rvalue = data3$estimate
          
          return(c(list(pvalue = pvalue), list(Rvalue = Rvalue)))
        }
        
        pRvalue = mclapply(trials, func_1, mc.cores = mc.cores)
        pRvalue = (unlist(pRvalue))
        keep = seq(1, length(pRvalue), 2)
        correlation_data$Pvalue = pRvalue[keep]
        
        keep = seq(2, length(pRvalue), 2)
        correlation_data$Rvalue = pRvalue[keep]
        
        if (a ==1) {
          correlation_data_all = correlation_data
        } else {
          correlation_data_all = rbind(correlation_data_all,correlation_data)
        }
      })
      
      data4 <- adjust.p(correlation_data_all$Pvalue, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data4 <- data4$adjp
      correlation_data_all$adj_p <- data4$adjusted.p
      correlation_data_all = correlation_data_all[order(correlation_data_all$Pvalue,decreasing = F),]
      write.csv(correlation_data_all,'Correlation_neonatal_taxa_factors_NS.csv', row.names = F)
    }
  }
  
  # correlation with maternal microbiome
  {
    data = data.frame(Shannon_index_NS = output, Shannon_index_M = input[,1])
    linearMod <- lm(Shannon_index_NS ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NS)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MV microbiome', 
           y = 'Shanon index of the NS microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NS_MV.pdf',width=3, height=3)
    
    data = data.frame(Shannon_index_NS = output, Shannon_index_M = input[,4])
    linearMod <- lm(Shannon_index_NS ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NS)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MR microbiome', 
           y = 'Shanon index of the NS microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NS_MR.pdf',width=3, height=3)
    
    data = data.frame(Shannon_index_NS = output, Shannon_index_M = input[,7])
    linearMod <- lm(Shannon_index_NS ~ Shannon_index_M, data=data)
    linearMod = summary(linearMod)
    r_squared <- linearMod$r.squared
    pvalue <- as.numeric(linearMod$coefficients[,4][2])
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    
    ggplot(data, aes(x = Shannon_index_M, y = Shannon_index_NS)) + 
      geom_point() +
      stat_smooth(method = "lm", col = "blue")+
      labs(x = 'Shanon index of the MB microbiome', 
           y = 'Shanon index of the NS microbiome')+ 
      ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+ theme_bw()
    
    ggsave('Shannon_NS_MB.pdf',width=3, height=3)
  }
  
  # random forest model
  {
    # cross validation
    ### 
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    write.csv(input,'mom_baby_shannon_NS.filted.csv')
    write.csv(output,'mom_baby_shannon_NS.output.csv')
    write.csv(p_cluster,'mom_baby_shannon_NS.p_cluster.csv')
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    ### below should be run on Fenn server ###
    $ conda activate ALDEx2
    #$ export PATH="/usr/global/R-4.0.2/bin:$PATH"
    $ R
    {
      library(pROC)
      library(doSNOW)
      library(parallel)
      library(matrixStats)
      library(caret)
      
      input = read.csv('mom_baby_shannon_NS.filted.csv', row.names = 1)
      input = input[,!str_detect(colnames(input),'baby')]
      output = read.csv('mom_baby_shannon_NS.output.csv', row.names = 1)
      output = output$x
      p_cluster = read.csv('mom_baby_shannon_NS.p_cluster.csv', row.names = 1)
      
      pvalue_RF = vector()
      r_squared_RF = vector()
      
      p_cluster$p_cluster = (as.character(p_cluster$p_cluster))
      p_cluster_list = unique(p_cluster$p_cluster)
      
      # RF
      {
        trials = c(1: nrow(input))
        func_1 = function(trial) {
          input_1 = input[-trial,]
          output_1 = output[-trial]
          
          input_2 = as.data.frame(rbind(input[trial,],input[trial,]))
          
          data = cbind(output_1,input_1)
          set.seed(2022)
          mtry_number = as.integer(ncol(input_3)/3)
          rfregFit <- train(output_1 ~ ., data = data, method = "ranger",importance="permutation",
                            tuneGrid = data.frame(mtry=mtry_number,
                                                  min.node.size = ncol(input) *2,
                                                  splitrule="variance"))
          
          predicted_value_RF_x = predict(rfregFit, input_2, na.action = na.pass)[1]
          
          return(predicted_value_RF_x)
        }
        predicted_value_RF_x = mclapply(trials, func_1, mc.cores = 50)
        predicted_value_RF_x = unlist(predicted_value_RF_x)
        
        data = data.frame(Shannon_index = output, Predicted_value = predicted_value_RF_x)
        linearMod <- lm(Predicted_value ~ Shannon_index, data=data)
        linearMod = summary(linearMod)
        r_squared <- linearMod$r.squared
        pvalue <- as.numeric(linearMod$coefficients[,4][2])
        pvalue
        
        ggplot(data, aes(x = Shannon_index, y = Predicted_value)) + 
          geom_point(size = 1) +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = 'Shanon index of the NS microbiome', 
               y = 'Predicted Shannon index \n of the NS microbiome')+ 
          ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        
        ggsave('Shannon_NS_RF.pdf',width=4, height=3)
      }
      
      input = read.csv('mom_baby_shannon_NS.filted.csv', row.names = 1)
      for (b in 1:length(p_cluster_list)) {
        input_3 = input
        n = which(p_cluster$p_cluster == p_cluster_list[b])
        n = row.names(p_cluster)[n]
        n = which(colnames(input) %in% n)
        input_3 = input_3[,n]
        
        trials = c(1: nrow(input_3))
        func_1 = function(trial) {
          set.seed(2022+123*trial)
          input_1 = input_3[-trial,]
          output_1 = output[-trial]
          
          input_2 = as.data.frame(rbind(input_3[trial,],input_3[trial,]))
          
          data_test = cbind(output_1,input_1)
          
          mtry_number = as.integer(ncol(input_3)/3)
          if (mtry_number <1) {
            mtry_number = 1
          }
          
          rfregFit <- train(output_1 ~ ., data = data_test, method = "ranger",importance="permutation",
                            tuneGrid = data.frame(mtry= mtry_number,
                                                  min.node.size = ncol(input_3) *2,
                                                  splitrule="variance"))
          y = predict(rfregFit, input_2)[1]
          
          return(y)
        }
        predicted_value_RF_x = mclapply(trials, func_1, mc.cores = 50)
        predicted_value_RF_x = unlist(predicted_value_RF_x)
        
        data = data.frame(Shannon_index = output, Predicted_value = predicted_value_RF_x)
        linearMod <- lm(Predicted_value ~ Shannon_index, data=data)
        linearMod = summary(linearMod)
        r_squared <- linearMod$r.squared
        pvalue <- as.numeric(linearMod$coefficients[,4][2])
        
        pvalue_RF[b] = pvalue
        r_squared_RF[b] = r_squared
      }
      
      data = data.frame(Variable_name = p_cluster_list, r_squared_RF = r_squared_RF, pvalue_RF = pvalue_RF)
      write.csv(data,'mom_baby_shannon_NS.importance.csv')
      data[order(data$pvalue_RF),]
      
    }
    
    ### back to local ###
    setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
    importance_RF = read.csv('mom_baby_shannon_NS.importance.csv')
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    
    importance_RF$X = NULL
    importance_RF$Importance = -log10(importance_RF$pvalue_RF)
    colnames(importance_RF)[1] = 'Factor_cluster'
    #    importance_RF$Missing_values  = missing_values_metadata_2
    write.csv(importance_RF,'mom_baby_shannon_NS_importance.csv', row.names = F)
    
    data = cbind(output,input)
    set.seed(2022)
    rfregFit_final_NS <- train(output ~ ., data = data, method = "ranger",importance="permutation",
                               tuneGrid = data.frame(mtry=as.integer(ncol(input)/3),
                                                     min.node.size = ncol(input) *2,
                                                     splitrule="variance"))
    
  }
  
  
  
  
}

# sample size
{
  {
    setwd('/Users/binzhu/Desktop/Mom-baby/results_multivariates/')
    NB_effect_size = read.csv('Adonis_NB.csv', row.names = 1)
    NB_effect_size$power = NA
    
    setwd('/Users/binzhu/Desktop/Mom-baby/results_multivariates/sample_size_power_adonis_test')
    {
      file_list = list.files()
      for (a in 1: nrow(NB_effect_size)) {
        factor_name = NB_effect_size$X[a]
        effect_size = NB_effect_size$R2[a]
        effect_size
        
        file_name = paste0('NB_power_size_model_lm_',factor_name,'.RData')
        
        if (file_name %in% file_list) {
          bp30_model = readRDS(file_name)
          
          size_power = 10^predict(bp30_model, newdata=data.frame(log_omega2=log10(effect_size)))
          size_power
          NB_effect_size$power[a] = size_power
        }
        
      }
    }
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    write.csv(NB_effect_size,'Adonis_NB.csv')
  }
  {
    setwd('/Users/binzhu/Desktop/Mom-baby/results_multivariates/')
    NR_effect_size = read.csv('Adonis_NR.csv', row.names = 1)
    NR_effect_size$power = NA
    
    setwd('/Users/binzhu/Desktop/Mom-baby/results_multivariates/sample_size_power_adonis_test')
    {
      file_list = list.files()
      for (a in 1: nrow(NR_effect_size)) {
        factor_name = NR_effect_size$X[a]
        effect_size = NR_effect_size$R2[a]
        effect_size
        
        file_name = paste0('NR_power_size_model_lm_',factor_name,'.RData')
        
        if (file_name %in% file_list) {
          bp30_model = readRDS(file_name)
          
          size_power = 10^predict(bp30_model, newdata=data.frame(log_omega2=log10(effect_size)))
          size_power
          NR_effect_size$power[a] = size_power
        }
        
      }
    }
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    write.csv(NR_effect_size,'Adonis_NR.csv')
  }
  {
    setwd('/Users/binzhu/Desktop/Mom-baby/results_multivariates/')
    NS_effect_size = read.csv('Adonis_NS.csv', row.names = 1)
    NS_effect_size$power = NA
    
    setwd('/Users/binzhu/Desktop/Mom-baby/results_multivariates/sample_size_power_adonis_test')
    {
      file_list = list.files()
      for (a in 1: nrow(NS_effect_size)) {
        factor_name = NS_effect_size$X[a]
        effect_size = NS_effect_size$R2[a]
        effect_size
        
        file_name = paste0('NS_power_size_model_lm_',factor_name,'.RData')
        
        if (file_name %in% file_list) {
          bp30_model = readRDS(file_name)
          
          size_power = 10^predict(bp30_model, newdata=data.frame(log_omega2=log10(effect_size)))
          size_power
          NS_effect_size$power[a] = size_power
        }
        
      }
    }
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    write.csv(NS_effect_size,'Adonis_NS.csv')
  }
}




















