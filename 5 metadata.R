setwd('/Users/binzhu/Desktop/Mom-baby/results/script_5/metadata')
color1=c('#E43B2D','#0070FF','#4DAF4A','#FEF840','#F27E33','#F180BF')
##### compare sample collection time of maternal microbiome and metadata #####
# get case-match
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
                            BS1D_birth_0_day = NA, BS1D_birth_1_days = NA,BS1D_birth_2_days = NA, 
                            cytokine = NA, lipid = NA)
  for (a in 1:nrow(metadata_2)) {
    x = which(output_table$participant == metadata_2$ParticipantID[a])
    y = which(colnames(output_table) == metadata_2$SampleType[a])
    output_table[x,y] = metadata_2$SampleID[a]
  }
  
  output_table_2 = data.frame(participant = unique(metadata_2$ParticipantID), MCKD= NA, MRCD=NA, MV1D = NA)
  for (a in 1:nrow(metadata_2)) {
    x = which(output_table_2$participant == metadata_2$ParticipantID[a])
    y = which(colnames(output_table_2) == metadata_2$SampleType[a])
    output_table_2[x,y] = metadata_2$weeks_pregnant[a]
  }
  
  for (a in 1:nrow(metadata)) {
    x = which(output_table$participant == metadata$ParticipantID[a])
    y = which(colnames(output_table) == metadata$Flag[a])
    output_table[x,y] = metadata$SampleID[a]
  }
  
  metadata_cytokine = metadata_16s[metadata_16s$SampleID %in% colnames(reads_table_cytokines),]
  metadata_cytokine = metadata_cytokine[order(metadata_cytokine$VisitNum, decreasing = T),]
  metadata_cytokine = metadata_cytokine[!duplicated(metadata_cytokine$ParticipantID),]
  
  for (a in 1:nrow(metadata_cytokine)) {
    x = which(output_table$participant == metadata_cytokine$ParticipantID[a])
    y = which(colnames(output_table) == 'cytokine')
    output_table[x,y] = metadata_cytokine$SampleID[a]
  }
  metadata_cytokine = metadata_cytokine[match(output_table$participant,metadata_cytokine$ParticipantID),]
  
  metadata_lipid = metadata_16s[metadata_16s$SampleID %in% colnames(reads_table_lipid),]
  metadata_lipid = metadata_lipid[order(metadata_lipid$VisitNum, decreasing = T),]
  metadata_lipid = metadata_lipid[!duplicated(metadata_lipid$ParticipantID),]
  
  for (a in 1:nrow(metadata_lipid)) {
    x = which(output_table$participant == metadata_lipid$ParticipantID[a])
    y = which(colnames(output_table) == 'lipid')
    output_table[x,y] = metadata_lipid$SampleID[a]
  }
  metadata_lipid = metadata_lipid[match(output_table$participant,metadata_lipid$ParticipantID),]
}
#
{
  data = output_table_2[,c(2,3,4)]
  data_2 = data; data_2[is.na(data_2)] =100
  data = cbind(data, apply(data_2,1,min))
  data$lipid = metadata_lipid$weeks_pregnant; data$cytokine = metadata_cytokine$weeks_pregnant;
  data = gather(data)
  data$Microbiome = c(rep("MB",nrow(output_table)),rep("MR",nrow(output_table)),
                      rep("MV",nrow(output_table)),rep("Metadata",nrow(output_table)),
                      rep("Lipid",nrow(output_table)),rep("Cytokine",nrow(output_table)))
  colnames(data)[2] = 'weeks_pregnant'
  data = data[!is.na(data$weeks_pregnant),]
  data = data[data$weeks_pregnant != 100,]
  
  data$Microbiome = factor(data$Microbiome, levels = c('MB','MR','MV','Metadata','Lipid','Cytokine'))
  data = data[data$Microbiome != 'Lipid',]
  ggplot(data, aes(x=Microbiome, y=weeks_pregnant)) +
    geom_violin(trim=T,aes(fill=Microbiome))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1)+ scale_fill_manual(values=c(color1[c(4,5,6)], '#FFACAC','#A1FF6B','#7471FF')) +theme_bw()+
    labs(x = NULL, y = "Weeks of pregnancy")+
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))
  ggsave('microbiome_weeks_pregnant_compare.pdf',width = 4, height = 2)
  
  data = output_table_2[,c(2,3,4)]
  data_2 = data; data_2[is.na(data_2)] =100
  data = cbind(data, apply(data_2,1,min))
  data$lipid = metadata_lipid$weeks_pregnant; data$cytokine = metadata_cytokine$weeks_pregnant;
  data = gather(data)
  data$Microbiome = c(rep("MB",nrow(output_table)),rep("MR",nrow(output_table)),
                      rep("MV",nrow(output_table)),rep("Metadata",nrow(output_table)),
                      rep("Lipid",nrow(output_table)),rep("Cytokine",nrow(output_table)))
  colnames(data)[2] = 'weeks_pregnant'
  data = data[!is.na(data$weeks_pregnant),]
  data = data[data$weeks_pregnant != 100,]
  
  data$Microbiome = factor(data$Microbiome, levels = c('MB','MR','MV','Metadata','Lipid','Cytokine'))
  data = data[data$Microbiome != 'Lipid',]
  ggplot(data, aes(x=Microbiome, y=weeks_pregnant)) +
    geom_violin(trim=T,aes(fill=Microbiome))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1)+ scale_fill_manual(values=c(color1[c(4,5,6)], '#FFACAC','#A1FF6B','#7471FF')) +theme_bw()+
    labs(x = NULL, y = "Weeks of pregnancy")+
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))
  ggsave('microbiome_weeks_pregnant_compare.pdf',width = 4, height = 2)
  
  p = vector()
  p[1] = paste0('MB MR: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MB'],data$weeks_pregnant[data$Microbiome == 'MR'])$p.value)
  p[2] = paste0('MB MV: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MB'],data$weeks_pregnant[data$Microbiome == 'MV'])$p.value)
  p[3] = paste0('MV MR: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MV'],data$weeks_pregnant[data$Microbiome == 'MR'])$p.value)
  p[4] = paste0('MB Metadata: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MB'],data$weeks_pregnant[data$Microbiome == 'Metadata'])$p.value)
  p[5] = paste0('MR Metadata: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MR'],data$weeks_pregnant[data$Microbiome == 'Metadata'])$p.value)
  p[6] = paste0('MV Metadata: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MV'],data$weeks_pregnant[data$Microbiome == 'Metadata'])$p.value)
  
  p[11] = paste0('Metadata Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'Metadata'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)
  p[12] = paste0('MB Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MB'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)
  p[13] = paste0('MV Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MV'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)
  p[14] = paste0('MR Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MR'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)

  write.table(p,'microbiome_weeks_pregnant_compare.txt', quote = F,row.names = F, col.names = F)
}

##### convert microbiome data to numeric data #####
# use t-SNE coordinate to represent the difference of main components in the maternal microbiome
pic_3 = pic[,c(1:3)] # pic is tested in script 4
{
  # MB
  {
    reads_table = as.matrix(t(reads_table_16s_MCKD))
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
    reads_table = as.data.frame(reads_table$otu.tab.rff)
    
    Shannon_MB <- as.numeric(as.character(diversity(reads_table)))
    Evenness_MB <- Shannon_MB/log(specnumber(reads_table))
    Number_of_taxa_MB <- as.numeric(as.character(colSums(t(reads_table) != 0)))
    
    metadata_16s$Shannon_MB = NA; metadata_16s$Evenness_MB = NA;metadata_16s$Number_of_taxa_MB = NA;
    metadata_16s$MB_t_SNE1 = NA;metadata_16s$MB_t_SNE2 = NA;
    for (a in 1: ncol(reads_table_16s_MCKD)) {
      n = which(metadata_16s$SampleID == colnames(reads_table_16s_MCKD)[a])
      metadata_16s$Shannon_MB[n] = Shannon_MB[a]
      metadata_16s$Evenness_MB[n] = Evenness_MB[a]
      metadata_16s$Number_of_taxa_MB[n] = Number_of_taxa_MB[a]
      
      m = which(pic_3$SampleID == colnames(reads_table_16s_MCKD)[a])
      metadata_16s$MB_t_SNE1[n] = pic_3$X1[m]
      metadata_16s$MB_t_SNE2[n] = pic_3$X2[m]
    }
    
    # differential abundance of MB_t_SNE1
    {
      reads_table = reads_table_16s_MCKD[,colnames(reads_table_16s_MCKD) %in% 
                                           metadata_16s$SampleID[!is.na(metadata_16s$MB_t_SNE1)]]
      
      metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
      reads_table = reads_table[metadata$SampleID]
      metadata = metadata$X1
      keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                                         paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
      output$p
      ggsave('MB_t_SNE1_difference.pdf', width = 5, height = 5)
    }
    # differential abundance of MB_t_SNE2
    {
      reads_table = reads_table_16s_MCKD[,colnames(reads_table_16s_MCKD) %in% 
                                           metadata_16s$SampleID[!is.na(metadata_16s$MB_t_SNE2)]]
      
      metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
      reads_table = reads_table[metadata$SampleID]
      metadata = metadata$X2
      keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                             paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
      output$p
      ggsave('MB_t_SNE2_difference.pdf', width = 5, height = 6)
    }
  }
  # MR
  {
    reads_table = as.matrix(t(reads_table_16s_MRCD))
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
    reads_table = as.data.frame(reads_table$otu.tab.rff)
    
    Shannon_MR <- as.numeric(as.character(diversity(reads_table)))
    Evenness_MR <- Shannon_MR/log(specnumber(reads_table))
    Number_of_taxa_MR <- as.numeric(as.character(colSums(t(reads_table) != 0)))
    
    metadata_16s$Shannon_MR = NA; metadata_16s$Evenness_MR = NA;metadata_16s$Number_of_taxa_MR = NA;
    metadata_16s$MR_t_SNE1 = NA;metadata_16s$MR_t_SNE2 = NA;
    for (a in 1: ncol(reads_table_16s_MRCD)) {
      n = which(metadata_16s$SampleID == colnames(reads_table_16s_MRCD)[a])
      metadata_16s$Shannon_MR[n] = Shannon_MR[a]
      metadata_16s$Evenness_MR[n] = Evenness_MR[a]
      metadata_16s$Number_of_taxa_MR[n] = Number_of_taxa_MR[a]
      
      m = which(pic_3$SampleID == colnames(reads_table_16s_MRCD)[a])
      metadata_16s$MR_t_SNE1[n] = pic_3$X1[m]
      metadata_16s$MR_t_SNE2[n] = pic_3$X2[m]
    }
    
    # differential abundance of MR_t_SNE1
    {
      reads_table = reads_table_16s_MRCD[,colnames(reads_table_16s_MRCD) %in% 
                                           metadata_16s$SampleID[!is.na(metadata_16s$MR_t_SNE1)]]
      
      metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
      reads_table = reads_table[metadata$SampleID]
      metadata = metadata$X1
      keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                             paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
      output$p
      ggsave('MR_t_SNE1_difference.pdf', width = 5, height = 10)
    }
    # differential abundance of MR_t_SNE2
    {
      reads_table = reads_table_16s_MRCD[,colnames(reads_table_16s_MRCD) %in% 
                                           metadata_16s$SampleID[!is.na(metadata_16s$MR_t_SNE2)]]
      
      metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
      reads_table = reads_table[metadata$SampleID]
      metadata = metadata$X2
      keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                             paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
      output$p
      ggsave('MR_t_SNE2_difference.pdf', width = 5, height = 15)
    }
  }
  # MV
  {
    reads_table = as.matrix(t(reads_table_16s_MV1D))
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
    reads_table = as.data.frame(reads_table$otu.tab.rff)
    
    Shannon_MV <- as.numeric(as.character(diversity(reads_table)))
    Evenness_MV <- Shannon_MV/log(specnumber(reads_table))
    Number_of_taxa_MV <- as.numeric(as.character(colSums(t(reads_table) != 0)))
    
    metadata_16s$Shannon_MV = NA; metadata_16s$Evenness_MV = NA;metadata_16s$Number_of_taxa_MV = NA;
    metadata_16s$MV_t_SNE1 = NA;metadata_16s$MV_t_SNE2 = NA;
    for (a in 1: ncol(reads_table_16s_MV1D)) {
      n = which(metadata_16s$SampleID == colnames(reads_table_16s_MV1D)[a])
      metadata_16s$Shannon_MV[n] = Shannon_MV[a]
      metadata_16s$Evenness_MV[n] = Evenness_MV[a]
      metadata_16s$Number_of_taxa_MV[n] = Number_of_taxa_MV[a]
      
      m = which(pic_3$SampleID == colnames(reads_table_16s_MV1D)[a])
      metadata_16s$MV_t_SNE1[n] = pic_3$X1[m]
      metadata_16s$MV_t_SNE2[n] = pic_3$X2[m]
    }
    
    # differential abundance of MV_t_SNE1
    {
      reads_table = reads_table_16s_MV1D[,colnames(reads_table_16s_MV1D) %in% 
                                           metadata_16s$SampleID[!is.na(metadata_16s$MV_t_SNE1)]]
      
      metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
      reads_table = reads_table[metadata$SampleID]
      metadata = metadata$X1
      keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                             paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
      output$p
      ggsave('MV_t_SNE1_difference.pdf', width = 5, height = 3)
    }
    # differential abundance of MV_t_SNE2
    {
      reads_table = reads_table_16s_MV1D[,colnames(reads_table_16s_MV1D) %in% 
                                           metadata_16s$SampleID[!is.na(metadata_16s$MV_t_SNE2)]]
      
      metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
      reads_table = reads_table[metadata$SampleID]
      metadata = metadata$X2
      keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                             paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
      output$p
      ggsave('MV_t_SNE2_difference.pdf', width = 5, height = 5)
    }
  }
  
}

{
  # NB
  {
    # day 0
    {
      reads_table_16s_BCKD_day_0 = reads_table_16s_BCKD[,colnames(reads_table_16s_BCKD) %in% output_table$BCKD_birth_0_day[!is.na(output_table$BCKD_birth_0_day)]]
      reads_table = as.matrix(t(reads_table_16s_BCKD_day_0))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table = as.data.frame(reads_table$otu.tab.rff)
      
      Shannon_NB <- as.numeric(as.character(diversity(reads_table)))
      Evenness_NB <- Shannon_NB/log(specnumber(reads_table))
      Number_of_taxa_NB <- as.numeric(as.character(colSums(t(reads_table) != 0)))
      
      metadata_16s$Shannon_NB_day_0 = NA; metadata_16s$Evenness_NB_day_0 = NA;metadata_16s$Number_of_taxa_NB_day_0 = NA;
      metadata_16s$NB_t_SNE1_day_0 = NA;metadata_16s$NB_t_SNE2_day_0 = NA;
      for (a in 1: ncol(reads_table_16s_BCKD_day_0)) {
        n = which(metadata_16s$SampleID == colnames(reads_table_16s_BCKD_day_0)[a])
        metadata_16s$Shannon_NB_day_0[n] = Shannon_NB[a]
        metadata_16s$Evenness_NB_day_0[n] = Evenness_NB[a]
        metadata_16s$Number_of_taxa_NB_day_0[n] = Number_of_taxa_NB[a]
        
        m = which(pic_3$SampleID == colnames(reads_table_16s_BCKD_day_0)[a])
        metadata_16s$NB_t_SNE1_day_0[n] = pic_3$X1[m]
        metadata_16s$NB_t_SNE2_day_0[n] = pic_3$X2[m]
      }
      
      # differential abundance of NB_t_SNE1_day_0
      {
        reads_table = reads_table_16s_BCKD_day_0[,colnames(reads_table_16s_BCKD_day_0) %in% 
                                             metadata_16s$SampleID[!is.na(metadata_16s$NB_t_SNE1_day_0)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X1
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        write.csv(output$data, 'NB_t_SNE1_day_0.csv')
        
      }
      # differential abundance of NB_t_SNE2_day_0
      {
        reads_table = reads_table_16s_BCKD_day_0[,colnames(reads_table_16s_BCKD_day_0) %in% 
                                                   metadata_16s$SampleID[!is.na(metadata_16s$NB_t_SNE1_day_0)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X2
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        output$p
        ggsave('NB_t_SNE2_day_0_difference.pdf', width = 5, height = 3.8)
      }      
      
    }
    # days 1 2
    {
      reads_table_16s_BCKD_days_12 = reads_table_16s_BCKD[,colnames(reads_table_16s_BCKD) %in% 
                                                            c(output_table$BCKD_birth_1_days[!is.na(output_table$BCKD_birth_1_days)],output_table$BCKD_birth_2_days[!is.na(output_table$BCKD_birth_2_days)])]
      reads_table = as.matrix(t(reads_table_16s_BCKD_days_12))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table = as.data.frame(reads_table$otu.tab.rff)
      
      Shannon_NB <- as.numeric(as.character(diversity(reads_table)))
      Evenness_NB <- Shannon_NB/log(specnumber(reads_table))
      Number_of_taxa_NB <- as.numeric(as.character(colSums(t(reads_table) != 0)))
      
      metadata_16s$Shannon_NB_days_12 = NA; metadata_16s$Evenness_NB_days_12 = NA;metadata_16s$Number_of_taxa_NB_days_12 = NA;
      metadata_16s$NB_t_SNE1_days_12 = NA;metadata_16s$NB_t_SNE2_days_12 = NA;
      for (a in 1: ncol(reads_table_16s_BCKD_days_12)) {
        n = which(metadata_16s$SampleID == colnames(reads_table_16s_BCKD_days_12)[a])
        metadata_16s$Shannon_NB_days_12[n] = Shannon_NB[a]
        metadata_16s$Evenness_NB_days_12[n] = Evenness_NB[a]
        metadata_16s$Number_of_taxa_NB_days_12[n] = Number_of_taxa_NB[a]
        
        m = which(pic_3$SampleID == colnames(reads_table_16s_BCKD_days_12)[a])
        metadata_16s$NB_t_SNE1_days_12[n] = pic_3$X1[m]
        metadata_16s$NB_t_SNE2_days_12[n] = pic_3$X2[m]
      }
      
      # differential abundance of NB_t_SNE1_days_12
      {
        reads_table = reads_table_16s_BCKD_days_12[,colnames(reads_table_16s_BCKD_days_12) %in% 
                                                   metadata_16s$SampleID[!is.na(metadata_16s$NB_t_SNE1_days_12)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X1
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        write.csv(output$data, 'NB_t_SNE1_day_0.csv')
        output$p
        ggsave('NB_t_SNE1_days_12_difference.pdf', width = 5, height = 2.3)
      }
      
      # differential abundance of NB_t_SNE2_days_12
      {
        reads_table = reads_table_16s_BCKD_days_12[,colnames(reads_table_16s_BCKD_days_12) %in% 
                                                     metadata_16s$SampleID[!is.na(metadata_16s$NB_t_SNE2_days_12)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X2
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        output$p
        ggsave('NB_t_SNE2_days_12_difference.pdf', width = 4, height = 1.8)
        ggsave('NB_t_SNE2_days_12_difference_2.pdf', width =4, height = 4)
      }
    }
    
  }
  # NR
  {
    # day 0
    {
      reads_table_16s_BRCD_day_0 = reads_table_16s_BRCD[,colnames(reads_table_16s_BRCD) %in% output_table$BRCD_birth_0_day[!is.na(output_table$BRCD_birth_0_day)]]
      reads_table = as.matrix(t(reads_table_16s_BRCD_day_0))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table = as.data.frame(reads_table$otu.tab.rff)
      
      Shannon_NR <- as.numeric(as.character(diversity(reads_table)))
      Evenness_NR <- Shannon_NR/log(specnumber(reads_table))
      Number_of_taxa_NR <- as.numeric(as.character(colSums(t(reads_table) != 0)))
      
      metadata_16s$Shannon_NR_day_0 = NA; metadata_16s$Evenness_NR_day_0 = NA;metadata_16s$Number_of_taxa_NR_day_0 = NA;
      metadata_16s$NR_t_SNE1_day_0 = NA;metadata_16s$NR_t_SNE2_day_0 = NA;
      for (a in 1: ncol(reads_table_16s_BRCD_day_0)) {
        n = which(metadata_16s$SampleID == colnames(reads_table_16s_BRCD_day_0)[a])
        metadata_16s$Shannon_NR_day_0[n] = Shannon_NR[a]
        metadata_16s$Evenness_NR_day_0[n] = Evenness_NR[a]
        metadata_16s$Number_of_taxa_NR_day_0[n] = Number_of_taxa_NR[a]
        
        m = which(pic_3$SampleID == colnames(reads_table_16s_BRCD_day_0)[a])
        metadata_16s$NR_t_SNE1_day_0[n] = pic_3$X1[m]
        metadata_16s$NR_t_SNE2_day_0[n] = pic_3$X2[m]
      }
      
      # differential abundance of NR_t_SNE1_day_0
      {
        reads_table = reads_table_16s_BRCD_day_0[,colnames(reads_table_16s_BRCD_day_0) %in% 
                                                   metadata_16s$SampleID[!is.na(metadata_16s$NR_t_SNE1_day_0)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X1
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        write.csv(output$data, 'NR_t_SNE1_day_0.csv')
        
      }
      # differential abundance of NR_t_SNE2_day_0
      {
        reads_table = reads_table_16s_BRCD_day_0[,colnames(reads_table_16s_BRCD_day_0) %in% 
                                                   metadata_16s$SampleID[!is.na(metadata_16s$NR_t_SNE1_day_0)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X2
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        output$p
        ggsave('NR_t_SNE2_day_0.pdf', width = 5, height = 2.5)
      }  
    }
    # days 1 and 2
    {
      reads_table_16s_BRCD_days_12 = reads_table_16s_BRCD[,colnames(reads_table_16s_BRCD) %in% 
                                                            c(output_table$BRCD_birth_1_days[!is.na(output_table$BRCD_birth_1_days)],output_table$BRCD_birth_2_days[!is.na(output_table$BRCD_birth_2_days)])]
      reads_table = as.matrix(t(reads_table_16s_BRCD_days_12))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table = as.data.frame(reads_table$otu.tab.rff)
      
      Shannon_NR <- as.numeric(as.character(diversity(reads_table)))
      Evenness_NR <- Shannon_NR/log(specnumber(reads_table))
      Number_of_taxa_NR <- as.numeric(as.character(colSums(t(reads_table) != 0)))
      
      metadata_16s$Shannon_NR_days_12 = NA; metadata_16s$Evenness_NR_days_12 = NA;metadata_16s$Number_of_taxa_NR_days_12 = NA;
      metadata_16s$NR_t_SNE1_days_12 = NA;metadata_16s$NR_t_SNE2_days_12 = NA;
      for (a in 1: ncol(reads_table_16s_BRCD_days_12)) {
        n = which(metadata_16s$SampleID == colnames(reads_table_16s_BRCD_days_12)[a])
        metadata_16s$Shannon_NR_days_12[n] = Shannon_NR[a]
        metadata_16s$Evenness_NR_days_12[n] = Evenness_NR[a]
        metadata_16s$Number_of_taxa_NR_days_12[n] = Number_of_taxa_NR[a]
        
        m = which(pic_3$SampleID == colnames(reads_table_16s_BRCD_days_12)[a])
        metadata_16s$NR_t_SNE1_days_12[n] = pic_3$X1[m]
        metadata_16s$NR_t_SNE2_days_12[n] = pic_3$X2[m]
      }
      
      # differential abundance of NR_t_SNE1_days_12
      {
        reads_table = reads_table_16s_BRCD_days_12[,colnames(reads_table_16s_BRCD_days_12) %in% 
                                                     metadata_16s$SampleID[!is.na(metadata_16s$NR_t_SNE1_days_12)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X1
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        output$p
        ggsave('NR_t_SNE1_days_12_difference.pdf', width = 5, height = 2.5)
      }
      
      # differential abundance of NR_t_SNE2_days_12
      {
        reads_table = reads_table_16s_BRCD_days_12[,colnames(reads_table_16s_BRCD_days_12) %in% 
                                                     metadata_16s$SampleID[!is.na(metadata_16s$NR_t_SNE1_days_12)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X2
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        output$p
        ggsave('NR_t_SNE2_days_12_difference.pdf', width = 5, height = 1.5)
      }
    }
    
  }
  # NS
  {
    # days 1 and 2
    {
      reads_table_16s_BS1D_days_12 = reads_table_16s_BS1D[,colnames(reads_table_16s_BS1D) %in% 
                                                            c(output_table$BS1D_birth_1_days[!is.na(output_table$BS1D_birth_1_days)],output_table$BS1D_birth_2_days[!is.na(output_table$BS1D_birth_2_days)])]
      reads_table = as.matrix(t(reads_table_16s_BS1D_days_12))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table = as.data.frame(reads_table$otu.tab.rff)
      
      Shannon_NS <- as.numeric(as.character(diversity(reads_table)))
      Evenness_NS <- Shannon_NS/log(specnumber(reads_table))
      Number_of_taxa_NS <- as.numeric(as.character(colSums(t(reads_table) != 0)))
      
      metadata_16s$Shannon_NS_days_12 = NA; metadata_16s$Evenness_NS_days_12 = NA;metadata_16s$Number_of_taxa_NS_days_12 = NA;
      metadata_16s$NS_t_SNE1_days_12 = NA;metadata_16s$NS_t_SNE2_days_12 = NA;
      for (a in 1: ncol(reads_table_16s_BS1D_days_12)) {
        n = which(metadata_16s$SampleID == colnames(reads_table_16s_BS1D_days_12)[a])
        metadata_16s$Shannon_NS_days_12[n] = Shannon_NS[a]
        metadata_16s$Evenness_NS_days_12[n] = Evenness_NS[a]
        metadata_16s$Number_of_taxa_NS_days_12[n] = Number_of_taxa_NS[a]
        
        m = which(pic_3$SampleID == colnames(reads_table_16s_BS1D_days_12)[a])
        metadata_16s$NS_t_SNE1_days_12[n] = pic_3$X1[m]
        metadata_16s$NS_t_SNE2_days_12[n] = pic_3$X2[m]
      }
      
      # differential abundance of NS_t_SNE1_days_12
      {
        reads_table = reads_table_16s_BS1D_days_12[,colnames(reads_table_16s_BS1D_days_12) %in% 
                                                     metadata_16s$SampleID[!is.na(metadata_16s$NS_t_SNE1_days_12)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X1
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        output$p
        ggsave('NS_t_SNE1_days_12_difference.pdf', width = 5, height = 2.5)
      }
      
      # differential abundance of NS_t_SNE2_days_12
      {
        reads_table = reads_table_16s_BS1D_days_12[,colnames(reads_table_16s_BS1D_days_12) %in% 
                                                     metadata_16s$SampleID[!is.na(metadata_16s$NS_t_SNE1_days_12)]]
        
        metadata = pic_3[pic_3$SampleID %in% colnames(reads_table),]
        reads_table = reads_table[metadata$SampleID]
        metadata = metadata$X2
        keep = metadata > mean(metadata); metadata[keep] = 'High';metadata[!keep] = 'Low'
        
        output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, 
                               paired_test = F, order_reverse = F, style = 1, order = c('Low','High'))
        output$p
        ggsave('NS_t_SNE2_days_12_difference.pdf', width = 5, height = 2)
      }
    }
    
  }
  
}

x = as.data.frame(metadata_16s[,c(1:5,(ncol(metadata_16s) -50):(ncol(metadata_16s)))])

{
  metadata_instruction = read.csv('/Users/binzhu/Desktop/Mom-baby/metadata_instruction copy 2.csv')
  metadata_16s = metadata_16s[,colnames(metadata_16s) %in% metadata_instruction$Metadata_subject_.Unless_otherwise_stated._listed_subjects_are_collected_from_mothers.]
  participant_list = unique(metadata_16s$ParticipantID)
  
  for (a in 1: length(participant_list)) {
    n = which(metadata_16s$ParticipantID == participant_list[a])
    metadata = metadata_16s[n,]
    for (b in 1: ncol(metadata)) {
      if (length(unique(metadata[,b])[!is.na(unique(metadata[,b]))]) == 1) {
        metadata[,b] = unique(metadata[,b])[!is.na(unique(metadata[,b]))]
      }
    }
    metadata_16s[n,] = metadata
  }
}



##### network analysis of the maternal factors ######
# prepare data
{
  metadata = metadata_16s[metadata_16s$Mombaby == 'Mom',]
  metadata = metadata[,!str_detect(colnames(metadata), 'baby|NB_|_NB|NR_|_NR|NS_|_NS')]
  metadata = metadata[order(metadata$VisitNum, decreasing = T),]
  metadata = metadata[!duplicated(metadata$ParticipantID),]
  metadata = metadata[match(output_table$participant,metadata$ParticipantID),]
  
  # cytokine data
  {
    reads_table = reads_table_cytokines
    reads_table = as.data.frame(t(reads_table))
    reads_table$ParticipantID = NA
    reads_table$ParticipantID = sapply(1:nrow(reads_table), function(j) (reads_table$ParticipantID[j] = metadata_16s$ParticipantID[which(metadata_16s$SampleID == row.names(reads_table)[j])]))
    reads_table = reads_table[match(metadata$ParticipantID, reads_table$ParticipantID),]
    reads_table$ParticipantID = NULL
    metadata = cbind(metadata, reads_table)
  }
  
  keep = sapply(1:ncol(metadata), function(j) (length(unique(metadata[,j])) >= 2 )); sum(keep)
  metadata = metadata[,keep]
  
  metadata = as.data.frame(metadata)
  metadata$VisitNum = NULL
  
  keep = colSums(is.na(metadata)) <= (nrow(metadata) / 2); sum(keep)
  metadata = metadata[,keep]
  
  for (a in 1: ncol(metadata)) {
    factor_level = as.data.frame(table((metadata[,a])))
    factor_level = factor_level[order(factor_level$Freq, decreasing = T),]
    
    if ((nrow(factor_level) <= 6 & factor_level[2,2] <=5) | nrow(factor_level) <2){
      metadata[,a] = NA
    }
  }
  
  keep = sapply(1:ncol(metadata), function(j) (length(unique(metadata[,j])[!is.na(unique(metadata[,j]))]) >= 2 )); sum(keep)
  metadata = metadata[,keep]
  write.csv(colnames(metadata),'involved_factors.csv', row.names = F)
  
  table(metadata$prenatal_care_start)
  metadata$prenatal_care_start[metadata$prenatal_care_start == 1] = 'A' 
  metadata$prenatal_care_start[metadata$prenatal_care_start == 2] = 'B' 
  metadata$prenatal_care_start[metadata$prenatal_care_start == 3] = 'C' 
  metadata$prenatal_care_start[metadata$prenatal_care_start == 4] = 'D' 
  
  metadata = metadata[,!str_detect(colnames(metadata), 'ParticipantID|SampleID|KitType|SampleType|BodySite|Flag')]
  
  metadata$changes_since_last_time = NULL
  # summary of the data
  {
    x = list()
    for(a in 1: ncol(metadata)) {
      if (length(unique(metadata[,a])) <=7) {
        data = table(metadata[,a])
        x[[a]] = data
      } else {
        data = quantile(metadata[!is.na(metadata[,a]),a])
        x[[a]] = data
      }
    }
    
    library(erer)
    names(x) = colnames(metadata)
    write.list(x, 'metadata_summary.csv')
  }
  
}

# correlation
{
  keep = sapply(1:ncol(metadata), function(j) (is.numeric(metadata[,j]))); sum(keep)
  metadata_1 = metadata[,keep]; metadata_2 = metadata[,!keep];  # metadata_1 is interval
  
  keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[!is.na(metadata_2[,j]),j])) ==2)); sum(keep)
  metadata_3 = metadata_2[,keep]; metadata_2 = metadata_2[,!keep];  # metadata_2 is ordinal  # metadata_3 is binary
  
  # The spearman's correlation for interval and ordinal
  {
    data =  t(cbind(metadata_1, metadata_2)); 
    
    output = newwork_rcorr(data, normalization_method = NA, type = 'spearman', pvalue = 0.05, 
                           cor_parameter= 0, style = 1, bar_max = 1, bar_min = -1, 
                           pheatmap_fontsize = 5, treeheight = 20, alpha = 0.05,FDR = T)
    
    correlated_variables_network_1 = data.frame(p_value = output$pvalue, adj_pvalue = output$adj_pvalue, R_value = output$Rvalue)
    correlated_variables_network_1 = correlated_variables_network_1[,c(3,5,4,1,2)]
    colnames(correlated_variables_network_1) = c('Variable 1','Variable 2', 'R-value','P-value','Adj-P-value')
    write.csv(correlated_variables_network_1,'correlation_interval_and_ordinal.csv',row.names = F)
  }
  
  
  # The Point-Biserial Correlation between binary and interval/ordinal
  {
    correlated_variables_network_2 = as.data.frame(matrix(data = NA, ncol = 5, nrow = 100000))
    colnames(correlated_variables_network_2) = c('Variable 1','Variable 2', 'R-value','P-value','Adj-P-value')
    
    data =  as.data.frame(cbind(metadata_1, metadata_2))
    n=0
    for (a in 1: ncol(data)) {
      for (b in 1: ncol(metadata_3)) {
        
        data_2 = data.frame(V1 = xtfrm(data[,a]), V2 = xtfrm(metadata_3[,b]))
        data_2 = data_2[!is.na(data_2$V1) & !is.na(data_2$V2),]
        
        if (nrow(data_2) < nrow(data) /10) {next}
        
        n=n+1
        correlated_variables_network_2$`Variable 1`[n] = colnames(data)[a]
        correlated_variables_network_2$`Variable 2`[n] = colnames(metadata_3)[b]
        
        x = cor.test(data_2$V1, data_2$V2)
        correlated_variables_network_2$`R-value`[n] = as.numeric(as.character(x$estimate))
        correlated_variables_network_2$`P-value`[n] = as.numeric(as.character(x$p.value))
      }
    }
    correlated_variables_network_2 = correlated_variables_network_2[!is.na(correlated_variables_network_2$`Variable 1`),]
    correlated_variables_network_2$`Adj-P-value` = fdr(correlated_variables_network_2$`P-value`)
    write.csv(correlated_variables_network_2,'correlation_binary_and_interval_or_ordinal.csv',row.names = F)
  }
  
  # The Phi Coefficient among binary
  {
    output = newwork_chi_squared(t(metadata_3), style = 1,pheatmap_fontsize = 5, treeheight = 50, alpha = 0.05,FDR = T)
    
    correlated_variables_network_3 = output$pvalue
    colnames(correlated_variables_network_3) = c('Variable 1','Variable 2', 'R-value','P-value','Adj-P-value')
    write.csv(correlated_variables_network_3,'correlation_binary.csv',row.names = F)
  }
  
  colnames(correlated_variables_network_3) = colnames(correlated_variables_network_1)
  correlated_variables_network = rbind(correlated_variables_network_1,correlated_variables_network_2,correlated_variables_network_3)
  correlated_variables_network = correlated_variables_network[correlated_variables_network$`Adj-P-value`<=0.05 & !is.na(correlated_variables_network$`Adj-P-value`),]
  correlated_variables_network = data.frame(Source = correlated_variables_network$`Variable 1`,
                                            Weight = abs(correlated_variables_network$`R-value`),
                                            Target = correlated_variables_network$`Variable 2`,
                                            Correlation= correlated_variables_network$`R-value`)
  
  write.csv(correlated_variables_network,'correlated_variables_network.csv', row.names = F)
  
  
}

# gephi modules
{
  modele_list = read.csv('gephi_node.csv')
  modele_weight = modele_list[,c(2,8)]
  modele_weight$betweenesscentrality[modele_weight$betweenesscentrality == 0] = 0.1
  
  modele_list = data.frame(Model = modele_list$modularity_class, Metadata = modele_list$Id)
  modele_list = modele_list[!str_detect(modele_list$Metadata, 'NB_|_NB|NR_|_NR|NS_|_NS|baby'),]
  
  modele_list_2 = data.frame(Model = max(modele_list$Model)+1, Metadata = colnames(metadata)[!str_remove_all(colnames(metadata),' ') %in% modele_list$Metadata])
  modele_list = rbind(modele_list,modele_list_2)
  #  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  #  myColor <- colorRampPalette(c("white", "red"))(50)
  
  #  metadata = as.data.frame(scale(metadata))
  module_sample_list = as.data.frame(matrix(data = NA, ncol = length(unique(modele_list$Model)), nrow = nrow(metadata)))
  colnames(module_sample_list) = paste0('Module_',c(0:(length(unique(modele_list$Model))-1)))
  
  library(StatMatch)
  # The normalization method for Ward.D clustering is to scale the data so that all the features have the same variance and mean.
  for (a in 0: (length(unique(modele_list$Model))-1)) {
    metadata_2 = metadata[,str_remove_all(colnames(metadata),' ') %in% modele_list$Metadata[modele_list$Model == a]]
    
    if (is.numeric(metadata_2[!is.na(metadata_2)])) {
      metadata_2[,1:ncol(metadata_2)] <- lapply(metadata_2[,1:ncol(metadata_2)],as.numeric)
      row.names(metadata_2) = c(1:nrow(metadata_2))
      
      for (b in 1: ncol(metadata_2)) {
        keep = is.nan(metadata_2[,b])
        metadata_2[keep,b] = NA
      }
      
      data = metadata_2
      for (b in 1: ncol(metadata_2)) {
        data[,b] = !is.na(metadata_2[,b])
      }
      
      keep = colSums(data == F & !is.na(data)) <= nrow(metadata_2) /1.1; sum(keep)
      metadata_2 = metadata_2[,keep]
      
      keep = rowSums(data == F & !is.na(data)) <= ncol(metadata_2) /1.1; sum(keep)
      metadata_2 = metadata_2[keep,]
      
      metadata_2 = as.data.frame(scale(metadata_2))
      
      # add weight
      
      
      #    for (b in 1: ncol(metadata_2)) {
      #      metadata_2[,b][!is.na(metadata_2[,b])] = range01(metadata_2[,b][!is.na(metadata_2[,b])])
      #    }
      
      x = pheatmap(t(metadata_2), show_rownames=T, show_colnames=F, fontsize = 6,
                   treeheight_row = 20, treeheight_col = 20,
                   clustering_method = "ward.D", clustering_distance_rows = 'minkowski', 
                   clustering_distance_cols = 'minkowski')
      
      x_cluster = sort(cutree(x$tree_col, k=2))
      keep = x_cluster == 1
      x_cluster[keep] = 'A'; x_cluster[!keep] = 'B'
      x_cluster = as.data.frame(x_cluster)
      colnames(x_cluster) = 'Cluster'
      x_cluster$X = 'X'
      
      x_cluster = x_cluster[match(row.names(metadata_2), row.names(x_cluster)),];x_cluster$X = NULL
      
      for (c in 1: nrow(module_sample_list)) {
        y = as.character(x_cluster$Cluster[which(row.names(x_cluster) == c)])
        if (length(y) == 0) {
          module_sample_list[c,a+1] = NA
        } else {module_sample_list[c,a+1] = y}
      }
      
      x_cluster$Cluster = factor(x_cluster$Cluster, levels = c('B','A'))
      x = pheatmap(t(metadata_2), show_rownames=T, show_colnames=F, 
                   fontsize = 6,treeheight_row = 20, treeheight_col = 20,
                   clustering_method = "ward.D", 
                   clustering_distance_rows = 'minkowski', clustering_distance_cols = 'minkowski',
                   annotation_col = x_cluster)
      
      save_heatmap_pdf(x, paste0('Module_',a,'.pdf'), width=6, height=0.55+0.09*ncol(metadata_2))
    } else {
      row.names(metadata_2) = c(1:nrow(metadata_2))
      
      for (b in 1: ncol(metadata_2)) {
        keep = is.nan(metadata_2[,b])
        metadata_2[keep,b] = NA
      }
      
      data = metadata_2
      for (b in 1: ncol(metadata_2)) {
        data[,b] = !is.na(metadata_2[,b])
      }
      
      keep = colSums(data == F & !is.na(data)) <= nrow(metadata_2) /1.1; sum(keep)
      metadata_2 = metadata_2[,keep]
      
      keep = rowSums(data == F & !is.na(data)) <= ncol(metadata_2) /1.1; sum(keep)
      metadata_2 = metadata_2[keep,]
      
      if (is.numeric(metadata_2[!is.na(metadata_2)])) {
        d_dist<-dist(metadata_2, method = "minkowski")
        d_dist[is.nan(d_dist) | is.na(d_dist)] = 1
        d_dist = hclust(d_dist, method="ward.D")
      } else {
        d_dist<-gower.dist(metadata_2, rngs=NULL, KR.corr=TRUE, var.weights = NULL, robcb=NULL)
        row.names(d_dist) = row.names(metadata_2); colnames(d_dist) = row.names(metadata_2)
        d_dist[is.nan(d_dist) | is.na(d_dist)] = 1
        d_dist<- as.dist(d_dist)
        d_dist = hclust(d_dist, method="ward.D")
      }
      
      x_cluster = sort(cutree(d_dist, k=2))
      keep = x_cluster == 1
      x_cluster[keep] = 'A'; x_cluster[!keep] = 'B'
      x_cluster = as.data.frame(x_cluster)
      colnames(x_cluster) = 'Cluster'
      x_cluster$X = row.names(x_cluster)
      
      metadata_2 = metadata_2[match(row.names(x_cluster),row.names(metadata_2)),]
      
      data = x_cluster
      
      if (nrow(data) < nrow(metadata)) {
        data_2 = data.frame(Cluster = NA, X = c(1:nrow(metadata))[!c(1:nrow(metadata)) %in% data$X]) 
        data = rbind(data,data_2)
      }
      
      for (c in 1: nrow(module_sample_list)) {
        y = as.character(x_cluster$Cluster[which(row.names(x_cluster) == c)])
        if (length(y) == 0) {
          module_sample_list[c,a+1] = NA
        } else {module_sample_list[c,a+1] = y}
      }
      
      x_cluster$Cluster = factor(x_cluster$Cluster, levels = c('A','B'))
      x_cluster$X = NULL
      
      for (b in 1:ncol(metadata_2)) {
        if (!is.numeric(metadata_2[,b])) {
          metadata_2[,b] = xtfrm(metadata_2[,b])
        }
      }
      
      metadata_2 = as.data.frame(scale(metadata_2))
      # add weight
      #    for (b in 1: ncol(metadata_2)) {
      #      metadata_2[,b][!is.na(metadata_2[,b])] = range01(metadata_2[,b][!is.na(metadata_2[,b])])
      #    }
      
      x_cluster$Cluster = factor(x_cluster$Cluster, levels = c('B','A'))
      x = pheatmap(t(metadata_2), show_rownames=T, show_colnames=F, cluster_cols = F,
                   fontsize = 6,treeheight_row = 20, treeheight_col = 20,
                   clustering_method = "ward.D", 
                   clustering_distance_rows = 'minkowski', 
                   annotation_col = x_cluster)
      
      save_heatmap_pdf(x, paste0('Module_',a,'.pdf'), width=6, height=0.55+0.09*ncol(metadata_2))
    }
    

  }
  
  write.csv(module_sample_list,'module_sample_list.csv',row.names = T)
}

# correlation among modules
{
  module_correlation = as.data.frame(matrix(data = NA, ncol = 5, nrow = 1000))
  colnames(module_correlation) = c('Module_1st','Module_2nd','X_squared','P_value','Adj_P_value')
  n=0
  for (a in 1:(ncol(module_sample_list)-1)) {
    for (b in (a+1) : ncol(module_sample_list)) {
      keep = !is.na(module_sample_list[,a]) & !is.na(module_sample_list[,b])
      x = chisq.test(module_sample_list[keep,a], module_sample_list[keep,b])
      
      n = n+1
      module_correlation$Module_1st[n] = colnames(module_sample_list)[a];module_correlation$Module_2nd[n] = colnames(module_sample_list)[b];
      module_correlation$X_squared[n] = as.numeric(x$statistic);module_correlation$P_value[n] = as.numeric(x$p.value);
    }
  }
  module_correlation = module_correlation[!is.na(module_correlation$Module_1st),]
  
  x <- adjust.p(module_correlation$P_value, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  x <- x$adjp
  module_correlation$Adj_P_value <- x$adjusted.p
  module_correlation = module_correlation[order(module_correlation$P_value),]
  
  write.csv(module_correlation,'module_correlation.csv',row.names = F)
}

###### association #####
# prepare metadata
{
  metadata_baby = metadata_16s[metadata_16s$Mombaby == 'Baby',]
  metadata_baby = metadata_baby[,str_detect(colnames(metadata_baby), 'baby|VisitNum|ParticipantID|NB_|_NB|NR_|_NR|NS_|_NS')]
  metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
  metadata_baby = metadata_baby[!duplicated(metadata_baby$ParticipantID),]
  metadata_baby = metadata_baby[match(output_table$participant, metadata_baby$ParticipantID),]
  
  sum(metadata$ParticipantID != metadata_baby$ParticipantID)

  metadata_baby = metadata_baby[,!str_detect(colnames(metadata_baby), 'ParticipantID|SampleID|KitType|SampleType|
                                             BodySite|Flag|VisitNum|Mombaby')]
  
  metadata_baby_nicu = metadata_baby$baby_nicu
  
  keep = sapply(1:ncol(metadata_baby), function(j) (is.numeric(metadata_baby[,j]))) ; sum(keep)
  metadata_baby_alpha = metadata_baby[,keep]
}

# test association with alpha
{
  reads_table1 = as.matrix(t(module_sample_list))
  reads_table2 = as.matrix(t(metadata_baby_alpha))
  
  correlation_data_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = 0))
  colnames(correlation_data_all) = c('Mom_module','Baby_alpha','Pvalue','Adj_P_value','Test')
  
  for (a in 1: nrow(reads_table2)) {
    correlation_data = as.data.frame(matrix(data = NA, ncol = 5, nrow = nrow(reads_table1)))
    colnames(correlation_data) = c('Mom_module','Baby_alpha','Pvalue','Adj_P_value','Test')
    correlation_data$Baby_alpha = row.names(reads_table2)[a]
    correlation_data$Test = a
    
    for (b in 1: nrow(reads_table1)) {
      correlation_data$Mom_module[b] = row.names(reads_table1)[b]
      
      data1 = as.numeric(reads_table2[a,])
      data2 = (reads_table1[b,])
      
      keep = (!is.na(data1)) & (!is.na(data2))
      data1 = data1[keep]
      data2 = data2[keep]
      data = data.frame(V1 = data1, V2 = data2)
      data$V2 = as.factor(data$V2)
      
      data3 = wilcox.test(V1~V2,data)
      correlation_data$Pvalue[b] = data3$p.value
      
    }
    
    data4 <- adjust.p(correlation_data$Pvalue, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
    data4 <- data4$adjp
    correlation_data$Adj_P_value <- data4$adjusted.p
    
    
    for (b in 1: nrow(reads_table1)) {
      correlation_data$Mom_module[b] = row.names(reads_table1)[b]
      
      data1 = as.numeric(reads_table2[a,])
      data2 = (reads_table1[b,])
      
      keep = (!is.na(data1)) & (!is.na(data2))
      data1 = data1[keep]
      data2 = data2[keep]
      data = data.frame(V1 = data1, V2 = data2)
      data$V2 = as.factor(data$V2)

      if (correlation_data$Pvalue[b] <= 0.05 & correlation_data$Adj_P_value[b] > 0.05) {
        ggplot(data=data, aes(x=V2, y=V1)) +geom_violin(trim=T,aes(fill=V2))+
          geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
          geom_jitter(size = 0.1)+theme_bw()+
          xlab(paste0('Module ',b-1,' cluster')) +
          ylab(correlation_data$Baby_alpha[b])+
          ggtitle(paste0('Case number: ','\n',sum(data$V2 == 'A'), ' (A); ',sum(data$V2 == 'B'),' (B)','\n',
                         'P = ', format(round(correlation_data$Pvalue[b], 3), nsmall = 3)))+theme_bw()+
          theme(axis.title = element_text(size = 6), 
                axis.text = element_text(size = 6), 
                plot.title = element_text(size = 6),
                legend.position = "none")
        ggsave(paste0('alpha_No_FDR_',correlation_data$Mom_module[b],'_',correlation_data$Baby_alpha[b],'.pdf'),width=1, height=1.8)
      } else if (correlation_data$Adj_P_value[b] <= 0.05) {
        ggplot(data=data, aes(x=V2, y=V1)) +geom_violin(trim=T,aes(fill=V2))+
          geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
          geom_jitter(size = 0.1)+theme_bw()+
          xlab(paste0('Module ',b-1,' cluster')) +
          ylab(correlation_data$Baby_alpha[b])+
          ggtitle(paste0('Case number: ','\n',sum(data$V2 == 'A'), ' (A); ',sum(data$V2 == 'B'),' (B)','\n',
                         'P = ', '\n',
                         'FDR = '))+theme_bw()+
          theme(axis.title = element_text(size = 6), 
                axis.text = element_text(size = 6), 
                plot.title = element_text(size = 6),
                legend.position = "none")
        ggsave(paste0('alpha__',correlation_data$Mom_module[b],'_',correlation_data$Baby_alpha[b],'.pdf'),width=1, height=1.8)
      }
    }
    
    
    correlation_data_all = rbind(correlation_data_all,correlation_data)
  }
  correlation_data_all = correlation_data_all[order(correlation_data_all$Adj_P_value),]
  write.csv(correlation_data_all,'correlation_alpha.csv')
  

}

# test association with NICU by Odds Ratio and 95% Confidence Interval
{
  library(epitools)
  reads_table1 = as.data.frame(t(module_sample_list))
  NICU = metadata_baby_nicu; NICU[NICU == 1] = 'No';NICU[NICU == '10'] = 'Yes'
  
  keep = !is.na(NICU); reads_table1 = reads_table1[,keep]; NICU = NICU[keep]
  
  correlation_data = as.data.frame(matrix(data = NA, ncol = 8, nrow = nrow(reads_table1)))
  colnames(correlation_data) = c('Mom_module','OR','Lower','Upper',
                                 'Pvalue_fisher_exact','Adj_P_value','NICU_in_cluster_A','NICU_in_cluster_B')
  
  for (b in 1: nrow(reads_table1)) {
    correlation_data$Mom_module[b] = row.names(reads_table1)[b]
    
    data1 = NICU; data2 = (reads_table1[b,])
    
    keep = (!is.na(data1)) & (!is.na(data2));data1 = data1[keep];data2 = data2[keep]
    data = data.frame(V1 = data1, V2 = data2)
    x= table(data$V2,data$V1); x[x ==0] = 0.6
    data3 = oddsratio(x)
    
    correlation_data$OR[b] = data3$measure[2,1]
    correlation_data$Lower[b] = data3$measure[2,2]; correlation_data$Upper[b] = data3$measure[2,3]; 
    
    correlation_data$Pvalue_fisher_exact[b] = data3$p.value[2,2]
    correlation_data$NICU_in_cluster_A[b] = paste0(sum(data$V2 == 'A'), ' (', sum(data$V2 == 'A' & data$V1 == 'Yes'), ')')
    correlation_data$NICU_in_cluster_B[b] = paste0(sum(data$V2 == 'B'), ' (', sum(data$V2 == 'B' & data$V1 == 'Yes'), ')')
  }
  x <- adjust.p(correlation_data$Pvalue_fisher_exact, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  x <- x$adjp
  correlation_data$Adj_P_value <- x$adjusted.p
  
  correlation_data = correlation_data[order(correlation_data$Pvalue),]
  write.csv(correlation_data,'correlation_NICU.csv')
  
  correlation_data$Mom_module = factor(correlation_data$Mom_module, levels = correlation_data$Mom_module) 
  correlation_data$Color = 'black'
  correlation_data$Color[correlation_data$Pvalue_fisher_exact <= 0.05] = 'red'
  
  ggplot(data = correlation_data, aes(x = OR, y = Mom_module)) + 
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), linewidth = 0.5, linetype = 1, height = 0.3,
                   colour = correlation_data$Color) +
    geom_point(colour = correlation_data$Color, size = 1) + xlim(-2,8)+
    geom_vline(aes(xintercept = 1), linetype = 2) + scale_x_continuous(breaks = seq(-1, 15, 2))+
    theme_bw() +
    theme(legend.position = "none")
  ggsave('OR_NICU_mom_modules.pdf',width = 8,height = 3)
  
  ggplot(data = correlation_data, aes(x = OR, y = Mom_module)) + 
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), linewidth = 0.5, linetype = 1, height = 0.3,
                   colour = correlation_data$Color) +
    geom_point(colour = correlation_data$Color, size = 1) + 
    geom_vline(aes(xintercept = 1), linetype = 2) + scale_x_continuous(breaks = seq(-1, 69, 20))+
    theme_bw() +
    theme(legend.position = "none")
  ggsave('OR_NICU_mom_modules_2.pdf',width = 2.0,height = 3)
}

# test association with composition of the microbiome
{
  output_1 = as.data.frame(matrix(data= NA, ncol = 5, nrow = ncol(module_sample_list)))
  colnames(output_1) = c('NB day 0','NB days 12','NR day 0','NR days 12','NS days 12')
  row.names(output_1) = paste0('Module_',c(0:(ncol(module_sample_list)-1)))
  
  output_2 = as.data.frame(matrix(data= NA, ncol = 5, nrow = ncol(module_sample_list)))
  colnames(output_2) = c('NB day 0','NB days 12','NR day 0','NR days 12','NS days 12')
  row.names(output_2) = paste0('Module_',c(0:(ncol(module_sample_list)-1)))
  
  # NB day 0
  {
    reads_table1 = as.data.frame(module_sample_list)
    reads_table1 = reads_table1[!is.na(output_table$BCKD_birth_0_day),]
    dim(reads_table1)
    
    reads_table2 = reads_table_16s_BCKD_day_0
    reads_table2 = reads_table2[output_table$BCKD_birth_0_day[!is.na(output_table$BCKD_birth_0_day)]]
    reads_table2 = as.data.frame(t(reads_table2))
    reads_table2 = Rarefy(reads_table2, depth = min(rowSums(reads_table2)))
    reads_table2 <- reads_table2$otu.tab.rff
    reads_table2 <- as.data.frame(reads_table2)
    
    pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(reads_table1)))
    row.names(pvalue_all) = colnames(reads_table1)
    
    for (a in 1:(nrow(pvalue_all))) {
      metadata = reads_table1[,a]
      keep = !is.na(metadata)
      metadata = as.data.frame(metadata[keep])
      reads_table = reads_table2[keep,]
      
      pvalue <- adonis2(reads_table ~ ., data = metadata, 
                        method = "bray", parallel = 8)
      pvalue_all[a,] = pvalue[1,]
    }
    colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
    output_1[,1] = pvalue_all$P_value
    
    pvalue_all = pvalue_all[order(pvalue_all$P_value, decreasing = F),]
    pvalue_all$Adj_P_value = NA
    x <- adjust.p(pvalue_all$P_value, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
    x <- x$adjp
    pvalue_all$Adj_P_value <- x$adjusted.p
    
    reads_table1 = reads_table1[row.names(pvalue_all)]
    
    metadata = as.data.frame(reads_table1)
    metadata[is.na(metadata)] = 'C'
    pvalue <- adonis2(reads_table2 ~ ., data = metadata,method = "bray", parallel = 8)
    x = as.data.frame(pvalue); x = x[c(1:(nrow(x) -2)),]; x = x[match(row.names(output_2),row.names(x)),]
    output_2[,1] = x$`Pr(>F)`
    
    # differential abundance analysis
    {
      metadata = reads_table1[,colnames(reads_table1) == 'Module_3']
      keep = !is.na(metadata)
      metadata = metadata[keep]
      reads_table = reads_table2[keep,]
      reads_table = as.data.frame(t(reads_table))
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, 
                             fold_change_th = 1, paired_test = F, 
                             order_reverse = F, style = 1, order = c('A','B'))
      x = output$data
      x = x[,c(-12,-13)]
      x = x[order(x$we.ep,decreasing = F),]
      write.csv(x,'abunance_change_BCKD_birth_0_day_Module_3.csv')
      
      metadata = reads_table1[,colnames(reads_table1) == 'Module_5']
      keep = !is.na(metadata)
      metadata = metadata[keep]
      reads_table = reads_table2[keep,]
      reads_table = as.data.frame(t(reads_table))
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, 
                             fold_change_th = 1, paired_test = F, 
                             order_reverse = F, style = 1, order = c('A','B'))
      x = output$data
      x = x[,c(-12,-13)]
      x = x[order(x$we.ep,decreasing = F),]
      write.csv(x,'abunance_change_BCKD_birth_0_day_Module_5.csv')
      
      metadata = reads_table1[,colnames(reads_table1) == 'Module_11']
      keep = !is.na(metadata)
      metadata = metadata[keep]
      reads_table = reads_table2[keep,]
      reads_table = as.data.frame(t(reads_table))
      
      output = dif_abundance(reads_table,metadata, pvalue_th = 0.05, 
                             fold_change_th = 1, paired_test = F, 
                             order_reverse = F, style = 1, order = c('A','B'))
      x = output$data
      x = x[,c(-12,-13)]
      x = x[order(x$we.ep,decreasing = F),]
      write.csv(x,'abunance_change_BCKD_birth_0_day_Module_11.csv')
    }
    
  }
  
  # NB day 1 2
  {
    reads_table1 = as.data.frame(module_sample_list)
    reads_table1 = reads_table1[!is.na(output_table$BCKD_birth_1_days) | !is.na(output_table$BCKD_birth_2_days),]
    dim(reads_table1)
    sample_list = vector(length= nrow(output_table))
    for (a in 1: length(sample_list)) {
      if (!is.na(output_table$BCKD_birth_1_days[a])) {
        sample_list[a] = output_table$BCKD_birth_1_days[a]
      } else if (!is.na(output_table$BCKD_birth_2_days[a])) {
        sample_list[a] = output_table$BCKD_birth_2_days[a]
      } else {
        sample_list[a] = NA
      }
    }
    sample_list = sample_list[!is.na(sample_list)]
    
    reads_table2 = reads_table_16s_BCKD_days_12
    reads_table2 = reads_table2[,colnames(reads_table2) %in% sample_list]
    reads_table2 = reads_table2[sample_list]
    reads_table2 = as.data.frame(t(reads_table2))
    reads_table2 = Rarefy(reads_table2, depth = min(rowSums(reads_table2)))
    reads_table2 <- reads_table2$otu.tab.rff
    reads_table2 <- as.data.frame(reads_table2)
    
    pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(reads_table1)))
    row.names(pvalue_all) = colnames(reads_table1)
    
    for (a in 1:(nrow(pvalue_all))) {
      metadata = reads_table1[,a]
      keep = !is.na(metadata)
      metadata = as.data.frame(metadata[keep])
      reads_table = reads_table2[keep,]
      
      pvalue <- adonis2(reads_table ~ ., data = metadata, 
                        method = "bray", parallel = 8)
      pvalue_all[a,] = pvalue[1,]
    }
    colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
    output_1[,2] = pvalue_all$P_value
    
    pvalue_all = pvalue_all[order(pvalue_all$P_value, decreasing = F),]

    
    reads_table1 = reads_table1[row.names(pvalue_all)]
    
    metadata = as.data.frame(reads_table1)
    metadata[is.na(metadata)] = 'C'
    pvalue <- adonis2(reads_table2 ~ ., data = metadata,method = "bray", parallel = 8)
    x = as.data.frame(pvalue); x = x[c(1:(nrow(x) -2)),]; x = x[match(row.names(output_2),row.names(x)),]
    output_2[,2] = x$`Pr(>F)`
    
  }
  
  # NR day 0
  {
    reads_table1 = as.data.frame(module_sample_list)
    reads_table1 = reads_table1[!is.na(output_table$BRCD_birth_0_day),]
    dim(reads_table1)
    
    reads_table2 = reads_table_16s_BRCD_day_0
    reads_table2 = reads_table2[output_table$BRCD_birth_0_day[!is.na(output_table$BRCD_birth_0_day)]]
    reads_table2 = as.data.frame(t(reads_table2))
    reads_table2 = Rarefy(reads_table2, depth = min(rowSums(reads_table2)))
    reads_table2 <- reads_table2$otu.tab.rff
    reads_table2 <- as.data.frame(reads_table2)
    
    pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(reads_table1)))
    row.names(pvalue_all) = colnames(reads_table1)
    
    for (a in 1:(nrow(pvalue_all))) {
      metadata = reads_table1[,a]
      keep = !is.na(metadata)
      metadata = as.data.frame(metadata[keep])
      reads_table = reads_table2[keep,]
      
      pvalue <- adonis2(reads_table ~ ., data = metadata, 
                        method = "bray", parallel = 8)
      pvalue_all[a,] = pvalue[1,]
    }
    colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
    output_1[,3] = pvalue_all$P_value
    
    reads_table1 = reads_table1[row.names(pvalue_all)]
    
    metadata = as.data.frame(reads_table1)
    metadata[is.na(metadata)] = 'C'
    pvalue <- adonis2(reads_table2 ~ ., data = metadata,method = "bray", parallel = 8)
    x = as.data.frame(pvalue); x = x[c(1:(nrow(x) -2)),]; x = x[match(row.names(output_2),row.names(x)),]
    output_2[,3] = x$`Pr(>F)`
    
    
  }
  
  # NR day 1 2
  {
    reads_table1 = as.data.frame(module_sample_list)
    reads_table1 = reads_table1[!is.na(output_table$BRCD_birth_1_days) | !is.na(output_table$BRCD_birth_2_days),]
    dim(reads_table1)
    sample_list = vector(length= nrow(output_table))
    for (a in 1: length(sample_list)) {
      if (!is.na(output_table$BRCD_birth_1_days[a])) {
        sample_list[a] = output_table$BRCD_birth_1_days[a]
      } else if (!is.na(output_table$BRCD_birth_2_days[a])) {
        sample_list[a] = output_table$BRCD_birth_2_days[a]
      } else {
        sample_list[a] = NA
      }
    }
    sample_list = sample_list[!is.na(sample_list)]
    
    reads_table2 = reads_table_16s_BRCD_days_12
    reads_table2 = reads_table2[,colnames(reads_table2) %in% sample_list]
    reads_table2 = reads_table2[sample_list]
    reads_table2 = as.data.frame(t(reads_table2))
    reads_table2 = Rarefy(reads_table2, depth = min(rowSums(reads_table2)))
    reads_table2 <- reads_table2$otu.tab.rff
    reads_table2 <- as.data.frame(reads_table2)
    
    pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(reads_table1)))
    row.names(pvalue_all) = colnames(reads_table1)
    
    for (a in 1:(nrow(pvalue_all))) {
      metadata = reads_table1[,a]
      keep = !is.na(metadata)
      metadata = as.data.frame(metadata[keep])
      reads_table = reads_table2[keep,]
      
      pvalue <- adonis2(reads_table ~ ., data = metadata, 
                        method = "bray", parallel = 8)
      pvalue_all[a,] = pvalue[1,]
    }
    colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
    output_1[,4] = pvalue_all$P_value
    
    
    reads_table1 = reads_table1[row.names(pvalue_all)]
    
    metadata = as.data.frame(reads_table1)
    metadata[is.na(metadata)] = 'C'
    pvalue <- adonis2(reads_table2 ~ ., data = metadata,method = "bray", parallel = 8)
    x = as.data.frame(pvalue); x = x[c(1:(nrow(x) -2)),]; x = x[match(row.names(output_2),row.names(x)),]
    output_2[,4] = x$`Pr(>F)`
    
  }
  
  # NS day 1 2
  {
    reads_table1 = as.data.frame(module_sample_list)
    reads_table1 = reads_table1[!is.na(output_table$BS1D_birth_1_days) | !is.na(output_table$BS1D_birth_2_days),]
    dim(reads_table1)
    sample_list = vector(length= nrow(output_table))
    for (a in 1: length(sample_list)) {
      if (!is.na(output_table$BS1D_birth_1_days[a])) {
        sample_list[a] = output_table$BS1D_birth_1_days[a]
      } else if (!is.na(output_table$BS1D_birth_2_days[a])) {
        sample_list[a] = output_table$BS1D_birth_2_days[a]
      } else {
        sample_list[a] = NA
      }
    }
    sample_list = sample_list[!is.na(sample_list)]
    
    reads_table2 = reads_table_16s_BS1D_days_12
    reads_table2 = reads_table2[,colnames(reads_table2) %in% sample_list]
    reads_table2 = reads_table2[sample_list]
    reads_table2 = as.data.frame(t(reads_table2))
    reads_table2 = Rarefy(reads_table2, depth = min(rowSums(reads_table2)))
    reads_table2 <- reads_table2$otu.tab.rff
    reads_table2 <- as.data.frame(reads_table2)
    
    pvalue_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = ncol(reads_table1)))
    row.names(pvalue_all) = colnames(reads_table1)
    
    for (a in 1:(nrow(pvalue_all))) {
      metadata = reads_table1[,a]
      keep = !is.na(metadata)
      metadata = as.data.frame(metadata[keep])
      reads_table = reads_table2[keep,]
      
      pvalue <- adonis2(reads_table ~ ., data = metadata, 
                        method = "bray", parallel = 8)
      pvalue_all[a,] = pvalue[1,]
    }
    colnames(pvalue_all) = c('Df','SumOfSqs','R2','F','P_value')
    output_1[,5] = pvalue_all$P_value
    
    pvalue_all = pvalue_all[order(pvalue_all$P_value, decreasing = F),]
    
    
    reads_table1 = reads_table1[row.names(pvalue_all)]
    
    metadata = as.data.frame(reads_table1)
    metadata[is.na(metadata)] = 'C'
    pvalue <- adonis2(reads_table2 ~ ., data = metadata,method = "bray", parallel = 8)
    x = as.data.frame(pvalue); x = x[c(1:(nrow(x) -2)),]; x = x[match(row.names(output_2),row.names(x)),]
    output_2[,5] = x$`Pr(>F)`
  }
  
  myColor <- colorRampPalette(c("#FF8F8F","white"))(50)

  x = output_1; x[x >0.05] = 1;
  p = pheatmap(as.matrix(x),cluster_rows = F,cluster_cols = F,color=myColor,
               display_numbers = as.matrix(output_1))
  save_heatmap_pdf(p, 'Beta_p.pdf', width=3, height=3.5)
  
  x = output_2; x[x >0.05] = 1;
  p = pheatmap(as.matrix(x),cluster_rows = F,cluster_cols = F,color=myColor,
               display_numbers = as.matrix(output_2))
  save_heatmap_pdf(p, 'Beta_p_multiple_variables.pdf', width=3, height=3.5)

}
