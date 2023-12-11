setwd('/Users/binzhu/Desktop/Mom-baby/results/script_6')
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
  
  data$Microbiome = factor(data$Microbiome, levels = c('Metadata','MB','MR','MV','Lipid','Cytokine'))
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
  
  p[7] = paste0('MV Lipid: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MV'],data$weeks_pregnant[data$Microbiome == 'Lipid'])$p.value)
  p[8] = paste0('MB Lipid: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MB'],data$weeks_pregnant[data$Microbiome == 'Lipid'])$p.value)
  p[9] = paste0('MR Lipid: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MR'],data$weeks_pregnant[data$Microbiome == 'Lipid'])$p.value)
  p[10] = paste0('Metadata Lipid: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'Metadata'],data$weeks_pregnant[data$Microbiome == 'Lipid'])$p.value)

  p[11] = paste0('Metadata Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'Metadata'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)
  p[12] = paste0('MB Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MB'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)
  p[13] = paste0('MV Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MV'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)
  p[14] = paste0('MR Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'MR'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)
  p[15] = paste0('Lipid Cytokine: ', wilcox.test(data$weeks_pregnant[data$Microbiome == 'Lipid'],data$weeks_pregnant[data$Microbiome == 'Cytokine'])$p.value)
  
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
  }
  # MB
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



##### ML ######
# prepare input maternal data
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
  
  # one-hot conversion
  metadata[metadata == 'Yes'] = 10; metadata[metadata == 'No'] = 1; 
  
  metadata[metadata == 'A'] = 1;metadata[metadata == 'B'] = 2;metadata[metadata == 'C'] = 3;
  metadata[metadata == 'D'] = 4;metadata[metadata == 'E'] = 5;metadata[metadata == 'F'] = 6;
  
  for (b in 1 : ncol(metadata)) {
    if (length(unique(metadata[!is.na(metadata[,b]),b])) <= 6 & length(unique(metadata[!is.na(metadata[,b]),b])) > 2) {
      data = as.numeric(as.character(metadata[,b]))
      keep = median(data[!is.na(data)])
      keep_1 = data > keep & !is.na(data); keep_2 = data <= keep & !is.na(data); 
      metadata[keep_1,b] = 10; metadata[keep_2,b] = 1; 
    }
  }
  
  data_all = metadata
  
  for (a in 1:ncol(data_all)) {
    data = data_all[,a]
    if (!is.numeric(data[!is.na(data)])) {
      data_all[,a] = xtfrm(data)
    }
  }
  
  data_all = as.data.frame(scale(data_all))
  
  # filter low quality metadata
  {
    # remove input with less than 2 levels
    keep = sapply(1:ncol(data_all), function(j) (length(unique(data_all[,j])) >= 2 ))
    data_all = data_all[,keep]
    
    keep = colSums(is.na(data_all)) <= nrow(data_all)/4
    sum(keep)
    data_all = data_all[,keep]
    
    keep = sapply(1:ncol(data_all), function(j) (length(unique(data_all[,j][!is.na(data_all[,j])])) >= 2 ))
    sum(keep)
    data_all = data_all[,keep]
  }
}

# prepare baby independent variables
{
  metadata_baby = metadata_16s[metadata_16s$Mombaby == 'Baby',]
  metadata_baby = metadata_baby[,str_detect(colnames(metadata_baby), 'baby|VisitNum|ParticipantID|NB_|_NB|NR_|_NR|NS_|_NS|t_SNE')]
  metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
  metadata_baby = metadata_baby[!duplicated(metadata_baby$ParticipantID),]
  metadata_baby = metadata_baby[match(output_table$participant, metadata_baby$ParticipantID),]
  
  sum(metadata$ParticipantID != metadata_baby$ParticipantID)
  
  keep = str_detect(colnames(metadata_baby),'NR_|NB_|NS_|nicu'); sum(keep)
  metadata_baby = metadata_baby[,keep]
  
  metadata_baby_nicu = metadata_baby$baby_nicu
  metadata_baby_alpha = metadata_baby[,-1]
  for (a in 1:ncol(metadata_baby_alpha)) {
    data = metadata_baby_alpha[,a] 
    keep = data > median(data[!is.na(data)])
    metadata_baby_alpha[keep & !is.na(keep),a] = 'Yes'
    metadata_baby_alpha[!keep & !is.na(keep),a] = 'No'
  }
  metadata_baby = cbind(metadata_baby_nicu,metadata_baby_alpha)
}

# modeling
{
  library(pROC)
  library(doSNOW)
  library(parallel)
  library(matrixStats)
  library(erer) # export list
  library(randomForest)
  
  importance_ML = list()
  myPvs = vector()
  myPvs2 = vector()
  myPvs3 = vector()
  myPvs4 = vector()
  
  output_all = as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  colnames(output_all) = c('Model_name','AUC_ROC_curve','p_value','Case_number')
  
  for (a in 1: ncol(metadata_baby)) {

    keep = !is.na(metadata_baby[,a]); sum(keep)
    output = metadata_baby[keep,a]
    output = as.factor(output)
    
    input = data_all[keep,]
    
    input = input[,colnames(input) != 'ga_at_delivery' & colnames(input) !='Preterm']
    
    set.seed(2023+123*a)

    input = FilterFeatures(input, output); colnames(input)
    
    if (sum(is.na(input)) >0) {
      input = apply(input, 2, as.numeric)
      output = as.factor(output)
      data = cbind(input,output)
      data_2 <- rfImpute(output ~ ., data)    # Impute missing values in predictor data using proximity from randomForest.
      input = data_2[,-1]
      input = as.matrix(input)
    } 
    
    data = xxx_get_result(input, output, ntree= ncol(input)) 
    
    data$p_roc
    ggsave(paste0('ML_',colnames(metadata_baby)[a],'.pdf'),width=3, height=3)
    
    myPvs[a]=data$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
    myPvs2[a]=data$pvalue   # test correlation between expected values and true values (PTB yes or no).
    myPvs3[a] = data$err_median
    importance_ML[[a]] = data$importance_ML
    myPvs4[a] = paste0(length(output),' (',sum(output == 'Yes'),')')
    myPvs[a]
    
    model_name = colnames(metadata_baby)[a]
    
    data = cbind(model_name,myPvs[a],myPvs2[a], myPvs4[a])
    data = as.data.frame(data)
    colnames(data) = c('Model_name','AUC_ROC_curve','p_value','Case_number')
    
    output_all = rbind(output_all,data)

  }
  
  write.csv(output_all,'ML_output_all.csv')
#  names(importance_ML) = colnames(metadata_baby)
  library(erer)
  write.list(importance_ML, 'importance_ML.csv')
}

# Test association among risk factors of NICU
{
  data = importance_ML[[1]]
  data = data$name[data$MeanDecreaseGini>0]
  
  data = data_all[,colnames(data_all) %in% data]
  data = as.data.frame(t(data))

  {
    reads_table = data; normalization_method = NA; type = 'spearman'; pvalue = 0.05; cor_parameter= 0; 
    style = 1; bar_max = 2; bar_min = -2; pheatmap_fontsize = 5; treeheight = 10; alpha = 0.05;
    FDR = T; ABS = F
    #  reads_table = reads_table_2
    if (is.na(normalization_method)) {
      reads_table = as.matrix(t(reads_table))
    } else if (normalization_method == 'clr') {
      reads_table = as.data.frame(t(reads_table))
      reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
      
    } else if (normalization_method == 'rarefy') {
      reads_table = sweep(reads_table,2,colSums(reads_table),"/")
      reads_table = as.matrix(t(reads_table))
    }
    
    # convert ordinal data to numeric as introduced in 'cor' function instruction
    for (a in 1:ncol(reads_table)) {
      if (!is.numeric(reads_table[,a])) {
        reads_table[,a] = xtfrm(reads_table[,a])
      }
    }
    
    otu.cor <- rcorr(reads_table, type= type)
    otu.pval <- forceSymmetric(otu.cor$P)
    otu.pval <- otu.pval@x
    pvalue_output = otu.pval
    
    if (FDR == T) {
      otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      otu.pval <- otu.pval$adjp
      otu.pval <- otu.pval$adjusted.p
    }
    
    p.yes <- otu.pval< pvalue  
    
    r.val = otu.cor$r # select all the correlation values 
    r.val_output = gather(as.data.frame(r.val)); r.val_output$key2 = rep(colnames(r.val), ncol(r.val))
    
    p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
    
    p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
    p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.
    
    p.yes.rr[is.na(p.yes.rr)] = 0
    
    keep = abs(colSums(p.yes.rr)) > 0
    p.yes.rr = p.yes.rr[,keep]
    p.yes.rr = p.yes.rr[keep,]
  }
  
  pdf("NICU_risk_factor_correlation_1.pdf", width=8, height=6)
  corrplot(p.yes.rr, type="upper", order="hclust", insig = "blank",
           tl.col = "black",na.label = " ", tl.srt = 45)
  dev.off()

  data = importance_ML[[1]]
  data = data$name[data$MeanDecreaseGini>0]
  
  data = data_all[,colnames(data_all) %in% data]
  data = data[,str_detect(colnames(data), 'SNE') | !str_detect(colnames(data), 'MB_|MR_|MV_')]
  data = as.data.frame(t(data))
  
  output = newwork_rcorr(data, normalization_method = NA, type = 'spearman', pvalue = 0.05, cor_parameter= 0, 
                         style = 1, bar_max = 2, bar_min = -2, pheatmap_fontsize = 5, treeheight = 20, alpha = 0.05,
                         FDR = T, ABS = T)
  
  save_heatmap_pdf(output$p, 'NICU_risk_factor_correlation_2.pdf', width=3, height=2.8)

}

# test p-values by NICU groups Mann-Whitely
{
  keep = !is.na(metadata_baby_nicu); sum(keep)
  data = data_all[keep,]
  metadata_baby_nicu_2 = metadata_baby_nicu[keep]
  
  output = data.frame(variable = colnames(data), P_value = NA, Adj_P_value = NA, Mean_value_high_in_NICU = NA)
  for (a in 1:ncol(data)) {
    data_2 = data.frame(NICU = as.factor(metadata_baby_nicu_2), variable = data[,a])
    output$P_value[a] = wilcox.test(variable~NICU,data_2)$p.value
    
    if (mean(data_2$variable[data_2$NICU == 'Yes' & !is.na(data_2$variable)]) > mean(data_2$variable[data_2$NICU == 'No'& !is.na(data_2$variable)])) {
      output$Mean_value_high_in_NICU[a] = 'True'
    } else if (mean(data_2$variable[data_2$NICU == 'Yes'& !is.na(data_2$variable)]) < mean(data_2$variable[data_2$NICU == 'No'& !is.na(data_2$variable)])) {
      output$Mean_value_high_in_NICU[a] = 'False'
    } else {
      output$Mean_value_high_in_NICU[a] = NA
    }
    
    if (output$P_value[a] <=0.05) {
      ggplot(data_2, aes(x=NICU, y=variable)) + 
        geom_boxplot(aes(fill=NICU) , color="black", outlier.shape=NA, width=0.1) +
        geom_jitter(size = 0.1)+ theme_bw()+ 
        labs(y = 'NICU', x = colnames(data)[a])+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))+ coord_flip()
      ggsave(paste0('NICU_', colnames(data)[a], '.pdf'), width = 3,height = 1.5)
    }
    
  }
  x = adjust.p(output$P_value, alpha=0.05, pz=0.05)
  x = x$adjp
  output$Adj_P_value = x$adjusted.p
  
  write.csv(output,'NICU_wilcox.csv',row.names = F)
}

# case match for NICU
{
  metadata = metadata_16s[metadata_16s$Mombaby == 'Mom' & !is.na(metadata_16s$Shannon_MB) & !is.na(metadata_16s$baby_nicu),]
  metadata = metadata[order(metadata$VisitNum, decreasing = T),]
  metadata = metadata[!duplicated(metadata$ParticipantID),]
  
  metadata_case = metadata[metadata$baby_nicu == "10" & !is.na(metadata$baby_nicu),]
  metadata_control = metadata_case; metadata_control[1:nrow(metadata_control),1:ncol(metadata_control)] = NA
  
  metadata_2 = metadata[metadata$baby_nicu == "1" & !is.na(metadata$baby_nicu),]
  
  for (a in 1: nrow(metadata_case)) {
    
    n =1
    keep = metadata_2$date_of_birth >= metadata_case$date_of_birth[a] -n & 
      metadata_2$date_of_birth <= metadata_case$date_of_birth[a] +n & 
      metadata_2$SampleID != metadata_case$SampleID[a] & 
      metadata_2$african_american == metadata_case$african_american[a] & 
      metadata_2$caucasian == metadata_case$caucasian[a] & 
      metadata_2$education == metadata_case$education[a] & 
      metadata_2$delivery_csection == metadata_case$delivery_csection[a] & 
      metadata_2$weeks_pregnant >= metadata_case$weeks_pregnant[a] -n & 
      metadata_2$weeks_pregnant <= metadata_case$weeks_pregnant[a] +n; sum(keep[!is.na(keep)])
    while (sum(keep[!is.na(keep)]) == 0 & n < 40) {
      n = n+1
      keep = metadata_2$date_of_birth >= metadata_case$date_of_birth[a] -n & 
        metadata_2$date_of_birth <= metadata_case$date_of_birth[a] +n & 
        metadata_2$SampleID != metadata_case$SampleID[a] & 
        metadata_2$african_american == metadata_case$african_american[a] & 
        metadata_2$caucasian == metadata_case$caucasian[a] & 
        metadata_2$education == metadata_case$education[a] & 
        metadata_2$delivery_csection == metadata_case$delivery_csection[a] & 
        metadata_2$weeks_pregnant >= metadata_case$weeks_pregnant[a] -n & 
        metadata_2$weeks_pregnant <= metadata_case$weeks_pregnant[a] +n; sum(keep[!is.na(keep)])
    }

    if (n < 40) {
      x = metadata_2[which(keep),]
      x = x[order(x$ga_at_delivery, decreasing = F),]
      metadata_control[a,] = x[1,]
      metadata_2 = metadata_2[-which(keep)[1],]
    }
    
  }
  
  sum(is.na(metadata_control$ParticipantID))
  
  keep = !is.na(metadata_control$SampleID); sum(keep)
  metadata_control = metadata_control[keep,]; metadata_case = metadata_case[keep,]
  
  data = data.frame(Shannon_MB = c(metadata_case$Shannon_MB,metadata_control$Shannon_MB), 
                    Preterm = c(metadata_case$Preterm,metadata_control$Preterm), 
                    NICU = c(rep('Yes',length(metadata_case$Shannon_MB)), rep('No',length(metadata_control$Shannon_MB))))
  
  table(data$Preterm,data$NICU)
  
  
  
  ggplot(data=data, aes(x=NICU, y=Shannon_MB)) +geom_violin(trim=T, aes(fill=NICU))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1)+theme_bw()+
    theme_classic()+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  ggsave('NICU_Shannon_MB.pdf',width=2, height=2)
  
  wilcox.test(Shannon_MB ~ NICU,data = data)
  
  
  # Two-Way ANOVA
  data$NICU[data$NICU == 'Yes'] =1;data$NICU[data$NICU == 'No'] =0
  res.aov2 <- aov(NICU ~ Preterm + Shannon_MB + Shannon_MB:Preterm, data = data)
  summary(res.aov2)
  

  data_1 = metadata_case[,colnames(metadata_case) %in% c('date_of_birth','african_american','caucasian','education',
                                     'delivery_csection','weeks_pregnant')]
  data_2 = metadata_control[,colnames(metadata_control) %in% c('date_of_birth','african_american','caucasian','education',
                                                         'delivery_csection','weeks_pregnant')]
  output = as.data.frame(matrix(data = NA, ncol = 4, nrow = 6))
  colnames(output) = c('Variable','Data_in_NICU; mean (sd)','Data_in_control; mean (sd)','P_value')
  output$Variable = colnames(data_1)[1:6]
  
  for (a in 1:6) {
    if (length(unique(data_1[,a])) > 2) {
      x = wilcox.test(data_1[,a], data_2[,a], paired = T)
      output$P_value[a] = paste0('Mann-Whitely U test P-value = ', x$p.value)
      output$`Data_in_NICU; mean (sd)`[a] = paste0(mean(data_1[,a]),' (', sd(data_1[,a]),')')
      output$`Data_in_control; mean (sd)`[a] = paste0(mean(data_2[,a]),' (', sd(data_2[,a]),')')
      
    } else {
      output$`Data_in_NICU; mean (sd)`[a] = sum(data_1[,a] == 'Yes') / nrow(data_1)
      output$`Data_in_control; mean (sd)`[a] = sum(data_2[,a] == 'Yes') / nrow(data_2)
      output$P_value[a] = 'The same in matched cases'
    }
    
  }
  write.csv(output, 'NICU_case_match.csv')
}




