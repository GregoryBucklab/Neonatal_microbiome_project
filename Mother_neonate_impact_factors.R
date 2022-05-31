##### correlation of mom's and baby's widely distributed taxa, paired mom-baby 16s samples (time closest between mom and baby) in a participantID #######
### find core taxa  (if not exclude other taxa that not widely distributed in these six microbiomes, 
# those uniqe taxa in some certain microbiome will impact the abundance at normalization and affect correlation analysis (no correlation will get with taxa not widely exist in these 6 microbiome))

### run this part (1676-1748) for 10 times and find taxa widely distributed in the 6 microbiomes
# randomly select 68 microbiomes from NR, RS, NB, MB, MR, MV microbiomes.
{
  # randomly select 68 microbiomes from NR, RS, NB, MB, MR, MV microbiomes.
  x = sort(sample(colnames((reads_table_16s_BCKD)), 68))
  x = c(x,sort(sample(colnames((reads_table_16s_BRCD)), 68)))
  x = c(x,sort(sample(colnames((reads_table_16s_BS1D)), 68)))
  x = c(x,sort(sample(colnames((reads_table_16s_MRCD)), 68)))
  x = c(x,sort(sample(colnames((reads_table_16s_MCKD)), 68)))
  x = c(x,sort(sample(colnames((reads_table_16s_MV1D)), 68)))
  
  keep = colnames(reads_table_16s) %in% x
  sum(keep)
  reads_table_16s_2 = reads_table_16s[,keep]
  metadata_16s_2 = metadata_16s[keep,]
  
  reads_table_16s_2 = prepare_reads_table(reads_table_16s_2, metadata_16s_2, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
  
  metadata_16s_2 = reads_table_16s_2$metadata
  reads_table_16s_2 = reads_table_16s_2$reads_table
  
  reads_table_16s_2_abundance = get_abundance_table(reads_table_16s_2, mc.cores = 8)
  
  # output widely distributed taxa and abundance
  reads_table = as.data.frame(matrix(data = NA, ncol = 6, nrow = nrow(reads_table_16s_2_abundance)))
  row.names(reads_table) = row.names(reads_table_16s_2_abundance)
  
  type_list = unique(metadata_16s_2$SampleType)
  
  for (a in 1: length(type_list)) {
    keep = metadata_16s_2$SampleType == type_list[a]
    reads_table_2 = reads_table_16s_2_abundance[,keep]
    reads_table[,a] = rowSums(reads_table_2) / ncol(reads_table_2)
    colnames(reads_table)[a] = type_list[a]
  }
  
  write.csv(reads_table,'Widely distributed taxa in 6 microbiomes.csv')
  
}

# the test of widely distributed taxa run 10 times, find the overlapped taxa in this 10 times
{
  data = read.csv('Widely distributed taxa in 6 microbiomes 10 times.csv', header = F)
  
  taxa_list = gather(data)
  taxa_list = unique(taxa_list$value)
  
  for (a in 1: length(taxa_list)) {
    n = 0
    for (b in 1:10) {
      if (taxa_list[a] %in% data[,b]) {
        n=n+1
      }
    }
    if (n < 5) {
      taxa_list[a] = NA
    }
  }
  
  taxa_list = taxa_list[!is.na(taxa_list)]
  taxa_list = taxa_list[taxa_list != ""]
  write.csv(taxa_list,'taxa_list_output.csv')
  # get reads table with only taxa in the taxa_list
  keep = row.names(reads_table_16s) %in% taxa_list
  sum(keep)
  reads_table_16s_2 = reads_table_16s[keep,]
  
  reads_table_16s_2 = prepare_reads_table(reads_table_16s_2, metadata_16s, total_reads_threshold = 5000, species_threshold = 0, mc.cores = 8)
  
  metadata_16s_2 = reads_table_16s_2$metadata
  reads_table_16s_2 = reads_table_16s_2$reads_table
  
}


### get correlation ###
{
  type1_all = c('BCKD','BRCD','BS1D')
  type2_all = c('MCKD','MRCD','MV1D')
  
  p = as.data.frame(matrix(data = 0, nrow = 18, ncol = 4))
  colnames(p) = c('Type','Number','Same_taxa','Number_samples')
  n=1
  
  for (a in 1: length(type1_all)) {
    for (b in 1: length(type2_all)) {
      
      type1 = type1_all[a]
      type2 = type2_all[b]
      
      # baby cheek and mom cheek
      keep = metadata_16s_2$SampleType == type1 | metadata_16s_2$SampleType == type2 
      metadata = metadata_16s_2[keep,]
      
      metadata_baby = metadata
      keep = metadata_baby$SampleType == type1
      metadata_baby = metadata_baby[keep,]
      
      metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
      metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]
      keep = duplicated(metadata_baby$ParticipantID)
      metadata_baby = metadata_baby[!keep,]
      
      p_number_1 = unique(metadata_baby$ParticipantID)
      length(p_number_1)
      keep = metadata$ParticipantID %in% p_number_1
      metadata = metadata[keep,]
      
      metadata_mom = metadata[metadata$SampleType == type2,]
      metadata_mom = metadata_mom[order(metadata_mom$VisitNum),]
      metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
      keep = duplicated(metadata_mom$ParticipantID)
      metadata_mom = metadata_mom[!keep,]
      
      keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
      metadata_baby = metadata_baby[keep,]
      
      reads_table_mom = reads_table_16s_2[metadata_mom$SampleID]
      reads_table_baby = reads_table_16s_2[metadata_baby$SampleID]
      
      sample_pair_name = cbind(colnames(reads_table_mom),colnames(reads_table_baby))
      write.csv(sample_pair_name,paste(type1, type2,'.csv'))
      
      row.names(reads_table_mom) = paste0('Mom_',row.names(reads_table_mom))
      row.names(reads_table_baby) = paste0('Baby_',row.names(reads_table_baby))
      colnames(reads_table_baby) = colnames(reads_table_mom) 
      
      if (ncol(reads_table_baby) == 0) {
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        p$Number[n*2-1] = 0
        p$Number[n*2] = 0
        
        n = n+1
        next
      }
      
      # normalize tables
      reads_table_mom = reads_table_mom[rowSums(reads_table_mom) > 0,]
      reads_table_baby = reads_table_baby[rowSums(reads_table_baby) > 0,]
      
      reads_table_mom = as.data.frame(t(reads_table_mom))
      reads_table_mom = reads_table_mom + 0.5
      reads_table_mom <- clr(reads_table_mom)      ### CLR normalization in rows
      reads_table_mom = as.data.frame(t(reads_table_mom))
      
      reads_table_baby = as.data.frame(t(reads_table_baby))
      reads_table_baby = reads_table_baby + 0.5
      reads_table_baby <- clr(reads_table_baby)      ### CLR normalization in rows
      reads_table_baby = as.data.frame(t(reads_table_baby))
      
      x=which(row.names(reads_table_baby) == 'Baby_Streptococcus_salivarius')
      if (b==1 & length(x) != 0 ) {
        reads_table_baby = reads_table_baby[-x,]
      }
      
      
      correlation = newwork_rcorr3(reads_table_mom, reads_table_baby, pvalue = 0.05, 
                                   cor_parameter= 0, style = 1, bar_max = 1, bar_min = -1, 
                                   pheatmap_fontsize = 5, alpha = 0.05,mc.cores =8)
      
      if (!is.na(correlation)) {
        write.csv(correlation$cor_matrix,paste0("network_",type1,'vs.',type2,'.csv'))
        
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        if (is.matrix(correlation)) {
          correlation = as.data.frame(as.matrix(correlation))
        } else {
          correlation = as.data.frame(correlation$cor_matrix)
        }
        
        c= gather(correlation)
        c$key2 = rep(row.names(correlation),ncol(correlation))
        c$taxa1 = str_replace_all(c$key, 'Baby_','')
        c$taxa2 = str_replace_all(c$key2, 'Mom_','')
        c$same =NA
        c$same[c$taxa1 == c$taxa2] = 'same'
        
        p$Number[n*2-1] = sum(c$taxa1 != c$taxa2 & c$value != 0)
        p$Number[n*2] = sum(c$taxa1 == c$taxa2 & c$value != 0)
        
        n = n+1
        
        rm(metadata_mom,metadata_baby,reads_table_mom,reads_table_baby)
      } else {
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        p$Number[n*2-1] = 0
        p$Number[n*2] = 0
        
        n = n+1
        rm(metadata_mom,metadata_baby,reads_table_mom,reads_table_baby)
      }
    }
  }
  
  p = p[p$Type != '0',]
  p$Same_taxa = factor(p$Same_taxa, levels= c('Y','N'))
  
  p2 = ggplot(data=p, aes(x=Type, y=Number, fill=Same_taxa)) +
    geom_bar(stat="identity")+  
    labs(x = NULL, y = "Number of taxa pairs\nwith significant correlation")+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
  p2
  ggsave(paste0('Correlation_baby_mom_1.pdf'),width=5, height=2.8)
  
  p2 = ggplot(data=p, aes(x=Type, y=Number, fill=Same_taxa)) +
    geom_bar(stat="identity")+  
    geom_text(data=p, 
              aes(Type, 0,  label= Number_samples), 
              position = position_dodge(width=0.9),
              size=2, angle=90, hjust=0, vjust = 0)+ 
    labs(x = NULL, y = "Number of taxa pairs\nwith significant correlation")+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
  p2
  ggsave(paste0('Correlation_baby_mom_2.pdf'),width=5, height=2.8)
  
  cor_number = as.data.frame(matrix(data = 0, nrow = 3, ncol = 3))
  colnames(cor_number) = type1_all
  row.names(cor_number) = type2_all
  
  n = 1
  for (a in 1: 3) {
    for (b in 1:3) {
      cor_number[b,a]= p$Number[n*2] + p$Number[n*2-1]
      n=n+1
    }
  }
  
  cor_number = as.matrix(cor_number)
  
  bar_max = 100
  bar_min = 0
  paletteLength <- 50
  myColor <- colorRampPalette(c("white", "red"))(paletteLength)
  myBreaks <- seq(0, bar_max,length.out=paletteLength)
  
  pdf("Correlation_baby_mom_3.pdf", width = 2, height = 1.6)
  p3 <- pheatmap(cor_number, color=myColor,breaks=myBreaks, fontsize = 10, 
                 display_numbers = cor_number, treeheight_row = 10, treeheight_col =10, 
                 cluster_rows = F,cluster_cols = F)
  dev.off()
}

# correlation for C section
{
  type1_all = c('BCKD','BRCD','BS1D')
  type2_all = c('MCKD','MRCD','MV1D')
  
  p = as.data.frame(matrix(data = 0, nrow = 18, ncol = 4))
  colnames(p) = c('Type','Number','Same_taxa','Number_samples')
  n=1
  
  for (a in 1: length(type1_all)) {
    for (b in 1: length(type2_all)) {
      
      type1 = type1_all[a]
      type2 = type2_all[b]
      
      # baby cheek and mom cheek
      keep = metadata_16s_2$SampleType == type1 | metadata_16s_2$SampleType == type2 
      metadata = metadata_16s_2[keep,]
      
      metadata_baby = metadata
      keep = metadata_baby$SampleType == type1 & metadata_baby$delivery_csection == "Yes" & 
        !is.na(metadata_baby$delivery_csection)
      metadata_baby = metadata_baby[keep,]
      
      metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
      metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]
      keep = duplicated(metadata_baby$ParticipantID)
      metadata_baby = metadata_baby[!keep,]
      
      p_number_1 = unique(metadata_baby$ParticipantID)
      length(p_number_1)
      keep = metadata$ParticipantID %in% p_number_1
      metadata = metadata[keep,]
      
      metadata_mom = metadata[metadata$SampleType == type2,]
      metadata_mom = metadata_mom[order(metadata_mom$VisitNum),]
      metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
      keep = duplicated(metadata_mom$ParticipantID)
      metadata_mom = metadata_mom[!keep,]
      
      keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
      metadata_baby = metadata_baby[keep,]
      
      reads_table_mom = reads_table_16s_2[metadata_mom$SampleID]
      reads_table_baby = reads_table_16s_2[metadata_baby$SampleID]
      
      sample_pair_name = cbind(colnames(reads_table_mom),colnames(reads_table_baby))
      write.csv(sample_pair_name,paste(type1, type2,'.csv'))
      
      row.names(reads_table_mom) = paste0('Mom_',row.names(reads_table_mom))
      row.names(reads_table_baby) = paste0('Baby_',row.names(reads_table_baby))
      colnames(reads_table_baby) = colnames(reads_table_mom) 
      
      if (ncol(reads_table_baby) == 0) {
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        p$Number[n*2-1] = 0
        p$Number[n*2] = 0
        
        n = n+1
        next
      }
      
      # normalize tables
      reads_table_mom = reads_table_mom[rowSums(reads_table_mom) > 0,]
      reads_table_baby = reads_table_baby[rowSums(reads_table_baby) > 0,]
      
      reads_table_mom = as.data.frame(t(reads_table_mom))
      reads_table_mom = reads_table_mom + 0.5
      reads_table_mom <- clr(reads_table_mom)      ### CLR normalization in rows
      reads_table_mom = as.data.frame(t(reads_table_mom))
      
      reads_table_baby = as.data.frame(t(reads_table_baby))
      reads_table_baby = reads_table_baby + 0.5
      reads_table_baby <- clr(reads_table_baby)      ### CLR normalization in rows
      reads_table_baby = as.data.frame(t(reads_table_baby))
      
      x=which(row.names(reads_table_baby) == 'Baby_Streptococcus_salivarius')
      if (b==1 & length(x) != 0 ) {
        reads_table_baby = reads_table_baby[-x,]
      }
      
      
      correlation = newwork_rcorr3(reads_table_mom, reads_table_baby, pvalue = 0.05, 
                                   cor_parameter= 0, style = 1, bar_max = 1, bar_min = -1, 
                                   pheatmap_fontsize = 5, alpha = 0.05,mc.cores =8)
      
      if (!is.na(correlation)) {
        write.csv(correlation$cor_matrix,paste0("network_",type1,'vs.',type2,'.csv'))
        
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        if (is.matrix(correlation)) {
          correlation = as.data.frame(as.matrix(correlation))
        } else {
          correlation = as.data.frame(correlation$cor_matrix)
        }
        
        c= gather(correlation)
        c$key2 = rep(row.names(correlation),ncol(correlation))
        c$taxa1 = str_replace_all(c$key, 'Baby_','')
        c$taxa2 = str_replace_all(c$key2, 'Mom_','')
        c$same =NA
        c$same[c$taxa1 == c$taxa2] = 'same'
        
        p$Number[n*2-1] = sum(c$taxa1 != c$taxa2 & c$value != 0)
        p$Number[n*2] = sum(c$taxa1 == c$taxa2 & c$value != 0)
        
        n = n+1
        
        rm(metadata_mom,metadata_baby,reads_table_mom,reads_table_baby)
      } else {
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        p$Number[n*2-1] = 0
        p$Number[n*2] = 0
        
        n = n+1
        rm(metadata_mom,metadata_baby,reads_table_mom,reads_table_baby)
      }
    }
  }
  
  p = p[p$Type != '0',]
  p$Same_taxa = factor(p$Same_taxa, levels= c('Y','N'))
  
  p2 = ggplot(data=p, aes(x=Type, y=Number, fill=Same_taxa)) +
    geom_bar(stat="identity")+  
    labs(x = NULL, y = "Number of taxa pairs\nwith significant correlation")+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
  p2
  ggsave(paste0('Correlation_baby_mom_1_csection.pdf'),width=5, height=2.8)
  
  p2 = ggplot(data=p, aes(x=Type, y=Number, fill=Same_taxa)) +
    geom_bar(stat="identity")+  
    geom_text(data=p, 
              aes(Type, 0,  label= Number_samples), 
              position = position_dodge(width=0.9),
              size=2, angle=90, hjust=0, vjust = 0)+ 
    labs(x = NULL, y = "Number of taxa pairs\nwith significant correlation")+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
  p2
  ggsave(paste0('Correlation_baby_mom_2_csection.pdf'),width=5, height=2.8)
  
  cor_number = as.data.frame(matrix(data = 0, nrow = 3, ncol = 3))
  colnames(cor_number) = type1_all
  row.names(cor_number) = type2_all
  
  n = 1
  for (a in 1: 3) {
    for (b in 1:3) {
      cor_number[b,a]= p$Number[n*2] + p$Number[n*2-1]
      n=n+1
    }
  }
  
  cor_number = as.matrix(cor_number)
  
  bar_max = 160
  bar_min = 0
  paletteLength <- 50
  myColor <- colorRampPalette(c("white", "red"))(paletteLength)
  myBreaks <- seq(0, bar_max,length.out=paletteLength)
  
  pdf("Correlation_baby_mom_3_csection.pdf", width = 2, height = 1.6)
  p3 <- pheatmap(cor_number, color=myColor,breaks=myBreaks, fontsize = 10, 
                 display_numbers = cor_number, treeheight_row = 10, treeheight_col =10, 
                 cluster_rows = F,cluster_cols = F)
  dev.off()
}

# correlation for vaginal delivery
{
  type1_all = c('BCKD','BRCD','BS1D')
  type2_all = c('MCKD','MRCD','MV1D')
  
  p = as.data.frame(matrix(data = 0, nrow = 18, ncol = 4))
  colnames(p) = c('Type','Number','Same_taxa','Number_samples')
  n=1
  
  for (a in 1: length(type1_all)) {
    for (b in 1: length(type2_all)) {
      
      type1 = type1_all[a]
      type2 = type2_all[b]
      
      # baby cheek and mom cheek
      keep = metadata_16s_2$SampleType == type1 | metadata_16s_2$SampleType == type2 
      metadata = metadata_16s_2[keep,]
      
      metadata_baby = metadata
      keep = metadata_baby$SampleType == type1 & metadata_baby$delivery_csection == "No" & 
        !is.na(metadata_baby$delivery_csection)
      metadata_baby = metadata_baby[keep,]
      
      metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
      metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]
      keep = duplicated(metadata_baby$ParticipantID)
      metadata_baby = metadata_baby[!keep,]
      
      p_number_1 = unique(metadata_baby$ParticipantID)
      length(p_number_1)
      keep = metadata$ParticipantID %in% p_number_1
      metadata = metadata[keep,]
      
      metadata_mom = metadata[metadata$SampleType == type2,]
      metadata_mom = metadata_mom[order(metadata_mom$VisitNum),]
      metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
      keep = duplicated(metadata_mom$ParticipantID)
      metadata_mom = metadata_mom[!keep,]
      
      keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
      metadata_baby = metadata_baby[keep,]
      
      reads_table_mom = reads_table_16s_2[metadata_mom$SampleID]
      reads_table_baby = reads_table_16s_2[metadata_baby$SampleID]
      
      sample_pair_name = cbind(colnames(reads_table_mom),colnames(reads_table_baby))
      write.csv(sample_pair_name,paste(type1, type2,'.csv'))
      
      row.names(reads_table_mom) = paste0('Mom_',row.names(reads_table_mom))
      row.names(reads_table_baby) = paste0('Baby_',row.names(reads_table_baby))
      colnames(reads_table_baby) = colnames(reads_table_mom) 
      
      if (ncol(reads_table_baby) == 0) {
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        p$Number[n*2-1] = 0
        p$Number[n*2] = 0
        
        n = n+1
        next
      }
      
      # normalize tables
      reads_table_mom = reads_table_mom[rowSums(reads_table_mom) > 0,]
      reads_table_baby = reads_table_baby[rowSums(reads_table_baby) > 0,]
      
      reads_table_mom = as.data.frame(t(reads_table_mom))
      reads_table_mom = reads_table_mom + 0.5
      reads_table_mom <- clr(reads_table_mom)      ### CLR normalization in rows
      reads_table_mom = as.data.frame(t(reads_table_mom))
      
      reads_table_baby = as.data.frame(t(reads_table_baby))
      reads_table_baby = reads_table_baby + 0.5
      reads_table_baby <- clr(reads_table_baby)      ### CLR normalization in rows
      reads_table_baby = as.data.frame(t(reads_table_baby))
      
      x=which(row.names(reads_table_baby) == 'Baby_Streptococcus_salivarius')
      if (b==1 & length(x) != 0 ) {
        reads_table_baby = reads_table_baby[-x,]
      }
      
      
      correlation = newwork_rcorr3(reads_table_mom, reads_table_baby, pvalue = 0.05, 
                                   cor_parameter= 0, style = 1, bar_max = 1, bar_min = -1, 
                                   pheatmap_fontsize = 5, alpha = 0.05,mc.cores =8)
      
      if (!is.na(correlation)) {
        write.csv(correlation$cor_matrix,paste0("network_",type1,'vs.',type2,'.csv'))
        
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        if (is.matrix(correlation)) {
          correlation = as.data.frame(as.matrix(correlation))
        } else {
          correlation = as.data.frame(correlation$cor_matrix)
        }
        
        c= gather(correlation)
        c$key2 = rep(row.names(correlation),ncol(correlation))
        c$taxa1 = str_replace_all(c$key, 'Baby_','')
        c$taxa2 = str_replace_all(c$key2, 'Mom_','')
        c$same =NA
        c$same[c$taxa1 == c$taxa2] = 'same'
        
        p$Number[n*2-1] = sum(c$taxa1 != c$taxa2 & c$value != 0)
        p$Number[n*2] = sum(c$taxa1 == c$taxa2 & c$value != 0)
        
        n = n+1
        
        rm(metadata_mom,metadata_baby,reads_table_mom,reads_table_baby)
      } else {
        p$Type[n*2-1]= paste0(type1,'vs.',type2)
        p$Type[n*2]= paste0(type1,'vs.',type2)
        
        p$Same_taxa[n*2-1]= 'N'
        p$Same_taxa[n*2]= 'Y'
        
        p$Number_samples[n*2-1] = ncol(reads_table_mom)
        p$Number_samples[n*2] = ncol(reads_table_mom)
        
        p$Number[n*2-1] = 0
        p$Number[n*2] = 0
        
        n = n+1
        rm(metadata_mom,metadata_baby,reads_table_mom,reads_table_baby)
      }
    }
  }
  
  p = p[p$Type != '0',]
  p$Same_taxa = factor(p$Same_taxa, levels= c('Y','N'))
  
  p2 = ggplot(data=p, aes(x=Type, y=Number, fill=Same_taxa)) +
    geom_bar(stat="identity")+  
    labs(x = NULL, y = "Number of taxa pairs\nwith significant correlation")+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
  p2
  ggsave(paste0('Correlation_baby_mom_1_vaginal.pdf'),width=5, height=2.8)
  
  p2 = ggplot(data=p, aes(x=Type, y=Number, fill=Same_taxa)) +
    geom_bar(stat="identity")+  
    geom_text(data=p, 
              aes(Type, 0,  label= Number_samples), 
              position = position_dodge(width=0.9),
              size=2, angle=90, hjust=0, vjust = 0)+ 
    labs(x = NULL, y = "Number of taxa pairs\nwith significant correlation")+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
  p2
  ggsave(paste0('Correlation_baby_mom_2_vaginal.pdf'),width=5, height=2.8)
  
  cor_number = as.data.frame(matrix(data = 0, nrow = 3, ncol = 3))
  colnames(cor_number) = type1_all
  row.names(cor_number) = type2_all
  
  n = 1
  for (a in 1: 3) {
    for (b in 1:3) {
      cor_number[b,a]= p$Number[n*2] + p$Number[n*2-1]
      n=n+1
    }
  }
  
  cor_number = as.matrix(cor_number)
  
  bar_max = 160
  bar_min = 0
  paletteLength <- 50
  myColor <- colorRampPalette(c("white", "red"))(paletteLength)
  myBreaks <- seq(0, bar_max,length.out=paletteLength)
  
  pdf("Correlation_baby_mom_3_vaginal.pdf", width = 2, height = 1.6)
  p3 <- pheatmap(cor_number, color=myColor,breaks=myBreaks, fontsize = 10, 
                 display_numbers = cor_number, treeheight_row = 10, treeheight_col =10, 
                 cluster_rows = F,cluster_cols = F)
  dev.off()
}

### get co-existence ###
{
  th = 1
  species_threshold = 0.001
  # get reads table with only taxa in the taxa_list
  keep = row.names(reads_table_16s) %in% taxa_list
  sum(keep)
  reads_table_16s_2 = reads_table_16s[keep,]
  reads_table_16s_2 = prepare_reads_table(reads_table_16s_2, metadata_16s, total_reads_threshold = 5000,
                                          species_threshold = species_threshold, mc.cores = 8)
  metadata_16s_2 = reads_table_16s_2$metadata
  reads_table_16s_2 = reads_table_16s_2$reads_table
  
  keep = metadata_16s_2$`delivery_vaginal ` == "Yes" & !is.na(metadata_16s_2$`delivery_vaginal `)
  metadata_16s_2 = metadata_16s_2[keep,]
  reads_table_16s_2 = reads_table_16s_2[,keep]
  
  # all days
  {
    type1_all = c('BCKD','BRCD','BS1D')
    type2_all = c('MCKD','MRCD','MV1D')
    
    p_pair = list()
    
    p_unpair = list()
    n=0
    
    for (a in 1: length(type1_all)) {
      for (b in 1: length(type2_all)) {
        n=n+1
        type1 = type1_all[a]
        type2 = type2_all[b]
        
        # baby cheek and mom cheek
        keep = metadata_16s_2$SampleType == type1 | metadata_16s_2$SampleType == type2 
        metadata = metadata_16s_2[keep,]
        
        metadata_baby = metadata
        keep = metadata_baby$SampleType == type1
        metadata_baby = metadata_baby[keep,]
        
        metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
        metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]
        keep = duplicated(metadata_baby$ParticipantID)
        metadata_baby = metadata_baby[!keep,]
        
        p_number_1 = unique(metadata_baby$ParticipantID)
        length(p_number_1)
        keep = metadata$ParticipantID %in% p_number_1
        metadata = metadata[keep,]
        
        metadata_mom = metadata[metadata$SampleType == type2,]
        metadata_mom = metadata_mom[order(metadata_mom$VisitNum),]
        metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
        keep = duplicated(metadata_mom$ParticipantID)
        metadata_mom = metadata_mom[!keep,]
        
        participant_list = intersect(metadata_mom$ParticipantID, metadata_baby$ParticipantID)
        
        keep = metadata_baby$ParticipantID %in% participant_list
        metadata_baby = metadata_baby[keep,]
        keep = metadata_mom$ParticipantID %in% participant_list
        metadata_mom = metadata_mom[keep,]
        
        reads_table_mom = reads_table_16s_2[metadata_mom$SampleID]
        reads_table_baby = reads_table_16s_2[metadata_baby$SampleID]
        
        reads_table_mom = reads_table_mom >= th
        reads_table_baby = reads_table_baby >= th
        colnames(reads_table_mom) = metadata_mom$ParticipantID
        colnames(reads_table_baby) = metadata_baby$ParticipantID
        
        co_exist=sapply(1:ncol(reads_table_mom), function(j) (reads_table_mom[,j]*reads_table_baby[,j]))
        co_exist = colSums(co_exist)
        
        p_pair[[n]] = co_exist
        names(p_pair)[n] = paste(type1,type2,sep='_')
        
        co_exist_all = vector()
        for (c in 1:ncol(reads_table_mom)) {
          x = seq(ncol(reads_table_mom))
          x = x[-c]
          co_exist=sapply(x, function(j) (reads_table_mom[,c]*reads_table_baby[,j]))
          co_exist = colSums(co_exist)
          co_exist_all = c(co_exist_all,co_exist)
        }
        
        p_unpair[[n]] = co_exist_all
        names(p_unpair)[n] = paste(type1,type2,sep='_')
      }
    }
    
    data = as.data.frame(matrix(data = NA, nrow=0, ncol = 3))
    pvalue = vector()
    for (a in 1:length(p_pair)) {
      data_2 = as.data.frame(p_pair[[a]])
      colnames(data_2) = 'V1'
      data_2$V2 = 'Paired'
      data_2$V3 = names(p_pair)[a]
      data = rbind(data,data_2)
      
      data_3 = as.data.frame(p_unpair[[a]])
      colnames(data_3) = 'V1'
      data_3$V2 = 'Unpaired'
      data_3$V3 = names(p_unpair)[a]
      data = rbind(data,data_3)
      
      pvalue[a] = wilcox.test(data_2$V1,data_3$V1)$p.value
      
    }
    
    ggboxplot(data, x = "V3", y = "V1", fill = "V2", outlier.shape = NA)+ 
      rotate_x_text(45)+ 
      coord_cartesian(ylim=c(0, 25))
    
    ggsave('co_existence_all_days.pdf',width=6, height=4)
    
    write.csv(pvalue,'co_existence_all_days.csv')
  }
  
  # day 0
  {
    type1_all = c('BCKD','BRCD')
    type2_all = c('MCKD','MRCD','MV1D')
    
    p_pair = list()
    
    p_unpair = list()
    n=0
    
    for (a in 1: length(type1_all)) {
      for (b in 1: length(type2_all)) {
        n=n+1
        type1 = type1_all[a]
        type2 = type2_all[b]
        
        # baby cheek and mom cheek
        keep = metadata_16s_2$SampleType == type1 | metadata_16s_2$SampleType == type2
        metadata = metadata_16s_2[keep,]
        
        metadata_baby = metadata
        keep = metadata_baby$SampleType == type1 & str_detect(metadata_baby$Flag, '0')
        metadata_baby = metadata_baby[keep,]
        
        metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
        metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]
        keep = duplicated(metadata_baby$ParticipantID)
        metadata_baby = metadata_baby[!keep,]
        
        p_number_1 = unique(metadata_baby$ParticipantID)
        length(p_number_1)
        keep = metadata$ParticipantID %in% p_number_1
        metadata = metadata[keep,]
        
        metadata_mom = metadata[metadata$SampleType == type2,]
        metadata_mom = metadata_mom[order(metadata_mom$VisitNum),]
        metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
        keep = duplicated(metadata_mom$ParticipantID)
        metadata_mom = metadata_mom[!keep,]
        
        participant_list = intersect(metadata_mom$ParticipantID, metadata_baby$ParticipantID)
        
        keep = metadata_baby$ParticipantID %in% participant_list
        metadata_baby = metadata_baby[keep,]
        keep = metadata_mom$ParticipantID %in% participant_list
        metadata_mom = metadata_mom[keep,]
        
        reads_table_mom = reads_table_16s_2[metadata_mom$SampleID]
        reads_table_baby = reads_table_16s_2[metadata_baby$SampleID]
        
        reads_table_mom = reads_table_mom >= th
        reads_table_baby = reads_table_baby >= th
        
        co_exist=sapply(1:ncol(reads_table_mom), function(j) (reads_table_mom[,j]*reads_table_baby[,j]))
        co_exist = colSums(co_exist)
        
        p_pair[[n]] = co_exist
        names(p_pair)[n] = paste(type1,type2,sep='_')
        
        co_exist_all = vector()
        for (c in 1:ncol(reads_table_mom)) {
          x = seq(ncol(reads_table_mom))
          x = x[-c]
          co_exist=sapply(x, function(j) (reads_table_mom[,c]*reads_table_baby[,j]))
          co_exist = colSums(co_exist)
          co_exist_all = c(co_exist_all,co_exist)
        }
        
        p_unpair[[n]] = co_exist_all
        names(p_unpair)[n] = paste(type1,type2,sep='_')
      }
    }
    
    data = as.data.frame(matrix(data = NA, nrow=0, ncol = 3))
    pvalue = vector()
    for (a in 1:length(p_pair)) {
      data_2 = as.data.frame(p_pair[[a]])
      colnames(data_2) = 'V1'
      data_2$V2 = 'Paired'
      data_2$V3 = names(p_pair)[a]
      data = rbind(data,data_2)
      
      data_3 = as.data.frame(p_unpair[[a]])
      colnames(data_3) = 'V1'
      data_3$V2 = 'Unpaired'
      data_3$V3 = names(p_unpair)[a]
      data = rbind(data,data_3)
      
      pvalue[a] = wilcox.test(data_2$V1,data_3$V1)$p.value
      
    }
    
    ggboxplot(data, x = "V3", y = "V1", fill = "V2", outlier.shape = NA)+ 
      rotate_x_text(45)+ 
      coord_cartesian(ylim=c(0, 50))
    
    ggsave('co_existence_day_0.pdf',width=4, height=4)
    write.csv(pvalue,'co_existence_day_0.csv')
  }
  
}



##### mom factors associated with Shannon index of baby's microbiome #####
setwd('/Users/binzhu/Desktop/Mom-baby/')
colnames(metadata_16s)[colnames(metadata_16s) == 'mom_delivery_method (C-section yes / vaginal no)'] = 'mom_delivery_method'

type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    Participant_list = metadata_baby$ParticipantID
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby = metadata_baby[keep,]
    reads_table_baby = reads_table_baby[,keep]
    
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]    # reads table
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    metadata_mom = metadata_mom[keep,]
    
    sample_pair_id = cbind(metadata_mom$ParticipantID ,metadata_mom$SampleID, metadata_baby$SampleID)
    colnames(sample_pair_id) = c('ParticipantID','Mom_SampleID','Baby_SampleID')
    write.csv(sample_pair_id,paste0('sample_pair_id_',type1_all[a],'.csv'))
    
    # get alpha diversity
    reads_table_baby_2 = t(reads_table_baby)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness <- as.matrix(alpha.shannon_diversity/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
  }
  
  # get numeric and character factors 
  {
    metadata_num = as.data.frame(matrix(data = NA, nrow = nrow(metadata_mom), ncol =0))
    metadata_cha = as.data.frame(matrix(data = NA, nrow = nrow(metadata_mom), ncol =0))
    n =1
    m =1
    metadata_mom[,1:ncol(metadata_mom)] <- lapply(metadata_mom[,1:ncol(metadata_mom)],as.character)
    
    c2 = c('0','1','2','3','4','5','6','7','8','9','.')
    
    for (a2 in 1:ncol(metadata_mom)) {
      b2 = metadata_mom[,a2]
      b2 = b2[!is.na(b2)]
      b2 = strsplit(b2,'*')
      b2 = unlist(b2)
      b2 = unique(b2)
      keep = b2 %in% c2
      
      if (sum(!keep) == 0) {
        metadata_mom[,a2] <- as.numeric(metadata_mom[,a2])
        metadata_num <- cbind(metadata_num,metadata_mom[,a2] )
        colnames(metadata_num)[n] = colnames(metadata_mom)[a2]
        n=n+1
      } else {
        metadata_cha <- cbind(metadata_cha,metadata_mom[,a2])
        colnames(metadata_cha)[m] = colnames(metadata_mom)[a2]
        m=m+1
      }
    }
    
  }
  
  # case numbers
  case_number_cha = as.data.frame(matrix(data = NA, nrow =7, ncol = ncol(metadata_cha)))
  colnames(case_number_cha) = colnames(metadata_cha)
  row.names(case_number_cha) = c("BCKD_y","BCKD_n", "BRCD_y", "BRCD_n","BS1D_y","BS1D_n",'Levels')
  
  case_number_num = as.data.frame(matrix(data = NA, nrow =3, ncol = ncol(metadata_num)))
  colnames(case_number_num) = colnames(metadata_num)
  row.names(case_number_num) = c("BCKD","BRCD","BS1D")
  ########### character factor ##########
  metadata = metadata_cha
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 8))
  colnames(correlation) = c('Type','Factor','shannon_pvalue','shannon_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)',
                            'evenness_pvalue','evenness_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)','otu_pvalue','otu_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)')
  
  correlation$Type = type1_all[a]
  
  
  for (b in 1: ncol(metadata)) {
    correlation$Factor[b] <- colnames(metadata)[b]
    
    keep = !is.na(metadata[,b])
    
    metadata_2 = metadata[keep,b]
    
    shannon_2 = alpha.shannon_diversity[keep]
    evenness_2 = alpha.evenness[keep]
    otu_2 = alpha.ovserved_OTU[keep]
    reads_table_5 = reads_table_baby[,keep]
    
#    days_rel2birth = vector()
#    for (c in 1: ncol(reads_table_5)) {
#      days_rel2birth[c] = metadata_16s$days_rel2birth[metadata_16s$SampleID == colnames(reads_table_5)[c]]
#    }
      
    ###
    factor_level = unique(metadata_2)
    factor_level = factor_level[order(factor_level, decreasing = T)]
    
    case_number_cha[7,b] = length(factor_level)
    
    factor_number = factor_level
    for (b3 in 1: length(factor_level)) {
      factor_number[b3] = sum(metadata_2 == factor_level[b3])
    }
    factor_number = as.numeric(factor_number)
    factor_number = factor_number[order(factor_number,decreasing = T)]
    
    # get case numbers
    case_number_cha[a*2-1,b] = factor_number[2]
    case_number_cha[a*2,b] = factor_number[1]
    
    if (sum(keep) <= 10) {
      next
    }
    
    if (length(factor_level) != 2) {
      next
    }
    
    if (factor_number[1] < 5 | factor_number[2] < 5) { 
      next
    }
    
    # calculate ratio of quantile
    data1 = shannon_2[metadata_2 == factor_level[1]]
    data2 = shannon_2[metadata_2 == factor_level[2]]
    quantile1 = quantile(data1)
    quantile2 = quantile(data2)
    quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
    quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
    quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
    quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
    quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
    quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
    correlation$`shannon_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
    
    data1 = evenness_2[metadata_2 == factor_level[1]]
    data2 = evenness_2[metadata_2 == factor_level[2]]
    quantile1 = quantile(data1)
    quantile2 = quantile(data2)
    quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
    quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
    quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
    quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
    quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
    quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
    correlation$`evenness_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
    
    data1 = otu_2[metadata_2 == factor_level[1]]
    data2 = otu_2[metadata_2 == factor_level[2]]
    quantile1 = quantile(data1)
    quantile2 = quantile(data2)
    quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
    quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
    quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
    quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
    quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
    quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
    correlation$`otu_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
    
    # 
    temp = as.data.frame(cbind(shannon_2,metadata_2))
    colnames(temp) = c('V1','V2')
    correlation$shannon_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
    
    if (correlation$shannon_pvalue[b] <= 0.05) {
      temp$V1 = as.numeric(as.character(temp$V1))
      ggplot(data=temp, aes(x=V2, y=V1)) +
        geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
        geom_jitter(shape=16, size = 0.5)+
        labs(x = correlation$Factor[b], y = paste0("Shannon index",'\n','of ',type1_all[a]))+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw()
      ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
      
      if (correlation$Factor[b] == 'mom_time_in_hospital') {
        ggplot(data=temp, aes(x=V2, y=V1)) +
          geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
          geom_jitter(shape=16, size = 0.5)+
          labs(x = correlation$Factor[b], y = paste0("Shannon index",'\n','of ',type1_all[a]))+ 
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
        ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=5)
      }
    }
    
    temp = as.data.frame(cbind(evenness_2,metadata_2))
    colnames(temp) = c('V1','V2')
    correlation$evenness_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
    
    if (correlation$evenness_pvalue[b] <= 0.05) {
      temp$V1 = as.numeric(as.character(temp$V1))
      ggplot(data=temp, aes(x=V2, y=V1)) +
        geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
        geom_jitter(shape=16, size = 0.5)+
        labs(x = correlation$Factor[b], y = "Evenness")+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
      ggsave(paste0('alpha_evenness_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
    }
    
    temp = as.data.frame(cbind(otu_2,metadata_2))
    colnames(temp) = c('V1','V2')
    correlation$otu_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
    
    if (correlation$otu_pvalue[b] <= 0.05) {
      temp$V1 = as.numeric(as.character(temp$V1))
      ggplot(data=temp, aes(x=V2, y=V1)) +
        geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
        geom_jitter(shape=16, size = 0.5)+
        labs(x = correlation$Factor[b], y = "Observed taxa")+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
      ggsave(paste0('alpha_otu_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
    }
    
    #      if (correlation$shannon_pvalue[b] <= 0.05 & length(factor_level) == 2) {
    #        das = dif_abundance3(reads_table_5,metadata_2) # differential analysis
    #        if (!is.na(das)) {
    #          write.csv(das$data,file = paste0('das_',type1_all[a],'_',colnames(metadata)[b],'.csv'))
    #      }
    
  }
  
  
  write.csv(case_number_cha,paste0('case_number_cha_',type1_all[a],'.csv'))
  
  
  
  if (a ==1) {
    correlation_all_cha = correlation
  } else {
    correlation_all_cha = rbind(correlation_all_cha,correlation)
  }
  
  
  
  
  
  
  
  
  
  
  
  ########### numeric factor ##########
  metadata = metadata_num
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 8))
  colnames(correlation) = c('Type','Factor','shannon_pvalue','shannon_lm_function',
                            'evenness_pvalue','evenness_lm_function','otu_pvalue','otu_lm_function')
  
  correlation$Type = type1_all[a]
  
  
  for (b in 1: ncol(metadata)) {
    correlation$Factor[b] <- colnames(metadata)[b]
    
    keep = !is.na(metadata[,b])
    
    # get case numbers
    case_number_num[a,b] = sum(keep)
    
    #
    if (sum(keep) <= 10) {
      next
    }
    
    metadata_2 = metadata[keep,b]
    shannon_2 = alpha.shannon_diversity[keep]
    evenness_2 = alpha.evenness[keep]
    otu_2 = alpha.ovserved_OTU[keep]
    reads_table_5 = reads_table_baby[,keep]
    
    ###
    factor_level = unique(metadata_2)
    
    if (length(factor_level) < 3) {
      next
    }
    
    # lm model to test shannon correlation
    temp = as.data.frame(cbind(shannon_2,metadata_2))
    colnames(temp) = c('V1','V2')
    temp$V1 = as.numeric(as.character(temp$V1))
    
    # remove outliers
    outliers_test = boxplot(temp$V2, plot=FALSE)$out
    
    if (length(outliers_test) > 0) {
      n = !(temp$V2 %in% outliers_test)
      temp = temp[n,]
    }
    
    if (length(unique(temp$V2)) > 3) {
      
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$shannon_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$shannon_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$shannon_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      if (correlation$shannon_pvalue[b] <= 0.05) {
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Shannon index of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        
        ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
    }
    
    
    
    # lm model to test evenness correlation
    temp = as.data.frame(cbind(evenness_2,metadata_2))
    colnames(temp) = c('V1','V2')
    
    temp$V1 = as.numeric(as.character(temp$V1))
    ##### remove outliers #####
    outliers_test = boxplot(temp$V2, plot=FALSE)$out
    
    if (length(outliers_test) > 0) {
      n = !(temp$V2 %in% outliers_test)
      temp = temp[n,]
    }
    
    if (length(unique(temp$V2)) > 3) {
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$evenness_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$evenness_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$evenness_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      
      if (correlation$evenness_pvalue[b] <= 0.05) {
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Evenness of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        ggsave(paste0('alpha_evenness_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
    }
    
    
    # lm model to test otu correlation
    temp = as.data.frame(cbind(otu_2,metadata_2))
    colnames(temp) = c('V1','V2')
    
    temp$V1 = as.numeric(as.character(temp$V1))
    ##### remove outliers #####
    outliers_test = boxplot(temp$V2, plot=FALSE)$out
    
    if (length(outliers_test) > 0) {
      n = !(temp$V2 %in% outliers_test)
      temp = temp[n,]
    }
    
    if (length(unique(temp$V2)) > 3) {
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$otu_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$otu_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$otu_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      if (correlation$otu_pvalue[b] <= 0.05) {
        
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Observed taxa of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        ggsave(paste0('alpha_otu_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
    }
    
  }
  
  
  if (a ==1) {
    correlation_all_num = correlation
  } else {
    correlation_all_num = rbind(correlation_all_num,correlation)
  }
  
  
  write.csv(case_number_num,paste0('case_number_num_',type1_all[a],'.csv'))
}
write.csv(correlation_all_num,'impact_of_factor_num_to_baby_microbiome_alpha.csv',  row.names = F)
write.csv(correlation_all_cha,'impact_of_factor_cha_to_baby_microbiome_alpha.csv',  row.names = F)








##### baby factors associated with Shannon index of baby's microbiome #####
setwd('/Users/binzhu/Desktop/Mom-baby/')
type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    Participant_list = metadata_baby$ParticipantID
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby = metadata_baby[keep,]
    reads_table_baby = reads_table_baby[,keep]
    
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]    # reads table
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    metadata_mom = metadata_mom[keep,]
    
    sample_pair_id = cbind(metadata_mom$ParticipantID ,metadata_mom$SampleID, metadata_baby$SampleID)
    colnames(sample_pair_id) = c('ParticipantID','Mom_SampleID','Baby_SampleID')
    write.csv(sample_pair_id,paste0('sample_pair_id_',type1_all[a],'.csv'))
    
    # get alpha diversity
    reads_table_baby_2 = t(reads_table_baby)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness <- as.matrix(alpha.shannon_diversity/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
  }
  
  # get numeric and character factors 
  {
    metadata_num = as.data.frame(matrix(data = NA, nrow = nrow(metadata_baby), ncol =0))
    metadata_cha = as.data.frame(matrix(data = NA, nrow = nrow(metadata_baby), ncol =0))
    n =1
    m =1
    metadata_baby[,1:ncol(metadata_baby)] <- lapply(metadata_baby[,1:ncol(metadata_baby)],as.character)
    
    c2 = c('0','1','2','3','4','5','6','7','8','9','.')
    
    for (a2 in 1:ncol(metadata_baby)) {
      b2 = metadata_baby[,a2]
      b2 = b2[!is.na(b2)]
      b2 = strsplit(b2,'*')
      b2 = unlist(b2)
      b2 = unique(b2)
      keep = b2 %in% c2
      
      if (sum(!keep) == 0) {
        metadata_baby[,a2] <- as.numeric(metadata_baby[,a2])
        metadata_num <- cbind(metadata_num,metadata_baby[,a2] )
        colnames(metadata_num)[n] = colnames(metadata_baby)[a2]
        n=n+1
      } else {
        metadata_cha <- cbind(metadata_cha,metadata_baby[,a2])
        colnames(metadata_cha)[m] = colnames(metadata_baby)[a2]
        m=m+1
      }
    }
    
  }
  
  ########### character factor ##########
  # case numbers
  case_number_cha = as.data.frame(matrix(data = NA, nrow =7, ncol = ncol(metadata_cha)))
  colnames(case_number_cha) = colnames(metadata_cha)
  row.names(case_number_cha) = c("BCKD_y","BCKD_n", "BRCD_y", "BRCD_n","BS1D_y","BS1D_n",'Levels')
  
  metadata = metadata_cha
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 8))
  colnames(correlation) = c('Type','Factor','shannon_pvalue','shannon_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)',
                            'evenness_pvalue','evenness_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)','otu_pvalue','otu_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)')
  
  correlation$Type = type1_all[a]
  
  try(
    for (b in 1: ncol(metadata)) {
      correlation$Factor[b] <- colnames(metadata)[b]
      
      keep = !is.na(metadata[,b])
      
      metadata_2 = metadata[keep,b]
      shannon_2 = alpha.shannon_diversity[keep]
      evenness_2 = alpha.evenness[keep]
      otu_2 = alpha.ovserved_OTU[keep]
      reads_table_5 = reads_table_baby[,keep]
      
      ###
      factor_level = unique(metadata_2)
      factor_level = factor_level[order(factor_level, decreasing = T)]
      
      case_number_cha[7,b] = length(factor_level)
      
      factor_number = factor_level
      for (b3 in 1: length(factor_level)) {
        factor_number[b3] = sum(metadata_2 == factor_level[b3])
      }
      factor_number = as.numeric(factor_number)
      factor_number = factor_number[order(factor_number,decreasing = T)]
      
      # get case numbers
      case_number_cha[a*2-1,b] = factor_number[2]
      case_number_cha[a*2,b] = factor_number[1]
      
      if (sum(keep) <= 10) {
        next
      }
      
      if (length(factor_level) != 2) {
        next
      }
      
      if (factor_number[1] < 5 | factor_number[2] < 5) { 
        next
      }
      
      # calculate ratio of quantile
      data1 = shannon_2[metadata_2 == factor_level[1]]
      data2 = shannon_2[metadata_2 == factor_level[2]]
      quantile1 = quantile(data1)
      quantile2 = quantile(data2)
      quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
      quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
      quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
      quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
      quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
      quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
      correlation$`shannon_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
      
      data1 = evenness_2[metadata_2 == factor_level[1]]
      data2 = evenness_2[metadata_2 == factor_level[2]]
      quantile1 = quantile(data1)
      quantile2 = quantile(data2)
      quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
      quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
      quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
      quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
      quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
      quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
      correlation$`evenness_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
      
      data1 = otu_2[metadata_2 == factor_level[1]]
      data2 = otu_2[metadata_2 == factor_level[2]]
      quantile1 = quantile(data1)
      quantile2 = quantile(data2)
      quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
      quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
      quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
      quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
      quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
      quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
      correlation$`otu_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
      
      # 
      temp = as.data.frame(cbind(shannon_2,metadata_2))
      colnames(temp) = c('V1','V2')
      correlation$shannon_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
      
      if (correlation$shannon_pvalue[b] <= 0.05) {
        temp$V1 = as.numeric(as.character(temp$V1))
        ggplot(data=temp, aes(x=V2, y=V1)) +
          geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
          geom_jitter(shape=16, size = 0.5)+
          labs(x = correlation$Factor[b], y = paste0("Shannon index",'\n','of ',type1_all[a]))+ 
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw()
        ggsave(paste0('alpha_shannon_baby_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
        
        if (correlation$Factor[b] == 'mom_time_in_hospital') {
          ggplot(data=temp, aes(x=V2, y=V1)) +
            geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
            geom_jitter(shape=16, size = 0.5)+
            labs(x = correlation$Factor[b], y = paste0("Shannon index",'\n','of ',type1_all[a]))+ 
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
          ggsave(paste0('alpha_shannon_baby_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=5)
        }
      }
      
      temp = as.data.frame(cbind(evenness_2,metadata_2))
      colnames(temp) = c('V1','V2')
      correlation$evenness_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
      
      if (correlation$evenness_pvalue[b] <= 0.05) {
        temp$V1 = as.numeric(as.character(temp$V1))
        ggplot(data=temp, aes(x=V2, y=V1)) +
          geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
          geom_jitter(shape=16, size = 0.5)+
          labs(x = correlation$Factor[b], y = "Evenness")+ 
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
        ggsave(paste0('alpha_evenness_baby_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
      }
      
      temp = as.data.frame(cbind(otu_2,metadata_2))
      colnames(temp) = c('V1','V2')
      correlation$otu_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
      
      if (correlation$otu_pvalue[b] <= 0.05) {
        temp$V1 = as.numeric(as.character(temp$V1))
        ggplot(data=temp, aes(x=V2, y=V1)) +
          geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
          geom_jitter(shape=16, size = 0.5)+
          labs(x = correlation$Factor[b], y = "Observed taxa")+ 
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
        ggsave(paste0('alpha_otu_baby_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
      }
      
      #      if (correlation$shannon_pvalue[b] <= 0.05 & length(factor_level) == 2) {
      #        das = dif_abundance3(reads_table_5,metadata_2) # differential analysis
      #        if (!is.na(das)) {
      #          write.csv(das$data,file = paste0('das_',type1_all[a],'_',colnames(metadata)[b],'.csv'))
      #      }
      
    }
  )
  
  write.csv(case_number_cha,paste0('case_number_baby_cha_',type1_all[a],'.csv'))
  
  
  
  if (a ==1) {
    correlation_all_cha = correlation
  } else {
    correlation_all_cha = rbind(correlation_all_cha,correlation)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  ########### numeric factor ##########
  metadata = as.data.frame(cbind(metadata_num$height, metadata_num$weight, metadata_num$BMI, metadata_num$pulse, metadata_num$days_rel2birth))
  colnames(metadata) = c('baby_height','baby_weight','baby_BMI','baby_pulse','Sample_collection_time')
  
  # case numbers
  case_number_num = as.data.frame(matrix(data = NA, nrow =3, ncol = ncol(metadata)))
  colnames(case_number_num) = colnames(metadata)
  row.names(case_number_num) = c("BCKD","BRCD","BS1D")
  
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 8))
  colnames(correlation) = c('Type','Factor','shannon_pvalue','shannon_lm_function',
                            'evenness_pvalue','evenness_lm_function','otu_pvalue','otu_lm_function')
  
  correlation$Type = type1_all[a]
  
  try(
    for (b in 1: ncol(metadata)) {
      correlation$Factor[b] <- colnames(metadata)[b]
      
      keep = !is.na(metadata[,b])
      
      # get case numbers
      case_number_num[a,b] = sum(keep)
      
      #
      if (sum(keep) <= 10) {
        next
      }
      
      metadata_2 = metadata[keep,b]
      shannon_2 = alpha.shannon_diversity[keep]
      evenness_2 = alpha.evenness[keep]
      otu_2 = alpha.ovserved_OTU[keep]
      reads_table_5 = reads_table_baby[,keep]
      
      ###
      factor_level = unique(metadata_2)
      
      if (length(factor_level) < 3) {
        next
      }
      
      # lm model to test shannon correlation
      temp = as.data.frame(cbind(shannon_2,metadata_2))
      colnames(temp) = c('V1','V2')
      
      temp$V1 = as.numeric(as.character(temp$V1))
      
      ##### remove outliers #####
      outliers_test = boxplot(temp$V2, plot=FALSE)$out
      
      if (length(outliers_test) > 0) {
        n = !(temp$V2 %in% outliers_test)
        temp = temp[n,]
      }
      
      #####
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$shannon_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$shannon_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$shannon_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      if (correlation$shannon_pvalue[b] <= 0.05) {
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Shannon index of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        
        ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
      
      
      
      # lm model to test evenness correlation
      temp = as.data.frame(cbind(evenness_2,metadata_2))
      colnames(temp) = c('V1','V2')
      temp$V1 = as.numeric(as.character(temp$V1))
      
      ##### remove outliers #####
      outliers_test = boxplot(temp$V2, plot=FALSE)$out
      
      if (length(outliers_test) > 0) {
        n = !(temp$V2 %in% outliers_test)
        temp = temp[n,]
      }
      
      #####
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$evenness_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$evenness_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$evenness_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      
      if (correlation$evenness_pvalue[b] <= 0.05) {
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Evenness of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        ggsave(paste0('alpha_evenness_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
      
      
      
      # lm model to test otu correlation
      temp = as.data.frame(cbind(otu_2,metadata_2))
      colnames(temp) = c('V1','V2')
      temp$V1 = as.numeric(as.character(temp$V1))
      
      ##### remove outliers #####
      outliers_test = boxplot(temp$V2, plot=FALSE)$out
      
      if (length(outliers_test) > 0) {
        n = !(temp$V2 %in% outliers_test)
        temp = temp[n,]
      }
      
      #####
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$otu_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$otu_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$otu_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      if (correlation$otu_pvalue[b] <= 0.05) {
        temp$V1 = as.numeric(as.character(temp$V1))
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Observed taxa of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        ggsave(paste0('alpha_otu_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
      
    }
  )
  
  if (a ==1) {
    correlation_all_num = correlation
  } else {
    correlation_all_num = rbind(correlation_all_num,correlation)
  }
  
  
  write.csv(case_number_num,paste0('case_number_baby_',type1_all[a],'.csv'))
}
write.csv(correlation_all_num,'impact_of_factor_baby_num_to_baby_microbiome_alpha.csv',  row.names = F)
write.csv(correlation_all_cha,'impact_of_factor_baby_cha_to_baby_microbiome_alpha.csv',  row.names = F)









##### Shannon vs. lipids/cytokines and mediation effect of lipids / cytokines on the influence of metadata on Shannon #####
library(bda) # for The Sobel mediation test
library(tidyverse)
library(knitr)
library(lavaan)  # for The sem mediation test
library(psych)
library(MBESS)

reads_table_lipid[reads_table_lipid == 0] = NA

type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]

correlation_all_cytokine = as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
colnames(correlation_all_cytokine) = c('Type','Factor','pvalue', 'Case_number')

correlation_all_lipid = as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
colnames(correlation_all_lipid) = c('Type','Factor','pvalue', 'Case_number')

sem_test = read.csv('correlation_shannon_lipid_input.csv',header = F)

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    Participant_list = metadata_baby$ParticipantID
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby = metadata_baby[keep,]
    reads_table_baby = reads_table_baby[,keep]
    
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]    # reads table
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    metadata_mom = metadata_mom[keep,]
    
    
    # get paired lipid data
    keep = colnames(reads_table_lipid) %in% metadata_mom$SampleID
    sum(keep)
    reads_table_lipid_mom =  reads_table_lipid[,keep]
    
    keep = metadata_mom$SampleID %in% colnames(reads_table_lipid_mom)
    sum(keep)
    metadata_mom_lipid = metadata_mom[keep,]  # lipid mom metadata
    metadata_baby_lipid = metadata_baby[keep,]  # lipid baby metadata
    reads_table_baby_lipid = reads_table_baby[,keep]
    reads_table_lipid_mom = reads_table_lipid_mom[metadata_mom_lipid$SampleID]  # lipid reads table
    reads_table_lipid_mom = log10(reads_table_lipid_mom)    # log10 transformation
    
    # get paired cytokine data
    keep = colnames(reads_table_cytokines) %in% metadata_mom$SampleID
    sum(keep)
    reads_table_cytokines_mom =  reads_table_cytokines[,keep]
    
    keep = metadata_mom$SampleID %in% colnames(reads_table_cytokines_mom)
    sum(keep)
    metadata_mom_cytokine = metadata_mom[keep,]  # cytokine mom metadata
    metadata_baby_cytokine = metadata_baby[keep,]  # cytokine baby metadata
    reads_table_baby_cytokine = reads_table_baby[,keep]
    reads_table_cytokines_mom = reads_table_cytokines_mom[metadata_mom_cytokine$SampleID]   # cytokine reads table
    reads_table_cytokines_mom = log10(reads_table_cytokines_mom)         # log10 transformation
    
    # get alpha diversity
    reads_table_baby_2 = t(reads_table_baby)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity <- as.numeric(diversity(reads_table_baby_2))   # lipid Shannon
    
    reads_table_baby_lipid = t(reads_table_baby_lipid)
    reads_table_baby_lipid = Rarefy(reads_table_baby_lipid, depth = min(rowSums(reads_table_baby_lipid)))
    reads_table_baby_lipid <- reads_table_baby_lipid$otu.tab.rff
    reads_table_baby_lipid <- as.data.frame(reads_table_baby_lipid)
    alpha.shannon_diversity_lipid <- as.numeric(diversity(reads_table_baby_lipid))   # lipid Shannon
    
    reads_table_baby_cytokine = t(reads_table_baby_cytokine)
    reads_table_baby_cytokine = Rarefy(reads_table_baby_cytokine, depth = min(rowSums(reads_table_baby_cytokine)))
    reads_table_baby_cytokine <- reads_table_baby_cytokine$otu.tab.rff
    reads_table_baby_cytokine <- as.data.frame(reads_table_baby_cytokine)
    alpha.shannon_diversity_cytokine <- as.numeric(diversity(reads_table_baby_cytokine))    # cytokine Shannon
    
    # output design
    case_table_cytokine = cbind(metadata_mom_cytokine$ParticipantID,colnames(reads_table_cytokines_mom), row.names(reads_table_baby_cytokine))
    case_table_lipid = cbind(metadata_mom_lipid$ParticipantID,colnames(reads_table_lipid_mom), row.names(reads_table_baby_lipid))
    write.csv(case_table_cytokine,paste0(type1_all[a],'_case_table_cytokine.csv'))
    write.csv(case_table_lipid,paste0(type1_all[a],'_case_table_lipid.csv'))
  }
  
  ### test correlation cytokine ###
  correlation_all = as.data.frame(matrix(data = NA, nrow = nrow(reads_table_cytokines_mom), ncol = 4))
  colnames(correlation_all) = c('Type','Factor','pvalue', 'Case_number')
  correlation_all[,1] = type1_all[a]
  correlation_all[,4] = ncol(reads_table_cytokines_mom)
  
  for (b in 1:nrow(reads_table_cytokines_mom)) {
    correlation_all[b,2] = row.names(reads_table_cytokines_mom)[b]
    
    pvalue = cor.test(alpha.shannon_diversity_cytokine, as.numeric(reads_table_cytokines_mom[b,]), method = "pearson")
    pvalue = pvalue$p.value
    correlation_all[b,3] = pvalue
    
    if (pvalue <= 0.05) {
      data = as.data.frame(matrix(data = NA, nrow = length(alpha.shannon_diversity_cytokine), ncol = 2))
      data$V1 = alpha.shannon_diversity_cytokine
      data$V2 = as.numeric(reads_table_cytokines_mom[b,])
      
      linearMod <- lm(V1 ~ V2, data=data)
      r_squared <- summary(linearMod)$r.squared
      pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])
      
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      ggplot(data, aes(x = V2, y = V1)) + 
        geom_point() +
        stat_smooth(method = "lm", col = "blue")+
        labs(x = paste0('Log10 ',row.names(reads_table_cytokines_mom)[b]), 
             y = paste0('Shannon index of ', type1_all[a]))+ 
        ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", pvalue))  +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))+ theme_bw()
      ggsave(paste0('Shannon_vs_cytokine_',type1_all[a],'_',row.names(reads_table_cytokines_mom)[b],'.pdf'),width=3, height=3)
      
      
      
    }
  }
  
  if (a == 1) {
    correlation_all_cytokine = correlation_all
  } else {
    correlation_all_cytokine = rbind(correlation_all_cytokine,correlation_all)
  }
  
  ### test correlation lipid ###
  correlation_all = as.data.frame(matrix(data = NA, nrow = nrow(reads_table_lipid_mom), ncol = 4))
  colnames(correlation_all) = c('Type','Factor','pvalue', 'Case_number')
  correlation_all[,1] = type1_all[a]
  correlation_all[,4] = ncol(reads_table_lipid_mom)
  
  for (b in 1:nrow(reads_table_lipid_mom)) {
    correlation_all[b,2] = row.names(reads_table_lipid_mom)[b]
    
    keep = !is.na(reads_table_lipid_mom[b,])
    if (sum(keep) < 8) {
      next
    }
    
    data = cbind(alpha.shannon_diversity_lipid,as.numeric(reads_table_lipid_mom[b,]))
    keep = !is.na(data[,2])
    data = data[keep,]
    
    pvalue = cor.test(alpha.shannon_diversity_lipid, as.numeric(reads_table_lipid_mom[b,]), method = "pearson")
    pvalue = pvalue$p.value
    correlation_all[b,3] = pvalue
    
    if (pvalue <= 0.05) {
      data = as.data.frame(matrix(data = NA, nrow = length(alpha.shannon_diversity_lipid), ncol = 2))
      data$V1 = alpha.shannon_diversity_lipid
      data$V2 = as.numeric(reads_table_lipid_mom[b,])
      
      keep = !is.na(reads_table_lipid_mom[b,])
      data = data[keep,]
      
      linearMod <- lm(V1 ~ V2, data=data)
      r_squared <- summary(linearMod)$r.squared
      pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])
      
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      ggplot(data, aes(x = V2, y = V1)) + 
        geom_point() +
        stat_smooth(method = "lm", col = "blue")+
        labs(x = paste0('Log10 ',row.names(reads_table_lipid_mom)[b]), 
             y = paste0('Shannon index of ', type1_all[a]))+ 
        ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", pvalue))  +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))+ theme_bw()
      ggsave(paste0('Shannon_vs_lipid_',type1_all[a],'_',row.names(reads_table_lipid_mom)[b],'.pdf'),width=3, height=3)
      
      
      
      # mediation analysis using sem   lavaan R package
      for (c in 1: ncol(metadata_baby_lipid)) {
        if (paste(type1_all[a],colnames(metadata_baby_lipid)[c], sep = '_') %in% sem_test$V3) {
          
          data = as.data.frame(matrix(data = NA, nrow = length(alpha.shannon_diversity_lipid), ncol = 3))
          data$V1 = alpha.shannon_diversity_lipid
          data$V2 = as.numeric(reads_table_lipid_mom[b,])
          data$V3 = metadata_baby_lipid[,c]
          data$V4 = metadata_baby_lipid[,c]
          data$V3 = as.numeric(factor(data$V3))
          
          keep = rowSums(is.na(data)) == 0 
          data = data[keep,]
          
          if (nrow(data) <8) {
            next
          }
          
          keep = length(unique(data$V1)) > 1 & length(unique(data$V2)) > 1 & length(unique(data$V3)) > 1
          if (!keep) {
            next
          }
          
          
          mod1 <- "# a path
        V2 ~ a * V3
        
        # b path
        V1 ~ b * V2
        
        # c prime path 
        V1 ~ cp * V3
        
        # indirect and total effects
        ab := a * b
        total := cp + ab"
          
          set.seed(1234)
          fsem1 <- sem(mod1, data = data, se = "bootstrap", bootstrap = 1000)
          
          pvalue_sem = parameterestimates(fsem1, boot.ci.type = "bca.simple", standardized = TRUE)
          pvalue_sem[pvalue_sem == 'V1' & !is.na(pvalue_sem)] = 'Shannon index'
          pvalue_sem[pvalue_sem == 'V2' & !is.na(pvalue_sem)] = row.names(reads_table_lipid_mom)[b]
          pvalue_sem[pvalue_sem == 'V3' & !is.na(pvalue_sem)] = colnames(metadata_baby_lipid)[c]
          
          if (pvalue_sem[7,8] <= 0.05 & pvalue_sem[1,8] <= 0.05 & pvalue_sem[2,8] <= 0.05 & pvalue_sem[3,8] <= 0.05) {
            write.csv(pvalue_sem,paste0('Mediation_',type1_all[a],'_',row.names(reads_table_lipid_mom)[b],'_',colnames(metadata_baby_lipid)[c],'.csv'))
            
            pdf(file = paste0('Mediation_',type1_all[a],'_',row.names(reads_table_lipid_mom)[b],'_',colnames(metadata_baby_lipid)[c],'.pdf'),width=5, height=4)
            with(data, mediation.effect.plot(x = V3, mediator = V2, dv = V1, ylab = paste0("Shannon index of ", type1_all[a]), xlab = paste0('Concentration of ', row.names(reads_table_lipid_mom)[b])))
            dev.off()
          }
        }
        
      }
      
      
      
      
      
      
    }
    
  }
  
  if (a == 1) {
    correlation_all_lipid = correlation_all
  } else {
    correlation_all_lipid = rbind(correlation_all_lipid,correlation_all)
  }
}

write.csv(correlation_all_lipid,'correlation_shannon_lipid.csv')
write.csv(correlation_all_cytokine,'correlation_shannon_cytokine.csv')











##### mom's factor associated with baby's beta diversity #####
setwd('/Users/binzhu/Desktop/Mom-baby/')
type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]
correlation_all = as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
colnames(correlation_all) = c('Type','Factor','Adonis_pvalue','Adj_pvalue')

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    Participant_list = metadata_baby$ParticipantID
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby = metadata_baby[keep,]
    reads_table_baby = reads_table_baby[,keep]
    
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]    # reads table
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    metadata_mom = metadata_mom[keep,]
    
  }
  
  metadata = metadata_mom
  reads_table = reads_table_baby
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 4))
  colnames(correlation) = c('Type','Factor','Adonis_pvalue','Adj_pvalue')
  
  correlation$Type = type1_all[a]
  
  for (b in 1: ncol(metadata)) {
    ##### prepare data #####
    # correlation by adonis2
    correlation$Factor[b] <- colnames(metadata)[b]
    
    keep = !is.na(metadata[,b])
    if (sum(keep) <= 10) {
      next
    }
    
    metadata_2 = metadata[keep,b]
    reads_table_2 = reads_table[,keep]
    reads_table_3 = reads_table[,keep]
    
    factor_level = unique(metadata_2)
    if (length(factor_level) < 2) {
      next
    }
    
    ##### convert numeric factor to charactor factor by grouping to two levels #####
    c2 = c('0','1','2','3','4','5','6','7','8','9','.')
    
    b2 = factor_level
    b2 = b2[!is.na(b2)]
    b2 = as.character(b2)
    b2 = strsplit(b2,'*')
    b2 = unlist(b2)
    b2 = unique(b2)
    keep = b2 %in% c2
    
    if (sum(!keep) == 0) {
      factor_level = as.numeric(as.character(factor_level))
      b2 = median(factor_level)
      keep = factor_level>=b2
      factor_level[keep] ='High'
      factor_level[!keep] ='Low'
      metadata_2 = as.numeric(metadata_2)
      keep = metadata_2>=b2
      metadata_2[keep] = 'High'
      metadata_2[!keep] = 'Low'
      metadata_2 = factor(metadata_2,levels = c('Low','High'))
      factor_level = unique(factor_level)
    } else {
      metadata_2 = factor(metadata_2,levels = c('No','Yes'))
    }
    
    if (sum(!is.na(metadata_2)) < 10) {
      next
    }
    
    factor_number = factor_level
    for (b3 in 1: length(factor_level)) {
      factor_number[b3] = sum(as.character(metadata_2) == factor_level[b3])
    }
    factor_number = as.numeric(factor_number)
    factor_number = factor_number[order(factor_number,decreasing = T)]
    
    if (sum(is.na(factor_number)) > 0) {
      next
    }
    
    if (factor_number[1] < 5 | factor_number[2] < 5) { 
      next
    } 
    
    
    ##### t-SNE #####
    reads_table_2 <- t(reads_table_2)
    reads_table_2 = Rarefy(reads_table_2, depth = min(rowSums(reads_table_2)))
    reads_table_2 <- reads_table_2$otu.tab.rff
    reads_table_2 <- as.data.frame(reads_table_2)
    
    metadata_2 = as.data.frame(metadata_2)
    pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")  
    correlation$Adonis_pvalue[b] <- pvalue[1,5]
    
    if (!is.na(pvalue[1,5])) {
      if (pvalue[1,5] <= 0.05) {
        
        # tsne
        tsne <- Rtsne(reads_table_2, dims = 2, perplexity=5, verbose=TRUE, max_iter = 5000)
        
        pic <- tsne$Y
        pic <- data.frame(pic,metadata_2$metadata_2)
        
        colnames(pic) <- c('X1','X2','Factor')
        
        ggplot(pic, aes(X1, X2,color = Factor))  +
          geom_point(size=0.5) +
          xlab(paste0("tsne1")) +
          ylab(paste0("tsne2"))+
          stat_ellipse(type = "t") + 
          coord_fixed()+ 
          ggtitle(colnames(metadata)[b])+
          theme(
            axis.title.x = element_text( size=7),
            axis.title.y = element_text( size=7),
            legend.text = element_text(size=7),
            legend.title = element_text(size=7),
            plot.title = element_text(hjust = 0.5, size = 7)
          ) + theme_bw() 
        ggsave(paste('tsne_',type1_all[a],colnames(metadata)[b],'.pdf', sep='_') , width=4, height=4)
        
        # difference analysis 
        keep <- metadata_2$metadata_2 == as.character(factor_level[1])
        if (length(factor_level) == 2 & sum(keep) >=5 & sum(!keep) >=5) {
          das = dif_abundance3(reads_table_3, metadata_2$metadata_2) # differential analysis
          if (!is.na(das)) {
            write.csv(das$data,file = paste('Abundance_',type1_all[a],colnames(metadata)[b],'.csv', sep='_'))
          }
        }
      }
    }
    
    correlation_all = rbind(correlation_all,correlation)
    
    
    
  }
  
  
}
correlation_all = correlation_all[!duplicated(correlation_all),]
write.csv(correlation_all,'impact_of_factor_to_baby_microbiome_beta.csv',  row.names = F)









##### baby's factor associated with baby's beta diversity #####
setwd('/Users/binzhu/Desktop/Mom-baby/')
type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]
correlation_all = as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
colnames(correlation_all) = c('Type','Factor','Adonis_pvalue','Adj_pvalue')

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    Participant_list = metadata_baby$ParticipantID
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby = metadata_baby[keep,]
    reads_table_baby = reads_table_baby[,keep]
    
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]    # reads table
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    metadata_mom = metadata_mom[keep,]
    
    
  }
  
  metadata = as.data.frame(cbind(metadata_baby$height, metadata_baby$weight, metadata_baby$BMI, metadata_baby$pulse, metadata_baby$days_rel2birth))
  colnames(metadata) = c('baby_height','baby_weight','baby_BMI','baby_pulse','Sample_collection_time')
  
  metadata = cbind(metadata,metadata_baby)
  
  reads_table = reads_table_baby
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 4))
  colnames(correlation) = c('Type','Factor','Adonis_pvalue','Adj_pvalue')
  
  correlation$Type = type1_all[a]
  
  for (b in 1: ncol(metadata)) {
    ##### prepare data #####
    # correlation by adonis2
    correlation$Factor[b] <- colnames(metadata)[b]
    
    keep = !is.na(metadata[,b])
    if (sum(keep) <= 10) {
      next
    }
    
    metadata_2 = metadata[keep,b]
    reads_table_2 = reads_table[,keep]
    reads_table_3 = reads_table[,keep]
    
    factor_level = unique(metadata_2)
    if (length(factor_level) < 2) {
      next
    }
    
    ##### convert numeric factor to charactor factor by grouping to two levels #####
    c2 = c('0','1','2','3','4','5','6','7','8','9','.')
    
    b2 = factor_level
    b2 = b2[!is.na(b2)]
    b2 = as.character(b2)
    b2 = strsplit(b2,'*')
    b2 = unlist(b2)
    b2 = unique(b2)
    keep = b2 %in% c2
    
    if (sum(!keep) == 0) {
      factor_level = as.numeric(as.character(factor_level))
      b2 = median(factor_level)
      keep = factor_level>=b2
      factor_level[keep] ='High'
      factor_level[!keep] ='Low'
      metadata_2 = as.numeric(metadata_2)
      keep = metadata_2>=b2
      metadata_2[keep] = 'High'
      metadata_2[!keep] = 'Low'
      metadata_2 = factor(metadata_2,levels = c('Low','High'))
      factor_level = unique(factor_level)
    } else {
      metadata_2 = factor(metadata_2,levels = c('No','Yes'))
    }
    
    if (sum(!is.na(metadata_2)) < 10) {
      next
    }
    
    factor_number = factor_level
    for (b3 in 1: length(factor_level)) {
      factor_number[b3] = sum(metadata_2 == factor_level[b3])
    }
    factor_number = as.numeric(factor_number)
    factor_number = factor_number[order(factor_number,decreasing = T)]
    
    if (sum(is.na(factor_number)) > 0) {
      next
    }
    
    if (factor_number[1] < 5 | factor_number[2] < 5) { 
      next
    } 
    
    
    ##### t-SNE #####
    reads_table_2 <- t(reads_table_2)
    reads_table_2 = Rarefy(reads_table_2, depth = min(rowSums(reads_table_2)))
    reads_table_2 <- reads_table_2$otu.tab.rff
    reads_table_2 <- as.data.frame(reads_table_2)
    
    metadata_2 = as.data.frame(metadata_2)
    pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")  
    correlation$Adonis_pvalue[b] <- pvalue[1,5]
    
    if (!is.na(pvalue[1,5])) {
      if (pvalue[1,5] <= 0.05) {
        
        # tsne
        tsne <- Rtsne(reads_table_2, dims = 2, perplexity=5, verbose=TRUE, max_iter = 5000)
        
        pic <- tsne$Y
        pic <- data.frame(pic,metadata_2$metadata_2)
        
        colnames(pic) <- c('X1','X2','Factor')
        
        ggplot(pic, aes(X1, X2,color = Factor))  +
          geom_point(size=0.5) +
          xlab(paste0("tsne1")) +
          ylab(paste0("tsne2"))+
          stat_ellipse(type = "t") + 
          coord_fixed()+ 
          ggtitle(colnames(metadata)[b])+
          theme(
            axis.title.x = element_text( size=7),
            axis.title.y = element_text( size=7),
            legend.text = element_text(size=7),
            legend.title = element_text(size=7),
            plot.title = element_text(hjust = 0.5, size = 7)
          ) + theme_bw() 
        ggsave(paste('tsne_',type1_all[a],colnames(metadata)[b],'.pdf', sep='_') , width=4, height=4)
        
        # difference analysis 
        keep <- metadata_2$metadata_2 == as.character(factor_level[1])
        if (length(factor_level) == 2 & sum(keep) >=5 & sum(!keep) >=5) {
          das = dif_abundance3(reads_table_3, metadata_2$metadata_2) # differential analysis
          if (!is.na(das)) {
            write.csv(das$data,file = paste('Abundance_',type1_all[a],colnames(metadata)[b],'.csv', sep='_'))
          }
        }
      }
    }
    
    correlation_all = rbind(correlation_all,correlation)
    
    
    
  }
  
  
}
correlation_all = correlation_all[!duplicated(correlation_all),]
write.csv(correlation_all,'impact_of_factor_baby_to_baby_microbiome_beta.csv',  row.names = F)










##### mom factors associated with Shannon index of baby's microbiome day 0 #####
setwd('/Users/binzhu/Desktop/Mom-baby/')
colnames(metadata_16s)[colnames(metadata_16s) == 'mom_delivery_method (C-section yes / vaginal no)'] = 'mom_delivery_method'

type1_all = unique(metadata_16s$SampleType)[c(1,5)]

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth) & (metadata_16s$days_rel2birth == 0)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    Participant_list = metadata_baby$ParticipantID
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby = metadata_baby[keep,]
    reads_table_baby = reads_table_baby[,keep]
    
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]    # reads table
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    metadata_mom = metadata_mom[keep,]
    
    sample_pair_id = cbind(metadata_mom$ParticipantID ,metadata_mom$SampleID, metadata_baby$SampleID)
    colnames(sample_pair_id) = c('ParticipantID','Mom_SampleID','Baby_SampleID')
    write.csv(sample_pair_id,paste0('sample_pair_id_',type1_all[a],'.csv'))
    
    # get alpha diversity
    reads_table_baby_2 = t(reads_table_baby)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness <- as.matrix(alpha.shannon_diversity/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
  }
  
  # get numeric and character factors 
  {
    metadata_num = as.data.frame(matrix(data = NA, nrow = nrow(metadata_mom), ncol =0))
    metadata_cha = as.data.frame(matrix(data = NA, nrow = nrow(metadata_mom), ncol =0))
    n =1
    m =1
    metadata_mom[,1:ncol(metadata_mom)] <- lapply(metadata_mom[,1:ncol(metadata_mom)],as.character)
    
    c2 = c('0','1','2','3','4','5','6','7','8','9','.')
    
    for (a2 in 1:ncol(metadata_mom)) {
      b2 = metadata_mom[,a2]
      b2 = b2[!is.na(b2)]
      b2 = strsplit(b2,'*')
      b2 = unlist(b2)
      b2 = unique(b2)
      keep = b2 %in% c2
      
      if (sum(!keep) == 0) {
        metadata_mom[,a2] <- as.numeric(metadata_mom[,a2])
        metadata_num <- cbind(metadata_num,metadata_mom[,a2] )
        colnames(metadata_num)[n] = colnames(metadata_mom)[a2]
        n=n+1
      } else {
        metadata_cha <- cbind(metadata_cha,metadata_mom[,a2])
        colnames(metadata_cha)[m] = colnames(metadata_mom)[a2]
        m=m+1
      }
    }
    
  }
  
  # case numbers
  case_number_cha = as.data.frame(matrix(data = NA, nrow =7, ncol = ncol(metadata_cha)))
  colnames(case_number_cha) = colnames(metadata_cha)
  row.names(case_number_cha) = c("BCKD_y","BCKD_n", "BRCD_y", "BRCD_n","BS1D_y","BS1D_n",'Levels')
  
  case_number_num = as.data.frame(matrix(data = NA, nrow =3, ncol = ncol(metadata_num)))
  colnames(case_number_num) = colnames(metadata_num)
  row.names(case_number_num) = c("BCKD","BRCD","BS1D")
  ########### character factor ##########
  metadata = metadata_cha
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 8))
  colnames(correlation) = c('Type','Factor','shannon_pvalue','shannon_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)',
                            'evenness_pvalue','evenness_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)','otu_pvalue','otu_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)')
  
  correlation$Type = type1_all[a]
  
  
  for (b in 1: ncol(metadata)) {
    correlation$Factor[b] <- colnames(metadata)[b]
    
    keep = !is.na(metadata[,b])
    
    metadata_2 = metadata[keep,b]
    shannon_2 = alpha.shannon_diversity[keep]
    evenness_2 = alpha.evenness[keep]
    otu_2 = alpha.ovserved_OTU[keep]
    reads_table_5 = reads_table_baby[,keep]
    
    ###
    factor_level = unique(metadata_2)
    factor_level = factor_level[order(factor_level, decreasing = T)]
    
    case_number_cha[7,b] = length(factor_level)
    
    factor_number = factor_level
    for (b3 in 1: length(factor_level)) {
      factor_number[b3] = sum(metadata_2 == factor_level[b3])
    }
    factor_number = as.numeric(factor_number)
    factor_number = factor_number[order(factor_number,decreasing = T)]
    
    # get case numbers
    case_number_cha[a*2-1,b] = factor_number[2]
    case_number_cha[a*2,b] = factor_number[1]
    
    if (sum(keep) <= 10) {
      next
    }
    
    if (length(factor_level) != 2) {
      next
    }
    
    if (factor_number[1] < 5 | factor_number[2] < 5) { 
      next
    }
    
    # calculate ratio of quantile
    data1 = shannon_2[metadata_2 == factor_level[1]]
    data2 = shannon_2[metadata_2 == factor_level[2]]
    
    quantile1 = quantile(data1)
    quantile2 = quantile(data2)
    quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
    quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
    quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
    quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
    quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
    quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
    correlation$`shannon_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
    
    data1 = evenness_2[metadata_2 == factor_level[1]]
    data2 = evenness_2[metadata_2 == factor_level[2]]
    quantile1 = quantile(data1)
    quantile2 = quantile(data2)
    quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
    quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
    quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
    quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
    quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
    quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
    correlation$`evenness_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
    
    data1 = otu_2[metadata_2 == factor_level[1]]
    data2 = otu_2[metadata_2 == factor_level[2]]
    quantile1 = quantile(data1)
    quantile2 = quantile(data2)
    quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
    quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
    quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
    quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
    quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
    quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
    correlation$`otu_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
    
    # 
    temp = as.data.frame(cbind(shannon_2,metadata_2))
    colnames(temp) = c('V1','V2')
    correlation$shannon_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
    
    if (correlation$shannon_pvalue[b] <= 0.05) {
      temp$V1 = as.numeric(as.character(temp$V1))
      ggplot(data=temp, aes(x=V2, y=V1)) +
        geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
        geom_jitter(shape=16, size = 0.5)+
        labs(x = correlation$Factor[b], y = paste0("Shannon index",'\n','of ',type1_all[a]))+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw()
      ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
      
      if (correlation$Factor[b] == 'mom_time_in_hospital') {
        ggplot(data=temp, aes(x=V2, y=V1)) +
          geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
          geom_jitter(shape=16, size = 0.5)+
          labs(x = correlation$Factor[b], y = paste0("Shannon index",'\n','of ',type1_all[a]))+ 
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
        ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=5)
      }
    }
    
    temp = as.data.frame(cbind(evenness_2,metadata_2))
    colnames(temp) = c('V1','V2')
    correlation$evenness_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
    
    if (correlation$evenness_pvalue[b] <= 0.05) {
      temp$V1 = as.numeric(as.character(temp$V1))
      ggplot(data=temp, aes(x=V2, y=V1)) +
        geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
        geom_jitter(shape=16, size = 0.5)+
        labs(x = correlation$Factor[b], y = "Evenness")+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
      ggsave(paste0('alpha_evenness_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
    }
    
    temp = as.data.frame(cbind(otu_2,metadata_2))
    colnames(temp) = c('V1','V2')
    correlation$otu_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
    
    if (correlation$otu_pvalue[b] <= 0.05) {
      temp$V1 = as.numeric(as.character(temp$V1))
      ggplot(data=temp, aes(x=V2, y=V1)) +
        geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
        geom_jitter(shape=16, size = 0.5)+
        labs(x = correlation$Factor[b], y = "Observed taxa")+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
      ggsave(paste0('alpha_otu_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
    }
    
    #      if (correlation$shannon_pvalue[b] <= 0.05 & length(factor_level) == 2) {
    #        das = dif_abundance3(reads_table_5,metadata_2) # differential analysis
    #        if (!is.na(das)) {
    #          write.csv(das$data,file = paste0('das_',type1_all[a],'_',colnames(metadata)[b],'.csv'))
    #      }
    
  }
  
  
  write.csv(case_number_cha,paste0('case_number_cha_',type1_all[a],'.csv'))
  
  
  
  if (a ==1) {
    correlation_all_cha = correlation
  } else {
    correlation_all_cha = rbind(correlation_all_cha,correlation)
  }
  
  
  
  
  
  
  
  
  
  
  
  ########### numeric factor ##########
  metadata = metadata_num
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 8))
  colnames(correlation) = c('Type','Factor','shannon_pvalue','shannon_lm_function',
                            'evenness_pvalue','evenness_lm_function','otu_pvalue','otu_lm_function')
  
  correlation$Type = type1_all[a]
  
  
  for (b in 1: ncol(metadata)) {
    correlation$Factor[b] <- colnames(metadata)[b]
    
    keep = !is.na(metadata[,b])
    
    # get case numbers
    case_number_num[a,b] = sum(keep)
    
    #
    if (sum(keep) <= 10) {
      next
    }
    
    metadata_2 = metadata[keep,b]
    shannon_2 = alpha.shannon_diversity[keep]
    evenness_2 = alpha.evenness[keep]
    otu_2 = alpha.ovserved_OTU[keep]
    reads_table_5 = reads_table_baby[,keep]
    
    ###
    factor_level = unique(metadata_2)
    
    if (length(factor_level) < 3) {
      next
    }
    
    # lm model to test shannon correlation
    temp = as.data.frame(cbind(shannon_2,metadata_2))
    colnames(temp) = c('V1','V2')
    temp$V1 = as.numeric(as.character(temp$V1))
    
    # remove outliers
    outliers_test = boxplot(temp$V2, plot=FALSE)$out
    
    if (length(outliers_test) > 0) {
      n = !(temp$V2 %in% outliers_test)
      temp = temp[n,]
    }
    
    if (length(unique(temp$V2)) > 3) {
      
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$shannon_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$shannon_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$shannon_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      if (correlation$shannon_pvalue[b] <= 0.05) {
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Shannon index of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        
        ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
    }
    
    
    
    # lm model to test evenness correlation
    temp = as.data.frame(cbind(evenness_2,metadata_2))
    colnames(temp) = c('V1','V2')
    
    temp$V1 = as.numeric(as.character(temp$V1))
    ##### remove outliers #####
    outliers_test = boxplot(temp$V2, plot=FALSE)$out
    
    if (length(outliers_test) > 0) {
      n = !(temp$V2 %in% outliers_test)
      temp = temp[n,]
    }
    
    if (length(unique(temp$V2)) > 3) {
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$evenness_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$evenness_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$evenness_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      
      if (correlation$evenness_pvalue[b] <= 0.05) {
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Evenness of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        ggsave(paste0('alpha_evenness_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
    }
    
    
    # lm model to test otu correlation
    temp = as.data.frame(cbind(otu_2,metadata_2))
    colnames(temp) = c('V1','V2')
    
    temp$V1 = as.numeric(as.character(temp$V1))
    ##### remove outliers #####
    outliers_test = boxplot(temp$V2, plot=FALSE)$out
    
    if (length(outliers_test) > 0) {
      n = !(temp$V2 %in% outliers_test)
      temp = temp[n,]
    }
    
    if (length(unique(temp$V2)) > 3) {
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$otu_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$otu_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$otu_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      if (correlation$otu_pvalue[b] <= 0.05) {
        
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Observed taxa of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        ggsave(paste0('alpha_otu_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
    }
    
  }
  
  
  if (a ==1) {
    correlation_all_num = correlation
  } else {
    correlation_all_num = rbind(correlation_all_num,correlation)
  }
  
  
  write.csv(case_number_num,paste0('case_number_num_',type1_all[a],'.csv'))
}
write.csv(correlation_all_num,'impact_of_factor_num_to_baby_microbiome_alpha.csv',  row.names = F)
write.csv(correlation_all_cha,'impact_of_factor_cha_to_baby_microbiome_alpha.csv',  row.names = F)








##### baby factors associated with Shannon index of baby's microbiome day 0 #####
setwd('/Users/binzhu/Desktop/Mom-baby/')
type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth) & (metadata_16s$days_rel2birth == 0)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    Participant_list = metadata_baby$ParticipantID
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby = metadata_baby[keep,]
    reads_table_baby = reads_table_baby[,keep]
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    metadata_mom = metadata_mom[keep,]
    
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    
    sample_pair_id = cbind(metadata_mom$ParticipantID ,metadata_mom$SampleID, metadata_baby$SampleID)
    colnames(sample_pair_id) = c('ParticipantID','Mom_SampleID','Baby_SampleID')
    write.csv(sample_pair_id,paste0('sample_pair_id_',type1_all[a],'.csv'))
    
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]    # reads table
    
    # get alpha diversity
    reads_table_baby_2 = t(reads_table_baby)
    reads_table_baby_2 = Rarefy(reads_table_baby_2, depth = min(rowSums(reads_table_baby_2)))
    reads_table_baby_2 <- reads_table_baby_2$otu.tab.rff
    reads_table_baby_2 <- as.data.frame(reads_table_baby_2)
    alpha.shannon_diversity <- as.numeric(diversity(reads_table_baby_2))   #  Shannon
    alpha.evenness <- as.matrix(alpha.shannon_diversity/log(specnumber(reads_table_baby_2)))
    alpha.ovserved_OTU <- as.matrix(data.frame(colSums(t(reads_table_baby_2) != 0)))
  }
  
  # get numeric and character factors 
  {
    metadata_num = as.data.frame(matrix(data = NA, nrow = nrow(metadata_baby), ncol =0))
    metadata_cha = as.data.frame(matrix(data = NA, nrow = nrow(metadata_baby), ncol =0))
    n =1
    m =1
    metadata_baby[,1:ncol(metadata_baby)] <- lapply(metadata_baby[,1:ncol(metadata_baby)],as.character)
    
    c2 = c('0','1','2','3','4','5','6','7','8','9','.')
    
    for (a2 in 1:ncol(metadata_baby)) {
      b2 = metadata_baby[,a2]
      b2 = b2[!is.na(b2)]
      b2 = strsplit(b2,'*')
      b2 = unlist(b2)
      b2 = unique(b2)
      keep = b2 %in% c2
      
      if (sum(!keep) == 0) {
        metadata_baby[,a2] <- as.numeric(metadata_baby[,a2])
        metadata_num <- cbind(metadata_num,metadata_baby[,a2] )
        colnames(metadata_num)[n] = colnames(metadata_baby)[a2]
        n=n+1
      } else {
        metadata_cha <- cbind(metadata_cha,metadata_baby[,a2])
        colnames(metadata_cha)[m] = colnames(metadata_baby)[a2]
        m=m+1
      }
    }
    
  }
  
  ########### character factor ##########
  # case numbers
  case_number_cha = as.data.frame(matrix(data = NA, nrow =7, ncol = ncol(metadata_cha)))
  colnames(case_number_cha) = colnames(metadata_cha)
  row.names(case_number_cha) = c("BCKD_y","BCKD_n", "BRCD_y", "BRCD_n","BS1D_y","BS1D_n",'Levels')
  
  metadata = metadata_cha
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 8))
  colnames(correlation) = c('Type','Factor','shannon_pvalue','shannon_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)',
                            'evenness_pvalue','evenness_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)','otu_pvalue','otu_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)')
  
  correlation$Type = type1_all[a]
  
  try(
    for (b in 1: ncol(metadata)) {
      correlation$Factor[b] <- colnames(metadata)[b]
      
      keep = !is.na(metadata[,b])
      
      metadata_2 = metadata[keep,b]
      shannon_2 = alpha.shannon_diversity[keep]
      evenness_2 = alpha.evenness[keep]
      otu_2 = alpha.ovserved_OTU[keep]
      reads_table_5 = reads_table_baby[,keep]
      
      ###
      factor_level = unique(metadata_2)
      factor_level = factor_level[order(factor_level, decreasing = T)]
      
      case_number_cha[7,b] = length(factor_level)
      
      factor_number = factor_level
      for (b3 in 1: length(factor_level)) {
        factor_number[b3] = sum(metadata_2 == factor_level[b3])
      }
      factor_number = as.numeric(factor_number)
      factor_number = factor_number[order(factor_number,decreasing = T)]
      
      # get case numbers
      case_number_cha[a*2-1,b] = factor_number[2]
      case_number_cha[a*2,b] = factor_number[1]
      
      if (sum(keep) <= 10) {
        next
      }
      
      if (length(factor_level) != 2) {
        next
      }
      
      if (factor_number[1] < 5 | factor_number[2] < 5) { 
        next
      }
      
      # calculate ratio of quantile
      data1 = shannon_2[metadata_2 == factor_level[1]]
      data2 = shannon_2[metadata_2 == factor_level[2]]
      quantile1 = quantile(data1)
      quantile2 = quantile(data2)
      quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
      quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
      quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
      quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
      quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
      quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
      correlation$`shannon_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
      
      data1 = evenness_2[metadata_2 == factor_level[1]]
      data2 = evenness_2[metadata_2 == factor_level[2]]
      quantile1 = quantile(data1)
      quantile2 = quantile(data2)
      quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
      quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
      quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
      quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
      quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
      quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
      correlation$`evenness_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
      
      data1 = otu_2[metadata_2 == factor_level[1]]
      data2 = otu_2[metadata_2 == factor_level[2]]
      quantile1 = quantile(data1)
      quantile2 = quantile(data2)
      quantile_25_ratio = as.numeric(quantile1[2]/quantile2[2])
      quantile_25_ratio = format(round(quantile_25_ratio, 3), nsmall = 3)
      quantile_50_ratio = as.numeric(quantile1[3]/quantile2[3])
      quantile_50_ratio = format(round(quantile_50_ratio, 3), nsmall = 3)
      quantile_75_ratio = as.numeric(quantile1[4]/quantile2[4])
      quantile_75_ratio = format(round(quantile_75_ratio, 3), nsmall = 3)
      correlation$`otu_Ratio of "yes" versus "No", median (25th empirical quartile, 75th empirical quartile)`[b] = paste0(quantile_50_ratio, ' (',quantile_25_ratio,', ', quantile_75_ratio,')')
      
      # 
      temp = as.data.frame(cbind(shannon_2,metadata_2))
      colnames(temp) = c('V1','V2')
      correlation$shannon_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
      
      if (correlation$shannon_pvalue[b] <= 0.05) {
        temp$V1 = as.numeric(as.character(temp$V1))
        ggplot(data=temp, aes(x=V2, y=V1)) +
          geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
          geom_jitter(shape=16, size = 0.5)+
          labs(x = correlation$Factor[b], y = paste0("Shannon index",'\n','of ',type1_all[a]))+ 
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw()
        ggsave(paste0('alpha_shannon_baby_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
        
        if (correlation$Factor[b] == 'mom_time_in_hospital') {
          ggplot(data=temp, aes(x=V2, y=V1)) +
            geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
            geom_jitter(shape=16, size = 0.5)+
            labs(x = correlation$Factor[b], y = paste0("Shannon index",'\n','of ',type1_all[a]))+ 
            theme(axis.title = element_text(size = 7), 
                  axis.text = element_text(size = 7), 
                  legend.text = element_text(size = 7), 
                  legend.title = element_text(size = 7),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
          ggsave(paste0('alpha_shannon_baby_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=5)
        }
      }
      
      temp = as.data.frame(cbind(evenness_2,metadata_2))
      colnames(temp) = c('V1','V2')
      correlation$evenness_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
      
      if (correlation$evenness_pvalue[b] <= 0.05) {
        temp$V1 = as.numeric(as.character(temp$V1))
        ggplot(data=temp, aes(x=V2, y=V1)) +
          geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
          geom_jitter(shape=16, size = 0.5)+
          labs(x = correlation$Factor[b], y = "Evenness")+ 
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
        ggsave(paste0('alpha_evenness_baby_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
      }
      
      temp = as.data.frame(cbind(otu_2,metadata_2))
      colnames(temp) = c('V1','V2')
      correlation$otu_pvalue[b] = kruskal.test(V1~V2, temp)$p.value
      
      if (correlation$otu_pvalue[b] <= 0.05) {
        temp$V1 = as.numeric(as.character(temp$V1))
        ggplot(data=temp, aes(x=V2, y=V1)) +
          geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
          geom_jitter(shape=16, size = 0.5)+
          labs(x = correlation$Factor[b], y = "Observed taxa")+ 
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
        ggsave(paste0('alpha_otu_baby_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=1.3, height=3)
      }
      
      #      if (correlation$shannon_pvalue[b] <= 0.05 & length(factor_level) == 2) {
      #        das = dif_abundance3(reads_table_5,metadata_2) # differential analysis
      #        if (!is.na(das)) {
      #          write.csv(das$data,file = paste0('das_',type1_all[a],'_',colnames(metadata)[b],'.csv'))
      #      }
      
    }
  )
  
  write.csv(case_number_cha,paste0('case_number_baby_cha_',type1_all[a],'.csv'))
  
  
  
  if (a ==1) {
    correlation_all_cha = correlation
  } else {
    correlation_all_cha = rbind(correlation_all_cha,correlation)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  ########### numeric factor ##########
  metadata = as.data.frame(cbind(metadata_num$height, metadata_num$weight, metadata_num$BMI, metadata_num$pulse, metadata_num$days_rel2birth))
  colnames(metadata) = c('baby_height','baby_weight','baby_BMI','baby_pulse','Sample_collection_time')
  
  # case numbers
  case_number_num = as.data.frame(matrix(data = NA, nrow =3, ncol = ncol(metadata)))
  colnames(case_number_num) = colnames(metadata)
  row.names(case_number_num) = c("BCKD","BRCD","BS1D")
  
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 8))
  colnames(correlation) = c('Type','Factor','shannon_pvalue','shannon_lm_function',
                            'evenness_pvalue','evenness_lm_function','otu_pvalue','otu_lm_function')
  
  correlation$Type = type1_all[a]
  
  try(
    for (b in 1: ncol(metadata)) {
      correlation$Factor[b] <- colnames(metadata)[b]
      
      keep = !is.na(metadata[,b])
      
      # get case numbers
      case_number_num[a,b] = sum(keep)
      
      #
      if (sum(keep) <= 10) {
        next
      }
      
      metadata_2 = metadata[keep,b]
      shannon_2 = alpha.shannon_diversity[keep]
      evenness_2 = alpha.evenness[keep]
      otu_2 = alpha.ovserved_OTU[keep]
      reads_table_5 = reads_table_baby[,keep]
      
      ###
      factor_level = unique(metadata_2)
      
      if (length(factor_level) < 3) {
        next
      }
      
      # lm model to test shannon correlation
      temp = as.data.frame(cbind(shannon_2,metadata_2))
      colnames(temp) = c('V1','V2')
      
      temp$V1 = as.numeric(as.character(temp$V1))
      
      ##### remove outliers #####
      outliers_test = boxplot(temp$V2, plot=FALSE)$out
      
      if (length(outliers_test) > 0) {
        n = !(temp$V2 %in% outliers_test)
        temp = temp[n,]
      }
      
      #####
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$shannon_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$shannon_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$shannon_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      if (correlation$shannon_pvalue[b] <= 0.05) {
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Shannon index of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        
        ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
      
      
      
      # lm model to test evenness correlation
      temp = as.data.frame(cbind(evenness_2,metadata_2))
      colnames(temp) = c('V1','V2')
      temp$V1 = as.numeric(as.character(temp$V1))
      
      ##### remove outliers #####
      outliers_test = boxplot(temp$V2, plot=FALSE)$out
      
      if (length(outliers_test) > 0) {
        n = !(temp$V2 %in% outliers_test)
        temp = temp[n,]
      }
      
      #####
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$evenness_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$evenness_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$evenness_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      
      if (correlation$evenness_pvalue[b] <= 0.05) {
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Evenness of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        ggsave(paste0('alpha_evenness_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
      
      
      
      # lm model to test otu correlation
      temp = as.data.frame(cbind(otu_2,metadata_2))
      colnames(temp) = c('V1','V2')
      temp$V1 = as.numeric(as.character(temp$V1))
      
      ##### remove outliers #####
      outliers_test = boxplot(temp$V2, plot=FALSE)$out
      
      if (length(outliers_test) > 0) {
        n = !(temp$V2 %in% outliers_test)
        temp = temp[n,]
      }
      
      #####
      linearMod <- lm(V1 ~ V2, data=temp)
      linearMod = summary(linearMod)
      r_squared <- linearMod$r.squared
      pvalue <- as.numeric(linearMod$coefficients[,4][2])
      r_squared = format(round(r_squared, 5), nsmall = 5)
      pvalue = format(round(pvalue, 5), nsmall = 5)
      
      correlation$otu_pvalue[b] = pvalue
      
      if (linearMod$coefficients[2,1]>0) {
        correlation$otu_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),'+',format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      } else {
        correlation$otu_lm_function[b] = paste0('y=',format(round(linearMod$coefficients[1,1], 3), nsmall = 3),format(round(linearMod$coefficients[2,1], 3), nsmall = 3),'x')
      }
      
      if (correlation$otu_pvalue[b] <= 0.05) {
        temp$V1 = as.numeric(as.character(temp$V1))
        
        ggplot(temp, aes(x = V2, y = V1)) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue")+
          labs(x = correlation$Factor[b], 
               y = paste0("Observed taxa of ",type1_all[a]))+ 
          ggtitle(paste0("R = ", r_squared, '\n',"pvalue = ", correlation$shannon_pvalue[b]))  +
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw()
        ggsave(paste0('alpha_otu_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      }
      
    }
  )
  
  if (a ==1) {
    correlation_all_num = correlation
  } else {
    correlation_all_num = rbind(correlation_all_num,correlation)
  }
  
  
  write.csv(case_number_num,paste0('case_number_baby_',type1_all[a],'.csv'))
}
write.csv(correlation_all_num,'impact_of_factor_baby_num_to_baby_microbiome_alpha.csv',  row.names = F)
write.csv(correlation_all_cha,'impact_of_factor_baby_cha_to_baby_microbiome_alpha.csv',  row.names = F)










##### mom's factor associated with baby's beta diversity day 0 #####
setwd('/Users/binzhu/Desktop/Mom-baby/')
type1_all = unique(metadata_16s$SampleType)[c(1,5)]
correlation_all = as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
colnames(correlation_all) = c('Type','Factor','Adonis_pvalue','Adj_pvalue')

for (a in 1: length(type1_all)) {
  # find paired mom and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth) & (metadata_16s$days_rel2birth == 0)
    sum(keep)
    
    metadata_baby = metadata_16s[keep,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]
    
    keep = colnames(reads_table_16s) %in% metadata_baby$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]
    
    Participant_list = metadata_baby$ParticipantID
    metadata_mom = metadata_16s[metadata_16s$ParticipantID %in% Participant_list & metadata_16s$SampleType == 'MV1D',]
    metadata_mom = metadata_mom[metadata_mom$Mombaby == 'Mom',]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    keep = metadata_baby$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby = metadata_baby[keep,]
    reads_table_baby = reads_table_baby[,keep]
    
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata
    reads_table_baby = reads_table_baby[metadata_baby$SampleID]    # reads table
    
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
    keep = metadata_baby$SampleID %in% colnames(reads_table_baby)
    metadata_baby = metadata_baby[keep,]
    metadata_mom = metadata_mom[keep,]
  }
  
  metadata = metadata_mom
  reads_table = reads_table_baby
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata), ncol = 4))
  colnames(correlation) = c('Type','Factor','Adonis_pvalue','Adj_pvalue')
  
  correlation$Type = type1_all[a]
  
  for (b in 1: ncol(metadata)) {
    ##### prepare data #####
    # correlation by adonis2
    correlation$Factor[b] <- colnames(metadata)[b]
    
    keep = !is.na(metadata[,b])
    if (sum(keep) <= 10) {
      next
    }
    
    metadata_2 = metadata[keep,b]
    reads_table_2 = reads_table[,keep]
    reads_table_3 = reads_table[,keep]
    
    factor_level = unique(metadata_2)
    if (length(factor_level) < 2) {
      next
    }
    
    ##### convert numeric factor to charactor factor by grouping to two levels #####
    c2 = c('0','1','2','3','4','5','6','7','8','9','.')
    
    b2 = factor_level
    b2 = b2[!is.na(b2)]
    b2 = as.character(b2)
    b2 = strsplit(b2,'*')
    b2 = unlist(b2)
    b2 = unique(b2)
    keep = b2 %in% c2
    
    if (sum(!keep) == 0) {
      factor_level = as.numeric(as.character(factor_level))
      b2 = median(factor_level)
      keep = factor_level>=b2
      factor_level[keep] ='High'
      factor_level[!keep] ='Low'
      metadata_2 = as.numeric(metadata_2)
      keep = metadata_2>=b2
      metadata_2[keep] = 'High'
      metadata_2[!keep] = 'Low'
      metadata_2 = factor(metadata_2,levels = c('Low','High'))
      factor_level = unique(factor_level)
    } else {
      metadata_2 = factor(metadata_2,levels = c('No','Yes'))
    }
    
    if (sum(!is.na(metadata_2)) < 10) {
      next
    }
    
    factor_number = factor_level
    for (b3 in 1: length(factor_level)) {
      factor_number[b3] = sum(metadata_2 == factor_level[b3])
    }
    factor_number = as.numeric(factor_number)
    factor_number = factor_number[order(factor_number,decreasing = T)]
    
    if (sum(is.na(factor_number)) > 0) {
      next
    }
    
    if (factor_number[1] < 5 | factor_number[2] < 5) { 
      next
    } 
    
    
    ##### t-SNE #####
    reads_table_2 <- t(reads_table_2)
    reads_table_2 = Rarefy(reads_table_2, depth = min(rowSums(reads_table_2)))
    reads_table_2 <- reads_table_2$otu.tab.rff
    reads_table_2 <- as.data.frame(reads_table_2)
    
    metadata_2 = as.data.frame(metadata_2)
    pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")  
    correlation$Adonis_pvalue[b] <- pvalue[1,5]
    
    if (!is.na(pvalue[1,5])) {
      if (pvalue[1,5] <= 0.05) {
        
        
        # difference analysis 
        keep <- metadata_2$metadata_2 == as.character(factor_level[1])
        if (length(factor_level) == 2 & sum(keep) >=5 & sum(!keep) >=5) {
          das = dif_abundance3(reads_table_3, metadata_2$metadata_2) # differential analysis
          if (!is.na(das)) {
            write.csv(das$data,file = paste('Abundance_',type1_all[a],colnames(metadata)[b],'.csv', sep='_'))
          }
        }
      }
    }
    
    correlation_all = rbind(correlation_all,correlation)
    
    
    
  }
  
  
}
correlation_all = correlation_all[!duplicated(correlation_all),]
write.csv(correlation_all,'impact_of_factor_to_baby_microbiome_beta.csv',  row.names = F)










##### clustering of taxa by relative abundance change #####
setwd('/Users/binzhu/Desktop/Mom-baby/beta_abundance')

file_list = list.files(path = ".", pattern = '*.csv')

for (a in 1:length(file_list)) {
  data = read.csv(file_list[a])
  
  if (a ==1) {
    taxa_list = data$X
  } else {
    taxa_list = unique(c(taxa_list, data$X))
  }
}
taxa_list = taxa_list[order(taxa_list)]

reads_table = as.data.frame(matrix(data = NA, ncol = length(file_list), nrow = length(taxa_list)))
row.names(reads_table) = taxa_list
colnames(reads_table) = str_remove_all(file_list,'_.csv')
colnames(reads_table) = str_remove_all(colnames(reads_table),'Abundance__')

reads_table_2 = as.data.frame(matrix(data = NA, ncol = length(file_list), nrow = length(taxa_list)))
row.names(reads_table_2) = taxa_list
colnames(reads_table_2) = str_remove_all(file_list,'_.csv')
colnames(reads_table_2) = str_remove_all(colnames(reads_table_2),'Abundance__')

for (a in 1:length(file_list)) {
  data = read.csv(file_list[a])
  n= !(taxa_list %in% data$X)
  
  data_2 = as.data.frame(matrix(data = NA, ncol = ncol(data), nrow = sum(n)))
  colnames(data_2) = colnames(data)
  data_2$X = taxa_list[n]
  
  data = rbind(data,data_2)
  
  data = data[order(data$X),]
  
  reads_table[,a] = data$diff.btw
  reads_table_2[,a] = data$wi.eBH
}

write.csv(reads_table, 'Abundance_change_beta_change.csv')
write.csv(reads_table_2, 'Abundance_change_beta_change_2.csv')

bar_max = 2
bar_min = -2
pheatmap_fontsize = 5
treeheight = 20
alpha = 0.05
paletteLength <- 50

p.yes.rr = reads_table
p.yes.rr[is.na(p.yes.rr)] = 0

p.yes.rr_2 = abs(p.yes.rr)
p.yes.rr = p.yes.rr[order(rowSums(p.yes.rr_2), decreasing = T),]
p.yes.rr = p.yes.rr[c(1:50),]

myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
              seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))

myBreaks <- unique(myBreaks)



pdf('Abundance_change_beta_change.pdf',width=4, height=6)
pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
         treeheight_row = treeheight, treeheight_col = treeheight, clustering_method = 'ward.D2')
dev.off()

setwd('/Users/binzhu/Desktop/Mom-baby')



##### ??? baby microboime and mom gene expression, paired mom-baby 16s samples (time closest between mom and baby) in a participantID ##########
# get human transcriptomic data from 'human_gene_expression_raw_reads.R'

type1_all = unique(metadata_16s$SampleType)[4:6]

for (type_num in 1: length(type1_all)) {
  ##### get metadata and paired reads tables (time closest) ####
  keep = metadata_16s$SampleType == type1_all[type_num] & 
    !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
  sum(keep)
  metadata_baby = metadata_16s[keep,]
  
  metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
  metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]
  keep = duplicated(metadata_baby$ParticipantID)
  metadata_baby = metadata_baby[!keep,]
  
  metadata_mom = metadata_16s[metadata_16s$SampleType == 'MV1D' & !is.na(metadata_16s$SampleType),]
  metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
  metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
  keep = duplicated(metadata_mom$ParticipantID)
  metadata_mom = metadata_mom[!keep,]
  
  keep = (metadata_mom$ParticipantID %in% metadata_baby$ParticipantID) & 
    (metadata_mom$SampleID %in% colnames(reads_table_human_MTG))
  sum(keep)
  
  metadata_mom = metadata_mom[keep,]
  reads_table_human_MTG_16s = reads_table_human_MTG[,colnames(reads_table_human_MTG) %in% metadata_mom$SampleID]
  reads_table_human_MTG_16s = reads_table_human_MTG_16s[metadata_mom$SampleID]
  
  metadata_baby = metadata_baby[metadata_baby$ParticipantID %in% metadata_mom$ParticipantID,]
  reads_table_16s_human_MTG = reads_table_16s[,colnames(reads_table_16s) %in% metadata_baby$SampleID]
  reads_table_16s_human_MTG = reads_table_16s_human_MTG[metadata_baby$SampleID]
  
  # normalization
  reads_table_16s_human_MTG = as.data.frame(t(reads_table_16s_human_MTG))
  reads_table_16s_human_MTG = reads_table_16s_human_MTG + 0.5
  reads_table_16s_human_MTG  <- clr(reads_table_16s_human_MTG)
  reads_table_16s_human_MTG = as.data.frame(t(reads_table_16s_human_MTG))
  
  reads_table_human_MTG_16s = as.matrix(reads_table_human_MTG_16s)
  reads_table_human_MTG_16s = varianceStabilizingTransformation(reads_table_human_MTG_16s) # sample in column    reads = FPKM*colSum/10e6
  reads_table_human_MTG_16s = as.data.frame(reads_table_human_MTG_16s)
  colnames(reads_table_human_MTG_16s) = colnames(reads_table_16s_human_MTG)
  
  ###### correlation get data #####
  cor_parameter= 0.8
  gephi_p.yes.rr_minimum_link_taxa_to_gene = 5
  
  pvalue = 0.05
  pheatmap_fontsize = 5
  alpha = 0.05
  
  
  reads_table_human_MTG_16s_2 = as.matrix(reads_table_human_MTG_16s)
  reads_table_16s_human_MTG = as.matrix(reads_table_16s_human_MTG)
  
  reads_table_human_MTG_16s_2_1 = reads_table_human_MTG_16s_2[1:5000,]
  reads_table_human_MTG_16s_2_2 = reads_table_human_MTG_16s_2[5001:10000,]
  reads_table_human_MTG_16s_2_3 = reads_table_human_MTG_16s_2[10001:15000,]
  reads_table_human_MTG_16s_2_4 = reads_table_human_MTG_16s_2[15001:21650,]
  rm(reads_table_human_MTG_16s_2)
  
  # loop1
  reads_table = rbind(reads_table_16s_human_MTG , reads_table_human_MTG_16s_2_1)
  reads_table = as.matrix(t(reads_table))
  
  otu.cor <- rcorr(reads_table, type="pearson")
  otu.pval <- forceSymmetric(otu.cor$P)
  otu.pval <- otu.pval@x
  otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  p.yes <- otu.pval< pvalue  
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
  p.yes.r <- p.yes.r*r.val # use logical vector for subscripting.
  
  gephi_p.yes.rr = as.matrix(p.yes.r)
  
  for (a in 1:nrow(p.yes.r)) {
    for (b in 1:ncol(p.yes.r)) {
      if (a >=b) {
        gephi_p.yes.rr[a,b] = NA
      }
    }
  }
  gephi_p.yes.rr <- as.data.frame(gephi_p.yes.rr)
  gephi_p.yes.rr = gather(gephi_p.yes.rr)
  gephi_p.yes.rr$Taxa = rep(row.names(p.yes.r),ncol(p.yes.r))
  gephi_p.yes.rr = gephi_p.yes.rr[!is.na(gephi_p.yes.rr[,2]),]
  gephi_p.yes.rr$value = as.numeric(as.character(gephi_p.yes.rr$value))
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr$value != 0,] 
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr[,2] >= cor_parameter | gephi_p.yes.rr[,2] <= -cor_parameter,]
  colnames(gephi_p.yes.rr) = c('Source','Weigth','Target')
  
  correlation_data = gephi_p.yes.rr
  
  # loop2
  reads_table = rbind(reads_table_16s_human_MTG , reads_table_human_MTG_16s_2_2)
  reads_table = as.matrix(t(reads_table))
  
  otu.cor <- rcorr(reads_table, type="pearson")
  otu.pval <- forceSymmetric(otu.cor$P)
  otu.pval <- otu.pval@x
  otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  p.yes <- otu.pval< pvalue  
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
  p.yes.r <- p.yes.r*r.val # use logical vector for subscripting.
  
  gephi_p.yes.rr = as.matrix(p.yes.r)
  
  for (a in 1:nrow(p.yes.r)) {
    for (b in 1:ncol(p.yes.r)) {
      if (a >=b) {
        gephi_p.yes.rr[a,b] = NA
      }
    }
  }
  gephi_p.yes.rr <- as.data.frame(gephi_p.yes.rr)
  gephi_p.yes.rr = gather(gephi_p.yes.rr)
  gephi_p.yes.rr$Taxa = rep(row.names(p.yes.r),ncol(p.yes.r))
  gephi_p.yes.rr = gephi_p.yes.rr[!is.na(gephi_p.yes.rr[,2]),]
  gephi_p.yes.rr$value = as.numeric(as.character(gephi_p.yes.rr$value))
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr$value != 0,] 
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr[,2] >= cor_parameter | gephi_p.yes.rr[,2] <= -cor_parameter,]
  colnames(gephi_p.yes.rr) = c('Source','Weigth','Target')
  
  correlation_data = rbind(correlation_data,gephi_p.yes.rr)
  
  # loop3
  reads_table = rbind(reads_table_16s_human_MTG , reads_table_human_MTG_16s_2_3)
  reads_table = as.matrix(t(reads_table))
  
  otu.cor <- rcorr(reads_table, type="pearson")
  otu.pval <- forceSymmetric(otu.cor$P)
  otu.pval <- otu.pval@x
  otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  p.yes <- otu.pval< pvalue  
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
  p.yes.r <- p.yes.r*r.val # use logical vector for subscripting.
  
  gephi_p.yes.rr = as.matrix(p.yes.r)
  
  for (a in 1:nrow(p.yes.r)) {
    for (b in 1:ncol(p.yes.r)) {
      if (a >=b) {
        gephi_p.yes.rr[a,b] = NA
      }
    }
  }
  gephi_p.yes.rr <- as.data.frame(gephi_p.yes.rr)
  gephi_p.yes.rr = gather(gephi_p.yes.rr)
  gephi_p.yes.rr$Taxa = rep(row.names(p.yes.r),ncol(p.yes.r))
  gephi_p.yes.rr = gephi_p.yes.rr[!is.na(gephi_p.yes.rr[,2]),]
  gephi_p.yes.rr$value = as.numeric(as.character(gephi_p.yes.rr$value))
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr$value != 0,] 
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr[,2] >= cor_parameter | gephi_p.yes.rr[,2] <= -cor_parameter,]
  colnames(gephi_p.yes.rr) = c('Source','Weigth','Target')
  
  correlation_data = rbind(correlation_data,gephi_p.yes.rr)
  
  # loop4
  reads_table = rbind(reads_table_16s_human_MTG , reads_table_human_MTG_16s_2_4)
  reads_table = as.matrix(t(reads_table))
  
  otu.cor <- rcorr(reads_table, type="pearson")
  otu.pval <- forceSymmetric(otu.cor$P)
  otu.pval <- otu.pval@x
  otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  p.yes <- otu.pval< pvalue  
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
  p.yes.r <- p.yes.r*r.val # use logical vector for subscripting.
  
  gephi_p.yes.rr = as.matrix(p.yes.r)
  
  for (a in 1:nrow(p.yes.r)) {
    for (b in 1:ncol(p.yes.r)) {
      if (a >=b) {
        gephi_p.yes.rr[a,b] = NA
      }
    }
  }
  gephi_p.yes.rr <- as.data.frame(gephi_p.yes.rr)
  gephi_p.yes.rr = gather(gephi_p.yes.rr)
  gephi_p.yes.rr$Taxa = rep(row.names(p.yes.r),ncol(p.yes.r))
  gephi_p.yes.rr = gephi_p.yes.rr[!is.na(gephi_p.yes.rr[,2]),]
  gephi_p.yes.rr$value = as.numeric(as.character(gephi_p.yes.rr$value))
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr$value != 0,] 
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr[,2] >= cor_parameter | gephi_p.yes.rr[,2] <= -cor_parameter,]
  colnames(gephi_p.yes.rr) = c('Source','Weigth','Target')
  
  correlation_data = rbind(correlation_data,gephi_p.yes.rr)
  
  # remove duplicate in correlation_data
  rm(p.yes.r,r.val,otu.cor,gephi_p.yes.rr,reads_table_human_MTG_16s_2_1,
     reads_table_human_MTG_16s_2_2,reads_table_human_MTG_16s_2_3,
     reads_table_human_MTG_16s_2_4)
  
  correlation_data$dup = paste0(correlation_data$Source,correlation_data$Target)
  keep = !duplicated(correlation_data$dup)
  sum(keep)
  correlation_data = correlation_data[keep,]
  correlation_data$dup = NULL
  
  correlation_data$Interaction = NA
  keep1 = correlation_data$Source %in% row.names(reads_table_16s_human_MTG)
  keep2 = correlation_data$Target %in% row.names(reads_table_16s_human_MTG)
  correlation_data$Interaction[keep1&keep2] = 'Taxa_Taxa'
  correlation_data$Interaction[(keep1&(!keep2)) | ((!keep1)&keep2)] = 'Taxa_Gene'
  correlation_data$Interaction[(!keep1)&(!keep2)] = 'Gene_Gene'
  
  #  keep = correlation_data$Interaction == 'Taxa_Gene'
  #  sum(keep)
  #  correlation_data = correlation_data[keep,]
  write.csv(correlation_data,paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_edge_all.csv'), row.names = F)
  
  # remove taxa with link less than gephi_p.yes.rr_minimum_link_taxa_to_gene
  unique_taxa = unique(correlation_data$Target)
  for (a in 1: length(unique_taxa)) {
    n = which(correlation_data$Target == unique_taxa[a])
    if (length(n) < gephi_p.yes.rr_minimum_link_taxa_to_gene) {
      correlation_data$Target[n] = NA
    }
  }
  correlation_data = correlation_data[!is.na(correlation_data$Target),]
  
  ###### correlation draw figures #####
  library(igraph)
  
  correlation_data_2 = correlation_data[correlation_data$Interaction == 'Taxa_Gene',]
  
  # get correlation matrix
  unique_gene = unique(correlation_data_2$Source)
  unique_taxa = unique(correlation_data_2$Target)
  
  p.yes.rr = as.data.frame(matrix(data = 0, nrow = length(unique_gene), ncol = length(unique_taxa) ))
  row.names(p.yes.rr) = unique_gene
  colnames(p.yes.rr) = unique_taxa
  
  for (a in 1: nrow(correlation_data_2)) {
    x = which(unique_gene == correlation_data_2$Source[a])
    y = which(unique_taxa == correlation_data_2$Target[a])
    
    p.yes.rr[x,y]= correlation_data_2$Weigth[a]
  }
  
  
  
  # heatmap
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  paletteLength <- 50
  
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
  p.yes.rr = as.data.frame(t(p.yes.rr))
  
  pdf(paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_ward.D.pdf'),width=6, height=4.5)
  cluster = pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, 
                     fontsize = 5,clustering_method = "ward.D")
  dev.off()
  
  # export taxa cluster
  taxa_cluster = sort(cutree(cluster$tree_row, k=4))
  write.csv(taxa_cluster,paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_taxa_cluster.csv'))
  
  detach("package:igraph", unload = TRUE)
  
  # output for gephi
  keep = correlation_data$Interaction == 'Taxa_Gene' | correlation_data$Interaction == 'Gene_Gene'
  correlation_data_3=correlation_data[keep,]
  
  keep1 = (correlation_data_3$Source %in% row.names(reads_table_16s_human_MTG)) |
    (correlation_data_3$Target %in% row.names(reads_table_16s_human_MTG))
  
  a= unique(correlation_data_2$Source)
  keep2 = (correlation_data_3$Source %in% a) & (correlation_data_3$Target %in% a)
  
  correlation_data_3=correlation_data_3[(keep1|keep2),]
  
  length(unique(correlation_data_2$Source))
  length(unique(correlation_data_3$Source))
  
  write.csv(correlation_data_3, paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_edge.csv'), row.names = F)
  
  node1 = as.data.frame(unique(correlation_data_3$Target))
  colnames(node1) = 'V1'
  node1$Type = NA
  node1$Size = NA
  node2 = as.data.frame(unique(correlation_data_3$Source))
  colnames(node2) = 'V1'
  node2$Type = NA
  node2$Size = NA
  node = rbind(node1,node2)
  colnames(node)[1] = 'Id'
  node = node[!duplicated(node$Id),]
  keep = node$Id %in% row.names(reads_table_16s_human_MTG)
  node$Type[keep] = 'Taxa'
  node$Type[!keep] = 'Gene'
  node$Size[keep] = 10
  node$Size[!keep] = 1
  
  write.csv(node,paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_node.csv'), row.names = F)
  
  # get gene cluster in network
  keep = correlation_data_3$Source %in% row.names(reads_table_16s) |
    correlation_data_3$Target %in% row.names(reads_table_16s)
  
  correlation_data_4 = correlation_data_3[!keep,]
  
  unique_gene = unique(c(correlation_data_4$Source,correlation_data_4$Target))
  
  p.yes.rr = as.data.frame(matrix(data = 0, nrow = length(unique_gene), ncol = length(unique_gene) ))
  row.names(p.yes.rr) = unique_gene
  colnames(p.yes.rr) = unique_gene
  
  for (a in 1: nrow(correlation_data_4)) {
    x = which(unique_gene == correlation_data_4$Source[a])
    y = which(unique_gene == correlation_data_4$Target[a])
    
    p.yes.rr[x,y]= correlation_data_4$Weigth[a]
    p.yes.rr[y,x]= correlation_data_4$Weigth[a]
  }
  
  
  # heatmap
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  paletteLength <- 50
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
  myBreaks = unique(myBreaks)
  
  if (sum(myBreaks < 0) == 0) {
    myColor <- colorRampPalette(c("white", "red"))(paletteLength)
    pdf(paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_gene_network.pdf'),width=18.3, height=18.3)
    cluster = pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, 
                       fontsize = 0.5,clustering_method = "ward.D")
    dev.off()
  } else if (sum(myBreaks > 0) == 0) {
    myColor <- colorRampPalette(c("blue", "white"))(paletteLength)
    pdf(paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_gene_network.pdf'),width=18.3, height=18.3)
    cluster = pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, 
                       fontsize = 0.5,clustering_method = "ward.D")
    dev.off()
  } else {
    pdf(paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_gene_network.pdf'),width=18.3, height=18.3)
    cluster = pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, 
                       fontsize = 0.5,clustering_method = "ward.D")
    dev.off()
  }
  
  
  # export gene cluster
  gene_cluster = sort(cutree(cluster$tree_col, k=7))
  write.csv(gene_cluster,paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_gene_cluster.csv'))
  
  # interaction within taxa 
  pdf(paste0('Mom_MTG_vs.baby_16s_',type1_all[type_num],'_taxa_network.pdf'),width=6, height=6)
  taxa_cor = newwork_rcorr(reads_table_16s, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 2, bar_min = -2, pheatmap_fontsize = 5, alpha = 0.05)
  dev.off()
  taxa_cor = taxa_cor$gephi_input
  taxa_cor = taxa_cor[taxa_cor$Source %in% node1$V1,]
  taxa_cor = taxa_cor[taxa_cor$Target %in% node1$V1,]
  
  
  write.csv(taxa_cor,'MTG_vs._16s_taxa_cor.csv', row.names = F)
  
  rm(correlation_data)
}


##### ??? baby microboime shannon index and mom gene expression, paired mom-baby 16s samples (time closest between mom and baby) in a participantID ##########
type1_all = unique(metadata_16s$SampleType)[4:6]

for (type_num in 1: length(type1_all)) {
  ##### get metadata and paired reads tables (time closest) ####
  keep = metadata_16s$SampleType == type1_all[type_num] & 
    !is.na(metadata_16s$SampleType) & !is.na(metadata_16s$days_rel2birth)
  sum(keep)
  metadata_baby = metadata_16s[keep,]
  
  metadata_baby = metadata_baby[order(metadata_baby$VisitNum, decreasing = F),]
  metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]
  keep = duplicated(metadata_baby$ParticipantID)
  metadata_baby = metadata_baby[!keep,]
  
  metadata_mom = metadata_16s[metadata_16s$SampleType == 'MV1D' & !is.na(metadata_16s$SampleType),]
  metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
  metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
  keep = duplicated(metadata_mom$ParticipantID)
  metadata_mom = metadata_mom[!keep,]
  
  keep = (metadata_mom$ParticipantID %in% metadata_baby$ParticipantID) & 
    (metadata_mom$SampleID %in% colnames(reads_table_human_MTG))
  sum(keep)
  
  metadata_mom = metadata_mom[keep,]
  reads_table_human_MTG_16s = reads_table_human_MTG[,colnames(reads_table_human_MTG) %in% metadata_mom$SampleID]
  reads_table_human_MTG_16s = reads_table_human_MTG_16s[metadata_mom$SampleID]
  
  metadata_baby = metadata_baby[metadata_baby$ParticipantID %in% metadata_mom$ParticipantID,]
  reads_table_16s_human_MTG = reads_table_16s[,colnames(reads_table_16s) %in% metadata_baby$SampleID]
  reads_table_16s_human_MTG = reads_table_16s_human_MTG[metadata_baby$SampleID]
  
  # normalization
  reads_table_16s_human_MTG = t(reads_table_16s_human_MTG)
  reads_table_16s_human_MTG = Rarefy(reads_table_16s_human_MTG, depth = min(rowSums(reads_table_16s_human_MTG)))
  reads_table_16s_human_MTG <- reads_table_16s_human_MTG$otu.tab.rff
  reads_table_16s_human_MTG <- as.data.frame(reads_table_16s_human_MTG)
  
  alpha.shannon_diversity <- data.frame(diversity(reads_table_16s_human_MTG))
  alpha.shannon_diversity <- as.data.frame(t(alpha.shannon_diversity))
  
  reads_table_human_MTG_16s = as.matrix(reads_table_human_MTG_16s)
  reads_table_human_MTG_16s = varianceStabilizingTransformation(reads_table_human_MTG_16s) # sample in column    reads = FPKM*colSum/10e6
  reads_table_human_MTG_16s = as.data.frame(reads_table_human_MTG_16s)
  colnames(reads_table_human_MTG_16s) = row.names(reads_table_16s_human_MTG)
  
  cor = newwork_rcorr3(reads_table_human_MTG_16s, alpha.shannon_diversity, pvalue = 0.05, 
                       cor_parameter= 0, style = 1, bar_max = 1, bar_min = -1, 
                       pheatmap_fontsize = 5 , alpha = 0.05, mc.cores =8, 
                       treeheight = 10)
  
  pdf(paste0("Paired_Mom_gene_",type1_all[a],'_shannon.pdf'), width=5, height=5)
  cor$p
  dev.off()
  
}
