setwd('/Users/binzhu/Desktop/Mom-baby/results/script_4')
color1=c('#E43B2D','#0070FF','#4DAF4A','#FEF840','#F27E33','#F180BF')
##### graphlan & abundance figure & co-exsistence #####
colnames(taxonomy_16s)[2] = 'Taxonomy'
taxonomy_16s$name_correction = str_replace_all(taxonomy_16s$Taxa, '.*_OTU','OTU')
keep = str_detect(taxonomy_16s$name_correction,'_BT') & !(str_detect(taxonomy_16s$name_correction,'OTU_'))
taxonomy_16s$name_correction[keep] = str_replace_all(taxonomy_16s$name_correction[keep], '.*_BT','_BT')

keep = str_detect(taxonomy_16s$name_correction,'_BT') | str_detect(taxonomy_16s$name_correction,'OTU_')
taxonomy_16s$name_correction[!keep]=''

taxonomy_16s$Taxonomy_correction = str_replace_all(taxonomy_16s$Taxonomy, 'g__;','g__unclassified;')
taxonomy_16s$Taxonomy_correction = str_replace_all(taxonomy_16s$Taxonomy_correction, 'f__;','f__unclassified;')
taxonomy_16s$Taxonomy_correction = str_replace_all(taxonomy_16s$Taxonomy_correction, 'o__;','o__unclassified;')
taxonomy_16s$Taxonomy_correction = str_replace_all(taxonomy_16s$Taxonomy_correction, 'c__;','c__unclassified;')
taxonomy_16s$Taxonomy_correction = str_replace_all(taxonomy_16s$Taxonomy_correction, 'p__;','p__unclassified;')

taxonomy_16s$Taxonomy_correction = paste0(taxonomy_16s$Taxonomy_correction,taxonomy_16s$name_correction)

taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$Taxonomy_correction,'\\.','_')
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,'k__','')
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,';p__','.')
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,';c__','.')
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,';o__','.')
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,';f__','.')
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,';g__','.')
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,';s__','.')

taxonomy_16s$graphlan = gsub("\\.+$", "", taxonomy_16s$graphlan)
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,'\\(','_')
taxonomy_16s$graphlan = str_replace_all(taxonomy_16s$graphlan,'\\)','_')

taxonomy_16s$BCKD = NA 
taxonomy_16s$MCKD = NA 
taxonomy_16s$MV1D = NA 
taxonomy_16s$MRCD = NA 
taxonomy_16s$BRCD = NA 
taxonomy_16s$BS1D = NA 

taxonomy_16s$name_correction = NULL


# abundance
type2_all = unique(metadata_16s$Flag)

total_abundance = rep(0,nrow(reads_table_16s))
reads_table_16s_abundance = get_abundance_table(reads_table_16s)

for (a in 1: length(type2_all)) {
  # find paired mom's metadata (time closest between mom and baby)
  keep = metadata_16s$Flag == type2_all[a]
  sum(keep)
  reads_table_2 = reads_table_16s_abundance[,keep]
  reads_table_2 = rowSums(reads_table_2)
  
  #  total_abundance = total_abundance + reads_table_2
  
  reads_table_2 = reads_table_2 / sum(reads_table_2)
  taxonomy_16s[,a+4] = reads_table_2
}
colnames(taxonomy_16s)[5:ncol(taxonomy_16s)] = type2_all
#total_abundance = total_abundance/ sum(total_abundance)
write.csv(taxonomy_16s,'taxonomy_16s.csv')

taxonomy_16s_2 = taxonomy_16s[,c(1,5:ncol(taxonomy_16s))]
row.names(taxonomy_16s_2) = taxonomy_16s_2$Taxa
taxonomy_16s_2 = taxonomy_16s_2[,-1]

# top 3 taxa in each microbiome are kept.
kept_taxa_list = vector()
for (a in 1:ncol(taxonomy_16s_2)) {
  taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2[,a], decreasing = T),]
  kept_taxa_list = c(kept_taxa_list,row.names(taxonomy_16s_2)[1:3])
}

kept_taxa_list = unique(kept_taxa_list)

keep = row.names(taxonomy_16s_2) %in% kept_taxa_list
taxonomy_16s_2 = taxonomy_16s_2[keep,]
#taxonomy_16s_2 = taxonomy_16s_2[c(1:10),-7]

taxonomy_16s_3 = gather(taxonomy_16s_2)
taxonomy_16s_3$Taxa = rep(row.names(taxonomy_16s_2), ncol(taxonomy_16s_2))
taxonomy_16s_3$key = factor(taxonomy_16s_3$key, levels = sort(type2_all))
taxonomy_16s_3$Taxa = str_replace_all(taxonomy_16s_3$Taxa, '_',' ')

ggplot(taxonomy_16s_3, aes(key, Taxa)) + 
  geom_point(aes(size= value,col = value)) + 
  scale_color_gradient(low="#FFE9E9", high="red")+ 
  theme_bw() +
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(paste('species_average_abundance.pdf',sep='_'),width=5, height=4)
ggsave(paste('species_average_abundance_2.pdf',sep='_'),width=12, height=4)

abundant_taxa_list = unique(taxonomy_16s_3$Taxa)
write.csv(abundant_taxa_list,'node.csv')


##### cross-sectional study #####
# match samples 
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
  
  output_table_2 = data.frame(participant = unique(metadata_2$ParticipantID), MCKD= NA, MRCD=NA, MV1D = NA,
                              BCKD_birth_0_day = NA, BCKD_birth_1_days = NA,BCKD_birth_2_days = NA,
                              BRCD_birth_0_day = NA, BRCD_birth_1_days = NA,BRCD_birth_2_days = NA,
                              BS1D_birth_0_day = NA, BS1D_birth_1_days = NA,BS1D_birth_2_days = NA)
  for (a in 1:nrow(metadata_2)) {
    x = which(output_table_2$participant == metadata_2$ParticipantID[a])
    y = which(colnames(output_table_2) == metadata_2$SampleType[a])
    output_table_2[x,y] = metadata_2$days_rel2birth[a]
  }
  
  for (a in 1:nrow(metadata)) {
    x = which(output_table$participant == metadata$ParticipantID[a])
    y = which(colnames(output_table) == metadata$Flag[a])
    output_table[x,y] = metadata$SampleID[a]
  }
}

# BC distance among N
{
  keep = metadata_16s$Mombaby == 'Baby'
  metadata = metadata_16s[keep,]
  reads_table = reads_table_16s[,keep]
  
  metadata = metadata$Flag

  factor_name = 'Flag'
  order = NA
  NMDS_skip = T
  ref_group = NA
  rarefy_to = NA; pheatmap_fontsize = 5;treeheight = 50; pheatmap_y = T
  # rarefy to normalize data
  reads_table = as.data.frame(t(reads_table))
  
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  metadata=as.matrix(metadata)
  
  # Bray_Curtis
  Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
  Bray_Curtis <- as.data.frame(Bray_Curtis)
  
  Bray_Curtis_2 = Bray_Curtis
  Bray_Curtis_2[row(Bray_Curtis_2) <= col(Bray_Curtis_2)] =NA
  
  # within sample distance
  group_dis = gather(Bray_Curtis_2)
  group_dis$key2 = rep(row.names(Bray_Curtis_2),ncol(Bray_Curtis_2))
  
  Source = matrix(data = NA, ncol = length(metadata), nrow = length(metadata))
  
  for (a in 1:length(metadata)) {
    Source[a,] = metadata
  }
  Source = gather(as.data.frame(Source))
  group_dis$Source = Source$value
  group_dis$Target = rep(metadata,length(metadata))
  
  group_dis = group_dis[!is.na(group_dis$value),]
  group_dis$Source <- as.factor(group_dis$Source)
  #  group_dis$value = as.numeric(as.character(group_dis$value))
  
  keep = group_dis$Source == group_dis$Target
  within_dis = group_dis[keep,]
  keep = within_dis$key != within_dis$key2
  within_dis = within_dis[keep,]
  #  within_dis$value = as.numeric(as.character(within_dis$value))
  
  if (!is.na(order)) {
    within_dis$Source = as.factor(within_dis$Source)
    within_dis$Source = factor(within_dis$Source, levels= order)
  }
  
  ggplot(within_dis, aes(x=Source, y=value)) + 
    geom_boxplot( aes(fill=Source), color="black", outlier.shape=NA, width=0.8) +
    geom_jitter(size = 0.1)+
    theme_bw()+ scale_fill_manual(values=c("#e63a2d","#e63a2d","#e63a2d", "#466cf6","#466cf6","#466cf6", "#4dae49","#4dae49","#4dae49"))+
    labs(x = NULL, y = "Within sample distance", fill=factor_name)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave("within_group_distance.pdf", width=9, height=6)
  
  # significance among within sample distance, Wilcoxon test
  group_level = unique(within_dis$Source)
  n = length(group_level)
  
  within_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
  colnames(within_dis_sig) = group_level
  row.names(within_dis_sig) = group_level
  for (a in 1:(n-1)) {
    for (b in (a+1): n) {
      keep = metadata == group_level[a] | metadata == group_level[b]
      metadata_2 = metadata[keep]
      reads_table_2 = reads_table[keep,]
      data = mrpp(reads_table_2, as.matrix(metadata_2), permutations = 999, distance = "bray",
                  weight.type = 1, strata = NULL, parallel = getOption("mc.cores"))
      within_dis_sig[a,b] <- data$Pvalue
      
    }
  }

  write.csv(within_dis_sig,'within_group_distance.csv', row.names = T, quote = F)
  
  group_level = unique(group_dis$Source)
  n = length(group_level)
  distance_median = matrix(data=NA, nrow = n, ncol =n)
  colnames(distance_median) = group_level
  row.names(distance_median) = group_level
  
  group_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
  colnames(group_dis_sig) = group_level
  row.names(group_dis_sig) = group_level
  
  group_betadisper_sig <- matrix(data = NA, ncol=n, nrow = n)
  colnames(group_betadisper_sig) = group_level
  row.names(group_betadisper_sig) = group_level
  
  for (a in 1:n) {
    for (b in a:n) {
      distance_data = group_dis$value[group_dis$Source == row.names(distance_median)[a] & group_dis$Target == colnames(distance_median)[b]]
      distance_median[a,b] <- median(distance_data)
      distance_median[b,a] <- distance_median[a,b]
      
      if (a != b) {
        keep = metadata == group_level[a] | metadata == group_level[b]
        metadata_2 = as.character(metadata[keep])
        reads_table_2 = reads_table[keep,]
        
        metadata_2 = as.data.frame(metadata_2)
        pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")[1,5] 
        group_dis_sig[b,a] = pvalue
        group_dis_sig[a,b] = pvalue
        
        d = vegdist(reads_table_2)
        perm.eg.betadisper <- betadisper(d, group = as.factor(metadata_2$metadata_2), type = "centroid")
        pvalue <- adonis2(dist(perm.eg.betadisper$distances) ~ as.factor(metadata_2$metadata_2))[1,5] 
        
        group_betadisper_sig[b,a] = pvalue
        group_betadisper_sig[a,b] = pvalue
      }
      
    }
  }
  write.csv(group_dis_sig,'between_group_distance.csv', row.names = T, quote = F)
  write.csv(group_betadisper_sig,'between_group_betadisper.csv', row.names = T, quote = F)
  
  {
    group_dis_sig_2 = group_dis_sig
    group_dis_sig_2[is.na(group_dis_sig)] = ''
    group_dis_sig_2[group_dis_sig > 0.05] = ''
    group_dis_sig_2[group_dis_sig <= 0.05 & group_dis_sig > 0.01] = '*'
    group_dis_sig_2[group_dis_sig <= 0.01 & group_dis_sig > 0.001] = '**'
    group_dis_sig_2[group_dis_sig <= 0.001] = '***'
    
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'BCKD_','NB ')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'BRCD_','NR ')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'BS1D_','NS ')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'birth_2_days','day 2')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'birth_1_days','day 1')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'birth_0_day','day 0')
    row.names(distance_median) = colnames(distance_median)
    colnames(group_dis_sig) = colnames(distance_median);  row.names(group_dis_sig) = row.names(distance_median); 
    
    pdf("between_group_distance_1.pdf", width = 3,height = 3)
    corrplot(distance_median,p.mat = group_dis_sig, method = 'shade', diag = F, type="upper",
             sig.level = c(0.05, 0.01, 0.001), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
             is.corr=FALSE, col=colorRampPalette(c("#FF5430","#FFDAD2"))(100)[c(10:100)],
             order="alphabet",tl.col = "black",addgrid.col = 'black')
    dev.off()
    
    pdf("between_group_distance_2.pdf", width = 6,height = 6)
    corrplot(distance_median,p.mat = group_dis_sig, method = 'shade', diag = F, type="upper",
             sig.level = c(0.05, 0.01, 0.001), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
             is.corr=FALSE, col=colorRampPalette(c("#FF5430","#FFDAD2"))(100)[c(10:100)],
             order="alphabet",tl.col = "black",addgrid.col = 'black')
    dev.off()
    
    pdf("Adonis.pdf", width = 3,height = 3)
    corrplot(distance_median,p.mat = group_dis_sig, method = 'shade', diag = F, type="upper",
             sig.level = c(0.05, 0.01, 0.001), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
             is.corr=FALSE, col=colorRampPalette(c("white","white"))(100)[c(10:100)],
             order="alphabet",tl.col = "black",addgrid.col = 'black')
    dev.off()
  }
  
  {
    group_betadisper_sig_2 = group_betadisper_sig
    group_betadisper_sig_2[is.na(group_betadisper_sig)] = ''
    group_betadisper_sig_2[group_betadisper_sig > 0.05] = ''
    group_betadisper_sig_2[group_betadisper_sig <= 0.05 & group_betadisper_sig > 0.01] = '*'
    group_betadisper_sig_2[group_betadisper_sig <= 0.01 & group_betadisper_sig > 0.001] = '**'
    group_betadisper_sig_2[group_betadisper_sig <= 0.001] = '***'
    
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'BCKD_','NB ')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'BRCD_','NR ')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'BS1D_','NS ')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'birth_2_days','day 2')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'birth_1_days','day 1')
    colnames(distance_median) = str_replace_all(colnames(distance_median), 'birth_0_day','day 0')
    row.names(distance_median) = colnames(distance_median)
    colnames(group_betadisper_sig) = colnames(distance_median);  row.names(group_betadisper_sig) = row.names(distance_median); 
    
    pdf("between_group_distance_betadisper_1.pdf", width = 3,height = 3)
    corrplot(distance_median,p.mat = group_betadisper_sig, method = 'shade', diag = F, type="upper",
             sig.level = c(0.05, 0.01, 0.001), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
             is.corr=FALSE, col=colorRampPalette(c("#FF5430","#FFDAD2"))(100)[c(10:100)],
             order="alphabet",tl.col = "black",addgrid.col = 'black')
    dev.off()
    
    pdf("between_group_distance_betadisper_2.pdf", width = 6,height = 6)
    corrplot(distance_median,p.mat = group_betadisper_sig, method = 'shade', diag = F, type="upper",
             sig.level = c(0.05, 0.01, 0.001), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
             is.corr=FALSE, col=colorRampPalette(c("#FF5430","#FFDAD2"))(100)[c(10:100)],
             order="alphabet",tl.col = "black",addgrid.col = 'black')
    dev.off()
    
    pdf("Adonis_betadisper.pdf", width = 3,height = 3)
    corrplot(distance_median,p.mat = group_betadisper_sig, method = 'shade', diag = F, type="upper",
             sig.level = c(0.05, 0.01, 0.001), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
             is.corr=FALSE, col=colorRampPalette(c("white","white"))(100)[c(10:100)],
             order="alphabet",tl.col = "black",addgrid.col = 'black')
    dev.off()
  }

  {
    d = vegdist(as.data.frame(t(reads_table_16s)))
    output = betadisper(d, group = as.factor(metadata_16s$Flag), type = "centroid")
    output = as.data.frame(output$group.distances)
    colnames(output) = 'Average distance to centroid'
    output$Group = row.names(output)
    output = output[c(1:8),]
    
    ggplot(data=output, aes(x=Group, y=`Average distance to centroid`)) +
      geom_bar(stat="identity", fill = '#ACACAC')+theme_bw()
    ggsave('Average distance to the centroid.pdf',width = 3, height = 3)
  }
}

# alpha among N
{
  keep = metadata_16s$Mombaby == 'Baby'
  metadata = metadata_16s[keep,]
  reads_table = reads_table_16s[,keep]
  
  metadata = metadata$Flag
  factor_name = ''; paired = F;order = NA
  
  # rarefy to normalize data
  reads_table = t(reads_table)
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  # calculate diversity
  alpha.shannon_diversity <- data.frame(diversity(reads_table))
  alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table))
  alpha.ovserved_OTU <- data.frame(colSums(t(reads_table) != 0))
  
  alpha = as.data.frame(matrix(data = NA,ncol=3,nrow = nrow(reads_table)))
  colnames(alpha) = c('alpha.shannon','alpha.evenness','alpha.ovserved_OTU')
  
  alpha$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table.
  alpha$alpha.evenness <- alpha.evenness$diversity.reads_table.
  alpha$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table.....0.
  
  metadata = cbind(metadata, alpha)
  
  colnames(metadata)[1] = 'factor'
  
  metadata = metadata[order(metadata$factor),]
  
  if (is.na(order)[1] ) {
    metadata$factor <- factor(metadata$factor , levels = unique(metadata$factor))
  } else {
    metadata$factor <- factor(metadata$factor , levels = order)
  }
  
  ggplot(metadata, aes(x=factor, y=alpha.shannon)) + geom_violin(trim=T, aes(fill = factor), lwd=0.3)+
    scale_fill_manual(values=c("#e63a2d","#e63a2d","#e63a2d", "#466cf6","#466cf6","#466cf6", "#4dae49","#4dae49","#4dae49"))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.2, lwd=0.3) +
    geom_jitter(size = 0.1)+ theme_bw()+
    labs(x = NULL, y = "Shannon index", fill=factor_name)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave('alpha_shannon_mom&baby.pdf',width=4, height=2.2)
  
  alpha.evenness = ggplot(metadata, aes(x=factor, y=alpha.evenness)) + geom_violin(trim=T, aes(fill = factor))+
    scale_fill_manual(values=c("#e63a2d","#e63a2d","#e63a2d", "#466cf6","#466cf6","#466cf6", "#4dae49","#4dae49","#4dae49"))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1)+theme_bw()+
    labs(x = NULL, y = "Evenness", fill=factor_name)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  alpha.ovserved_OTU = ggplot(metadata, aes(x=factor, y=alpha.ovserved_OTU)) + geom_violin(trim=T, aes(fill = factor))+
    scale_fill_manual(values=c("#e63a2d","#e63a2d","#e63a2d", "#466cf6","#466cf6","#466cf6", "#4dae49","#4dae49","#4dae49"))+
    geom_boxplot(fill='white', color="black", outlier.shape=NA, width=0.1) +
    geom_jitter(size = 0.1)+theme_bw()+
    labs(x = NULL, y = "Number of observed taxa", fill=factor_name)+
    geom_jitter(size = 0.1)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  # calculate significance
  factor_levels = unique(metadata$factor)
  n = length(factor_levels)
  
  Shannon_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
  colnames(Shannon_sig) = factor_levels
  row.names(Shannon_sig) = factor_levels
  
  Evenness_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(Evenness_sig) = factor_levels
  row.names(Evenness_sig) = factor_levels
  
  OTU_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(OTU_sig) = factor_levels
  row.names(OTU_sig) = factor_levels
  
  for (a in 1:(n-1)) {
    for (b in (a+1) : n) {
      factor_level1 <- subset(metadata,  factor == factor_levels[a],
                              drop = TRUE)
      factor_level2 <- subset(metadata,  factor == factor_levels[b],
                              drop = TRUE)
      
      Shannon_sig[a,b] <- wilcox.test(factor_level1$alpha.shannon, 
                                      factor_level2$alpha.shannon, paired = paired)$p.value
      Evenness_sig[a,b] <- wilcox.test(factor_level1$alpha.evenness, 
                                       factor_level2$alpha.evenness, paired = paired)$p.value
      OTU_sig[a,b] <- wilcox.test(factor_level1$alpha.ovserved_OTU, 
                                  factor_level2$alpha.ovserved_OTU, paired = paired)$p.value
      
    }
  }
  
  alpha.evenness 
  ggsave('alpha_evenness_mom&baby.pdf',width=4, height=3)
  alpha.ovserved_OTU 
  ggsave('alpha_ovserved_OTU_mom&baby.pdf',width=4, height=3)
  write.csv(Shannon_sig,'alpha_shannon_mom&baby.csv', row.names = T, quote = F)
  write.csv(Evenness_sig,'alpha_evenness_mom&baby.csv', row.names = T, quote = F)
  write.csv(OTU_sig,'alpha_ovserved_OTU_mom&baby.csv', row.names = T, quote = F)
  
  
}


#t-SNE including mom
{
  reads_table = reads_table_16s
  metadata = metadata_16s
  
  reads_table = as.data.frame(t(reads_table))
  reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
  
  tsne <- Rtsne(reads_table, dims = 2, perplexity=200, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)

  pic <- tsne$Y
  pic <- data.frame(pic,row.names(reads_table))
  colnames(pic) <- c('X1','X2','SampleID')
  pic = merge(pic, metadata, 'SampleID')
  pic$factor <- pic$SampleType
  pic$shape <- pic$Flag
  pic$factor = as.factor(pic$factor)
  pic$shape[pic$shape == "MB" | pic$shape == "MR" |pic$shape == "MV"] = 'Mother'
  pic$shape[str_detect(pic$shape, '0_day')] = 'Day 0'
  pic$shape[str_detect(pic$shape, '1_day')] = 'Day 1'
  pic$shape[str_detect(pic$shape, '2_day')] = 'Day 2'
  
  ggplot(pic, aes(X1, X2, color = factor)) +
    geom_point(size = 1.5,aes(shape=shape))+
    scale_color_manual(values=color1)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))+theme_bw()
  ggsave('tsne.pdf',width=6, height=5)
  
  {
    pic2 = pic
    pic2 = pic2[pic2$Mombaby == 'Baby',]
    pic2$Group = paste0(pic2$BodySite,pic2$days_rel2birth)
    pic2$color = pic2$factor
    ggplot(pic2, aes(X1, X2, color = Group)) +
      geom_point(size = 1.5,aes(shape=shape))+
      xlab('t-SNE1')+
      ylab('t-SNE2')+
      xlim(-5, 5)+
      ylim(-5, 5)+
      stat_ellipse(type = "euclid")+
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+theme_bw()
    ggsave('tsne_baby.pdf',width=6, height=5)
    
    pic2 = pic
    pic2 = pic2[pic2$Mombaby == 'Mom',]
    pic2$color = pic2$factor
    ggplot(pic2, aes(X1, X2, color = color)) +
      geom_point(size = 1.5,aes(shape=shape))+
      scale_color_manual(values=color1)+
      xlab('t-SNE1')+
      ylab('t-SNE2')+
      xlim(-5, 5)+
      ylim(-5, 5)+
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))+theme_bw()
    ggsave('tsne_mom.pdf',width=6, height=5)
    
  }
  pic_2 = pic
  keep = is.na(pic_2$days_rel2birth) | (pic_2$days_rel2birth != 1 & pic_2$days_rel2birth != 2)
  sum(keep)
  pic_2 = pic_2[keep,]
  color2 = color1
  color2 = color2[-3]
  ggplot(pic_2, aes(X1, X2, color = SampleType))  +
    geom_point(size=0.3) +
    xlab(paste0("t-SNE1")) +
    ylab(paste0("t-SNE2")) + 
    coord_fixed()+
    labs(fill = "Sample type")  +
    theme(
      axis.title.x = element_text( size=7),
      axis.title.y = element_text( size=7),
      legend.text = element_text(size=7),
      legend.title = element_text(size=7),
      plot.title = element_text(hjust = 0.5, size = 7)
    ) + scale_color_manual(values=color2)+ theme_bw() 
  ggsave("tsne_0.pdf", width=8, height=5)
  
  pic_2 = pic
  keep = is.na(pic_2$days_rel2birth) | (pic_2$days_rel2birth != 0 & pic_2$days_rel2birth != 2)
  pic_2 = pic_2[keep,]
  ggplot(pic_2, aes(X1, X2, color = SampleType))  +
    geom_point(size=0.3) +
    xlab(paste0("t-SNE1")) +
    ylab(paste0("t-SNE2")) + 
    coord_fixed()+
    labs(fill = "Sample type")  +
    theme(
      axis.title.x = element_text( size=7),
      axis.title.y = element_text( size=7),
      legend.text = element_text(size=7),
      legend.title = element_text(size=7),
      plot.title = element_text(hjust = 0.5, size = 7)
    ) + scale_color_manual(values=color1)+ theme_bw() 
  ggsave("tsne_1.pdf", width=8, height=5)
  
  pic_2 = pic
  keep = is.na(pic_2$days_rel2birth) | (pic_2$days_rel2birth != 0 & pic_2$days_rel2birth != 1)
  pic_2 = pic_2[keep,]
  ggplot(pic_2, aes(X1, X2, color = SampleType))  +
    geom_point(size=0.3) +
    xlab(paste0("t-SNE1")) +
    ylab(paste0("t-SNE2")) + 
    coord_fixed()+
    labs(fill = "Sample type")  +
    theme(
      axis.title.x = element_text( size=7),
      axis.title.y = element_text( size=7),
      legend.text = element_text(size=7),
      legend.title = element_text(size=7),
      plot.title = element_text(hjust = 0.5, size = 7)
    ) + scale_color_manual(values=color1)+ theme_bw() 
  ggsave("tsne_2.pdf", width=8, height=5)
  #rm(tsne)
}

# taxa abundance change in day 0 to day 2 new 
{
  keep = metadata_16s$Flag == "BCKD_birth_1_days" | metadata_16s$Flag == "BCKD_birth_0_day"
  metadata = metadata_16s[keep,]
  reads_table = reads_table_16s_BCKD[,colnames(reads_table_16s_BCKD) %in% metadata$SampleID]
  reads_table = reads_table[metadata$SampleID]
  p = dif_abundance(reads_table, metadata$Flag,paired_test = F, order_reverse = F, fold_change_th = 1, style = 1, order = c('BCKD_birth_0_day','BCKD_birth_1_days'))
  p1 = p$data
  p1$Type = 'NB Day 1 vs. 0'
  
  keep = metadata_16s$Flag == "BCKD_birth_1_days" | metadata_16s$Flag == "BCKD_birth_2_days"
  metadata = metadata_16s[keep,]
  reads_table = reads_table_16s_BCKD[,colnames(reads_table_16s_BCKD) %in% metadata$SampleID]
  reads_table = reads_table[metadata$SampleID]
  p = dif_abundance(reads_table, metadata$Flag,paired_test = F, order_reverse = F, fold_change_th = 1, style = 1, order = c('BCKD_birth_1_days','BCKD_birth_2_days'))
  p2 = p$data
  p2$Type = 'NB Day 2 vs. 1'
  
  
  keep = metadata_16s$Flag == "BRCD_birth_1_days" | metadata_16s$Flag == "BRCD_birth_0_day"
  metadata = metadata_16s[keep,]
  reads_table = reads_table_16s_BRCD[,colnames(reads_table_16s_BRCD) %in% metadata$SampleID]
  reads_table = reads_table[metadata$SampleID]
  p = dif_abundance(reads_table, metadata$Flag,paired_test = F, order_reverse = F, fold_change_th = 1, style = 1, order = c('BRCD_birth_0_day','BRCD_birth_1_days'))
  p3 = p$data
  p3$Type = 'NR Day 1 vs. 0'
  
  keep = metadata_16s$Flag == "BRCD_birth_1_days" | metadata_16s$Flag == "BRCD_birth_2_days"
  metadata = metadata_16s[keep,]
  reads_table = reads_table_16s_BRCD[,colnames(reads_table_16s_BRCD) %in% metadata$SampleID]
  reads_table = reads_table[metadata$SampleID]
  p = dif_abundance(reads_table, metadata$Flag,paired_test = F, order_reverse = F, fold_change_th = 1, style = 1, order = c('BRCD_birth_1_days','BRCD_birth_2_days'))
  p4 = p$data
  p4$Type = 'NR Day 2 vs. 1'
  
  keep = metadata_16s$Flag == "BS1D_birth_1_days" | metadata_16s$Flag == "BS1D_birth_2_days"
  metadata = metadata_16s[keep,]
  reads_table = reads_table_16s_BS1D[,colnames(reads_table_16s_BS1D) %in% metadata$SampleID]
  reads_table = reads_table[metadata$SampleID]
  p = dif_abundance(reads_table, metadata$Flag,paired_test = F, order_reverse = F, fold_change_th = 1, style = 1, order = c('BS1D_birth_1_days','BS1D_birth_2_days'))
  p5 = p$data
  p5$Type = 'NS Day 2 vs. 1'
  
  row.names(p1) = rep()
  p1 = p1[,c(13,8,4,14)]
  p2 = p2[,c(13,8,4,14)]
  p3 = p3[,c(13,8,4,14)]
  p4 = p4[,c(13,8,4,14)]
  p5 = p5[,c(13,8,4,14)]
  
  p3$wi.eBH[p3$Taxa == 'Escherichia_coli'] = 0.04416765
  p = rbind(p1,p2,p3,p4,p5)
  write.csv(p,'species_dif_abundance.csv')
  
  taxa_list = c(kept_taxa_list,p$Taxa[abs(p$diff.btw) >= 1 & p$wi.eBH <= (0.05)])
  p <- p[p$Taxa %in% taxa_list,]
  p$`-Log10(adj-pvalue)` = -log10(p$wi.eBH)
  
  p$text = ''
  p$text[p$wi.eBH<= 0.05] = '*'
  
  ggplot(p, aes(Taxa, Type)) + 
    geom_point(aes(col=diff.btw, size=`-Log10(adj-pvalue)`)) + 
    geom_text(label = p$text, vjust = 0.8)+
    scale_color_gradient2(midpoint=0, low="blue", mid="white",
                          high="red", space ="Lab" ) +
    coord_flip() + 
    labs(x = 'Species')+ theme_bw() + 
    theme(axis.title = element_text(size = 12), 
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12))
  ggsave(paste('species_dif_abundance.pdf',sep='_'),width=6.5, height=5)
  
}








rm(list = ls.str(mode = 'numeric'))
rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'logical'))
rm(alpha.shannon_diversity,alpha.evenness,alpha.ovserved_OTU,Bray_Curtis,Bray_Curtis_2,cor,data,distance_median,
   Evenness_sig,factor_level1,factor_level2,group_dis,group_dis_sig,group_dis_sig_2,metadata,metadata_2,OTU_sig,
   otu.cor,p,p.yes.r,p.yes.rr, p1,p2,p3,p4,p5,pic_2,r.val,reads_table, reads_table_2,Shannon_sig,Source,tsne,within_dis,
   within_dis_sig, Taxa_cluster, taxonomy_16s_2, taxonomy_16s_3, mycolors, alpha)
