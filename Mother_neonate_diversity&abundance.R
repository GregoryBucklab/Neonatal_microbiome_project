##### graphlan & abundance figure & co-exsistence #####
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
type2_all = unique(metadata_16s$SampleType)
total_abundance = rep(0,nrow(reads_table_16s))

for (a in 1: length(type2_all)) {
  # find paired mom's metadata (time closest between mom and baby)
  keep = metadata_16s$SampleType == type2_all[a]
  sum(keep)
  reads_table_2 = reads_table_16s[,keep]
  reads_table_2 = rowSums(reads_table_2)
  
  #  total_abundance = total_abundance + reads_table_2
  
  reads_table_2 = reads_table_2 / sum(reads_table_2)
  taxonomy_16s[,a+4] = reads_table_2
}
#total_abundance = total_abundance/ sum(total_abundance)
write.csv(taxonomy_16s,'taxonomy_16s.csv')

taxonomy_16s_2 = taxonomy_16s[,c(1,5:10)]
row.names(taxonomy_16s_2) = taxonomy_16s_2$Taxa
taxonomy_16s_2 = taxonomy_16s_2[,-1]
#taxonomy_16s_2$Total_abundance = total_abundance

# top 3 taxa in each microbiome are kept.
taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$BCKD, decreasing = T),]
kept_taxa_list = row.names(taxonomy_16s_2)[1:3]
taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$MCKD, decreasing = T),]
kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$MV1D, decreasing = T),]
kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$MRCD, decreasing = T),]
kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$BRCD, decreasing = T),]
kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$BS1D, decreasing = T),]
kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
kept_taxa_list = unique(kept_taxa_list)

keep = row.names(taxonomy_16s_2) %in% kept_taxa_list
taxonomy_16s_2 = taxonomy_16s_2[keep,]
#taxonomy_16s_2 = taxonomy_16s_2[c(1:10),-7]

taxonomy_16s_3 = gather(taxonomy_16s_2)
taxonomy_16s_3$Taxa = rep(row.names(taxonomy_16s_2), 6)
taxonomy_16s_3$key = factor(taxonomy_16s_3$key, levels = type2_all)

ggplot(taxonomy_16s_3, aes(key, Taxa)) + 
  geom_point(aes(size= value), col = '#888888') + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7)) + theme_bw() 
ggsave(paste('species_average_abundance.pdf',sep='_'),width=4, height=2.8)
ggsave(paste('species_average_abundance_2.pdf',sep='_'),width=4, height=4)

### day 0 only
{
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
  type2_all = unique(metadata_16s$SampleType)
  
  keep  = metadata_16s$Mombaby == 'Baby' & metadata_16s$days_rel2birth == 0
  metadata_16s_day0 = metadata_16s[keep,]
  reads_table_16s_day0 = reads_table_16s[,keep]
  
  total_abundance = rep(0,nrow(reads_table_16s_day0))
  
  for (a in 1: length(type2_all)) {
    # find paired mom's metadata (time closest between mom and baby)
    keep = metadata_16s_day0$SampleType == type2_all[a]
    sum(keep)
    reads_table_2 = reads_table_16s_day0[,keep]
    reads_table_2 = rowSums(reads_table_2)
    
    #  total_abundance = total_abundance + reads_table_2
    
    reads_table_2 = reads_table_2 / sum(reads_table_2)
    taxonomy_16s[,a+4] = reads_table_2
  }
  
  
  taxonomy_16s_2 = taxonomy_16s[,c(1,5:10)]
  row.names(taxonomy_16s_2) = taxonomy_16s_2$Taxa
  taxonomy_16s_2 = taxonomy_16s_2[,-1]
  #taxonomy_16s_2$Total_abundance = total_abundance
  
  # top 3 taxa in each microbiome are kept.
  taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$BCKD, decreasing = T),]
  kept_taxa_list = row.names(taxonomy_16s_2)[1:3]
  taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$MCKD, decreasing = T),]
  kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
  taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$MV1D, decreasing = T),]
  kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
  taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$MRCD, decreasing = T),]
  kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
  taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$BRCD, decreasing = T),]
  kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
  taxonomy_16s_2 = taxonomy_16s_2[order(taxonomy_16s_2$BS1D, decreasing = T),]
  kept_taxa_list = c(kept_taxa_list, row.names(taxonomy_16s_2)[1:3])
  kept_taxa_list = unique(kept_taxa_list)
  
  keep = row.names(taxonomy_16s_2) %in% kept_taxa_list
  taxonomy_16s_2 = taxonomy_16s_2[keep,]
  #taxonomy_16s_2 = taxonomy_16s_2[c(1:10),-7]
  
  taxonomy_16s_3 = gather(taxonomy_16s_2)
  taxonomy_16s_3$Taxa = rep(row.names(taxonomy_16s_2), 6)
  taxonomy_16s_3$key = factor(taxonomy_16s_3$key, levels = type2_all)
  
  taxonomy_16s_3 = taxonomy_16s_3[!is.na(taxonomy_16s_3$value),]
  ggplot(taxonomy_16s_3, aes(key, Taxa)) + 
    geom_point(aes(size= value), col = '#888888') + 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7)) + theme_bw() 
  ggsave(paste('species_average_abundance_day0.pdf',sep='_'),width=3.5, height=1.5)
  ggsave(paste('species_average_abundance_2_day0.pdf',sep='_'),width=3.5, height=3)
  
  
  
  
  
}




##### baby main taxa network #####
cor_parameter = 0

reads_table = reads_table_16s_BCKD
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BCKD.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BCKD.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BCKD_taxa_cluster.csv')

keep = metadata_16s$Flag == 'BCKD_birth_0_day'
sum(keep)
keep = colnames(reads_table_16s_BCKD) %in% metadata_16s$SampleID[keep]
reads_table = reads_table_16s_BCKD[,keep]
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BCKD_0.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BCKD_0.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BCKD_taxa_cluster_0.csv')

keep = metadata_16s$Flag == 'BCKD_birth_1_days'
sum(keep)
keep = colnames(reads_table_16s_BCKD) %in% metadata_16s$SampleID[keep]
reads_table = reads_table_16s_BCKD[,keep]
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BCKD_1.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BCKD_1.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BCKD_taxa_cluster_1.csv')

keep = metadata_16s$Flag == 'BCKD_birth_2_days'
sum(keep)
keep = colnames(reads_table_16s_BCKD) %in% metadata_16s$SampleID[keep]
reads_table = reads_table_16s_BCKD[,keep]
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BCKD_2.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BCKD_2.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BCKD_taxa_cluster_2.csv')




reads_table = reads_table_16s_BRCD
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BRCD.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BRCD.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BRCD_taxa_cluster.csv')

reads_table = reads_table_16s_BS1D
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BS1D.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BS1D.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BS1D_taxa_cluster.csv')


keep = metadata_16s$Flag == 'BRCD_birth_0_day'
sum(keep)
keep = colnames(reads_table_16s_BRCD) %in% metadata_16s$SampleID[keep]
reads_table = reads_table_16s_BRCD[,keep]
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BRCD_0.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BRCD_0.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BRCD_taxa_cluster_0.csv')

keep = metadata_16s$Flag == 'BRCD_birth_1_days'
sum(keep)
keep = colnames(reads_table_16s_BRCD) %in% metadata_16s$SampleID[keep]
reads_table = reads_table_16s_BRCD[,keep]
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BRCD_1.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BRCD_1.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BRCD_taxa_cluster_1.csv')

keep = metadata_16s$Flag == 'BRCD_birth_2_days'
sum(keep)
keep = colnames(reads_table_16s_BRCD) %in% metadata_16s$SampleID[keep]
reads_table = reads_table_16s_BRCD[,keep]
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BRCD_2.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BRCD_2.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BRCD_taxa_cluster_2.csv')




# NR without day 0
reads_table = reads_table_16s_BRCD
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]

keep = metadata_16s$SampleID %in% colnames(reads_table)
metadata = metadata_16s[keep,]
reads_table = reads_table[metadata$SampleID]
keep = metadata$Flag != "BRCD_birth_0_day"
reads_table = reads_table[,keep]

pdf("network_BRCD_no_day_0.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BRCD_no_day_0.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BRCD_taxa_cluster_no_day_0.csv')



keep = metadata_16s$Flag == 'BS1D_birth_1_days'
sum(keep)
keep = colnames(reads_table_16s_BS1D) %in% metadata_16s$SampleID[keep]
reads_table = reads_table_16s_BS1D[,keep]
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BS1D_1.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BS1D_1.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BS1D_taxa_cluster_1.csv')

keep = metadata_16s$Flag == 'BS1D_birth_2_days'
sum(keep)
keep = colnames(reads_table_16s_BS1D) %in% metadata_16s$SampleID[keep]
reads_table = reads_table_16s_BS1D[,keep]
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BS1D_2.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BS1D_2.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BS1D_taxa_cluster_2.csv')



##### diversity #####
# impact of participants on the diversity
type1_all = c('BCKD','BRCD','BS1D')

reads_table = reads_table_16s_BCKD
keep = metadata_16s$SampleID %in% colnames(reads_table)
sum(keep)
metadata = metadata_16s[keep,]
reads_table = reads_table[metadata$SampleID]

reads_table <- t(reads_table)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff

adonis2(reads_table ~ ParticipantID, data = metadata, method = "bray")

alpha.shannon_diversity <- data.frame(diversity(reads_table))
data = cbind(metadata$ParticipantID, alpha.shannon_diversity$diversity.reads_table.)
kruskal.test(alpha.shannon_diversity$diversity.reads_table. ~ metadata$ParticipantID)$p.value


reads_table = reads_table_16s_BRCD
keep = metadata_16s$SampleID %in% colnames(reads_table)
sum(keep)
metadata = metadata_16s[keep,]
reads_table = reads_table[metadata$SampleID]

reads_table <- t(reads_table)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff

adonis2(reads_table ~ ParticipantID, data = metadata, method = "bray")

alpha.shannon_diversity <- data.frame(diversity(reads_table))
data = cbind(metadata$ParticipantID, alpha.shannon_diversity$diversity.reads_table.)
kruskal.test(alpha.shannon_diversity$diversity.reads_table. ~ metadata$ParticipantID)$p.value


reads_table = reads_table_16s_BS1D
keep = metadata_16s$SampleID %in% colnames(reads_table)
sum(keep)
metadata = metadata_16s[keep,]
reads_table = reads_table[metadata$SampleID]

reads_table <- t(reads_table)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff

adonis2(reads_table ~ ParticipantID, data = metadata, method = "bray")

alpha.shannon_diversity <- data.frame(diversity(reads_table))
data = cbind(metadata$ParticipantID, alpha.shannon_diversity$diversity.reads_table.)
kruskal.test(alpha.shannon_diversity$diversity.reads_table. ~ metadata$ParticipantID)$p.value


# select maternal microbiomes at the last visit
# fig 1b
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

keep = is.na(output_table$BCKD_birth_0_day) & is.na(output_table$BCKD_birth_1_days) &is.na(output_table$BCKD_birth_2_days) &is.na(output_table$BRCD_birth_0_day) & is.na(output_table$BRCD_birth_1_days) &is.na(output_table$BRCD_birth_2_days) &is.na(output_table$BS1D_birth_0_day) & is.na(output_table$BS1D_birth_1_days) &is.na(output_table$BS1D_birth_2_days)
sum(keep)
write.csv(output_table , 'Sample_table.csv')

keep = gather(output_table)
keep = keep[!is.na(keep$value),]
reads_table_16s = reads_table_16s[,colnames(reads_table_16s) %in% keep$value]
reads_table_16s_MCKD = reads_table_16s_MCKD[,colnames(reads_table_16s_MCKD) %in% keep$value]
reads_table_16s_MRCD = reads_table_16s_MRCD[,colnames(reads_table_16s_MRCD) %in% keep$value]
reads_table_16s_MV1D = reads_table_16s_MV1D[,colnames(reads_table_16s_MV1D) %in% keep$value]
metadata_16s = metadata_16s[metadata_16s$SampleID %in% keep$value,]
metadata_16s$Flag[metadata_16s$SampleType == "MRCD"] = "MRCD"
metadata_16s$Flag[metadata_16s$SampleType == "MV1D"] = "MV1D"
metadata_16s$Flag[metadata_16s$SampleType == "MCKD"] = "MCKD"

# diversity analysis
reads_table = reads_table_16s
metadata = metadata_16s

beta_result = beta_diversity(reads_table, metadata = metadata$Flag, factor_name = 'Flag', treeheight = 50,
                             order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA,pheatmap_y = T)
#metadata = metadata$Flag; factor_name = 'Flag'; treeheight = 50;
#order = NA; NMDS_skip = T;ref_group = NA; rarefy_to = NA;pheatmap_y = T
pdf("beta_diversity_group_dis.pdf",width=3, height=2.7)
beta_result$group_dis_2_p
dev.off()
write.csv(beta_result$group_dis_sig,'group_dis_sig.csv', row.names = T)



reads_table_2 = reads_table[,metadata$Mombaby == 'Baby']
metadata_2 = metadata[metadata$Mombaby == 'Baby',]
beta2 = beta_diversity(reads_table_2, metadata = metadata_2$Flag, factor_name = 'Flag', order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA)
beta2$within_dis_p + theme_bw() 
ggsave('beta_within_dis.pdf',width=4, height=3)
write.csv(beta2$within_dis_sig,'beta_within_dis.csv', row.names = T)

# alpha diversity
alpha = alpha_diversity(reads_table, metadata = metadata$SampleType, factor_name = 'SampleType', paired = F,order = NA)
alpha$shannon + theme_bw() 
ggsave('alpha_shannon_mom&baby.pdf',width=4, height=2)
alpha$evenness + theme_bw() 
ggsave('alpha_evenness_mom&baby.pdf',width=4, height=2)
alpha$ovserved_OTU + theme_bw() 
ggsave('alpha_ovserved_OTU_mom&baby.pdf',width=4, height=2)
write.csv(alpha$sig_Shannon,'alpha_shannon_mom&baby.csv', row.names = T, quote = F)
write.csv(alpha$sig_Evenness,'alpha_evenness_mom&baby.csv', row.names = T, quote = F)
write.csv(alpha$sig_OTU,'alpha_ovserved_OTU_mom&baby.csv', row.names = T, quote = F)


alpha = alpha_diversity(reads_table_2, metadata = metadata_2$Flag, factor_name = 'Flag', paired = F,order = NA)
alpha$shannon + theme_bw() 
ggsave('alpha_shannon.pdf',width=4, height=2)
alpha$evenness + theme_bw() 
ggsave('alpha_evenness.pdf',width=4, height=2)
alpha$ovserved_OTU + theme_bw() 
ggsave('alpha_ovserved_OTU.pdf',width=4, height=2)
write.csv(alpha$sig_Shannon,'alpha_shannon.csv', row.names = T, quote = F)
write.csv(alpha$sig_Evenness,'alpha_evenness.csv', row.names = T, quote = F)
write.csv(alpha$sig_OTU,'alpha_ovserved_OTU.csv', row.names = T, quote = F)

# significance using Kruskal–Wallis test
reads_table_2 = t(reads_table_2)
reads_table_2 = Rarefy(reads_table_2, depth = min(rowSums(reads_table_2)))

reads_table_2 <- reads_table_2$otu.tab.rff
reads_table_2 <- as.data.frame(reads_table_2)

# calculate diversity
alpha.shannon_diversity <- data.frame(diversity(reads_table_2))
alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table_2))
alpha.ovserved_OTU <- data.frame(colSums(t(reads_table_2) != 0))

alpha = as.data.frame(matrix(data = NA,ncol=3,nrow = nrow(reads_table_2)))
colnames(alpha) = c('alpha.shannon','alpha.evenness','alpha.ovserved_OTU')

alpha$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table_2.
alpha$alpha.evenness <- alpha.evenness$diversity.reads_table_2.
alpha$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table_2.....0.

metadata_2 = cbind(metadata_2$Flag, metadata_2$SampleType,alpha)
kruskal.test(metadata_2$alpha.shannon[metadata_2$`metadata_2$SampleType` == 'BCKD'] ~ metadata_2$`metadata_2$Flag`[metadata_2$`metadata_2$SampleType` == 'BCKD'])$p.value
kruskal.test(metadata_2$alpha.shannon[metadata_2$`metadata_2$SampleType` == 'BRCD'] ~ metadata_2$`metadata_2$Flag`[metadata_2$`metadata_2$SampleType` == 'BRCD'])$p.value
kruskal.test(metadata_2$alpha.shannon[metadata_2$`metadata_2$SampleType` == 'BS1D'] ~ metadata_2$`metadata_2$Flag`[metadata_2$`metadata_2$SampleType` == 'BS1D'])$p.value

rm(alpha)

#t-SNE
reads_table = reads_table_16s
metadata = metadata_16s

reads_table = as.data.frame(t(reads_table))
reads_table = reads_table + 0.5

reads_table <- as.data.frame(clr(reads_table))      ### CLR normalization in rows
#reads_table = data.matrix(t(reads_table))

tsne <- Rtsne(reads_table, dims = 2, perplexity=200, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)

pic <- tsne$Y
pic <- data.frame(pic,row.names(reads_table))
colnames(pic) <- c('X1','X2','SampleID')
pic <- merge(pic,metadata,'SampleID')
pic <- pic[order(pic$SampleType, decreasing = T),]

ggplot(pic, aes(X1, X2, color = SampleType))  +
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
#  stat_ellipse(type = "norm")   ### If the question concerns the entire population as it is distributed, then the normal distribution should be used. If the question concerns the mean of the population then the t statistic may be used. 

ggsave("tsne.pdf", width=5, height=3)

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




# median distance among 6 microbiomes
keep = !is.na(metadata_16s$Flag)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

metadata_Flag_2 = str_replace_all(metadata$Flag, '_.*','')
beta_2 = beta_diversity(reads_table, metadata = metadata_Flag_2, factor_name = 'Flag', order = NA,
                        NMDS_skip = T, ref_group = NA, rarefy_to = NA,treeheight = 5,pheatmap_y = T)

pdf("beta_group_dis_2.pdf",width=1.8, height=1.6)
beta_2$group_dis_2_p
dev.off()
write.csv(beta_2$group_dis_sig,'group_dis_sig.csv', row.names = T)

# median distance among 6 microbiomes
keep = !is.na(metadata_16s$Flag)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

keep = str_detect(metadata$Flag, 'M')
metadata_Flag_2 = metadata$Flag
metadata_Flag_2[keep] = str_replace_all(metadata$Flag[keep], '_.*','')
beta_2 = beta_diversity(reads_table, metadata = metadata_Flag_2, factor_name = 'Flag', order = NA,
                        NMDS_skip = T, ref_group = NA, rarefy_to = NA,treeheight = 5,pheatmap_y = T)



# distance between three neonatal microbiomes
keep = !is.na(metadata_16s$Flag)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

beta = beta_diversity(reads_table, metadata = metadata$Flag, factor_name = 'Flag', order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA)
Bray_Curtis_2 = beta$Bray_Curtis
group_dis = gather(Bray_Curtis_2)
group_dis$key2 = rep(colnames(Bray_Curtis_2),ncol(Bray_Curtis_2))
keep = !(group_dis$key == group_dis$key2)
group_dis <- group_dis[keep,]
group_dis$Sample_type_1 = NA
group_dis$Sample_type_2 = NA

Sample_1_list = unique(group_dis$key)
for (a in 1: length(Sample_1_list)) {
  n = which(group_dis$key == Sample_1_list[a])
  m = which(metadata$SampleID == Sample_1_list[a])
  
  group_dis$Sample_type_1[n] = metadata$Flag[m]
}

Sample_2_list = unique(group_dis$key2)
for (a in 1: length(Sample_2_list)) {
  n = which(group_dis$key2 == Sample_2_list[a])
  m = which(metadata$SampleID == Sample_2_list[a])
  
  group_dis$Sample_type_2[n] = metadata$Flag[m]
}
group_dis$Sample = paste(group_dis$Sample_type_1,group_dis$Sample_type_2,sep = '_')

keep = group_dis$Sample == 'BCKD_birth_0_day_BRCD_birth_0_day' | 
  group_dis$Sample == 'BCKD_birth_1_days_BRCD_birth_1_days' |
  group_dis$Sample == 'BCKD_birth_2_days_BRCD_birth_2_days'
group_dis_2 <- group_dis[keep,]

ggplot(data=group_dis_2, aes(x=Sample, y=value)) +
  geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
  labs(x = 'Sample type', y = "Bray-Curtis distance")+ ylim(0.90,1) + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
ggsave('Bray-Curtis_NB_NR.pdf',width=2, height=3)

data = group_dis_2[,c(2,4)]
pvalue =  kruskal.test(value~Sample, data)$p.value
pvalue

keep = group_dis$Sample == 'BCKD_birth_1_days_BS1D_birth_1_days' |
  group_dis$Sample == 'BCKD_birth_2_days_BS1D_birth_2_days'
sum(keep)
group_dis_2 <- group_dis[keep,]

ggplot(data=group_dis_2, aes(x=Sample, y=value)) +
  geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
  labs(x = 'Sample type', y = "Bray-Curtis distance")+ ylim(0.90,1) + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
ggsave('Bray-Curtis_NB_NS.pdf',width=2, height=3)

data = group_dis_2[,c(2,4)]
pvalue =  kruskal.test(value~Sample, data)$p.value
pvalue

keep = group_dis$Sample == 'BRCD_birth_1_days_BS1D_birth_1_days' |
  group_dis$Sample == 'BRCD_birth_2_days_BS1D_birth_2_days'
sum(keep)
group_dis_2 <- group_dis[keep,]

ggplot(data=group_dis_2, aes(x=Sample, y=value)) +
  geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
  labs(x = 'Sample type', y = "Bray-Curtis distance")+ ylim(0.90,1) + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
ggsave('Bray-Curtis_NR_NS.pdf',width=2, height=3)

data = group_dis_2[,c(2,4)]
pvalue =  kruskal.test(value~Sample, data)$p.value
pvalue

# adonis to test the impact of day after delivery on the N microbiomes
keep = metadata_16s$SampleID %in% colnames(reads_table_16s_BCKD)
sum(keep)
metadata = metadata_16s[keep,]
reads_table = reads_table_16s_BCKD[metadata$SampleID]
metadata = metadata$Flag

reads_table <- t(reads_table)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)
metadata = as.data.frame(metadata)
pvalue <- adonis2(reads_table ~ ., data = metadata, method = "bray")  
pvalue


keep = metadata_16s$SampleID %in% colnames(reads_table_16s_BRCD)
sum(keep)
metadata = metadata_16s[keep,]
reads_table = reads_table_16s_BRCD[metadata$SampleID]
metadata = metadata$Flag

reads_table <- t(reads_table)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)
metadata = as.data.frame(metadata)
pvalue <- adonis2(reads_table ~ ., data = metadata, method = "bray")  
pvalue


keep = metadata_16s$SampleID %in% colnames(reads_table_16s_BS1D)
sum(keep)
metadata = metadata_16s[keep,]
reads_table = reads_table_16s_BS1D[metadata$SampleID]
metadata = metadata$Flag

reads_table <- t(reads_table)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)
metadata = as.data.frame(metadata)
pvalue <- adonis2(reads_table ~ ., data = metadata, method = "bray")  
pvalue

##### adonis multivariables #####
keep = metadata_16s$Mombaby == 'Baby'
sum(keep)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]
metadata_2 = metadata

#metadata_2 = cbind(metadata$days_rel2birth,metadata$respiratory_3rdtri,metadata$little_or_no_energy,metadata$yeast,metadata$birthcontrol_condoms,
#                 metadata$physical_activity_light,metadata$hpv_history,metadata$yogurt,metadata$pulse,metadata$student,metadata$SampleType)
#metadata_2 = as.data.frame(metadata_2)
#colnames(metadata_2) = c("Time_after_delivery","respiratory_3rdtri","little_or_no_energy","yeast","birthcontrol_condoms",
#                         "physical_activity_light","hpv_vaccination","yogurt","pulse","student",'Body_site')

metadata_2[metadata_2 == 'No'] = 0
metadata_2[metadata_2 == 'Yes'] = 1
metadata_2[metadata_2 == 'BCKD'] = 1
metadata_2[metadata_2 == 'BRCD'] = 2
metadata_2[metadata_2 == 'BS1D'] = 3
colnames(metadata_2)[colnames(metadata_2) == "days_rel2birth"] = "time_after_birth"

metadata_2 = apply(metadata_2, 2, as.numeric)

keep = colSums(is.na(metadata_2)) <= nrow(metadata_2)/3
metadata_2 = metadata_2[,keep]

# remove input with less than 2 levels
keep = sapply(1:ncol(metadata_2), function(j) (length(unique(metadata_2[,j])) >= 2 ))
metadata_2 = metadata_2[,keep]

metadata_2 = metadata_2[, !str_detect(colnames(metadata_2), 'average|baby|ga_at_delivery|event
                                     |pregnancy|last_weight|last_pulse
                                     |last_ga_at_delivery|mom_feeding_method|last_weight
                                     |hispanic_or_latino|vaginal_penetration_frequency
                                     |hepatitis_c|weight_change|3rdtri
                                     |poor_appetite|vaginal_penetration_frequency|2ndtri') ]

involved_factors = read.csv('involved_factors.csv', header = F)
keep = c(involved_factors$V1, "SampleType")

metadata_2 = metadata_2[,(colnames(metadata_2) %in% keep)]

# Imputation by mice
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
metadata_2 = data

# reads table
reads_table = as.data.frame(t(reads_table))
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)
pvalue <- adonis2(reads_table ~ ., data = metadata_2, method = "bray")
pvalue

row.names(metadata_2) = metadata$SampleID
write.csv(metadata_2,'multi_variable_analysis_input.csv')
write.csv(pvalue,'multi_variable_analysis_pvalue.csv')
##### compare Bray-Curtis between mom-baby cheek and mom-baby rectum #####
keep = !is.na(metadata_16s$Flag)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

beta = beta_diversity(reads_table, metadata = metadata$Flag, factor_name = 'Flag', order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA)
Bray_Curtis_2 = beta$Bray_Curtis
group_dis = gather(Bray_Curtis_2)
group_dis$key2 = rep(colnames(Bray_Curtis_2),ncol(Bray_Curtis_2))
keep = !(group_dis$key == group_dis$key2)
group_dis <- group_dis[keep,]

group_dis$key_2 = NA
group_dis$key2_2 = NA
ID_list = unique(group_dis$key)
for (a in 1:length(ID_list)) {
  n = which(group_dis$key == ID_list[a])
  x = metadata$SampleType[metadata$SampleID == ID_list[a]]
  group_dis$key_2[n] = x
}
ID_list = unique(group_dis$key2)
for (a in 1:length(ID_list)) {
  n = which(group_dis$key2 == ID_list[a])
  x = metadata$SampleType[metadata$SampleID == ID_list[a]]
  group_dis$key2_2[n] = x
}

group_dis$Sample = paste(group_dis$key_2,group_dis$key2_2, sep = '_')

keep = group_dis$Sample == 'BCKD_MCKD' | group_dis$Sample == 'BRCD_MRCD' 
group_dis_2 <- group_dis[keep,]

ggplot(data=group_dis_2, aes(x=Sample, y=value)) +
  geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
  labs(x = 'Sample type', y = "Bray-Curtis distance")+ ylim(0.5,1) + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
ggsave('Bray-Curtis_cheek_rectum.pdf',width=1.2, height=2.5)

group_dis_2$Sample = as.character(group_dis_2$Sample)
wilcox.test(value~Sample, group_dis_2) # where y1 is numeric and A is a factor

# compare Bray-Curtis between baby-paired mom and baby-unpaired mom 
# baby cheek to mom cheek
keep = group_dis$Sample == 'MCKD_BCKD'
group_dis_2 = group_dis[keep,]

group_dis_2$Participant = 'N'
group_dis_2$Participant_m = NA
group_dis_2$Participant_b = NA

Participant_list = unique(metadata$ParticipantID)
for (a in 1: length(Participant_list)) {
  n=which(metadata$ParticipantID == Participant_list[a])
  Sample_list = metadata$SampleID[n]
  
  m = group_dis_2$key %in% Sample_list
  group_dis_2$Participant_m[m] = Participant_list[a]
  
  m = group_dis_2$key2 %in% Sample_list
  group_dis_2$Participant_b[m] = Participant_list[a]
}

keep = group_dis_2$Participant_b == group_dis_2$Participant_m
group_dis_2$Participant[keep] = 'Y'


ggplot(data=group_dis_2, aes(x=Participant, y=value)) +
  geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
  labs(x = 'Paired mom', y = "Bray-Curtis distance\nfrom baby to mom")+ ylim(0.5,1) + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7))+ theme_bw() 
ggsave('Bray-Curtis_baby-cheek_paired_unpaired_mom_cheek.pdf',width=1.2, height=2.5)

group_dis_2$Participant = as.character(group_dis_2$Participant)
wilcox.test(value~Participant, group_dis_2) # where y1 is numeric and A is a factor

# baby cheek to mom vagina
keep = group_dis$Sample == 'MV1D_BCKD'
group_dis_2 = group_dis[keep,]

group_dis_2$Participant = 'N'
group_dis_2$Participant_m = NA
group_dis_2$Participant_b = NA

Participant_list = unique(metadata$ParticipantID)
for (a in 1: length(Participant_list)) {
  n=which(metadata$ParticipantID == Participant_list[a])
  Sample_list = metadata$SampleID[n]
  
  m = group_dis_2$key %in% Sample_list
  group_dis_2$Participant_m[m] = Participant_list[a]
  
  m = group_dis_2$key2 %in% Sample_list
  group_dis_2$Participant_b[m] = Participant_list[a]
}

keep = group_dis_2$Participant_b == group_dis_2$Participant_m
group_dis_2$Participant[keep] = 'Y'


ggplot(data=group_dis_2, aes(x=Participant, y=value)) +
  geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
  labs(x = 'Paired mom', y = "Bray-Curtis distance\nfrom baby to mom")+ ylim(0.996,1) + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7))+ theme_bw() 
ggsave('Bray-Curtis_baby-cheek_paired_unpaired_mom_vagina.pdf',width=1.2, height=2.5)

group_dis_2$Participant = as.character(group_dis_2$Participant)
wilcox.test(value~Participant, group_dis_2) # where y1 is numeric and A is a factor

# baby cheek to mom rectum
keep = group_dis$Sample == 'MRCD_BCKD'
group_dis_2 = group_dis[keep,]

group_dis_2$Participant = 'N'
group_dis_2$Participant_m = NA
group_dis_2$Participant_b = NA

Participant_list = unique(metadata$ParticipantID)
for (a in 1: length(Participant_list)) {
  n=which(metadata$ParticipantID == Participant_list[a])
  Sample_list = metadata$SampleID[n]
  
  m = group_dis_2$key %in% Sample_list
  group_dis_2$Participant_m[m] = Participant_list[a]
  
  m = group_dis_2$key2 %in% Sample_list
  group_dis_2$Participant_b[m] = Participant_list[a]
}

keep = group_dis_2$Participant_b == group_dis_2$Participant_m
group_dis_2$Participant[keep] = 'Y'


ggplot(data=group_dis_2, aes(x=Participant, y=value)) +
  geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
  labs(x = 'Paired mom', y = "Bray-Curtis distance\nfrom baby to mom")+ ylim(0.996,1) + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7))+ theme_bw() 
ggsave('Bray-Curtis_baby-cheek_paired_unpaired_mom_rectum.pdf',width=1.2, height=2.5)

group_dis_2$Participant = as.character(group_dis_2$Participant)
wilcox.test(value~Participant, group_dis_2) # where y1 is numeric and A is a factor

##### factors that impact the distance from baby cheek to mom cheek ##################
keep = !is.na(metadata_16s$Flag)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

# beta_diversity
beta_result = beta_diversity(reads_table, metadata = metadata$Flag, factor_name = 'Flag', order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA)

metadata[metadata == 'Nelson 6th Floor'] = NA

Bray_Curtis_2 = beta_result$Bray_Curtis
group_dis = gather(Bray_Curtis_2)
group_dis$key2 = rep(colnames(Bray_Curtis_2),ncol(Bray_Curtis_2))

group_dis$Sample = paste(str_sub(group_dis$key,1,4),str_sub(group_dis$key2,1,4), sep = '_')


type1_all = c('BCKD_MCKD','BRCD_MRCD')
for (type_num in 1: length(type1_all)) {
  keep = group_dis$Sample == type1_all[type_num]
  group_dis_2 <- group_dis[keep,]
  
  group_dis_2$Participant1 = NA
  group_dis_2$Participant2 = NA
  
  b = unique(c(group_dis_2$key, group_dis_2$key2))
  
  for (a in 1: length(b)) {
    n = which(metadata$SampleID == b[a])
    m = which(group_dis_2$key == b[a])
    
    group_dis_2$Participant1[m] = metadata$ParticipantID[n]
    
    m = which(group_dis_2$key2 == b[a])
    group_dis_2$Participant2[m] = metadata$ParticipantID[n]
  }
  
  keep = group_dis_2$Participant1 == group_dis_2$Participant2
  group_dis_2 = group_dis_2[keep,]
  
  metadata_2 = metadata[((metadata$days_rel2birth == '1') | (metadata$days_rel2birth == '2') | (metadata$days_rel2birth == '0')) & 
                          !is.na(metadata$days_rel2birth) & 
                          metadata$SampleType == str_replace_all(type1_all[type_num],'_.*',''),]
  
  group_dis_2$Mom_visit = NA
  
  for (a in 1: nrow(group_dis_2)) {
    n = which(metadata$SampleID == group_dis_2$key2[a]) 
    group_dis_2$Mom_visit[a] = metadata$VisitNum[n]
  }
  group_dis_2 = group_dis_2[order(group_dis_2$Mom_visit),]
  group_dis_2 = group_dis_2[!duplicated(group_dis_2$key),]
  group_dis_2 = group_dis_2[group_dis_2$key %in% metadata_2$SampleID,]
  
  metadata_2 = metadata_2[metadata_2$SampleID %in% group_dis_2$key,]
  metadata_2 = metadata_2[order(metadata_2$SampleID),]
  group_dis_2 = group_dis_2[order(group_dis_2$key),]
  
  metadata_2_cha = as.data.frame(matrix(data = NA, nrow = nrow(metadata_2), ncol =0))
  metadata_2_num = as.data.frame(matrix(data = NA, nrow = nrow(metadata_2), ncol =0))
  
  n =1
  m =1
  metadata_2[,1:ncol(metadata_2)] <- lapply(metadata_2[,1:ncol(metadata_2)],as.character)
  
  c = c('0','1','2','3','4','5','6','7','8','9','.')
  
  for (a in 1:ncol(metadata_2)) {
    b = metadata_2[,a]
    b = b[!is.na(b)]
    b = strsplit(b,'*')
    b = unlist(b)
    b = unique(b)
    keep = b %in% c
    
    if (sum(!keep) == 0) {
      metadata_2[,a] <- as.numeric(metadata_2[,a])
      metadata_2_num <- cbind(metadata_2_num,metadata_2[,a] )
      colnames(metadata_2_num)[n] = colnames(metadata_2)[a]
      n=n+1
    } else {
      metadata_2_cha <- cbind(metadata_2_cha,metadata_2[,a])
      colnames(metadata_2_cha)[m] = colnames(metadata_2)[a]
      m=m+1
    }
  }
  
  # calculate numeric factor's correlation
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata_2_num) , ncol = 4))
  colnames(correlation) = c('Factor','pearson_pvalue','pearson_rvalue','pearson_adj-pvalue')
  correlation$Factor = colnames(metadata_2_num)
  
  for (a in 1: nrow(correlation)) {
    data1 = group_dis_2$value
    data2 = metadata_2_num[,a]
    
    keep = (!is.na(data1)) & (!is.na(data2))
    data1 = data1[keep]
    data2 = data2[keep]
    
    if (length(data1) <5) {
      next
    }
    
    data3 = cor.test(data1, data2, method = "pearson")
    correlation$pearson_pvalue[a] = data3$p.value
    correlation$pearson_rvalue[a] = data3$estimate
    
    # a=9 for rectum 'days_rel2birth'
    if (a == 9) {
      data = as.data.frame(cbind(data1,data2))
      colnames(data) = c('BC_distance','factor')
      data$factor = as.character(data$factor)
      ggplot(data=data, aes(x=factor, y=BC_distance))+
        geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
        geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5)+
        xlab(paste0(correlation$Factor[a])) +
        ylab(paste0("Bray-Curtis distance of\nbaby cheek to mom cheek")) +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
      ggsave(paste0('BC_distance_baby_cheek_2_mom_cheek_baby_',correlation$Factor[a],'.pdf'), width=2, height=2.5)
    }
    
    if (!is.na(correlation$pearson_pvalue[a])) {
      if (correlation$pearson_pvalue[a] <= 0.05) {
        data = as.data.frame(cbind(data1,data2))
        colnames(data) = c('BC_distance','factor')
        
        ggplot(data, aes(factor, BC_distance))  +
          geom_point(size=0.5) +
          xlab(paste0(correlation$Factor[a])) +
          ylab(paste0("Bray-Curtis distance of\nbaby cheek to mom cheek")) +
          geom_smooth(method='lm', se=T)+
          theme(axis.title = element_text(size = 7), 
                axis.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))+ theme_bw() 
        ggsave(paste0('BC_distance_baby_cheek_2_mom_cheek_baby_',correlation$Factor[a],'.pdf'), width=2, height=2.5)
      }
    }
    
  }
  
  
  correlation = correlation[!is.na(correlation$pearson_pvalue),]
  otu.pval <- adjust.p(correlation$pearson_pvalue, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  correlation$`pearson_adj-pvalue` <- otu.pval$adjusted.p
  
  
  write.csv(correlation,paste0('BC_distance_baby_cheek_2_mom_cheek_day0-2_',type1_all[type_num],'_baby_num.csv'))
  
  
  # calculate character factor's correlation
  correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata_2_cha) , ncol = 3))
  colnames(correlation) = c('Factor','kruskal.test-pvalue','kruskal.test-adj-p')
  correlation$Factor = colnames(metadata_2_cha)
  for (a in 1: nrow(metadata_2_cha)) {
    V1 = group_dis_2$value
    V2 = metadata_2_cha[,a]
    
    data = as.data.frame(cbind(V1,V2))
    data = data[!is.na(data$V2),]
    
    unique_factor = unique(data$V2)
    for (b in 1: length(unique_factor)) {
      unique_factor[b] = sum(data$V2 == unique_factor[b])
    }
    unique_factor = unique_factor[order(unique_factor, decreasing = T)]
    unique_factor = as.numeric(unique_factor)
    if (length(unique(data$V2)) < 2 | unique_factor[1] < 1 | unique_factor[2] < 1) {
      next
    }
    data$V1 = as.numeric(as.character(data$V1))
    correlation$`kruskal.test-pvalue`[a] =   kruskal.test(V1~V2, data)$p.value
    
    if ((correlation$`kruskal.test-pvalue`[a] <= 0.05 & !is.na(correlation$`kruskal.test-pvalue`[a])) | 
        colnames(metadata_2_cha)[a] == 'delivery_vaginal ' | colnames(metadata_2_cha)[a] == 'delivery_csection') {
      colnames(data) = c('BC_distance','factor')
      
      ggplot(data=data, aes(x=factor, y=BC_distance))+
        geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
        geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5)+
        labs(x = colnames(metadata_2_cha)[a], y = 
               paste0('Bray-Curtis distance of\nbaby ',str_replace_all(type1_all[type_num],'_.*',''),' to mom ',
                      str_replace_all(type1_all[type_num],'_.*','')))+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw() 
      ggsave(paste0('BC_distance_baby_cheek_2_mom_cheek_day0-2_',type1_all[type_num],'_',
                    colnames(metadata_2_cha)[a],'.pdf'), width=2, height=2.5)
      
    }
  }
  
  correlation = correlation[!is.na(correlation$`kruskal.test-pvalue`),]
  otu.pval <- adjust.p(correlation$`kruskal.test-pvalue`, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  correlation$`kruskal.test-adj-p` <- otu.pval$adjusted.p
  
  
  write.csv(correlation,paste0('BC_distance_baby_cheek_2_mom_cheek_day0-2_',type1_all[type_num],'_baby_cha.csv'))
  
  
}

##### two-way ANOVA to test impact of another factor on the interaction of 'days_rel2birth' and distance between NB and MB #####
type_num =1

keep = group_dis$Sample == type1_all[type_num]
group_dis_2 <- group_dis[keep,]

group_dis_2$Participant1 = NA
group_dis_2$Participant2 = NA

b = unique(c(group_dis_2$key, group_dis_2$key2))

for (a in 1: length(b)) {
  n = which(metadata$SampleID == b[a])
  m = which(group_dis_2$key == b[a])
  
  group_dis_2$Participant1[m] = metadata$ParticipantID[n]
  
  m = which(group_dis_2$key2 == b[a])
  group_dis_2$Participant2[m] = metadata$ParticipantID[n]
}

keep = group_dis_2$Participant1 == group_dis_2$Participant2
group_dis_2 = group_dis_2[keep,]

metadata_2 = metadata[((metadata$days_rel2birth == '1') | (metadata$days_rel2birth == '2') | (metadata$days_rel2birth == '0')) & 
                        !is.na(metadata$days_rel2birth) & 
                        metadata$SampleType == str_replace_all(type1_all[type_num],'_.*',''),]

group_dis_2$Mom_visit = NA

for (a in 1: nrow(group_dis_2)) {
  n = which(metadata$SampleID == group_dis_2$key2[a]) 
  group_dis_2$Mom_visit[a] = metadata$VisitNum[n]
}
group_dis_2 = group_dis_2[order(group_dis_2$Mom_visit),]
group_dis_2 = group_dis_2[!duplicated(group_dis_2$key),]
group_dis_2 = group_dis_2[group_dis_2$key %in% metadata_2$SampleID,]

metadata_2 = metadata_2[metadata_2$SampleID %in% group_dis_2$key,]
metadata_2 = metadata_2[order(metadata_2$SampleID),]
group_dis_2 = group_dis_2[order(group_dis_2$key),]

# metadata_2 weight and ga_at_delivery to character form
metadata_2$weight_2 = 'Medium'
metadata_2$weight_2[metadata_2$weight <= 96] = 'Low'
metadata_2$weight_2[metadata_2$weight >= 116] = 'High'

metadata_2$ga_at_delivery_2 = 'Medium'
metadata_2$ga_at_delivery_2[metadata_2$ga_at_delivery <= 268] = 'Low'
metadata_2$ga_at_delivery_2[metadata_2$ga_at_delivery >= 279] = 'High'

metadata_2_num = as.data.frame(matrix(data = NA, nrow = nrow(metadata_2), ncol =0))
metadata_2_cha = as.data.frame(matrix(data = NA, nrow = nrow(metadata_2), ncol =0))
n =1
m =1
metadata_2[,1:ncol(metadata_2)] <- lapply(metadata_2[,1:ncol(metadata_2)],as.character)

c = c('0','1','2','3','4','5','6','7','8','9','.')

for (a in 1:ncol(metadata_2)) {
  b = metadata_2[,a]
  b = b[!is.na(b)]
  b = strsplit(b,'*')
  b = unlist(b)
  b = unique(b)
  keep = b %in% c
  
  if (sum(!keep) == 0) {
    metadata_2[,a] <- as.numeric(metadata_2[,a])
    metadata_2_num <- cbind(metadata_2_num,metadata_2[,a] )
    colnames(metadata_2_num)[n] = colnames(metadata_2)[a]
    n=n+1
  } else {
    metadata_2_cha <- cbind(metadata_2_cha,metadata_2[,a])
    colnames(metadata_2_cha)[m] = colnames(metadata_2)[a]
    m=m+1
  }
}

### calculate numeric factor's correlation ###
library('dplyr')
library("ggpubr")

correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata_2_num) , ncol = 4))
colnames(correlation) = c('Factor','Two-way_anova_first_factor_pvalue','Two-way_anova_second_factor_pvalue','Two-way_anova_interaction_pvalue')
correlation$Factor = colnames(metadata_2_num)

for (a in 1: nrow(correlation)) {
  data1 = group_dis_2$value
  data2 = metadata_2_num[,a]
  data_days_rel2birth = metadata_2_num[,9]  # 'days_rel2birth'
  
  keep = (!is.na(data1)) & (!is.na(data2)) & (!is.na(data_days_rel2birth))
  data1 = data1[keep]
  data2 = data2[keep]
  data_days_rel2birth = data_days_rel2birth[keep]
  
  if (length(data1) <5) {
    next
  }
  
  data3 = cbind(data1, data2, data_days_rel2birth)
  data3 = as.data.frame(data3)
  res.aov2 <- aov(data1 ~ data2 + data_days_rel2birth  + data2:data_days_rel2birth, data = data3)
  res.aov2 = summary(res.aov2)
  res.aov2 = res.aov2[[1]]
  
  correlation$`Two-way_anova_first_factor_pvalue`[a] = res.aov2$`Pr(>F)`[2]
  correlation$`Two-way_anova_second_factor_pvalue`[a] = res.aov2$`Pr(>F)`[1]
  correlation$`Two-way_anova_interaction_pvalue`[a] = res.aov2$`Pr(>F)`[3]
  
  if (is.na(res.aov2$`Pr(>F)`[3])) {
    next
  }
  
  if (res.aov2$`Pr(>F)`[3] <= 0.05) {
    data3$data_days_rel2birth = as.factor(data3$data_days_rel2birth)
    ggplot(data3, aes(x = data2, y = data1, color = data_days_rel2birth))  +
      geom_point(size=0.3) + theme_bw() 
    ggsave(paste0('Distance_NB-MB_two_way_ANOVA_with_days_rel2birth_',colnames(metadata_2_num)[a],'.pdf'), width=3, height=3)
    
  }
  
}
write.csv(correlation,'Distance_NB-MB_two_way_ANOVA_with_days_rel2birth_num.csv')

### calculate character factor's correlation ##
correlation = as.data.frame(matrix(data = NA, nrow = ncol(metadata_2_cha) , ncol = 4))
colnames(correlation) = c('Factor','Two-way_anova_first_factor_pvalue','Two-way_anova_second_factor_pvalue','Two-way_anova_interaction_pvalue')
correlation$Factor = colnames(metadata_2_cha)



for (a in 1: ncol(metadata_2_cha)) {
  V1 = group_dis_2$value
  V2 = metadata_2_cha[,a]
  #  V3 = as.character(metadata_2_num[,9])
  V3 = metadata_2_num[,9]
  
  data = as.data.frame(cbind(V1,V2, V3))
  data = data[!is.na(data$V2) & !is.na(data$V3),]
  
  data$factor_type = paste(data$V2,data$V3)
  
  unique_factor = unique(data$factor_type)
  
  for (b in 1: length(unique_factor)) {
    unique_factor[b] = sum(data$factor_type == unique_factor[b])
  }
  unique_factor = unique_factor[order(unique_factor, decreasing = T)]
  unique_factor = as.numeric(unique_factor)
  if (length(unique(data$V2)) < 2 | sum(unique_factor < 5) > 0) {    ### any factor with case number < 5 will not be compared
    next
  }
  
  data$V1 = as.numeric(as.character(data$V1))
  
  res.aov2 <- aov(V1 ~ V2 + V3 + V2:V3, data = data)
  res.aov2 = summary(res.aov2)
  res.aov2
  res.aov2 = res.aov2[[1]]
  
  correlation$`Two-way_anova_first_factor_pvalue`[a] = res.aov2$`Pr(>F)`[2]
  correlation$`Two-way_anova_second_factor_pvalue`[a] = res.aov2$`Pr(>F)`[1]
  correlation$`Two-way_anova_interaction_pvalue`[a] = res.aov2$`Pr(>F)`[3]
  
  if (is.na(res.aov2$`Pr(>F)`[3])) {
    next
  }
  
  if (res.aov2$`Pr(>F)`[3] <= 0.05) {
    data$V3 = as.factor(data$V3)
    
    ggboxplot(data, x = "V3", y = "V1", color = "V2", outlier.shape = NA)
    ggsave(paste0('Distance_NB-MB_two_way_ANOVA_with_days_rel2birth_',colnames(metadata_2_cha)[a],'.pdf'), width=3, height=3)
    
  }
  
  if (colnames(metadata_2_cha)[a] == 'vdischarg_pregnancy') {  # convert no_vdischarge to vdischarg_pregnancy
    data$V3 = as.factor(data$V3)
    
    ggboxplot(data, x = "V3", y = "V1", color = "V2", outlier.shape = NA)
    ggsave('Distance_NB-MB_two_way_ANOVA_with_days_rel2birth_vdischarge_pregnancy.pdf', width=3, height=3)
    
  }
  
}
write.csv(correlation,'Distance_NB-MB_two_way_ANOVA_with_days_rel2birth_cha.csv')

# abundance difference of miscarriage_or_stillbirth, no_vdischarge, vodor_1sttri
keep = !is.na(metadata_16s$miscarriage_or_stillbirth)
sum(keep)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

keep = metadata$SampleType == "MV1D"
sum(keep)
metadata = metadata[keep,]
reads_table = reads_table[,keep]

metadata = metadata[order(metadata$VisitNum, decreasing = T),]
metadata = metadata[order(metadata$ParticipantID),]
reads_table = reads_table[metadata$SampleID]

keep = duplicated(metadata$ParticipantID)
metadata = metadata[!keep,]
reads_table = reads_table[,!keep]
metadata = metadata$miscarriage_or_stillbirth

p = dif_abundance3(reads_table,metadata)
write.csv(p,'Abundance_miscarriage_or_stillbirth.csv')

reads_table <- t(reads_table)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)
metadata = as.data.frame(metadata)
pvalue <- adonis2(reads_table ~ ., data = metadata, method = "bray")  
pvalue

#
keep = !is.na(metadata_16s$vdischarge_pregnancy)
sum(keep)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

keep = metadata$SampleType == "MV1D"
sum(keep)
metadata = metadata[keep,]
reads_table = reads_table[,keep]

metadata = metadata[order(metadata$VisitNum, decreasing = T),]
metadata = metadata[order(metadata$ParticipantID),]
reads_table = reads_table[metadata$SampleID]

keep = duplicated(metadata$ParticipantID)
metadata = metadata[!keep,]
reads_table = reads_table[,!keep]
metadata = metadata$vdischarge_pregnancy

p = dif_abundance3(reads_table,metadata)
write.csv(p,'Abundance_vdischarge_pregnancy.csv')

reads_table <- t(reads_table)
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)
metadata = as.data.frame(metadata)
pvalue <- adonis2(reads_table ~ ., data = metadata, method = "bray")  
pvalue

#
keep = !is.na(metadata_16s$vodor_1sttri)
sum(keep)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

keep = metadata$SampleType == "MV1D"
sum(keep)
metadata = metadata[keep,]
reads_table = reads_table[,keep]

metadata = metadata[order(metadata$VisitNum, decreasing = T),]
metadata = metadata[order(metadata$ParticipantID),]
reads_table = reads_table[metadata$SampleID]

keep = duplicated(metadata$ParticipantID)
metadata = metadata[!keep,]
reads_table = reads_table[,!keep]
metadata = metadata$vodor_1sttri

p = dif_abundance3(reads_table,metadata)
write.csv(p,'Abundance_vodor_1sttri.csv')


##### abundance change on day 0 by miscarriage_or_stillbirth and vdischarge_pregnancy #####
keep = metadata_16s$SampleID %in% colnames(reads_table_16s_BCKD)
metadata = metadata_16s[keep,]
reads_table = reads_table_16s_BCKD[metadata$SampleID]

for (a in 1:3) {
  if (a == 1) {
    keep = metadata$Flag == "BCKD_birth_0_day" 
  } else if (a ==2) {
    keep = metadata$Flag == "BCKD_birth_1_days" 
  } else {
    keep = metadata$Flag == "BCKD_birth_2_days" 
  }
  reads_table_2 = reads_table[,keep]
  metadata_2 = metadata[keep,]
  
  keep = !is.na(metadata_2$miscarriage_or_stillbirth)
  reads_table_3 = reads_table_2[,keep]
  metadata_3 = metadata_2[keep,]
  p = dif_abundance3(reads_table_3,metadata_3$miscarriage_or_stillbirth)
  write.csv(p$data,paste0('species_dif_abundance_miscarriage_or_stillbirth_day_',a,'.csv'))
  
  
  keep = !is.na(metadata_2$vdischarge_pregnancy)
  reads_table_3 = reads_table_2[,keep]
  metadata_3 = metadata_2[keep,]
  p = dif_abundance3(reads_table_3,metadata_3$vdischarge_pregnancy)
  write.csv(p$data,paste0('species_dif_abundance_vdischarge_pregnancy_day_',a,'.csv'))
  
}


# impact on MV
keep = metadata_16s$SampleID %in% colnames(reads_table_16s_MV1D)
metadata = metadata_16s[keep,]
reads_table = reads_table_16s_MV1D[metadata$SampleID]

keep = !is.na(metadata$miscarriage_or_stillbirth)
reads_table_2 = reads_table[,keep]
metadata_2 = metadata[keep,]

p = dif_abundance3(reads_table_2,metadata_2$miscarriage_or_stillbirth)
write.csv(p$data,'species_dif_abundance_MV1D_miscarriage_or_stillbirth.csv')

metadata_2 = metadata_2$miscarriage_or_stillbirth
metadata_2 = as.data.frame(metadata_2)

reads_table_2 <- t(reads_table_2)
reads_table_2 = Rarefy(reads_table_2, depth = min(rowSums(reads_table_2)))
reads_table_2 <- reads_table_2$otu.tab.rff
reads_table_2 <- as.data.frame(reads_table_2)

pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")  
pvalue

keep = !is.na(metadata$vdischarge_pregnancy)
reads_table_2 = reads_table[,keep]
metadata_2 = metadata[keep,]

p = dif_abundance3(reads_table_2,metadata_2$vdischarge_pregnancy)
write.csv(p$data,'species_dif_abundance_MV1D_vdischarge_pregnancy.csv')

metadata_2 = metadata_2$vdischarge_pregnancy
metadata_2 = as.data.frame(metadata_2)

reads_table_2 <- t(reads_table_2)
reads_table_2 = Rarefy(reads_table_2, depth = min(rowSums(reads_table_2)))
reads_table_2 <- reads_table_2$otu.tab.rff
reads_table_2 <- as.data.frame(reads_table_2)

pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")  
pvalue


##### taxa abundance change in day 0 to day 2 #####
keep = metadata_16s$SampleType == "BCKD"
metadata = metadata_16s[keep,]
metadata_2 = metadata
metadata_2$Flag[metadata_2$Flag == 'BCKD_birth_1_days' | metadata_2$Flag == 'BCKD_birth_2_days'] = 'BCKD_birth_1-2_days'

p = dif_abundance(reads_table_16s_BCKD, metadata_2$Flag,paired_test = F, order_reverse = F, fold_change_th = 1, style = 1, order = c('BCKD_birth_1-2_day','BCKD_birth_0_days'))
p1 = p$data
p1$Type = 'BCKD'
p$p
ggsave(paste('species_dif_abundance_BCKD.pdf',sep='_'),width=5, height=1.5)
ggsave(paste('species_dif_abundance_BCKD_2.pdf',sep='_'),width=5, height=4)
write.csv(p$data,'species_dif_abundance_BCKD.csv')


keep = metadata_16s$SampleType == "BRCD"
metadata = metadata_16s[keep,]
metadata_2 = metadata
metadata_2$Flag[metadata_2$Flag == 'BRCD_birth_1_days' | metadata_2$Flag == 'BRCD_birth_2_days'] = 'BRCD_birth_1-2_days'

p = dif_abundance(reads_table_16s_BRCD,metadata_2$Flag,paired_test = F, order_reverse = F, fold_change_th = 1, style = 1, order = c('BCKD_birth_0_day','BCKD_birth_1-2_days'))
p2 = p$data
p2$Type = 'BRCD'
p$p
ggsave(paste('species_dif_abundance_BRCD.pdf',sep='_'),width=4.5, height=1.1)
ggsave(paste('species_dif_abundance_BRCD_2.pdf',sep='_'),width=4.5, height=1.7)

write.csv(p$data,'species_dif_abundance_BRCD.csv')

p1 = p1[,c(9,13,14)]
p2 = p2[,c(9,13,14)]
p = rbind(p1,p2)
p <- p[abs(p$diff.btw) >= 1 & p$`-Log10(adj-pvalue)` >= -log10(0.05),]
p$taxa = row.names(p)

ggplot(p, aes(taxa, Type)) + 
  geom_point(aes(col=diff.btw, size=`-Log10(adj-pvalue)`)) + 
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab" ) +
  coord_flip() +          # convert x y axis
  labs(x = 'Species')+ 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) + theme_bw() 
ggsave(paste('species_dif_abundance_BRCD&BCKD.pdf',sep='_'),width=7, height=4)


keep = metadata_16s$SampleType == "BS1D"
metadata = metadata_16s[keep,]
metadata_2 = metadata
p = dif_abundance(reads_table_16s_BS1D,metadata_2$Flag,paired_test = F, order_reverse = F, fold_change_th = 1, style = 1, order = c('BCKD_birth_1_days','BCKD_birth_2_days'))
p$p

p = dif_abundance3(reads_table_16s_BS1D,metadata_2$Flag)
write.csv(p$data,'species_dif_abundance_BS1D.csv')









##### MR and PTB #####
keep = !is.na(metadata_16s$`Preterm delivery`) & metadata_16s$SampleType == 'MRCD'
metadata = metadata_16s[keep,]
reads_table = reads_table_16s[,keep]
metadata = cbind(metadata$`Preterm delivery`, metadata$weeks_pregnant) 
metadata = as.data.frame(metadata)
colnames(metadata) = c("Preterm_delivery",'Weeks_pregnancy')
row.names(metadata) = colnames(reads_table)

for (a in 2:7) {
    keep = metadata$Weeks_pregnancy >= a*5 + 1 & metadata$Weeks_pregnancy < (a+1)*5 + 1
    metadata_2 = metadata[keep,]
    reads_table_2 = reads_table[,keep]
    
    data = dif_abundance(reads_table_2,metadata_2$Preterm_delivery, pvalue_th = 0.05, fold_change_th = 1, 
                  paired_test = F, order_reverse = F, style = 1, order = c('No','Yes'))
    
}


library(nlme)
library(tidyverse)
library(compositions)
source("programs/ancom.R")


# Random intercept model adjusted for other covariates
metadata$ID = row.names(metadata)
prepro = feature_table_pre_process(reads_table, metadata, 'ID', 
                                   group_var = 'Preterm_delivery', out_cut = 0.05, 
                                   zero_cut = 0.90, lib_cut = 0, neg_lb = T)

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "Preterm_delivery"; p_adj_method = "BH"; alpha = 0.05
adj_formula = "Weeks_pregnancy"; rand_formula = NULL; list(maxIter = 100, msMaxIter = 100, opt = "optim")

res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula, lme_control)

write_csv(res$out, "ANCOM_Preterm_delivery_adjusted_by_Weeks_pregnancy.csv")


# Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa - 1), 0.8 * (n_taxa - 1), 0.7 * (n_taxa - 1), 0.6 * (n_taxa - 1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

figure_annotation = res$fig$plot_env$dat_fig$taxa_id
keep = res$out$detected_0.7 & !res$fig$plot_env$num_struc_zero

figure_annotation[!keep] = NA
figure_annotation = str_replace_all(figure_annotation, 'cluster', ' cluster')
figure_annotation = str_replace_all(figure_annotation, '_', ' ')
figure_annotation = str_replace_all(figure_annotation, 'Clostridium cluster 2 cluster 3', 'Clostridium cluster 3')

library(ggrepel)
res$fig + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = -18, color = "orange", parse = TRUE)+
  geom_text_repel(aes(label = figure_annotation),size = 3.5, colour = 'black', 
                  min.segment.length = unit(0, 'lines'), 
                  nudge_y = .5) 
ggsave('ANCOM_Preterm_delivery_adjusted_by_Weeks_pregnancy.pdf',width=6, height=6)

# Run '16s_analysis_pipeline_BZ 03242022.R'
write.csv(reads_table,'reads_table.csv')
write.csv(metadata,'metadata.csv')
