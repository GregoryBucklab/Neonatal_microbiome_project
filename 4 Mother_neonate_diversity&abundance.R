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
type2_all = unique(metadata_16s$Flag)

total_abundance = rep(0,nrow(reads_table_16s))
reads_table_16s_abundance = get_abundance_table(reads_table_16s,  mc.cores = 8)

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

ggplot(taxonomy_16s_3, aes(key, Taxa)) + 
  geom_point(aes(size= value), col = '#888888') + 
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7)) + theme_bw() 
ggsave(paste('species_average_abundance.pdf',sep='_'),width=5, height=3)
ggsave(paste('species_average_abundance_2.pdf',sep='_'),width=12, height=4)

ggplot(taxonomy_16s_3, aes(key, Taxa)) + 
  geom_point(aes(size= value,col = value)) + 
  scale_color_gradient(low="#FFE9E9", high="red")+
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7)) + theme_bw() 
ggsave(paste('species_average_abundance.pdf',sep='_'),width=5, height=3)
ggsave(paste('species_average_abundance_2.pdf',sep='_'),width=12, height=4)

abundant_taxa_list = unique(taxonomy_16s_3$Taxa)
write.csv(abundant_taxa_list,'node.csv')

##### baby main taxa network #####
cor_parameter = 0

reads_table = reads_table_16s_BCKD
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BCKD.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BCKD.csv')
data = cor$gephi_input
data = data[data$Source %in% abundant_taxa_list & data$Target %in% abundant_taxa_list,]
write.csv(data,'baby_taxa_network_BCKD_main_taxa.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BCKD_taxa_cluster.csv')

{
  pvalue = 0.05; cor_parameter= 0; style = 1; bar_max = 2; bar_min = -2; pheatmap_fontsize = 5; treeheight = 50; alpha = 0.05
  reads_table = reads_table_16s_BCKD
  keep = rowSums(reads_table) > 0
  reads_table = reads_table[keep,]
  reads_table = as.data.frame(t(reads_table))
  reads_table = reads_table + 0.5
  reads_table <- clr(reads_table)
  reads_table = as.matrix((reads_table))
  
  otu.cor <- rcorr(reads_table, type="spearman")
  
  otu.pval <- forceSymmetric(otu.cor$P)
  otu.pval <- otu.pval@x
  otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  
  p.yes <- otu.pval< pvalue  
  
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  
  p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
  p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.
  
  p.yes.rr[is.na(p.yes.rr)] = 0
  
  keep = abs(colSums(p.yes.rr)) > 0
  p.yes.rr = p.yes.rr[,keep]
  p.yes.rr = p.yes.rr[keep,]
  
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  paletteLength <- 50
  myColor <- colorRampPalette(c("#4D4DFD", "white", "red"))(paletteLength)
  myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
  
  myBreaks <- unique(myBreaks)
  p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                treeheight_row = treeheight, treeheight_col = treeheight)
  
  Taxa_cluster = sort(cutree(p$tree_col, k=3))
  Taxa_cluster = as.data.frame(Taxa_cluster)
  Taxa_cluster = as.data.frame(Taxa_cluster[match(row.names(p.yes.rr),row.names(Taxa_cluster)),])
  row.names(Taxa_cluster) = row.names(p.yes.rr)
  colnames(Taxa_cluster) = 'Taxa_cluster'

  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==1] = 'a'
  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==2] = 'b'
  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==3] = 'c'
  
  mycolors = list(Taxa_cluster = c("a"='blue',"b"='red',"c"='green'))
                  
  pdf("network_BCKD.pdf", width=8, height=8)
  print(pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                 treeheight_row = treeheight, treeheight_col = treeheight,
                 annotation_col = Taxa_cluster,annotation_colors =mycolors,
                 annotation_row = Taxa_cluster))
  dev.off()
}







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
data = cor$gephi_input
data = data[data$Source %in% abundant_taxa_list & data$Target %in% abundant_taxa_list,]
write.csv(data,'baby_taxa_network_BRCD_main_taxa.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BRCD_taxa_cluster.csv')

{
  pvalue = 0.05; cor_parameter= 0; style = 1; bar_max = 2; bar_min = -2; pheatmap_fontsize = 5; treeheight = 50; alpha = 0.05
  reads_table = reads_table_16s_BRCD
  keep = rowSums(reads_table) > 0
  reads_table = reads_table[keep,]
  reads_table = as.data.frame(t(reads_table))
  reads_table = reads_table + 0.5
  reads_table <- clr(reads_table)
  reads_table = as.matrix((reads_table))
  
  otu.cor <- rcorr(reads_table, type="spearman")
  
  otu.pval <- forceSymmetric(otu.cor$P)
  otu.pval <- otu.pval@x
  otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  
  p.yes <- otu.pval< pvalue  
  
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  
  p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
  p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.
  
  p.yes.rr[is.na(p.yes.rr)] = 0
  
  keep = abs(colSums(p.yes.rr)) > 0
  p.yes.rr = p.yes.rr[,keep]
  p.yes.rr = p.yes.rr[keep,]
  
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  paletteLength <- 50
  myColor <- colorRampPalette(c("#4D4DFD", "white", "red"))(paletteLength)
  myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
  
  myBreaks <- unique(myBreaks)
  p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                treeheight_row = treeheight, treeheight_col = treeheight)
  
  Taxa_cluster = sort(cutree(p$tree_col, k=3))
  Taxa_cluster = as.data.frame(Taxa_cluster)
  Taxa_cluster = as.data.frame(Taxa_cluster[match(row.names(p.yes.rr),row.names(Taxa_cluster)),])
  row.names(Taxa_cluster) = row.names(p.yes.rr)
  colnames(Taxa_cluster) = 'Taxa_cluster'
  
  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==1] = 'a'
  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==2] = 'b'
  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==3] = 'c'
  
  mycolors = list(Taxa_cluster = c("a"='green',"b"='blue',"c"='red'))
  
  pdf("network_BRCD.pdf", width=8, height=8)
  pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                 treeheight_row = treeheight, treeheight_col = treeheight,
                 annotation_col = Taxa_cluster,annotation_colors =mycolors,
                 annotation_row = Taxa_cluster)
  dev.off()
}

reads_table = reads_table_16s_BS1D
keep = rowSums(reads_table) > 0
reads_table = reads_table[keep,]
pdf("network_BS1D.pdf",width = 8, height = 8)
cor = newwork_rcorr(reads_table, pvalue = 0.05, cor_parameter= cor_parameter, style = 1, bar_max = 1, bar_min = -1, pheatmap_fontsize = 5, alpha = 0.05)
dev.off()
write.csv(cor$gephi_input,'baby_taxa_network_BS1D.csv')
data = cor$gephi_input
data = data[data$Source %in% abundant_taxa_list & data$Target %in% abundant_taxa_list,]
write.csv(data,'baby_taxa_network_BS1D_main_taxa.csv')
taxa_cluster = sort(cutree(cor$p$tree_row, k=3))
write.csv(taxa_cluster,'baby_taxa_network_BS1D_taxa_cluster.csv')

{
  pvalue = 0.05; cor_parameter= 0; style = 1; bar_max = 2; bar_min = -2; pheatmap_fontsize = 5; treeheight = 50; alpha = 0.05
  reads_table = reads_table_16s_BS1D
  keep = rowSums(reads_table) > 0
  reads_table = reads_table[keep,]
  reads_table = as.data.frame(t(reads_table))
  reads_table = reads_table + 0.5
  reads_table <- clr(reads_table)
  reads_table = as.matrix((reads_table))
  
  otu.cor <- rcorr(reads_table, type="spearman")
  
  otu.pval <- forceSymmetric(otu.cor$P)
  otu.pval <- otu.pval@x
  otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  
  p.yes <- otu.pval< pvalue  
  
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  
  p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
  p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.
  
  p.yes.rr[is.na(p.yes.rr)] = 0
  
  keep = abs(colSums(p.yes.rr)) > 0
  p.yes.rr = p.yes.rr[,keep]
  p.yes.rr = p.yes.rr[keep,]
  
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  paletteLength <- 50
  myColor <- colorRampPalette(c("#4D4DFD", "white", "red"))(paletteLength)
  myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
  
  myBreaks <- unique(myBreaks)
  p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                treeheight_row = treeheight, treeheight_col = treeheight)
  
  Taxa_cluster = sort(cutree(p$tree_col, k=3))
  Taxa_cluster = as.data.frame(Taxa_cluster)
  Taxa_cluster = as.data.frame(Taxa_cluster[match(row.names(p.yes.rr),row.names(Taxa_cluster)),])
  row.names(Taxa_cluster) = row.names(p.yes.rr)
  colnames(Taxa_cluster) = 'Taxa_cluster'
  
  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==1] = 'a'
  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==2] = 'b'
  Taxa_cluster$Taxa_cluster[Taxa_cluster$Taxa_cluster ==3] = 'c'
  
  mycolors = list(Taxa_cluster = c("a"='green',"b"='red',"c"='blue'))
  
  pdf("network_BS1D.pdf", width=8, height=8)
  print(pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                 treeheight_row = treeheight, treeheight_col = treeheight,
                 annotation_col = Taxa_cluster,annotation_colors =mycolors,
                 annotation_row = Taxa_cluster))
  dev.off()
}

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

### MGS sample collection ### 
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
sum((!is.na(output_table$BCKD_birth_0_day) | !is.na(output_table$BCKD_birth_1_days)| !is.na(output_table$BCKD_birth_2_days)) & 
  (!is.na(output_table$BRCD_birth_0_day) | !is.na(output_table$BRCD_birth_1_days)| !is.na(output_table$BRCD_birth_2_days)) & 
  (!is.na(output_table$BS1D_birth_1_days)| !is.na(output_table$BS1D_birth_2_days)))

setwd('/Users/binzhu/Desktop/Mom-baby/')
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
                             order = NA, NMDS_skip = F, ref_group = NA, rarefy_to = NA,pheatmap_y = T)
pdf("beta_diversity_group_dis.pdf",width=3, height=2.7)
beta_result$group_dis_2_p
dev.off()
write.csv(beta_result$group_dis_sig,'group_dis_sig.csv', row.names = T)

# fig 1b new
{
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
  
  group_level = unique(group_dis$Source)
  n = length(group_level)
  distance_median = matrix(data=NA, nrow = n, ncol =n)
  colnames(distance_median) = group_level
  row.names(distance_median) = group_level
  
  group_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
  colnames(group_dis_sig) = group_level
  row.names(group_dis_sig) = group_level
  
  for (a in 1:(n-1)) {
    for (b in (a+1):n) {
      distance_data = group_dis$value[group_dis$Source == row.names(distance_median)[a] & group_dis$Target == colnames(distance_median)[b]]
      distance_median[a,b] <- mean(distance_data)
      distance_median[b,a] <- distance_median[a,b]
      
      keep = metadata == group_level[a] | metadata == group_level[b]
      metadata_2 = as.character(metadata[keep])
      reads_table_2 = reads_table[keep,]
      
      metadata_2 = as.data.frame(metadata_2)
      pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")[1,5] 
      group_dis_sig[b,a] = pvalue
      group_dis_sig[a,b] = pvalue
    }
  }
  
  group_dis_sig_2 = group_dis_sig
  group_dis_sig_2[is.na(group_dis_sig)] = ''
  group_dis_sig_2[group_dis_sig > 0.05] = ''
  group_dis_sig_2[group_dis_sig <= 0.05 & group_dis_sig > 0.01] = '*'
  group_dis_sig_2[group_dis_sig <= 0.01 & group_dis_sig > 0.001] = '**'
  group_dis_sig_2[group_dis_sig <= 0.001] = '***'
  
  #     distance_median[lower.tri(distance_median)] <- NA
  save_pdf <- function(x, filename, width=7, height=7) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  distance_median_cp = distance_median
  group_dis_sig_2_cp = group_dis_sig_2
  
  distance_median = distance_median_cp

  treeheight = 25
  pheatmap_fontsize =6
  group_dis_2_p = pheatmap(distance_median, cluster_rows=TRUE, show_rownames=TRUE, 
                           cluster_cols=T, show_colnames=T, 
                           color=colorRampPalette(c("red","white"))(100),
                           fontsize = pheatmap_fontsize, display_numbers = group_dis_sig_2, 
                           treeheight_row = treeheight, treeheight_col = treeheight)
  
  save_pdf(group_dis_2_p, 'pheatmap.pdf',width=3, height=2.7)
}
  
  


# NMDS
reads_table = as.data.frame(t(reads_table))
reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
reads_table <- reads_table$otu.tab.rff
reads_table <- as.data.frame(reads_table)
Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
Bray_Curtis <- as.data.frame(Bray_Curtis)
NMDS <-
  metaMDS(Bray_Curtis,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 20,
          wascores = TRUE)
mds_data <- as.data.frame(NMDS$points)
mds_data$factor <- metadata$SampleType
mds_data$shape <- metadata$Flag
mds_data$factor = as.factor(mds_data$factor)
mds_data$shape[mds_data$shape == "MCKD" | mds_data$shape == "MV1D" |mds_data$shape == "MRCD"] = 'Mother'
mds_data$shape[str_detect(mds_data$shape, '0_day')] = 'Day 0'
mds_data$shape[str_detect(mds_data$shape, '1_day')] = 'Day 1'
mds_data$shape[str_detect(mds_data$shape, '2_day')] = 'Day 2'

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = factor)) +
  geom_point(size = 1,aes(shape=shape))+
  scale_color_manual(values=color1)+
  stat_ellipse(type = "t")+
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7))+theme_bw()
ggsave('beta_NMDS.pdf',width=6, height=6)



keep = !is.na(metadata_16s$Flag)
reads_table = reads_table_16s[,keep]
metadata = metadata_16s[keep,]

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


pic <- tsne$Y
pic <- data.frame(pic,row.names(reads_table))
colnames(pic) <- c('X1','X2','SampleID')
pic = merge(pic, metadata, 'SampleID')
pic$factor <- metadata$SampleType
pic$shape <- metadata$Flag
pic$factor = as.factor(pic$factor)
pic$shape[pic$shape == "MCKD" | pic$shape == "MV1D" |pic$shape == "MRCD"] = 'Mother'
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



##### taxa abundance change in day 0 to day 2 new #####
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
  labs(x = 'Species')+ 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) + theme_bw() 
ggsave(paste('species_dif_abundance.pdf',sep='_'),width=5.5, height=4)








