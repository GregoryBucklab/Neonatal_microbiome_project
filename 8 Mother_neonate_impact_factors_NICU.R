setwd('/Users/binzhu/Desktop/Mom-baby/')
colnames(metadata_16s)[colnames(metadata_16s) == 'mom_delivery_method (C-section yes / vaginal no)'] = 'mom_delivery_method'

##### get data ######
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

# metadata_mom_MV1D, reads_table_mom_MV1D, metadata_2, metadata_baby$ParticipantID

##### baby metadata and mom metadata in the last visit ##### 
{
  # metadata last visit prepare data
  metadata_mom_NICU = as.data.frame(metadata_2)
  NICU = metadata_mom_NICU$baby_nicu
  
  keep = !is.na(NICU)
  sum(keep)
  metadata_mom_NICU = metadata_mom_NICU[keep,]
  NICU = NICU[keep]
  metadata_mom_NICU$baby_nicu = NULL
  NICU[NICU == 0] = 'No'
  NICU[NICU == 1] = 'Yes'
  
  NICU_test = data.frame(factor_name = colnames(metadata_mom_NICU), P_value = NA, adj_P_value = NA, case_number = NA, With_host_characteristic = NA, Without_host_characteristic = NA)
  for (a in 1:ncol(metadata_mom_NICU)) {
    keep = !is.na(metadata_mom_NICU[,a])
    if (sum(keep) == 0) {
      next
    }
    data_1 = metadata_mom_NICU[keep,a]
    NICU_1 = NICU[keep]
    data = data.frame(V1 = NICU_1, V2 = data_1)
    data$V2 = as.numeric(as.character(data$V2))
    
    if (sum(is.na(data$V2)) > 50) {
      next
    }

    temp = table(data)
    if (dim(temp)[2] > 2) {
      NICU_test$case_number[a] = paste0(sum(keep), ' (',sum(temp[2,]),')')
    } else {
      NICU_test$case_number[a] = paste0(sum(keep), ' (',temp[2,1]+temp[2,2],')')
    }
    
    if (dim(temp)[1] == 2 & dim(temp)[2] ==2) {
        NICU_test$Without_host_characteristic[a] = paste0(temp[1,1] + temp[2,1], ' (',temp[2,1],')')

        NICU_test$With_host_characteristic[a] = paste0(temp[1,2] + temp[2,2], ' (',temp[2,2],')')
    }
    
    if (nrow(data) <case_number_th) {
      next
    }
    if (nrow(data) <case_number_th) {
      next
    }
    if (length(unique(data$V2)) < 2 | length(unique(data$V1)) < 2) {
      next
    }
    data_2 = as.data.frame(table(data))
    if (nrow(data_2) <2) {
      next
    }
    
    NICU_test[a,2] =  wilcox.test(V2~V1, data)$p.value
    if (NICU_test[a,2] <= 0.05) {
      data$V2 = as.numeric(as.character(data$V2))
      ggplot(data=data, aes(x=V1, y=V2)) + geom_violin(trim=T)+
        geom_boxplot(fill='gray', color="black", outlier.shape=NA, width=0.3) +
        geom_point(shape=16, size = 0.5)+
        labs(x = 'NICU', y = NICU_test[a,1])+ 
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_bw()
      ggsave(paste0('NICU_',NICU_test[a,1],'.pdf'),width=1.3, height=3)
    }
    
  }
  otu.pval <- adjust.p(NICU_test[,2], pi0.method="bky", alpha = 0.05,pz = 0.05)
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  otu.pval
  NICU_test[,3] = otu.pval
  NICU_test = NICU_test[order(NICU_test$P_value,decreasing = F),]
  NICU_test = NICU_test[!str_detect(NICU_test$factor_name, 'worri|baby'),]
  write.csv(NICU_test,'NICU_mom_metadata_last_visit.csv',  row.names = F)
}



##### NICU ML model last visit #####
{
  # metadata last visit prepare data
  metadata_mom_NICU = as.data.frame(metadata_2)
  NICU = metadata_mom_NICU$baby_nicu
  
  keep = !is.na(NICU)
  sum(keep)
  metadata_mom_NICU = metadata_mom_NICU[keep,]
  NICU = NICU[keep]
  metadata_mom_NICU$baby_nicu = NULL
  NICU[NICU == 0] = 'No'
  NICU[NICU == 1] = 'Yes'

}
{
  # imputation
  {
    library(mice)
    data = as.data.frame(metadata_mom_NICU)
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
        m = which((row.names(data) == row.names(x.impmi_3)[c]))
        data[m,n]= x.impmi_3$result[c]
      }
    }
    
    metadata_mom_imputation = data
    
    for (a in 1: ncol(metadata_mom_imputation)) {
      metadata_mom_imputation[,a][is.na(metadata_mom_imputation[,a])] <-  median(metadata_mom_imputation[,a][!is.na(metadata_mom_imputation[,a])])
    }
    
    metadata_mom_imputation = as.data.frame(scale(metadata_mom_imputation))
    
    metadata_mom_imputation = as.data.frame(metadata_mom_imputation)
    
  }
  
  # modeling
  {
    library(pROC)
    library(doSNOW)
    library(parallel)
    library(matrixStats)
    library(erer) # export list
    
    importance_ML = list()
    myPvs = vector()
    myPvs2 = vector()
    myPvs3 = vector()
    myPvs4 = vector()
    
    n_ML = 0
    
    n_ML = n_ML+1
    output = NICU
    output = as.factor(output)
    input = metadata_mom_imputation
    
    set.seed(2022)
    
    input = FilterFeatures(input, output)
    
    input = input[, !str_detect(colnames(input), 'baby') ]
    
    data = xxx_get_result(input, output, ntree= ncol(input)*2) 
    
    data$p_roc
    ggsave('icu_last_visit_M_meta_&_M_microbiome_Sensitivities_Specificities.pdf',width=3, height=3)
    
    myPvs[n_ML]=data$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
    myPvs2[n_ML]=data$pvalue   # test correlation between expected values and true values (PTB yes or no).
    myPvs3[n_ML] = data$err_median
    importance_ML[[n_ML]] = data$importance_ML
    myPvs4[n_ML] = paste0(length(output),' (',sum(output == 'Yes'),')')
    myPvs[n_ML]
    
    model_name ='metadata_1st_visit'
    
    data = cbind(model_name,myPvs,myPvs2, myPvs4)
    data = as.data.frame(data)
    colnames(data) = c('Model_name','AUC_ROC_curve','p_value','Case_number')
    write.csv(data,'icu_M_meta_&_M_microbiome.csv')
    names(importance_ML) = 'metadata_1st_visit'
    library(erer)
    write.list(importance_ML, 'importance_ML_icu_last_visit_M_meta_&_M_microbiome_.csv')
  }
}
##### baby metadata and microbiome ##### 
sample_list = vector()
# NICU baby
{
  # output a file for lefse
  taxonomy_16s$temp = taxonomy_16s$Taxonomy
  taxonomy_16s$temp = str_remove_all(taxonomy_16s$temp, ';s__.*')
  taxonomy_16s$lefse = str_replace_all(taxonomy_16s$temp,';p__|;c__|;o__|;f__|;g__','|')
  taxonomy_16s$lefse = str_replace_all(taxonomy_16s$lefse,'k__','')
  taxonomy_16s$lefse = paste0(taxonomy_16s$lefse,'|',taxonomy_16s$Taxa)
  taxonomy_16s$lefse = str_replace_all(taxonomy_16s$lefse,'\\|\\|\\|\\||\\|\\|\\||\\|\\|','\\|')
}
type1_all = unique(metadata_16s$SampleType)[c(1,5,6)]
for (a in 1: length(type1_all)) {

  # find paired baby microbiome and baby's metadata (time closest between mom and baby) and prepare reads tables
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
    
    sample_list = c(sample_list, colnames(reads_table_baby))
  }
  
  # nicu
  {
    keep = !is.na(metadata_baby$baby_nicu)
    NICU = metadata_baby$baby_nicu[keep]
    reads_table_baby = reads_table_baby[,keep]
    NICU = data.frame(Sample_ID = colnames(reads_table_baby), NICU = NICU)

  
    sum(NICU == 'Yes')
    sum(NICU == 'No')
    nrow(NICU)
  }
  
  setwd('/Users/binzhu/Desktop/Mom-baby/16s_analysis_pipeline_BZ-main')
  write.csv(NICU, 'metadata.csv',row.names = F)
  write.csv(reads_table_baby, 'reads_table.csv')

  # lefse output
  {
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 0, species_threshold = 0.001, mc.cores = 8)

    data = vector()
    
    for (b in 1:nrow(reads_table_baby)) {
      data[b] = taxonomy_16s$lefse[which(taxonomy_16s$Taxa == row.names(reads_table_baby)[b])]
    }
    reads_table_baby = cbind(data,reads_table_baby)
    
    
    output_lefse = rbind(c(NA,NICU$NICU),c(NA,colnames(reads_table_baby)), reads_table_baby)
    output_lefse[1,1] = 'NICU'; output_lefse[2,1] = 'id'
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    write.table(output_lefse, paste0(type1_all[a],'_reads_table.txt'),quote = F,col.names = F,row.names = F, sep = '\t')
  }
  
  
  ### GO TO '16s_analysis_pipeline_BinZhu_09282022.R' ###
}
# NICU mom last visit
type1_all = unique(metadata_16s$SampleType)[c(2,3,4)]
for (a in 1: length(type1_all)) {
  
  # find paired baby microbiome and baby's metadata (time closest between mom and baby) and prepare reads tables
  {
    keep = metadata_16s$SampleType == type1_all[a] & 
      !is.na(metadata_16s$SampleType)
    sum(keep)
    metadata_mom = metadata_16s[keep,]
    metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = T),]
    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
    keep = !duplicated(metadata_mom$ParticipantID)
    metadata_mom = metadata_mom[keep,]      
    
    Participant_list = metadata_mom$ParticipantID
    keep = !is.na(metadata_16s$SampleType) & metadata_16s$Mombaby == "Baby"
    sum(keep)
    metadata_baby = metadata_16s[keep & metadata_16s$ParticipantID %in% Participant_list,]
    metadata_baby = metadata_baby[order(metadata_baby$days_rel2birth),]
    keep = duplicated(metadata_baby$ParticipantID)
    metadata_baby =metadata_baby[!keep,]

    metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]   # mom's metadata 
    metadata_baby = metadata_baby[order(metadata_baby$ParticipantID),]   # baby's metadata

    keep = colnames(reads_table_16s) %in% metadata_mom$SampleID
    reads_table_baby = reads_table_16s[,keep]
    reads_table_baby = reads_table_baby[metadata_mom$SampleID]
    
    sample_list = c(sample_list, colnames(reads_table_baby))
  }
  
  # nicu
  {
    keep = !is.na(metadata_baby$baby_nicu)
    NICU = metadata_baby$baby_nicu[keep]
    reads_table_baby = reads_table_baby[,keep]
    NICU = data.frame(Sample_ID = colnames(reads_table_baby), NICU = NICU)
  }
  
  sum(NICU == 'Yes')
  sum(NICU == 'No')
  nrow(NICU)
  
  setwd('/Users/binzhu/Desktop/Mom-baby/16s_analysis_pipeline_BZ-main')
  write.csv(NICU, 'metadata.csv',row.names = F)
  write.csv(reads_table_baby, 'reads_table.csv')
  
  # lefse output
  {
    reads_table_baby = prepare_reads_table_2(reads_table_baby, total_reads_threshold = 0, species_threshold = 0.001, mc.cores = 8)
    data = vector()
    
    for (b in 1:nrow(reads_table_baby)) {
      data[b] = taxonomy_16s$lefse[which(taxonomy_16s$Taxa == row.names(reads_table_baby)[b])]
    }
    reads_table_baby = cbind(data,reads_table_baby)
    
    
    output_lefse = rbind(c(NA,NICU$NICU),c(NA,colnames(reads_table_baby)), reads_table_baby)
    output_lefse[1,1] = 'NICU'; output_lefse[2,1] = 'id'
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    write.table(output_lefse, paste0(type1_all[a],'_reads_table.txt'),quote = F,col.names = F,row.names = F, sep = '\t')
  }
  ### GO TO '16s_analysis_pipeline_BinZhu_09282022.R' ###
}


 



##### lefse visualization #####
setwd('/Users/binzhu/Desktop/Mom-baby/lefse')
file_list = list.files(path = ".", pattern = '.lefse_internal_res')

for (a in 1:length(file_list)) {
  data_1 = read.delim(file_list[a], header = F)
  data_1 = data_1[,c(1,2,5,3)]
  colnames(data_1) = c('Taxa','LDA_SCORE','pvalue','NICU')
  data_1$Microbiome = str_remove(file_list[a], '.lefse_internal_res')
  
  if (a == 1) {
    data = data_1
  } else {
    data = rbind(data , data_1)
  }
}
setwd('/Users/binzhu/Desktop/Mom-baby/')

data$pvalue[data$pvalue == '-'] = 1
data$pvalue = as.numeric(as.character(data$pvalue)); data$LDA_SCORE = as.numeric(as.character(data$LDA_SCORE))
write.csv(data,'lefse.csv')

data_1 = data[data$pvalue <= 0.001 & data$LDA_SCORE >= 3,]
data_1$LDA_SCORE[data_1$NICU == 'No'] = -data_1$LDA_SCORE[data_1$NICU == 'No']
data_1$pvalue_2 = -log10(data_1$pvalue)
data_1$Taxa_2 = str_replace_all(data_1$Taxa,'.*\\.','')

ggplot(data_1, aes(Taxa_2, Microbiome)) + 
  geom_point(aes(col=LDA_SCORE, size=pvalue_2)) + 
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab" ) +
  coord_flip() +theme_bw()+          # convert x y axis
  labs(x = 'Taxa')+ 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12))
ggsave(paste('NICU_lefse.pdf',sep='_'),width=12, height=10)

data_2 = data_1[data_1$Microbiome == 'MB',]
data_2 = data_2[order(data_2$LDA_SCORE, decreasing = F),]
data_2$Taxa_2 = factor(data_2$Taxa_2, levels = data_2$Taxa_2)

ggplot(data_2, aes(x = Taxa_2, y = LDA_SCORE,, fill = NICU))+
  geom_bar(stat="identity", width=0.5)+
  labs(x = '', y = "LDA score (Log 10)")+ 
  coord_flip()
ggsave(paste('NICU_lefse_MB.pdf',sep='_'),width=5, height=1.3)

data_2 = data_1[data_1$Microbiome == 'NB',]
data_2 = data_2[order(data_2$LDA_SCORE, decreasing = F),]
data_2$Taxa_2 = factor(data_2$Taxa_2, levels = data_2$Taxa_2)

data_2$NICU = factor(data_2$NICU, levels = c('Yes','No'))
ggplot(data_2, aes(x = Taxa_2, y = LDA_SCORE, fill = NICU))+
  geom_bar(stat="identity", width=0.5)+
  labs(x = '', y = "LDA score (Log 10)")+ 
  coord_flip()
ggsave(paste('NICU_lefse_NB.pdf',sep='_'),width=6.5, height=4)



















