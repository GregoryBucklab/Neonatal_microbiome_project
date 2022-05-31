library(pROC)
library(doSNOW)
library(parallel)
library(matrixStats)
library(erer) # export list


##### modeling #####
cl <- makeSOCKcluster(detectCores())
registerDoSNOW(cl)

importance_ML = list()
myPvs = vector()
myPvs2 = vector()
myPvs3 = vector()
myPvs4 = vector()

n_ML = 0

baby_weight_low = 88.1849 # Low birthweight is a term used to describe babies who are born weighing less than 2,500 grams (5 pounds, 8 ounces).

# microbiome
for (sample_type in c("MV1D","MRCD","MCKD")) {

  
  #prepare dataset
  { 
    keep = metadata_16s$Mombaby == "Baby"
    metadata_baby_all = metadata_16s[keep,]
    keep = !is.na(metadata_baby_all$weight )
    metadata_baby_all = metadata_baby_all[keep,]
    
    metadata_baby_all = metadata_baby_all[order(metadata_baby_all$VisitNum, decreasing = F),]
    metadata_baby_all = metadata_baby_all[order(metadata_baby_all$ParticipantID),]
    keep = !duplicated(metadata_baby_all$ParticipantID)
    metadata_baby_all = metadata_baby_all[keep,]
    
    keep = metadata_baby_all$weight <= baby_weight_low
    metadata_baby_all$weight[keep] = 'Yes'
    metadata_baby_all$weight[!keep] = 'No'
    
    keep = metadata_16s$Mombaby == "Mom"
    metadata_mom = metadata_16s[keep,]
    
    keep = metadata_mom$SampleType == sample_type
    metadata_mom =metadata_mom[keep,] 
    metadata_mom$baby_weight = NA
    
    keep = metadata_baby_all$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby_all = metadata_baby_all[keep,]
    keep = metadata_mom$ParticipantID %in% metadata_baby_all$ParticipantID
    metadata_mom = metadata_mom[keep,]
    
    metadata_mom$baby_weight = sapply(1:nrow(metadata_mom), 
                                      function(j) (metadata_mom$baby_weight[j] = 
                                                     metadata_baby_all$weight[which(metadata_baby_all$ParticipantID 
                                                                                    == metadata_mom$ParticipantID[j])] ))

    Baby_Weight = metadata_mom$baby_weight
#    metadata_mom$baby_weight = NULL
      
    data = data.frame(Baby_Weight=Baby_Weight, week_pregnancy=metadata_mom$weeks_pregnant)
    week_unique = unique(data$week_pregnancy)
    data_2 = as.data.frame(matrix(data = NA, nrow = length(week_unique)*2, ncol = 3))
    colnames(data_2) = c('week_pregnancy', 'Baby_Weight', 'Case_number')
    data_2$week_pregnancy = rep(week_unique,2)
    data_2$Baby_Weight = c(rep('Yes', length(week_unique)), rep('No', length(week_unique)))
    
    for (a in seq(length(week_unique))) {
      data_2$Case_number[a] = sum(data$Baby_Weight == 'Yes' & data$week_pregnancy == week_unique[a])
      data_2$Case_number[a+length(week_unique)] = sum(data$Baby_Weight == 'No' & data$week_pregnancy == week_unique[a])
    }
    data_2 = data_2[order(data_2$week_pregnancy),]
    ggplot(data_2, aes(week_pregnancy, Case_number, fill=Baby_Weight)) + 
      geom_bar(stat="identity", position="stack", width=1) 
    
    ggsave(paste0(sample_type,'_Baby_Weight_case_number.pdf'), width=8, height=3)
    
  }
  
  for (a in 2:7) {
    # prepare dataset
    {
      keep = metadata_mom$weeks_pregnant >= a*5 + 1 & metadata_mom$weeks_pregnant < (a+1)*5 + 1
      metadata = metadata_mom[keep,]
      
      keep = colnames(reads_table_16s) %in% metadata$SampleID
      sum(keep)
      reads_table_16s_mom_mv = reads_table_16s[,keep]
      
      keep = metadata$SampleID %in% colnames(reads_table_16s_mom_mv)
      metadata = metadata[keep,]
      reads_table_16s_mom_mv = reads_table_16s_mom_mv[metadata$SampleID]
      
      reads_table_16s_mom_mv = prepare_reads_table(reads_table_16s_mom_mv, metadata, total_reads_threshold = 0, species_threshold = 0.001, mc.cores = 8)
      metadata = reads_table_16s_mom_mv$metadata
      reads_table_16s_mom_mv = reads_table_16s_mom_mv$reads_table
      Baby_Weight_MV = metadata$baby_weight
    }
    
    # modeling
    {
      n_ML = n_ML+1
      output = Baby_Weight_MV
      
      reads_table = reads_table_16s_mom_mv
      reads_table = t(reads_table)
      reads_table = reads_table + 0.5
      reads_table <- clr(reads_table) 
      
      input = as.matrix((reads_table))
      
      set.seed(2020)
      input = FilterFeatures(input, output)
      if (is.vector(input)) {
        myPvs[n_ML]=NA
        myPvs2[n_ML]=NA
        myPvs3[n_ML] = NA
        importance_ML[[n_ML]] = NA
        myPvs4[n_ML] = NA
        
        next
      }
      
      if (ncol(input) == 0) {
        myPvs[n_ML]=NA
        myPvs2[n_ML]=NA
        myPvs3[n_ML] = NA
        importance_ML[[n_ML]] = NA
        myPvs4[n_ML] = NA
        next
      }
      data = xxx_get_result(input, output, ntree= ncol(input)*2) 
      
      myPvs[n_ML]=data$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
      myPvs2[n_ML]=data$pvalue   # test correlation between expected values and true values (PTB yes or no).
      myPvs3[n_ML] = data$err_median
      importance_ML[[n_ML]] = data$importance_ML
      myPvs4[n_ML] = paste0(length(output),' (',sum(output == 'Yes'),')')
      myPvs[n_ML]
    }
  }
  
}

#  metadata 1st visit
{
  setwd('/Users/binzhu/Desktop/Mom-baby/')
  keep = metadata_16s$Mombaby == "Baby"
  metadata_baby_all = metadata_16s[keep,]
  keep = !is.na(metadata_baby_all$weight )
  metadata_baby_all = metadata_baby_all[keep,]
  
  metadata_baby_all = metadata_baby_all[order(metadata_baby_all$VisitNum, decreasing = F),]
  metadata_baby_all = metadata_baby_all[order(metadata_baby_all$ParticipantID),]
  keep = !duplicated(metadata_baby_all$ParticipantID)
  metadata_baby_all = metadata_baby_all[keep,]
  
  keep = metadata_baby_all$weight <= baby_weight_low
  metadata_baby_all$weight[keep] = 'Yes'
  metadata_baby_all$weight[!keep] = 'No'
  
  keep = metadata_16s$Mombaby == "Mom"
  metadata_mom = metadata_16s[keep,]
  
  keep = metadata_mom$SampleType == sample_type
  metadata_mom =metadata_mom[keep,] 
  metadata_mom$baby_weight = NA
  
  keep = metadata_baby_all$ParticipantID %in% metadata_mom$ParticipantID
  metadata_baby_all = metadata_baby_all[keep,]
  keep = metadata_mom$ParticipantID %in% metadata_baby_all$ParticipantID
  metadata_mom = metadata_mom[keep,]
  
  metadata_mom$baby_weight = sapply(1:nrow(metadata_mom), 
                                    function(j) (metadata_mom$baby_weight[j] = 
                                                   metadata_baby_all$weight[which(metadata_baby_all$ParticipantID 
                                                                                  == metadata_mom$ParticipantID[j])] ))
  
  metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = F),]
  metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
  keep = !duplicated(metadata_mom$ParticipantID)
  metadata_mom = metadata_mom[keep,]      
  
  Baby_Weight = metadata_mom$baby_weight
  Baby_Weight = as.factor(Baby_Weight)
  metadata_mom$baby_weight = NULL
  
  metadata_mom = as.matrix(metadata_mom)
  metadata_mom[metadata_mom=='Yes'] = 1
  metadata_mom[metadata_mom=='No'] = 0
  metadata_mom = apply(metadata_mom, 2, as.numeric)
  metadata_mom = as.matrix(metadata_mom)
  
  keep = colSums(is.na(metadata_mom)) <= nrow(metadata_mom)/3
  metadata_mom = metadata_mom[,keep]

}
{
  n_ML = n_ML+1
  
  output = Baby_Weight
  output = as.factor(output)
  input = metadata_mom
  
  set.seed(2018)
  
  input = FilterFeatures(input, output)
  
  input = input[, !str_detect(colnames(input), 'average|baby|ga_at_delivery|event
                                     |pregnancy|last_weight|last_pulse
                                     |last_ga_at_delivery|mom_feeding_method|last_weight
                                     |hispanic_or_latino|vaginal_penetration_frequency
                                     |hepatitis_c|weight_change|3rdtri
                                     |poor_appetite|vaginal_penetration_frequency|2ndtri') ]
  
  data = xxx_get_result(input, output, ntree= ncol(input)*2) 
  
  data$p_roc
  ggsave('Baby_Weight_M_meta_Sensitivities_Specificities.pdf',width=3, height=3)
  
  myPvs[n_ML]=data$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
  myPvs2[n_ML]=data$pvalue   # test correlation between expected values and true values (PTB yes or no).
  myPvs3[n_ML] = data$err_median
  importance_ML[[n_ML]] = data$importance_ML
  myPvs4[n_ML] = paste0(length(output),' (',sum(output == 'Yes'),')')
  myPvs[n_ML]
  
}

# early metadata modeling, weeks 6 - 11 
{ 
  keep = metadata_16s$Mombaby == "Baby"
  metadata_baby_all = metadata_16s[keep,]
  keep = !is.na(metadata_baby_all$weight )
  metadata_baby_all = metadata_baby_all[keep,]
  
  metadata_baby_all = metadata_baby_all[order(metadata_baby_all$VisitNum, decreasing = F),]
  metadata_baby_all = metadata_baby_all[order(metadata_baby_all$ParticipantID),]
  keep = !duplicated(metadata_baby_all$ParticipantID)
  metadata_baby_all = metadata_baby_all[keep,]
  
  keep = metadata_baby_all$weight <= baby_weight_low
  metadata_baby_all$weight[keep] = 'Yes'
  metadata_baby_all$weight[!keep] = 'No'
  
  keep = metadata_16s$Mombaby == "Mom"
  metadata_mom = metadata_16s[keep,]
  
  keep = metadata_mom$SampleType == sample_type
  metadata_mom =metadata_mom[keep,] 
  metadata_mom$baby_weight = NA
  
  keep = metadata_baby_all$ParticipantID %in% metadata_mom$ParticipantID
  metadata_baby_all = metadata_baby_all[keep,]
  keep = metadata_mom$ParticipantID %in% metadata_baby_all$ParticipantID
  metadata_mom = metadata_mom[keep,]
  
  metadata_mom$baby_weight = sapply(1:nrow(metadata_mom), 
                                    function(j) (metadata_mom$baby_weight[j] = 
                                                   metadata_baby_all$weight[which(metadata_baby_all$ParticipantID 
                                                                                  == metadata_mom$ParticipantID[j])] ))
  
  a = 2
  keep = metadata_mom$weeks_pregnant >= a*5 + 1 & metadata_mom$weeks_pregnant < (a+1)*5 + 1
  metadata_mom = metadata_mom[keep,]
  
  metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = F),]
  metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
  keep = !duplicated(metadata_mom$ParticipantID)
  metadata_mom = metadata_mom[keep,]      
  
  Baby_Weight = metadata_mom$baby_weight
  Baby_Weight = as.factor(Baby_Weight)
  metadata_mom$baby_weight = NULL
  
  metadata_mom = as.matrix(metadata_mom)
  metadata_mom[metadata_mom=='Yes'] = 1
  metadata_mom[metadata_mom=='No'] = 0
  metadata_mom = apply(metadata_mom, 2, as.numeric)
  metadata_mom = as.matrix(metadata_mom)
  
  keep = colSums(is.na(metadata_mom)) <= nrow(metadata_mom)/3
  metadata_mom = metadata_mom[,keep]
  
  
}
{
  n_ML = n_ML+1
  
  output = Baby_Weight
  output = as.factor(output)
  input = metadata_mom
  
  input = input[, !str_detect(colnames(input), 'average|baby|ga_at_delivery|event
                                     |pregnancy|last_weight|last_pulse
                                     |last_ga_at_delivery|mom_feeding_method|last_weight
                                     |hispanic_or_latino|vaginal_penetration_frequency
                                     |hepatitis_c|weight_change|3rdtri
                                     |poor_appetite|vaginal_penetration_frequency|2ndtri') ]
  set.seed(2018)
  
  data = xxx_get_result(input, output, ntree= ncol(input)*2) 
  
  data$p_roc
  ggsave('Baby_Weight_M_meta_week6_11_Sensitivities_Specificities.pdf',width=3, height=3)
  
  myPvs[n_ML]=data$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
  myPvs2[n_ML]=data$pvalue   # test correlation between expected values and true values (PTB yes or no).
  myPvs3[n_ML] = data$err_median
  importance_ML[[n_ML]] = data$importance_ML
  myPvs4[n_ML] = paste0(length(output),' (',sum(output == 'Yes'),')')
  myPvs[n_ML]
  
}

# <= 21 weeks 
{ 
  keep = metadata_16s$Mombaby == "Baby"
  metadata_baby_all = metadata_16s[keep,]
  keep = !is.na(metadata_baby_all$weight )
  metadata_baby_all = metadata_baby_all[keep,]
  
  metadata_baby_all = metadata_baby_all[order(metadata_baby_all$VisitNum, decreasing = F),]
  metadata_baby_all = metadata_baby_all[order(metadata_baby_all$ParticipantID),]
  keep = !duplicated(metadata_baby_all$ParticipantID)
  metadata_baby_all = metadata_baby_all[keep,]
  
  keep = metadata_baby_all$weight <= baby_weight_low
  metadata_baby_all$weight[keep] = 'Yes'
  metadata_baby_all$weight[!keep] = 'No'
  
  keep = metadata_16s$Mombaby == "Mom"
  metadata_mom = metadata_16s[keep,]
  
  keep = metadata_mom$SampleType == sample_type
  metadata_mom =metadata_mom[keep,] 
  metadata_mom$baby_weight = NA
  
  keep = metadata_baby_all$ParticipantID %in% metadata_mom$ParticipantID
  metadata_baby_all = metadata_baby_all[keep,]
  keep = metadata_mom$ParticipantID %in% metadata_baby_all$ParticipantID
  metadata_mom = metadata_mom[keep,]
  
  metadata_mom$baby_weight = sapply(1:nrow(metadata_mom), 
                                    function(j) (metadata_mom$baby_weight[j] = 
                                                   metadata_baby_all$weight[which(metadata_baby_all$ParticipantID 
                                                                                  == metadata_mom$ParticipantID[j])] ))
  
  keep = metadata_mom$weeks_pregnant <= 21
  metadata_mom = metadata_mom[keep,]
  
  metadata_mom = metadata_mom[order(metadata_mom$VisitNum, decreasing = F),]
  metadata_mom = metadata_mom[order(metadata_mom$ParticipantID),]
  keep = !duplicated(metadata_mom$ParticipantID)
  metadata_mom = metadata_mom[keep,]      
  
  Baby_Weight = metadata_mom$baby_weight
  Baby_Weight = as.factor(Baby_Weight)
  metadata_mom$baby_weight = NULL
  
  metadata_mom = as.matrix(metadata_mom)
  metadata_mom[metadata_mom=='Yes'] = 1
  metadata_mom[metadata_mom=='No'] = 0
  metadata_mom = apply(metadata_mom, 2, as.numeric)
  metadata_mom = as.matrix(metadata_mom)
  
  keep = colSums(is.na(metadata_mom)) <= nrow(metadata_mom)/3
  metadata_mom = metadata_mom[,keep]
  
  
}
{
  n_ML = n_ML+1
  
  output = Baby_Weight
  output = as.factor(output)
  input = metadata_mom
  input = input[, !str_detect(colnames(input), 'average|baby|ga_at_delivery|event
                                     |pregnancy|last_weight|last_pulse
                                     |last_ga_at_delivery|mom_feeding_method|last_weight
                                     |hispanic_or_latino|vaginal_penetration_frequency
                                     |hepatitis_c|weight_change|3rdtri
                                     |poor_appetite|vaginal_penetration_frequency|2ndtri') ]
  set.seed(2018)
  
  data = xxx_get_result(input, output, ntree= ncol(input)*2) 
  
  data$p_roc
  ggsave('Baby_Weight_M_meta_half_Sensitivities_Specificities.pdf',width=3, height=3)
  
  myPvs[n_ML]=data$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
  myPvs2[n_ML]=data$pvalue   # test correlation between expected values and true values (PTB yes or no).
  myPvs3[n_ML] = data$err_median
  importance_ML[[n_ML]] = data$importance_ML
  myPvs4[n_ML] = paste0(length(output),' (',sum(output == 'Yes'),')')
  myPvs[n_ML]
  
}
for (sample_type in c("MV1D","MRCD","MCKD")) {
  
  
  #prepare dataset
  { 
    keep = metadata_16s$Mombaby == "Baby"
    metadata_baby_all = metadata_16s[keep,]
    keep = !is.na(metadata_baby_all$weight )
    metadata_baby_all = metadata_baby_all[keep,]
    
    metadata_baby_all = metadata_baby_all[order(metadata_baby_all$VisitNum, decreasing = F),]
    metadata_baby_all = metadata_baby_all[order(metadata_baby_all$ParticipantID),]
    keep = !duplicated(metadata_baby_all$ParticipantID)
    metadata_baby_all = metadata_baby_all[keep,]
    
    keep = metadata_baby_all$weight <= baby_weight_low
    metadata_baby_all$weight[keep] = 'Yes'
    metadata_baby_all$weight[!keep] = 'No'
    
    keep = metadata_16s$Mombaby == "Mom"
    metadata_mom = metadata_16s[keep,]
    
    keep = metadata_mom$SampleType == sample_type
    metadata_mom =metadata_mom[keep,] 
    metadata_mom$baby_weight = NA
    
    keep = metadata_baby_all$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby_all = metadata_baby_all[keep,]
    keep = metadata_mom$ParticipantID %in% metadata_baby_all$ParticipantID
    metadata_mom = metadata_mom[keep,]
    
    metadata_mom$baby_weight = sapply(1:nrow(metadata_mom), 
                                      function(j) (metadata_mom$baby_weight[j] = 
                                                     metadata_baby_all$weight[which(metadata_baby_all$ParticipantID 
                                                                                    == metadata_mom$ParticipantID[j])] ))
    
    Baby_Weight = metadata_mom$baby_weight
    #    metadata_mom$baby_weight = NULL
    
    data = data.frame(Baby_Weight=Baby_Weight, week_pregnancy=metadata_mom$weeks_pregnant)
    week_unique = unique(data$week_pregnancy)
    data_2 = as.data.frame(matrix(data = NA, nrow = length(week_unique)*2, ncol = 3))
    colnames(data_2) = c('week_pregnancy', 'Baby_Weight', 'Case_number')
    data_2$week_pregnancy = rep(week_unique,2)
    data_2$Baby_Weight = c(rep('Yes', length(week_unique)), rep('No', length(week_unique)))
    
    for (a in seq(length(week_unique))) {
      data_2$Case_number[a] = sum(data$Baby_Weight == 'Yes' & data$week_pregnancy == week_unique[a])
      data_2$Case_number[a+length(week_unique)] = sum(data$Baby_Weight == 'No' & data$week_pregnancy == week_unique[a])
    }
    data_2 = data_2[order(data_2$week_pregnancy),]
    ggplot(data_2, aes(week_pregnancy, Case_number, fill=Baby_Weight)) + 
      geom_bar(stat="identity", position="stack", width=1) 
    
    ggsave(paste0(sample_type,'_Baby_Weight_case_number.pdf'), width=8, height=3)
    
  }
  
  {
    # prepare dataset
    {
      keep = metadata_mom$weeks_pregnant <= 21
      metadata = metadata_mom[keep,]
      
      keep = colnames(reads_table_16s) %in% metadata$SampleID
      sum(keep)
      reads_table_16s_mom_mv = reads_table_16s[,keep]
      
      keep = metadata$SampleID %in% colnames(reads_table_16s_mom_mv)
      metadata = metadata[keep,]
      reads_table_16s_mom_mv = reads_table_16s_mom_mv[metadata$SampleID]
      
      reads_table_16s_mom_mv = prepare_reads_table(reads_table_16s_mom_mv, metadata, total_reads_threshold = 0, species_threshold = 0.001, mc.cores = 8)
      metadata = reads_table_16s_mom_mv$metadata
      reads_table_16s_mom_mv = reads_table_16s_mom_mv$reads_table
      Baby_Weight_MV = metadata$baby_weight
    }
    
    # modeling
    {
      n_ML = n_ML+1
      output = Baby_Weight_MV
      
      reads_table = reads_table_16s_mom_mv
      reads_table = t(reads_table)
      reads_table = reads_table + 0.5
      reads_table <- clr(reads_table) 
      
      input = as.matrix((reads_table))
      
      set.seed(2020)
      input = FilterFeatures(input, output)
      if (is.vector(input)) {
        myPvs[n_ML]=NA
        myPvs2[n_ML]=NA
        myPvs3[n_ML] = NA
        importance_ML[[n_ML]] = NA
        myPvs4[n_ML] = NA
        
        next
      }
      
      if (ncol(input) == 0) {
        myPvs[n_ML]=NA
        myPvs2[n_ML]=NA
        myPvs3[n_ML] = NA
        importance_ML[[n_ML]] = NA
        myPvs4[n_ML] = NA
        next
      }
      data = xxx_get_result(input, output, ntree= ncol(input)*2) 
      
      myPvs[n_ML]=data$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
      myPvs2[n_ML]=data$pvalue   # test correlation between expected values and true values (PTB yes or no).
      myPvs3[n_ML] = data$err_median
      importance_ML[[n_ML]] = data$importance_ML
      myPvs4[n_ML] = paste0(length(output),' (',sum(output == 'Yes'),')')
      myPvs[n_ML]
    }
  }
  
}

# combined model <= 21 weeks + MV+MB+MR
{
  # prepare data
  {
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    keep = metadata_16s$Mombaby == "Baby"
    metadata_baby_all = metadata_16s[keep,]
    keep = !is.na(metadata_baby_all$weight )
    metadata_baby_all = metadata_baby_all[keep,]
    
    metadata_baby_all = metadata_baby_all[order(metadata_baby_all$VisitNum, decreasing = F),]
    metadata_baby_all = metadata_baby_all[order(metadata_baby_all$ParticipantID),]
    keep = !duplicated(metadata_baby_all$ParticipantID)
    metadata_baby_all = metadata_baby_all[keep,]
    
    keep = metadata_baby_all$weight <= baby_weight_low
    metadata_baby_all$weight[keep] = 'Yes'
    metadata_baby_all$weight[!keep] = 'No'
    
    keep = metadata_16s$Mombaby == "Mom"
    metadata_mom = metadata_16s[keep,]
    
    metadata_mom$baby_weight = NA
    
    keep = metadata_baby_all$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby_all = metadata_baby_all[keep,]
    keep = metadata_mom$ParticipantID %in% metadata_baby_all$ParticipantID
    metadata_mom = metadata_mom[keep,]
    
    metadata_mom$baby_weight = sapply(1:nrow(metadata_mom), 
                                      function(j) (metadata_mom$baby_weight[j] = 
                                                     metadata_baby_all$weight[which(metadata_baby_all$ParticipantID 
                                                                                    == metadata_mom$ParticipantID[j])] ))
    
    keep = metadata_mom$weeks_pregnant <=21
    metadata_mom_all = metadata_mom[keep,]
    
    
    # get paired samples
    keep = metadata_mom_all$SampleType == "MV1D"
    metadata_Baby_Weight_MV = metadata_mom_all[keep,]
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[order(metadata_Baby_Weight_MV$VisitNum, decreasing = F),]
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[order(metadata_Baby_Weight_MV$ParticipantID),]
    keep = !duplicated(metadata_Baby_Weight_MV$ParticipantID)
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[keep,]
    
    keep = metadata_mom_all$SampleType == "MRCD"
    metadata_Baby_Weight_MR = metadata_mom_all[keep,]
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[order(metadata_Baby_Weight_MR$VisitNum, decreasing = F),]
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[order(metadata_Baby_Weight_MR$ParticipantID),]
    keep = !duplicated(metadata_Baby_Weight_MR$ParticipantID)
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[keep,] 
    
    keep = metadata_mom_all$SampleType == "MCKD"
    metadata_Baby_Weight_MB = metadata_mom_all[keep,]
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[order(metadata_Baby_Weight_MB$VisitNum, decreasing = F),]
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[order(metadata_Baby_Weight_MB$ParticipantID),]
    keep = !duplicated(metadata_Baby_Weight_MB$ParticipantID)
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[keep,] 
    
    # metadata 
    keep = metadata_mom_all$ParticipantID
    keep_list =   keep[keep %in% metadata_Baby_Weight_MV$ParticipantID & 
                         keep %in% metadata_Baby_Weight_MB$ParticipantID& 
                         keep %in% metadata_Baby_Weight_MR$ParticipantID& 
                         keep %in% metadata_Baby_Weight_MV$ParticipantID]
    
    keep = metadata_Baby_Weight_MV$ParticipantID %in% keep_list
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[keep,]
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[order(metadata_Baby_Weight_MV$ParticipantID),]
    
    keep = metadata_Baby_Weight_MR$ParticipantID %in% keep_list
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[keep,]
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[order(metadata_Baby_Weight_MR$ParticipantID),]
    
    keep = metadata_Baby_Weight_MB$ParticipantID %in% keep_list
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[keep,]
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[order(metadata_Baby_Weight_MB$ParticipantID),]
    
    keep = metadata_mom_all$ParticipantID %in% keep_list
    metadata_mom_all = metadata_mom_all[keep,]
    
    Participant_list = metadata_mom_all$ParticipantID
    metadata_mom_all = metadata_mom_all[order(metadata_mom_all$VisitNum, decreasing = F),]
    keep = !duplicated(metadata_mom_all$ParticipantID)
    metadata_mom_all = metadata_mom_all[keep,]
    
    metadata_mom_all = metadata_mom_all[order(metadata_mom_all$ParticipantID),]
    
    # get paired microbiomes
    get_data_for_Baby_Weight_1st <- function(metadata, reads_table) {
      
      keep = colnames(reads_table) %in% metadata$SampleID
      sum(keep)
      reads_table_output = reads_table[,keep]
      reads_table_output = reads_table_output[metadata$SampleID]
      
      reads_table_output = prepare_reads_table_2(reads_table_output, total_reads_threshold = 0, species_threshold = 0.001, mc.cores = 8)
      keep = metadata$SampleID %in% colnames(reads_table_output)
      metadata = metadata[keep,]
      Baby_Weight_output = metadata$baby_weight
      
      reads_table_output = t(reads_table_output)
      reads_table_output = reads_table_output + 0.5
      reads_table_output <- clr(reads_table_output) 
      reads_table_output = as.matrix(reads_table_output)
      
      return(c(reads_table_output = list(reads_table_output), 
               Baby_Weight_output = list(Baby_Weight_output), 
               metadata = list(metadata)))
    }
    
    data = get_data_for_Baby_Weight_1st(metadata_Baby_Weight_MV, reads_table_16s)
    reads_table_16s_MV = data$reads_table_output
    metadata_Baby_Weight_MV = data$metadata
    
    data = get_data_for_Baby_Weight_1st(metadata_Baby_Weight_MB, reads_table_16s)
    reads_table_16s_MB = data$reads_table_output
    metadata_Baby_Weight_MB = data$metadata
    
    data = get_data_for_Baby_Weight_1st(metadata_Baby_Weight_MR, reads_table_16s)
    reads_table_16s_MR = data$reads_table_output
    metadata_Baby_Weight_MR = data$metadata
    
    Baby_Weight_2 = metadata_mom_all$baby_weight
    
    metadata_mom_all[metadata_mom_all=='Yes'] = 1
    metadata_mom_all[metadata_mom_all=='No'] = 0
    
    names = metadata_mom_all$SampleID
    metadata = apply(metadata_mom_all, 2, as.numeric)
    
    keep = colSums(is.na(metadata)) < nrow(metadata)/3
    metadata = metadata[,keep]
    row.names(metadata) = names
    
    metadata = FilterFeatures(metadata, Baby_Weight_2)
    metadata = metadata[,!str_detect(colnames(metadata),'average|baby|ga_at_delivery|event
                                     |pregnancy|last_weight|last_pulse
                                     |last_ga_at_delivery|weeks_pregnant
                                     |mom_feeding_method|last_weight
                                     |hispanic_or_latino|vaginal_penetration_frequency
                                     |hepatitis_c|weight_change|contractions_2ndtri
                                     |poor_appetite|vaginal_penetration_frequency|2ndtri')]
  }
  # two layer modeling
  {
    output = as.factor(Baby_Weight_2)
    input = reads_table_16s_MV
    set.seed(2022)
    input = FilterFeatures(input, output)
    prd_MV_1st=foreach(i=seq(nrow(input))) %dopar% twoLayer(input, output, i, ntree= ncol(input)*2)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
    ppp_MV_1st=vector()
    for(i in seq(nrow(input)))
      ppp_MV_1st[i] = prd_MV_1st[[i]]$Layer1$p2   # all 81 expected values in the first layer model

    output = as.factor(Baby_Weight_2)
    input = reads_table_16s_MR
    set.seed(2022)
    input = FilterFeatures(input, output)
    prd_MR_1st=foreach(i=seq(nrow(input))) %dopar% twoLayer(input, output, i, ntree= ncol(input)*2)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
    ppp_MR_1st=vector()
    for(i in seq(nrow(input)))
      ppp_MR_1st[i] = prd_MR_1st[[i]]$Layer1$p2   # all 81 expected values in the first layer model
    
    output = as.factor(Baby_Weight_2)
    input = reads_table_16s_MB
    set.seed(2022)
    input = FilterFeatures(input, output)
    prd_MB_1st=foreach(i=seq(nrow(input))) %dopar% twoLayer(input, output, i, ntree= ncol(input)*2)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
    ppp_MB_1st=vector()
    for(i in seq(nrow(input)))
      ppp_MB_1st[i] = prd_MB_1st[[i]]$Layer1$p2   # all 81 expected values in the first layer model
    
    
    output = as.factor(Baby_Weight_2)
    input = metadata
    set.seed(2022)
    input = FilterFeatures(input, output)
    prd_metadata_1st=foreach(i=seq(nrow(input))) %dopar% twoLayer(input, output, i, ntree= ncol(input)*2)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
    ppp_metadata_1st=vector()
    for(i in seq(nrow(input)))
      ppp_metadata_1st[i] = prd_metadata_1st[[i]]$Layer1$p2   # all 81 expected values in the first layer model
    
  }
  
  # calculate weight
  {
    n_ML = n_ML+1
    
    CMBW=vector()
    prf_all = list()
    prf_all_sum = vector()
    for(i in seq(nrow(metadata)) )
    {
      prf=vector()
      
      ppp_M_meta_combine=vector()   # second layer prediction expected values
      ppp_MV_1st_combine=vector()   # second layer prediction expected values
      ppp_MB_1st_combine=vector()   # second layer prediction expected values
      ppp_MR_1st_combine=vector()   # second layer prediction expected values
      
      for(j in seq((nrow(metadata)-1)))
      {
        ppp_M_meta_combine[j]=prd_metadata_1st[[i]]$Layer2[[j]]$p2
        ppp_MV_1st_combine[j]=prd_MV_1st[[i]]$Layer2[[j]]$p2
        ppp_MB_1st_combine[j]=prd_MB_1st[[i]]$Layer2[[j]]$p2
        ppp_MR_1st_combine[j]=prd_MR_1st[[i]]$Layer2[[j]]$p2
      }
      
      prf[1]=wilcox.test(ppp_M_meta_combine , as.numeric(as.factor(Baby_Weight_2[-i])))$p.value
      prf[2]=wilcox.test(ppp_MV_1st_combine , as.numeric(as.factor(Baby_Weight_2[-i])))$p.value
      prf[3]=wilcox.test(ppp_MB_1st_combine , as.numeric(as.factor(Baby_Weight_2[-i])))$p.value
      prf[4]=wilcox.test(ppp_MR_1st_combine , as.numeric(as.factor(Baby_Weight_2[-i])))$p.value
      
      prf=exp(-log10(prf))
      prf_all[[i]] = prf
      prf_all_sum[i] = sum(prf)
#      CMBW[i]=(prf[1]*ppp_metadata_1st[i]+prf[2]*ppp_MV_1st[i])+prf[3]*ppp_MR_1st[i]/sum(prf)   # first layer expected value * -log(second layer pvalue) / sum(-log(second layer pvalue))
      CMBW[i]=(prf[1]*ppp_metadata_1st[i]+prf[2]*ppp_MV_1st[i]+prf[3]*ppp_MB_1st[i]+prf[4]*ppp_MR_1st[i])/sum(prf)   # first layer expected value * -log(second layer pvalue) / sum(-log(second layer pvalue))
    }
    
#    CMBW = CMBW/1000000000
    myPvs[n_ML]=roc(Baby_Weight_2, CMBW)$auc
    myPvs2[n_ML]=-log10(wilcox.test(CMBW ~ as.factor(Baby_Weight_2))$p.value)
    myPvs4[n_ML] = paste0(length(Baby_Weight_2),' (',sum(Baby_Weight_2 == 'Yes'),')')
    myPvs[n_ML]
  }
  
}





##### visualize results #####
{
  roc = roc(Baby_Weight_2, CMBW)
  data = cbind(roc$sensitivities,roc$specificities)
  data = as.data.frame(data)
  colnames(data) = c('Sensitivities','Specificities')
  data = data[order(data$Specificities, decreasing = F),]
  
  ggplot(data, aes(x = Sensitivities, y = Specificities)) +
    geom_point(size = 1)+geom_step()+scale_x_reverse()+
    theme(axis.title = element_text(size = 6), 
          axis.text = element_text(size = 6), 
          legend.text = element_text(size = 6), 
          legend.title = element_text(size = 6)) + theme_bw() 
  ggsave('Sensitivities_Specificities_models_integrated.pdf',width=3, height=3)
  
  
  
  data = data.frame(Real_value=Baby_Weight_2, Predicted_value= CMBW)
  
  ggplot(data, aes(x=Real_value, y=Predicted_value,fill=Real_value)) +
    geom_boxplot() + geom_point()
  ggsave('model_integrated.pdf',width=2.5, height=3)
  
  
  
  model_name =c('MV 11-15 weeks','MV 16-20 weeks','MV 21-25 weeks',
                'MV 26-30 weeks','MV 31-35 weeks','MV 36-40 weeks',
                'MR 11-15 weeks','MR 16-20 weeks','MR 21-25 weeks',
                'MR 26-30 weeks','MR 31-35 weeks','MR 36-40 weeks',
                'MB 11-15 weeks','MB 16-20 weeks','MB 21-25 weeks',
                'MB 26-30 weeks','MB 31-35 weeks','MB 36-40 weeks',
                'Metadata 1st visit','Metadata 11-15 weeks','Metadata <= 21 weeks',
                'MV <= 21 weeks','MR <= 21 weeks','MB <= 21 weeks','Integrated model')
  
  data = cbind(model_name,myPvs,myPvs2, myPvs4)
  data = as.data.frame(data)
  colnames(data) = c('Model_name','AUC_ROC_curve','p_value','Case_number')
  
  data$Model_name <- factor(data$Model_name, levels=c('MV 11-15 weeks','MV 16-20 weeks','MV 21-25 weeks',
                                                      'MV 26-30 weeks','MV 31-35 weeks','MV 36-40 weeks',
                                                      'MR 11-15 weeks','MR 16-20 weeks','MR 21-25 weeks',
                                                      'MR 26-30 weeks','MR 31-35 weeks','MR 36-40 weeks',
                                                      'MB 11-15 weeks','MB 16-20 weeks','MB 21-25 weeks',
                                                      'MB 26-30 weeks','MB 31-35 weeks','MB 36-40 weeks',
                                                      'Metadata 1st visit','Metadata 11-15 weeks','Metadata <= 21 weeks',
                                                      'MV <= 21 weeks','MR <= 21 weeks','MB <= 21 weeks','Integrated model'))
  
  data$AUC_ROC_curve = as.numeric(data$AUC_ROC_curve)
  data$p_value = as.numeric(data$p_value)
  
  ggplot(data, aes(Model_name, AUC_ROC_curve)) + 
    geom_bar(stat="identity", width = 0.7,color = 'black', fill = 'gray') + scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    geom_text(aes(label=Case_number), position=position_dodge(width=0.9), vjust=0.5,hjust=2,size = 2,angle=90)+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 7), 
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6), 
          legend.title = element_text(size = 7))
  ggsave('AUC_ROC_curve_integrated.pdf',width=4, height=3)
  
  
  ggplot(data, aes(Model_name, p_value)) + 
    geom_bar(stat="identity", width = 0.7,color = 'black', fill = 'gray') + scale_y_continuous(breaks = seq(0, 10, 1)) +
    geom_text(aes(label=Case_number), position=position_dodge(width=0.9), vjust=0.5,hjust=200,size = 2,angle=90)+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title = element_text(size = 7), 
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6), 
          legend.title = element_text(size = 7))+ 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red", size=1)
  ggsave('p_value_models_integrated.pdf',width=4, height=3)
  
  names(importance_ML) = c('MV 11-15 weeks','MV 16-20 weeks','MV 21-25 weeks',
                           'MV 26-30 weeks','MV 31-35 weeks','MV 36-40 weeks',
                           'MR 11-15 weeks','MR 16-20 weeks','MR 21-25 weeks',
                           'MR 26-30 weeks','MR 31-35 weeks','MR 36-40 weeks',
                           'MB 11-15 weeks','MB 16-20 weeks','MB 21-25 weeks',
                           'MB 26-30 weeks','MB 31-35 weeks','MB 36-40 weeks',
                           'Metadata 1st visit','Metadata 11-15 weeks','Metadata <= 21 weeks',
                           'MV <= 21 weeks','MR <= 21 weeks','MB <= 21 weeks')
  library(erer)
  write.list(importance_ML, 'importance_ML_metadata.csv')
}




### final model ###
{
  library(randomForest)
  n = which(prf_all_sum == max(prf_all_sum))[1]
  
  set.seed(2022)
  
  imputation_input <- function(input, output) {
    if (sum(is.na(input)) >0) {
      input = apply(input, 2, as.numeric)
      output = as.factor(output)
      data = cbind(input,output)
      data_2 <- rfImpute(output ~ ., data)    # Impute missing values in predictor data using proximity from randomForest.
      input = data_2[,-1]
      input = as.matrix(input)
    } 
    return(input)
  }
  
  output = as.factor(Baby_Weight_2)
  
  input_1 = metadata
  input_1 = imputation_input(input_1, output)
  keep = colnames(input_1) %in% row.names(prd_metadata_1st[[n]]$Layer1$RF[["importance"]])
  input_1 = input_1[,keep]
  
  input_2 = reads_table_16s_MV
  input_2 = imputation_input(input_2, output)
  keep = colnames(input_2) %in% row.names(prd_MV_1st[[n]]$Layer1$RF[["importance"]])
  input_2 = input_2[,keep]
  
  input_3 = reads_table_16s_MR
  input_3 = imputation_input(input_3, output)
  keep = colnames(input_3) %in% row.names(prd_MR_1st[[n]]$Layer1$RF[["importance"]])
  input_3 = input_3[,keep]
  
  prediction_value = vector()
  for (i in seq(length(output))) {
    prediction_value[i]= (prf_all[[n]][1] *predict(prd_metadata_1st[[n]]$Layer1$RF, rbind(input_1[i,], input_1[i,]), type='prob')[1] +
                            prf_all[[n]][2] *predict(prd_MV_1st[[n]]$Layer1$RF, rbind(input_2[i,], input_2[i,]), type='prob')[1] +
                            prf_all[[n]][3] *predict(prd_MR_1st[[n]]$Layer1$RF, rbind(input_3[i,], input_3[i,]), type='prob')[1])/prf_all_sum[n]
    
  }
  
}


##### combined model 1st visit not good #####
{
  # prepare data
  {
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    keep = metadata_16s$Mombaby == "Baby"
    metadata_baby_all = metadata_16s[keep,]
    keep = !is.na(metadata_baby_all$weight )
    metadata_baby_all = metadata_baby_all[keep,]
    
    metadata_baby_all = metadata_baby_all[order(metadata_baby_all$VisitNum, decreasing = F),]
    metadata_baby_all = metadata_baby_all[order(metadata_baby_all$ParticipantID),]
    keep = !duplicated(metadata_baby_all$ParticipantID)
    metadata_baby_all = metadata_baby_all[keep,]
    
    keep = metadata_baby_all$weight <= baby_weight_low
    metadata_baby_all$weight[keep] = 'Yes'
    metadata_baby_all$weight[!keep] = 'No'
    
    keep = metadata_16s$Mombaby == "Mom"
    metadata_mom = metadata_16s[keep,]
    
    metadata_mom$baby_weight = NA
    
    keep = metadata_baby_all$ParticipantID %in% metadata_mom$ParticipantID
    metadata_baby_all = metadata_baby_all[keep,]
    keep = metadata_mom$ParticipantID %in% metadata_baby_all$ParticipantID
    metadata_mom = metadata_mom[keep,]
    
    metadata_mom$baby_weight = sapply(1:nrow(metadata_mom), 
                                      function(j) (metadata_mom$baby_weight[j] = 
                                                     metadata_baby_all$weight[which(metadata_baby_all$ParticipantID 
                                                                                    == metadata_mom$ParticipantID[j])] ))
    
    
    metadata_mom_all = metadata_mom
    
    
    # get paired samples
    keep = metadata_mom_all$SampleType == "MV1D"
    metadata_Baby_Weight_MV = metadata_mom_all[keep,]
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[order(metadata_Baby_Weight_MV$VisitNum, decreasing = F),]
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[order(metadata_Baby_Weight_MV$ParticipantID),]
    keep = !duplicated(metadata_Baby_Weight_MV$ParticipantID)
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[keep,]
    
    keep = metadata_mom_all$SampleType == "MRCD"
    metadata_Baby_Weight_MR = metadata_mom_all[keep,]
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[order(metadata_Baby_Weight_MR$VisitNum, decreasing = F),]
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[order(metadata_Baby_Weight_MR$ParticipantID),]
    keep = !duplicated(metadata_Baby_Weight_MR$ParticipantID)
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[keep,] 
    
    keep = metadata_mom_all$SampleType == "MCKD"
    metadata_Baby_Weight_MB = metadata_mom_all[keep,]
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[order(metadata_Baby_Weight_MB$VisitNum, decreasing = F),]
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[order(metadata_Baby_Weight_MB$ParticipantID),]
    keep = !duplicated(metadata_Baby_Weight_MB$ParticipantID)
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[keep,] 
    
    # metadata 
    keep = metadata_mom_all$ParticipantID
    keep_list =   keep[keep %in% metadata_Baby_Weight_MV$ParticipantID & 
                         keep %in% metadata_Baby_Weight_MB$ParticipantID& 
                         keep %in% metadata_Baby_Weight_MR$ParticipantID& 
                         keep %in% metadata_Baby_Weight_MV$ParticipantID]
    
    keep = metadata_Baby_Weight_MV$ParticipantID %in% keep_list
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[keep,]
    metadata_Baby_Weight_MV = metadata_Baby_Weight_MV[order(metadata_Baby_Weight_MV$ParticipantID),]
    
    keep = metadata_Baby_Weight_MR$ParticipantID %in% keep_list
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[keep,]
    metadata_Baby_Weight_MR = metadata_Baby_Weight_MR[order(metadata_Baby_Weight_MR$ParticipantID),]
    
    keep = metadata_Baby_Weight_MB$ParticipantID %in% keep_list
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[keep,]
    metadata_Baby_Weight_MB = metadata_Baby_Weight_MB[order(metadata_Baby_Weight_MB$ParticipantID),]
    
    keep = metadata_mom_all$ParticipantID %in% keep_list
    metadata_mom_all = metadata_mom_all[keep,]
    
    Participant_list = metadata_mom_all$ParticipantID
    metadata_mom_all = metadata_mom_all[order(metadata_mom_all$VisitNum, decreasing = F),]
    keep = !duplicated(metadata_mom_all$ParticipantID)
    metadata_mom_all = metadata_mom_all[keep,]
    
    metadata_mom_all = metadata_mom_all[order(metadata_mom_all$ParticipantID),]
    
    # get paired microbiomes
    get_data_for_Baby_Weight_1st <- function(metadata, reads_table) {
      
      keep = colnames(reads_table) %in% metadata$SampleID
      sum(keep)
      reads_table_output = reads_table[,keep]
      reads_table_output = reads_table_output[metadata$SampleID]
      
      reads_table_output = prepare_reads_table_2(reads_table_output, total_reads_threshold = 0, species_threshold = 0.001, mc.cores = 8)
      keep = metadata$SampleID %in% colnames(reads_table_output)
      metadata = metadata[keep,]
      Baby_Weight_output = metadata$baby_weight
      
      reads_table_output = t(reads_table_output)
      reads_table_output = reads_table_output + 0.5
      reads_table_output <- clr(reads_table_output) 
      reads_table_output = as.matrix(reads_table_output)
      
      return(c(reads_table_output = list(reads_table_output), 
               Baby_Weight_output = list(Baby_Weight_output), 
               metadata = list(metadata)))
    }
    
    data = get_data_for_Baby_Weight_1st(metadata_Baby_Weight_MV, reads_table_16s)
    reads_table_16s_MV = data$reads_table_output
    metadata_Baby_Weight_MV = data$metadata
    
    data = get_data_for_Baby_Weight_1st(metadata_Baby_Weight_MB, reads_table_16s)
    reads_table_16s_MB = data$reads_table_output
    metadata_Baby_Weight_MB = data$metadata
    
    data = get_data_for_Baby_Weight_1st(metadata_Baby_Weight_MR, reads_table_16s)
    reads_table_16s_MR = data$reads_table_output
    metadata_Baby_Weight_MR = data$metadata
    
    Baby_Weight_2 = metadata_mom_all$baby_weight
    
    metadata_mom_all[metadata_mom_all=='Yes'] = 1
    metadata_mom_all[metadata_mom_all=='No'] = 0
    
    names = metadata_mom_all$SampleID
    metadata = apply(metadata_mom_all, 2, as.numeric)
    
    keep = colSums(is.na(metadata)) < nrow(metadata)/3
    metadata = metadata[,keep]
    row.names(metadata) = names
    
    metadata = FilterFeatures(metadata, Baby_Weight_2)
    metadata = metadata[,!str_detect(colnames(metadata),'average|baby|ga_at_delivery|event
                                     |pregnancy|last_weight|last_pulse
                                     |last_ga_at_delivery|mom_feeding_method|last_weight
                                     |hispanic_or_latino|vaginal_penetration_frequency
                                     |hepatitis_c|weight_change|3rdtri
                                     |poor_appetite|vaginal_penetration_frequency|2ndtri')]
  }
  # two layer modeling
  {
    output = as.factor(Baby_Weight_2)
    input = reads_table_16s_MV
    set.seed(2022)
    input = FilterFeatures(input, output)
    prd_MV_1st=foreach(i=seq(nrow(input))) %dopar% twoLayer(input, output, i, ntree= ncol(input)*2)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
    ppp_MV_1st=vector()
    for(i in seq(nrow(input)))
      ppp_MV_1st[i] = prd_MV_1st[[i]]$Layer1$p2   # all 81 expected values in the first layer model
    
    output = as.factor(Baby_Weight_2)
    input = reads_table_16s_MR
    set.seed(2022)
    input = FilterFeatures(input, output)
    prd_MR_1st=foreach(i=seq(nrow(input))) %dopar% twoLayer(input, output, i, ntree= ncol(input)*2)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
    ppp_MR_1st=vector()
    for(i in seq(nrow(input)))
      ppp_MR_1st[i] = prd_MR_1st[[i]]$Layer1$p2   # all 81 expected values in the first layer model
    
    output = as.factor(Baby_Weight_2)
    input = reads_table_16s_MB
    set.seed(2022)
    input = FilterFeatures(input, output)
    prd_MB_1st=foreach(i=seq(nrow(input))) %dopar% twoLayer(input, output, i, ntree= ncol(input)*2)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
    ppp_MB_1st=vector()
    for(i in seq(nrow(input)))
      ppp_MB_1st[i] = prd_MB_1st[[i]]$Layer1$p2   # all 81 expected values in the first layer model
    
    
    output = as.factor(Baby_Weight_2)
    input = metadata
    set.seed(2022)
    input = FilterFeatures(input, output)
    prd_metadata_1st=foreach(i=seq(nrow(input))) %dopar% twoLayer(input, output, i, ntree= ncol(input)*2)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
    ppp_metadata_1st=vector()
    for(i in seq(nrow(input)))
      ppp_metadata_1st[i] = prd_metadata_1st[[i]]$Layer1$p2   # all 81 expected values in the first layer model
    
  }
  
  # calculate weight
  {
    n_ML = n_ML+1
    
    CMBW=vector()
    prf_all = list()
    prf_all_sum = vector()
    for(i in seq(nrow(metadata)) )
    {
      prf=vector()
      
      ppp_M_meta_combine=vector()   # second layer prediction expected values
      ppp_MV_1st_combine=vector()   # second layer prediction expected values
      ppp_MB_1st_combine=vector()   # second layer prediction expected values
      #      ppp_MR_1st_combine=vector()   # second layer prediction expected values
      
      for(j in seq((nrow(metadata)-1)))
      {
        ppp_M_meta_combine[j]=prd_metadata_1st[[i]]$Layer2[[j]]$p2
        ppp_MV_1st_combine[j]=prd_MV_1st[[i]]$Layer2[[j]]$p2
        ppp_MB_1st_combine[j]=prd_MB_1st[[i]]$Layer2[[j]]$p2
        #        ppp_MR_1st_combine[j]=prd_MR_1st[[i]]$Layer2[[j]]$p2
      }
      
      prf[1]=wilcox.test(ppp_M_meta_combine , as.numeric(as.factor(Baby_Weight_2[-i])))$p.value
      prf[2]=wilcox.test(ppp_MV_1st_combine , as.numeric(as.factor(Baby_Weight_2[-i])))$p.value
      prf[3]=wilcox.test(ppp_MB_1st_combine , as.numeric(as.factor(Baby_Weight_2[-i])))$p.value
      #      prf[3]=wilcox.test(ppp_MR_1st_combine , as.numeric(as.factor(Baby_Weight_2[-i])))$p.value
      
      prf=exp(-log10(prf))
      prf_all[[i]] = prf
      prf_all_sum[i] = sum(prf)
      CMBW[i]=(prf[1]*ppp_metadata_1st[i]+prf[2]*ppp_MV_1st[i]+prf[3]*ppp_MB_1st[i])/sum(prf)   # first layer expected value * -log(second layer pvalue) / sum(-log(second layer pvalue))
      
      #      CMBW[i]=(prf[1]*ppp_metadata_1st[i]+prf[2]*ppp_MV_1st[i]+prf[3]*ppp_MB_1st[i]+prf[4]*ppp_MR_1st[i])/sum(prf)   # first layer expected value * -log(second layer pvalue) / sum(-log(second layer pvalue))
    }
    
    roc(Baby_Weight_2, CMBW)$auc
    
    myPvs[n_ML]=roc(Baby_Weight_2, CMBW)$auc
    myPvs2[n_ML]=-log10(wilcox.test(CMBW ~ as.factor(Baby_Weight_2))$p.value)
    myPvs4[n_ML] = paste0(length(Baby_Weight_2),' (',sum(Baby_Weight_2 == 'Yes'),')')
    myPvs[n_ML]
  }
  
}