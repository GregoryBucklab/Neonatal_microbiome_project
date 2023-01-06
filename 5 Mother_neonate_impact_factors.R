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
  correlation$size_power = NA
  
  
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
    
    if (factor_number[1] < 3 | factor_number[2] < 3) { 
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
    temp_2 = temp
    temp_2$V2[temp_2$V2 == 'No'] =0
    temp_2$V2[temp_2$V2 == 'Yes'] = 1
    temp_2$V2 = as.numeric(as.character(temp_2$V2))
    temp_2$V1 = as.numeric(as.character(temp_2$V1))
    correlation$shannon_pvalue[b] = wilcox.test(temp_2$V2,temp_2$V1)$p.value
    
    # calculate sample size tested by Mann-Whitely U test
    library(fitdistrplus)
    pdf(paste0('distribution_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=7, height=7)
    distribution = descdist(shannon_2, discrete=F)
    dev.off()
    
    library(rstatix)
    data = wilcox_effsize(temp_2,V2~V1)
    data = data$effsize
    data = data[!is.nan(data)]
    
    library(wmwpow)
    output = shiehpow(n = factor_number[1], m = factor_number[2], p = median(data), alpha = 0.05, dist = "norm", sides = "two.sided")
    correlation$size_power[b] = output$power
      
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
  correlation$size_power = NA
  
  
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
    
    if (str_detect(colnames(metadata)[b], 'worries')) {
      temp = as.data.frame(cbind(shannon_2,metadata_2))
      colnames(temp) = c('V1','V2')
      temp$V1 = as.numeric(as.character(temp$V1))
      temp$V2 = (as.factor(temp$V2))
      
      correlation$shannon_lm_function[b] = kruskal.test(temp$V1, temp$V2)$p.value
      
      ggplot(temp, aes(x = V2, y = V1)) + 
        geom_boxplot() +
        geom_jitter() +
        labs(x = correlation$Factor[b], 
             y = paste0("Shannon index of ",type1_all[a]))+ 
        ggtitle(paste0("pvalue = ", correlation$shannon_pvalue[b]))  +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))+ theme_bw()
      
      ggsave(paste0('alpha_shannon_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=3, height=3)
      
    }
    
    
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
      
      # sample size power
      output = wp.regression(n = nrow(temp), p1 = 1, p2 = 0, f2 = as.numeric(linearMod$fstatistic[1]), alpha = 0.05,
                                                power = NULL, type="regular")
      correlation$size_power[b] = output$power
      
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
  correlation$size_power = NA
  
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
      
      # calculate sample size tested by Mann-Whitely U test
      library(fitdistrplus)
      pdf(paste0('distribution_',type1_all[a],'_',colnames(metadata)[b],'.pdf'),width=7, height=7)
      distribution = descdist(shannon_2, discrete=F)
      dev.off()
      
      library(wmwpow)
      output = shiehpow(n = factor_number[1], m = factor_number[2], p = 0.80, alpha = 0.05, dist = "norm", sides = "two.sided")
      correlation$size_power[b] = output$power
      
      
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
  correlation$size_power = NA
  
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
      
      # sample size power
      output = wp.regression(n = nrow(temp), p1 = 1, p2 = 0, f2 = as.numeric(linearMod$fstatistic[1]), alpha = 0.05,
                             power = NULL, type="regular")
      correlation$size_power[b] = output$power
      
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

correlation_all_cytokine = as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(correlation_all_cytokine) = c('Type','Factor','pvalue', 'Case_number','Rvalue','Slope')

correlation_all_lipid = as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(correlation_all_lipid) = c('Type','Factor','pvalue', 'Case_number','Rvalue','Slope')

involved_factors = read.csv('involved_factors.csv', header = F)
involved_factors = involved_factors$V1
mediation_factors = read.csv('Mediation.csv', header = F)

n_case_number_list = 0
n_abundance_analysis = 0
case_number_list = vector()
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
  correlation_all = as.data.frame(matrix(data = NA, nrow = nrow(reads_table_cytokines_mom), ncol = 6))
  colnames(correlation_all) = c('Type','Factor','pvalue', 'Case_number','Rvalue','Slope')
  correlation_all[,1] = type1_all[a]
  correlation_all[,4] = ncol(reads_table_cytokines_mom)
  correlation_all$size_power = NA
  
  for (b in 1:nrow(reads_table_cytokines_mom)) {
    correlation_all[b,2] = row.names(reads_table_cytokines_mom)[b]
    
    #    pvalue = cor.test(alpha.shannon_diversity_cytokine, as.numeric(reads_table_cytokines_mom[b,]), method = "pearson")
    #    pvalue = pvalue$p.value
    data = as.data.frame(matrix(data = NA, nrow = length(alpha.shannon_diversity_cytokine), ncol = 2))
    data$V1 = alpha.shannon_diversity_cytokine
    data$V2 = as.numeric(reads_table_cytokines_mom[b,])
    
    linearMod <- lm(V1 ~ V2, data=data)
    r_squared <- summary(linearMod)$r.squared
    pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])
    slope <- as.numeric(linearMod$coefficients[2])  
    
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    slope = format(round(slope, 5), nsmall = 5)
    
    correlation_all[b,3] = pvalue
    correlation_all[b,5] = r_squared
    correlation_all[b,6] = slope
    
    # sample size power
    output = wp.regression(n = nrow(data), p1 = 1, p2 = 0, f2 = as.numeric(summary(linearMod)$fstatistic[1]), alpha = 0.05,
                           power = NULL, type="regular")
    correlation_all$size_power[b] = output$power
    
    
    if (pvalue <= 0.05) {
      ggplot(data, aes(x = V2, y = V1)) + 
        geom_point() +
        stat_smooth(method = "lm", col = "blue")+
        labs(x = paste0('Log10 ',row.names(reads_table_cytokines_mom)[b]), 
             y = paste0('Shannon index of ', type1_all[a]))+ 
        ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue, '\n',"Slope = ", slope))  +
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
  
  ### mediation cytokines ###
  {
    # mediation analysis using sem   lavaan R package
    for (b in 1:nrow(reads_table_cytokines_mom)) {
      if (row.names(reads_table_cytokines_mom)[b] %in% mediation_factors$V2 == F) {
        next
      }
      
      {
        n = which(mediation_factors$V2 == row.names(reads_table_cytokines_mom)[b])
        if (type1_all[a] %in% mediation_factors$V1[n] == F) {
          next
        }
      }
      
      
      #      metadata_baby_cytokine = metadata_baby_cytokine[,colnames(metadata_baby_cytokine) %in% involved_factors]
      
      
      reads_table_cytokines_mom_2 = as.numeric(reads_table_cytokines_mom[b,])
      
      trials = c(1: ncol(metadata_baby_cytokine))
      func_1 = function(trial) {
        data = as.data.frame(matrix(data = NA, nrow = length(alpha.shannon_diversity_cytokine), ncol = 3))
        data$V1 = alpha.shannon_diversity_cytokine
        
        data$V2 = reads_table_cytokines_mom_2
        
        data$V3 = metadata_baby_cytokine[,trial]
        data$V4 = metadata_baby_cytokine[,trial]
        data$V3 = as.numeric(factor(data$V3))
        
        keep = rowSums(is.na(data)) == 0 
        data = data[keep,]
        
        if (nrow(data) <8) {
          return(NA)
        }
        
        keep = length(unique(data$V1)) > 1 & length(unique(data$V2)) > 1 & length(unique(data$V3)) > 1
        if (!keep) {
          return(NA)
        }
        
        unique_factor = unique(data$V4)
        if (length(unique_factor) == 2) {
          n1 = sum(data$V4 == unique_factor[1])
          n2 = sum(data$V4 == unique_factor[2])
          
          if (n1 < 3 | n2 < 3) {
            return(NA)
          }
          
          case_number = paste0(n1+n2, " (",sum(data$V4 == 'Yes'),")")
        } else (case_number = nrow(data))
        
        
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
        pvalue_sem[pvalue_sem == 'V2' & !is.na(pvalue_sem)] = row.names(reads_table_cytokines_mom)[b]
        pvalue_sem[pvalue_sem == 'V3' & !is.na(pvalue_sem)] = colnames(metadata_baby_cytokine)[trial]
        
        return(c(list(data = data),list(pvalue_sem = pvalue_sem),list(case_number = case_number)))
      }
      data_sem = mclapply(trials, func_1, mc.cores = 8)
      
      
      for (c in 1: ncol(metadata_baby_cytokine)) {
        if (length(data_sem[[c]])==1) {
          next
        }
        
        pvalue_sem = data_sem[[c]]$pvalue_sem
        data = data_sem[[c]]$data
        
        case_number = data_sem[[c]]$case_number
        n_case_number_list = n_case_number_list +1
        case_number_list[n_case_number_list] = paste0('Mediation_',colnames(metadata_baby_cytokine)[c],'_',row.names(reads_table_cytokines_mom)[b],'_',type1_all[a],'_',case_number)
        
        if (pvalue_sem[7,8] <= 0.05 & pvalue_sem[1,8] <= 0.05 & pvalue_sem[2,8] <= 0.05 & pvalue_sem[3,8] <= 0.05) {
          write.csv(pvalue_sem,paste0('Mediation_',colnames(metadata_baby_cytokine)[c],'_',row.names(reads_table_cytokines_mom)[b],'_',type1_all[a],'.csv'))
          
          pdf(file = paste0('Mediation_',colnames(metadata_baby_cytokine)[c],'_',row.names(reads_table_cytokines_mom)[b],'_',type1_all[a],'.pdf'),width=5, height=4)
          with(data, mediation.effect.plot(x = V3, mediator = V2, dv = V1, ylab = paste0("Shannon index of ", type1_all[a]), xlab = row.names(reads_table_cytokines_mom)[b]))
          dev.off()
        }
        
        
      }
    }
    
    
    
  }
  
  ### test correlation lipid ###
  correlation_all = as.data.frame(matrix(data = NA, nrow = nrow(reads_table_lipid_mom), ncol = 6))
  colnames(correlation_all) = c('Type','Factor','pvalue', 'Case_number','Rvalue','Slope')
  correlation_all[,1] = type1_all[a]
  correlation_all[,4] = ncol(reads_table_lipid_mom)
  correlation_all$size_power = NA
  
  for (b in 1:nrow(reads_table_lipid_mom)) {
    correlation_all[b,2] = row.names(reads_table_lipid_mom)[b]
    
    keep = !is.na(reads_table_lipid_mom[b,])
    if (sum(keep) < 8) {
      next
    }
    
    data = cbind(alpha.shannon_diversity_lipid,as.numeric(reads_table_lipid_mom[b,]))
    keep = !is.na(data[,2])
    data = data[keep,]
    
    #    pvalue = cor.test(alpha.shannon_diversity_lipid, as.numeric(reads_table_lipid_mom[b,]), method = "pearson")
    #    pvalue = pvalue$p.value
    #    correlation_all[b,3] = pvalue
    
    data = as.data.frame(matrix(data = NA, nrow = length(alpha.shannon_diversity_lipid), ncol = 2))
    data$V1 = alpha.shannon_diversity_lipid
    data$V2 = as.numeric(reads_table_lipid_mom[b,])
    
    keep = !is.na(reads_table_lipid_mom[b,])
    data = data[keep,]
    
    linearMod <- lm(V1 ~ V2, data=data)
    r_squared <- summary(linearMod)$r.squared
    pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])
    slope <- as.numeric(linearMod$coefficients[2])  
    
    r_squared = format(round(r_squared, 5), nsmall = 5)
    pvalue = format(round(pvalue, 5), nsmall = 5)
    slope = format(round(slope, 5), nsmall = 5)
    
    correlation_all[b,3] = pvalue
    correlation_all[b,5] = r_squared
    correlation_all[b,6] = slope
    
    # sample size power
    output = wp.regression(n = nrow(data), p1 = 1, p2 = 0, f2 = as.numeric(summary(linearMod)$fstatistic[1]), alpha = 0.05,
                           power = NULL, type="regular")
    correlation_all$size_power[b] = output$power
    
    if (pvalue <= 0.05) {
      ggplot(data, aes(x = V2, y = V1)) + 
        geom_point() +
        stat_smooth(method = "lm", col = "blue")+
        labs(x = paste0('Log10 ',row.names(reads_table_lipid_mom)[b]), 
             y = paste0('Shannon index of ', type1_all[a]))+ 
        ggtitle(paste0("R-value = ", r_squared, '\n',"P-value = ", pvalue, '\n',"Slope = ", slope))  +
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))+ theme_bw()
      ggsave(paste0('Shannon_vs_lipid_',type1_all[a],'_',row.names(reads_table_lipid_mom)[b],'.pdf'),width=3, height=3)
      
    }
    
    
    
    
    
    
  }
  
  if (a == 1) {
    correlation_all_lipid = correlation_all
  } else {
    correlation_all_lipid = rbind(correlation_all_lipid,correlation_all)
  }
  
  ### mediation lipid using sem  lavaan R package ###
  {
    for (b in 47:nrow(reads_table_lipid_mom)) {
      if (row.names(reads_table_lipid_mom)[b] %in% mediation_factors$V2 == F) {
        next
      }
      
      {
        n = which(mediation_factors$V2 == row.names(reads_table_lipid_mom)[b])
        if (type1_all[a] %in% mediation_factors$V1[n] == F) {
          next
        }
      }
      
      #      metadata_baby_lipid = metadata_baby_lipid[,colnames(metadata_baby_lipid) %in% involved_factors]
      
      
      reads_table_lipid_mom_2 = as.numeric(reads_table_lipid_mom[b,])
      
      trials = c(1: ncol(metadata_baby_lipid))
      func_1 = function(trial) {
        data = as.data.frame(matrix(data = NA, nrow = length(alpha.shannon_diversity_lipid), ncol = 3))
        data$V1 = alpha.shannon_diversity_lipid
        
        data$V2 = reads_table_lipid_mom_2
        
        data$V3 = metadata_baby_lipid[,trial]
        data$V4 = metadata_baby_lipid[,trial]
        data$V3 = as.numeric(factor(data$V3))
        
        keep = rowSums(is.na(data)) == 0 
        data = data[keep,]
        
        if (nrow(data) <8) {
          return(NA)
        }
        
        keep = length(unique(data$V1)) > 1 & length(unique(data$V2)) > 1 & length(unique(data$V3)) > 1
        if (!keep) {
          return(NA)
        }
        
        unique_factor = unique(data$V4)
        if (length(unique_factor) == 2) {
          n1 = sum(data$V4 == unique_factor[1])
          n2 = sum(data$V4 == unique_factor[2])
          
          if (n1 < 3 | n2 < 3) {
            return(NA)
          }
          
          case_number = paste0(n1+n2, " (",sum(data$V4 == 'Yes'),")")
        } else (case_number = nrow(data))
        
        
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
        pvalue_sem[pvalue_sem == 'V3' & !is.na(pvalue_sem)] = colnames(metadata_baby_lipid)[trial]
        
        return(c(list(data = data),list(pvalue_sem = pvalue_sem)))
      }
      data_sem = mclapply(trials, func_1, mc.cores = 8)
      
      for (c in 1: ncol(metadata_baby_lipid)) {
        if (length(data_sem[[c]])==1) {
          next
        }
        
        pvalue_sem = data_sem[[c]]$pvalue_sem
        data = data_sem[[c]]$data
        
        case_number = data_sem[[c]]$case_number
        n_case_number_list = n_case_number_list +1
        case_number_list[n_case_number_list] = paste0('Mediation_',colnames(metadata_baby_lipid)[c],'_',row.names(reads_table_lipid_mom)[b],'_',type1_all[a],'_',case_number)
        
        if (pvalue_sem[7,8] <= 0.05 & pvalue_sem[1,8] <= 0.05 & pvalue_sem[2,8] <= 0.05 & pvalue_sem[3,8] <= 0.05) {
          write.csv(pvalue_sem,paste0('Mediation_',colnames(metadata_baby_lipid)[c],'_',row.names(reads_table_lipid_mom)[b],'_',type1_all[a],'.csv'))
          
          pdf(file = paste0('Mediation_',colnames(metadata_baby_lipid)[c],'_',row.names(reads_table_lipid_mom)[b],'_',type1_all[a],'.pdf'),width=5, height=4)
          with(data, mediation.effect.plot(x = V3, mediator = V2, dv = V1, ylab = paste0("Shannon index of ", type1_all[a]), xlab = row.names(reads_table_lipid_mom)[b]))
          dev.off()
        }
        
        
      }
    }
  }
}

write.csv(correlation_all_lipid,'correlation_shannon_lipid.csv')
write.csv(correlation_all_cytokine,'correlation_shannon_cytokine.csv')
write.csv(case_number_list,'Mediation_analysis_case_number_list.csv')
write.csv(abundance_analysis_all,'abundance_analysis_all.csv')









