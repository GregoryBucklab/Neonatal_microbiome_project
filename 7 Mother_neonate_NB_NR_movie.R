setwd('/Users/binzhu/Desktop/Mom-baby/')
# N # follow script 4
{
  # NR
  pic_2 = pic
  pic_2$days_rel2birth[pic_2$SampleType == "MV1D"] = 0
  pic_2$days_rel2birth[pic_2$SampleType == "MRCD"] = 0
  pic_2$days_rel2birth[pic_2$SampleType == "MCKD"] = 0
  keep = !is.na(pic_2$days_rel2birth) & (pic_2$days_rel2birth ==0 |pic_2$days_rel2birth == 1 | pic_2$days_rel2birth == 2) & 
    (pic_2$SampleType == "BRCD" | pic_2$SampleType == "MV1D" | pic_2$SampleType == "MRCD" | pic_2$SampleType == "MCKD")
  sum(keep)
  pic_2 = pic_2[keep,]
  pic_2$ParticipantID_2 = pic_2$ParticipantID
  pic_2$ParticipantID_2[pic_2$SampleType == "MV1D"] = "MV"
  pic_2$ParticipantID_2[pic_2$SampleType == "MRCD"] = "MR"
  pic_2$ParticipantID_2[pic_2$SampleType == "MCKD"] = "MB"
  
  
  keep = pic_2$ParticipantID_2 != "MV" & pic_2$ParticipantID_2 != "MR" & pic_2$ParticipantID_2 != "MB"
  pic_5 = pic_2[keep,]
  pic_4 = pic_2[!keep,]
  pic_5$ParticipantID_2 = paste0('NR_',pic_5$ParticipantID)
  
  pid_list = unique(pic_5$ParticipantID)
  data = as.data.frame(matrix(data = 0, ncol = 3, nrow = length(pid_list)))
  row.names(data) = pid_list
  for (a in 1: nrow(pic_5)) {
    x = which(pid_list == pic_5$ParticipantID[a])
    y = pic_5$days_rel2birth[a] +1
    
    data[x,y]=data[x,y]+1
  }
  data = data[rowSums(data) == 2 & data$V1 != 0,]
  
  pic_4 = pic_4[pic_4$ParticipantID %in% c(row.names(data)),] 
  pic_5 = pic_5[pic_5$ParticipantID %in% c(row.names(data)),] 
  
  
  
  keep = pic_5$ParticipantID %in% row.names(data)[data$V2 == 0] & pic_5$days_rel2birth ==0
  sum(keep)
  pic_5_plus = pic_5[keep,]
  pic_5_plus$days_rel2birth = 1
  pic_5 = rbind(pic_5,pic_5_plus)
  
  keep = pic_5$ParticipantID %in% row.names(data)[data$V3 == 0] & pic_5$days_rel2birth ==1
  sum(keep)
  pic_5_plus = pic_5[keep,]
  pic_5_plus$days_rel2birth = 2
  pic_5 = rbind(pic_5,pic_5_plus)
  
  pid_list = unique(pic_5$ParticipantID)
  data = as.data.frame(matrix(data = 0, ncol = 3, nrow = length(pid_list)))
  row.names(data) = pid_list
  for (a in 1: nrow(pic_5)) {
    x = which(pid_list == pic_5$ParticipantID[a])
    y = pic_5$days_rel2birth[a] +1
    
    data[x,y]=data[x,y]+1
  }
  
  pic_2$days_rel2birth = as.numeric(as.character(pic_2$days_rel2birth))
  
  
  
  
  keep = pic_4$SampleType == "MV1D" | pic_4$SampleType == "MCKD" | pic_4$SampleType == "MRCD"
  pic_3 = pic_4[keep,]
  pic_3$days_rel2birth = 1
  pic_6 = pic_4[keep,]
  pic_6$days_rel2birth = 2
  
  pic_4 = rbind(pic_4,pic_3,pic_6)
  
  pic_2_NR = rbind(pic_4,pic_5)
  
  
  
  # NB
  pic_2 = pic
  pic_2$days_rel2birth[pic_2$SampleType == "MV1D"] = 0
  pic_2$days_rel2birth[pic_2$SampleType == "MRCD"] = 0
  pic_2$days_rel2birth[pic_2$SampleType == "MCKD"] = 0
  keep = !is.na(pic_2$days_rel2birth) & (pic_2$days_rel2birth ==0 |pic_2$days_rel2birth == 1 | pic_2$days_rel2birth == 2) & 
    (pic_2$SampleType == "BCKD" | pic_2$SampleType == "MV1D" | pic_2$SampleType == "MRCD" | pic_2$SampleType == "MCKD")
  sum(keep)
  pic_2 = pic_2[keep,]
  pic_2$ParticipantID_2 = pic_2$ParticipantID
  pic_2$ParticipantID_2[pic_2$SampleType == "MV1D"] = "MV"
  pic_2$ParticipantID_2[pic_2$SampleType == "MRCD"] = "MR"
  pic_2$ParticipantID_2[pic_2$SampleType == "MCKD"] = "MB"
  
  
  keep = pic_2$ParticipantID_2 != "MV" & pic_2$ParticipantID_2 != "MR" & pic_2$ParticipantID_2 != "MB"
  pic_5 = pic_2[keep,]
  pic_4 = pic_2[!keep,]
  pic_5$ParticipantID_2 = paste0('NB_',pic_5$ParticipantID)
  
  pid_list = unique(pic_5$ParticipantID)
  data = as.data.frame(matrix(data = 0, ncol = 3, nrow = length(pid_list)))
  row.names(data) = pid_list
  for (a in 1: nrow(pic_5)) {
    x = which(pid_list == pic_5$ParticipantID[a])
    y = pic_5$days_rel2birth[a] +1
    
    data[x,y]=data[x,y]+1
  }
  data = data[rowSums(data) == 2 & data$V1 != 0,]
  
  pic_4 = pic_4[pic_4$ParticipantID %in% c(row.names(data)),] 
  pic_5 = pic_5[pic_5$ParticipantID %in% c(row.names(data)),] 
  
  
  
  keep = pic_5$ParticipantID %in% row.names(data)[data$V2 == 0] & pic_5$days_rel2birth ==0
  sum(keep)
  pic_5_plus = pic_5[keep,]
  pic_5_plus$days_rel2birth = 1
  pic_5 = rbind(pic_5,pic_5_plus)
  
  keep = pic_5$ParticipantID %in% row.names(data)[data$V3 == 0] & pic_5$days_rel2birth ==1
  sum(keep)
  pic_5_plus = pic_5[keep,]
  pic_5_plus$days_rel2birth = 2
  pic_5 = rbind(pic_5,pic_5_plus)
  
  pid_list = unique(pic_5$ParticipantID)
  data = as.data.frame(matrix(data = 0, ncol = 3, nrow = length(pid_list)))
  row.names(data) = pid_list
  for (a in 1: nrow(pic_5)) {
    x = which(pid_list == pic_5$ParticipantID[a])
    y = pic_5$days_rel2birth[a] +1
    
    data[x,y]=data[x,y]+1
  }
  
  pic_2$days_rel2birth = as.numeric(as.character(pic_2$days_rel2birth))
  
  
  
  
  keep = pic_4$SampleType == "MV1D" | pic_4$SampleType == "MCKD" | pic_4$SampleType == "MRCD"
  pic_3 = pic_4[keep,]
  pic_3$days_rel2birth = 1
  pic_6 = pic_4[keep,]
  pic_6$days_rel2birth = 2
  
  pic_4 = rbind(pic_4,pic_3,pic_6)
  
  pic_2_NB = rbind(pic_4,pic_5)
  
  
  
  # NR + NB
  pic_2 = rbind(pic_2_NR,pic_2_NB)
  pic_2$SampleType <- pic_2$SampleType
  pic_2$shape <- pic_2$Flag
  pic_2$SampleType = as.factor(pic_2$SampleType)
  
  ggplot(pic_2, aes(X1, X2, color = SampleType)) +
    geom_point(size = 1.5)+
    scale_color_manual(values=color1[c(1,2,4,5,6)])+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))+theme_bw()
  ggsave('tsne_111.pdf',width=6, height=5)
  
  setwd('/Users/binzhu/Desktop/Mom-baby/')
  library(gganimate)
  
  pic_2 = pic_2[order(pic_2$ParticipantID_2),]
  pic_2$ParticipantID_2 = as.factor(pic_2$ParticipantID_2)
  p <- ggplot(pic_2, 
              aes(x = X1, y=X2, colour = ParticipantID_2)) +
    geom_point(show.legend = F, alpha = 0.7) + 
    scale_color_manual(values=c("#FEF840", "#F27E33", "#F180BF", rep("#E43B2D",25), rep("#0070FF",13)))+
    xlab(paste0("t-SNE1")) +
    ylab(paste0("t-SNE2"))+theme_bw()
  
  p + transition_time(days_rel2birth) +
    labs(title = "Days after birth: {frame_time}")
  
  
  pic_3 = pic_2
  pic_3 = pic_3[pic_3$SampleType %in% c('BCKD','MCKD','MV1D','MRCD'),]
  
  pic_3$ParticipantID_2 = as.factor(pic_3$ParticipantID_2)
  p <- ggplot(pic_3, 
              aes(x = X1, y=X2, colour = ParticipantID_2)) +
    geom_point(show.legend = F, alpha = 0.7) + 
    scale_color_manual(values=c("#FEF840","#F27E33", "#F180BF",  rep("#E43B2D",25)))+
    xlab(paste0("t-SNE1")) +
    ylab(paste0("t-SNE2"))+theme_bw()
  
  p + transition_time(days_rel2birth) +
    labs(title = "Days after birth: {frame_time}")
  
  
  pic_3 = pic_2
  pic_3 = pic_3[pic_3$SampleType %in% c('BRCD','MCKD','MV1D','MRCD'),]
  
  pic_3$ParticipantID_2 = as.factor(pic_3$ParticipantID_2)
  p <- ggplot(pic_3, 
              aes(x = X1, y=X2, colour = ParticipantID_2)) +
    geom_point(show.legend = F, alpha = 0.7) + 
    scale_color_manual(values=c("#FEF840","#F27E33", "#F180BF",  rep("#0070FF",13)))+
    xlab(paste0("t-SNE1")) +
    ylab(paste0("t-SNE2"))+theme_bw()
  
  p + transition_time(days_rel2birth) +
    labs(title = "Days after birth: {frame_time}")
}

# M # follow script 1
{
  keep = metadata_16s$Mombaby == 'Mom'
  metadata_16s = metadata_16s[keep,]
  reads_table_16s = reads_table_16s[,keep]
  
  data = prepare_reads_table(reads_table_16s, metadata_16s, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores =8)
  
  reads_table = data$reads_table
  metadata = data$metadata
  
  reads_table = as.data.frame(t(reads_table))
  reads_table = as.matrix(decostand(reads_table,method = 'rclr'))
  
  tsne <- Rtsne(reads_table, dims = 2, perplexity=200, verbose=TRUE, max_iter = 2000,check_duplicates = FALSE)
  
  pic <- tsne$Y
  pic <- data.frame(pic,row.names(reads_table))
  colnames(pic) <- c('X1','X2','SampleID')
  pic = merge(pic, metadata, 'SampleID')
  pic$factor <- pic$SampleType
  pic$factor = as.factor(pic$factor)
  
  color1=c('#FEF840','#F27E33','#F180BF')
  ggplot(pic, aes(X1, X2, color = factor)) +
    geom_point(size = 1.5)+
    scale_color_manual(values=color1)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))+theme_bw()
  ggsave('tsne_mom.pdf',width=6, height=5)
  
  {
    {
      pic_3 = pic
      pic_3 = pic_3[pic_3$SampleType == 'MCKD',]
      pic_3$days_rel2birth = as.numeric(as.character(pic_3$days_rel2birth))
      pic_3 = pic_3[!is.na(pic_3$days_rel2birth),]
      pic_3 = pic_3[pic_3$days_rel2birth <=0,]
      pic_3 = pic_3[,c(1,2,3,4,7,21)]
      
      time_list = unique(pic_3$days_rel2birth)
      time_list = time_list[order(time_list)]
      participant_list = unique(pic_3$ParticipantID)
      
      for (a in 1: length(participant_list)) {
        n = which(pic_3$ParticipantID == participant_list[a])
        pic_4 = pic_3[n,]
        pic_4 = pic_4[order(pic_4$days_rel2birth),]
        
        x = which(time_list %in% (pic_4$days_rel2birth))
        
        pic_5 = as.data.frame(matrix(data = NA, ncol = ncol(pic_3), nrow = length(time_list)))
        colnames(pic_5) = colnames(pic_3)
        
        m = 0
        for (b in 1: length(time_list)) {
          if (time_list[b] %in% time_list[x]) {
            m = m+1
            pic_5[b,] = pic_4[m,]
            
          } else if (m == 0) {next} else {
            pic_5[b,] = pic_4[m,]
            pic_5$days_rel2birth[b] = time_list[b]
          }
          
        }
        
        if (a == 1) {
          pic_MCKD = pic_5
        } else {
          pic_MCKD = rbind(pic_MCKD, pic_5)
        }
      }
      pic_MCKD = pic_MCKD[!is.na(pic_MCKD$SampleID),]
    }
    {
      pic_3 = pic
      pic_3 = pic_3[pic_3$SampleType == 'MV1D',]
      pic_3$days_rel2birth = as.numeric(as.character(pic_3$days_rel2birth))
      pic_3 = pic_3[!is.na(pic_3$days_rel2birth),]
      pic_3 = pic_3[pic_3$days_rel2birth <=0,]
      pic_3 = pic_3[,c(1,2,3,4,7,21)]
      
      time_list = unique(pic_3$days_rel2birth)
      time_list = time_list[order(time_list)]
      participant_list = unique(pic_3$ParticipantID)
      
      for (a in 1: length(participant_list)) {
        n = which(pic_3$ParticipantID == participant_list[a])
        pic_4 = pic_3[n,]
        pic_4 = pic_4[order(pic_4$days_rel2birth),]
        
        x = which(time_list %in% (pic_4$days_rel2birth))
        
        pic_5 = as.data.frame(matrix(data = NA, ncol = ncol(pic_3), nrow = length(time_list)))
        colnames(pic_5) = colnames(pic_3)
        
        m = 0
        for (b in 1: length(time_list)) {
          if (time_list[b] %in% time_list[x]) {
            m = m+1
            pic_5[b,] = pic_4[m,]
            
          } else if (m == 0) {next} else {
            pic_5[b,] = pic_4[m,]
            pic_5$days_rel2birth[b] = time_list[b]
          }
          
        }
        
        if (a == 1) {
          pic_MV1D = pic_5
        } else {
          pic_MV1D = rbind(pic_MV1D, pic_5)
        }
      }
      pic_MV1D = pic_MV1D[!is.na(pic_MV1D$SampleID),]
    }
    {
      pic_3 = pic
      pic_3 = pic_3[pic_3$SampleType == 'MRCD',]
      pic_3$days_rel2birth = as.numeric(as.character(pic_3$days_rel2birth))
      pic_3 = pic_3[!is.na(pic_3$days_rel2birth),]
      pic_3 = pic_3[pic_3$days_rel2birth <=0,]
      pic_3 = pic_3[,c(1,2,3,4,7,21)]
      
      time_list = unique(pic_3$days_rel2birth)
      time_list = time_list[order(time_list)]
      participant_list = unique(pic_3$ParticipantID)
      
      for (a in 1: length(participant_list)) {
        n = which(pic_3$ParticipantID == participant_list[a])
        pic_4 = pic_3[n,]
        pic_4 = pic_4[order(pic_4$days_rel2birth),]
        
        x = which(time_list %in% (pic_4$days_rel2birth))
        
        pic_5 = as.data.frame(matrix(data = NA, ncol = ncol(pic_3), nrow = length(time_list)))
        colnames(pic_5) = colnames(pic_3)
        
        m = 0
        for (b in 1: length(time_list)) {
          if (time_list[b] %in% time_list[x]) {
            m = m+1
            pic_5[b,] = pic_4[m,]
            
          } else if (m == 0) {next} else {
            pic_5[b,] = pic_4[m,]
            pic_5$days_rel2birth[b] = time_list[b]
          }
          
        }
        
        if (a == 1) {
          pic_MRCD = pic_5
        } else {
          pic_MRCD = rbind(pic_MRCD, pic_5)
        }
      }
      pic_MRCD = pic_MRCD[!is.na(pic_MRCD$SampleID),]
    }
    
    pic_all = rbind(pic_MCKD, pic_MRCD,pic_MV1D)
    pic_all$SampleType = as.factor(pic_all$SampleType)
    
    p <- ggplot(pic_all, 
                aes(x = X1, y=X2, colour = SampleType)) +
      geom_point(show.legend = F, alpha = 0.7) + 
      scale_color_manual(values=c("#FEF840","#F27E33", "#F180BF"))+
      xlab(paste0("t-SNE1")) +
      ylab(paste0("t-SNE2"))+theme_bw()
    
    p <- p + transition_time(days_rel2birth) +
      labs(title = "Days after birth: {frame_time}")
    setwd('/Users/binzhu/Desktop/Mom-baby/1')
    animate(p, nframes = 5000, fps = 1,  end_pause = 0)
  }


  
  
}


