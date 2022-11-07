setwd('/Users/binzhu/Desktop/Mom-baby/')
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
pic_2$SampleType = as.factor(pic$SampleType)

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
  scale_color_manual(values=c("#FEF840", "#F27E33", "#F180BF", rep("#E43B2D",24), rep("#0070FF",13)))+
  xlab(paste0("t-SNE1")) +
  ylab(paste0("t-SNE2"))+theme_bw()

p + transition_time(days_rel2birth) +
  labs(title = "Days after birth: {frame_time}")

