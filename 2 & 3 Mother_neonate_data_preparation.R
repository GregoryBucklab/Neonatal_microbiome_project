#library('dplyr')
#library("ggpubr")

##### get metadata and 16s, cytokine reads table #####
# run 'get_pop2_ptb47_metadata.R' 'get_lipidomics.R'
setwd('/Users/binzhu/Desktop/Mom-baby/')

##### data preparation #####
# modify metadata
{
  metadata_16s$BMI = metadata_16s$weight * 0.0283495 / ((metadata_16s$height * 0.0254)^2)
  
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_vdischarge'] = 'vdischarge_pregnancy'
  metadata_16s$vdischarge_pregnancy[metadata_16s$vdischarge_pregnancy == 'No'] = 'xx'
  metadata_16s$vdischarge_pregnancy[metadata_16s$vdischarge_pregnancy == 'Yes'] = 'No'
  metadata_16s$vdischarge_pregnancy[metadata_16s$vdischarge_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_vodor'] = 'vodor_pregnancy'
  metadata_16s$vodor_pregnancy[metadata_16s$vodor_pregnancy == 'No'] = 'xx'
  metadata_16s$vodor_pregnancy[metadata_16s$vodor_pregnancy == 'Yes'] = 'No'
  metadata_16s$vodor_pregnancy[metadata_16s$vodor_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_diarrhea '] = 'diarrhea_pregnancy'
  metadata_16s$diarrhea_pregnancy[metadata_16s$diarrhea_pregnancy == 'No'] = 'xx'
  metadata_16s$diarrhea_pregnancy[metadata_16s$diarrhea_pregnancy == 'Yes'] = 'No'
  metadata_16s$diarrhea_pregnancy[metadata_16s$diarrhea_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_respiratory'] = 'respiratory_pregnancy'
  metadata_16s$respiratory_pregnancy[metadata_16s$respiratory_pregnancy == 'No'] = 'xx'
  metadata_16s$respiratory_pregnancy[metadata_16s$respiratory_pregnancy == 'Yes'] = 'No'
  metadata_16s$respiratory_pregnancy[metadata_16s$respiratory_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_contractions'] = 'contractions_pregnancy'
  metadata_16s$contractions_pregnancy[metadata_16s$contractions_pregnancy == 'No'] = 'xx'
  metadata_16s$contractions_pregnancy[metadata_16s$contractions_pregnancy == 'Yes'] = 'No'
  metadata_16s$contractions_pregnancy[metadata_16s$contractions_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_vitching'] = 'vitching_pregnancy'
  metadata_16s$vitching_pregnancy[metadata_16s$vitching_pregnancy == 'No'] = 'xx'
  metadata_16s$vitching_pregnancy[metadata_16s$vitching_pregnancy == 'Yes'] = 'No'
  metadata_16s$vitching_pregnancy[metadata_16s$vitching_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_false_labor'] = 'false_labor_pregnancy'
  metadata_16s$false_labor_pregnancy[metadata_16s$false_labor_pregnancy == 'No'] = 'xx'
  metadata_16s$false_labor_pregnancy[metadata_16s$false_labor_pregnancy == 'Yes'] = 'No'
  metadata_16s$false_labor_pregnancy[metadata_16s$false_labor_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_high_bp'] = 'high_bp_pregnancy'
  metadata_16s$high_bp_pregnancy[metadata_16s$high_bp_pregnancy == 'No'] = 'xx'
  metadata_16s$high_bp_pregnancy[metadata_16s$high_bp_pregnancy == 'Yes'] = 'No'
  metadata_16s$high_bp_pregnancy[metadata_16s$high_bp_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_pets_owner'] = 'pets_owner_pregnancy'
  metadata_16s$pets_owner_pregnancy[metadata_16s$pets_owner_pregnancy == 'No'] = 'xx'
  metadata_16s$pets_owner_pregnancy[metadata_16s$pets_owner_pregnancy == 'Yes'] = 'No'
  metadata_16s$pets_owner_pregnancy[metadata_16s$pets_owner_pregnancy == 'xx'] = 'Yes'
  colnames(metadata_16s)[colnames(metadata_16s) == 'no_vbleeding'] = 'vbleeding_pregnancy'
  metadata_16s$vbleeding_pregnancy[metadata_16s$vbleeding_pregnancy == 'No'] = 'xx'
  metadata_16s$vbleeding_pregnancy[metadata_16s$vbleeding_pregnancy == 'Yes'] = 'No'
  metadata_16s$vbleeding_pregnancy[metadata_16s$vbleeding_pregnancy == 'xx'] = 'Yes'
  
  colnames(metadata_16s)[colnames(metadata_16s) == 'mom_delivery_method'] = 'C_section'
  colnames(metadata_16s)[colnames(metadata_16s) == 'Delivery'] = 'Preterm'
  colnames(metadata_16s)[colnames(metadata_16s) == 'baby_icu'] = 'baby_nicu'
  colnames(metadata_16s)[colnames(metadata_16s) == 'baby_problems_none'] = 'baby_problems'
  
  metadata_16s[metadata_16s == 'PTB' | metadata_16s == 'Yes_I_plan_to_breast_feed' | 
                 metadata_16s == 'Yes_my_baby_was_put_in_the_ICU' |metadata_16s == 'Yes_I_have_been_diagnosed_with_hyperemesis' |
                 metadata_16s == 'C_Section' |metadata_16s == 'Cesarean_section' |  
                 metadata_16s == 'Female'| metadata_16s == 'Yes_98_100'| metadata_16s == 'Yes_101_103'| metadata_16s == 'Yes_membranes_ruptured_1_18_hours_before_labor'| metadata_16s == 'Yes_membranes_ruptured_before_37_weeks_of_pregnancy'] = 'Yes'
  
  metadata_16s[metadata_16s == 'TB' | metadata_16s == 'No_I_have_not_been_exposed_to_any_chemicals' | metadata_16s == 'No_I_have_not_been_to_the_hospital' | metadata_16s == 'No_I_have_not_been_diagnosed_with_hyperemesis_and_I_do_not_have_symptoms' |
                 metadata_16s == 'No_my_baby_was_not_put_in_the_ICU' | metadata_16s == 'No_I_do_not_plan_to_breast_feed' | metadata_16s == 'No_I_am_not_taking_any_new_medications' | metadata_16s == 'No_I_have_not_been_diagnosed_with_hyperemesis_but_I_have_symptoms' |
                 metadata_16s == 'Vaginal' | metadata_16s == 'Vaginal_delivery' |  metadata_16s == 'No_I_have_not_had_a_fever_since_my_last_visit' |
                 metadata_16s == 'Male' | metadata_16s == 'none'| metadata_16s == 'No_there_has_been_no_premature_rupture_of_membranes'] = 'No'
  
  metadata_16s[metadata_16s == 'Not_Sure' | 
                 metadata_16s == 'Not_sure' | 
                 metadata_16s == 'Uncertain' |
                 metadata_16s == 'Do_not_know'|
                 metadata_16s == 'Nelson 6th Floor'|
                 metadata_16s == 'Baby_NICU' | metadata_16s == 'NaN' | metadata_16s == 'Don_t_know'] = NA
  
  # write.csv(metadata_16s,'metadata_16s_ori.csv') ### run when output contains sample ID or participant ID
  
  # smoker_current
  metadata_16s[metadata_16s == 'Yes (trying to quit)'] = 'Yes'
  metadata_16s[metadata_16s == 'Yes (not trying to quit)'] = 'Yes'
  
  
  
  # smoker_second_hand
  metadata_16s[metadata_16s == 'Never'] = 0
  metadata_16s[metadata_16s == 'Rarely'] = 1
  metadata_16s[metadata_16s == 'Almost_every_day'] = 2
  metadata_16s[metadata_16s == 'Every_day'] = 3
  
  # various worries
  metadata_16s[metadata_16s == 'Much_Less'] = 1
  metadata_16s[metadata_16s == 'A_Bit_Less'] = 2
  metadata_16s[metadata_16s == 'About_Average'] = 3
  metadata_16s[metadata_16s == 'A_Bit_More'] = 4
  metadata_16s[metadata_16s == 'Much_More'] = 5
  
  metadata_16s[metadata_16s == 'Never'] = 0
  metadata_16s[metadata_16s == 'Almost_Never'] = 1
  metadata_16s[metadata_16s == 'Sometimes'] = 2
  metadata_16s[metadata_16s == 'Fairly_Often'] = 3
  metadata_16s[metadata_16s == 'Very_Often'] = 4
  
  # received_oral_frequency, douche_frequency
  metadata_16s$received_oral_frequency[metadata_16s$received_oral_frequency == 'No'] = 0
  metadata_16s$douche_frequency[metadata_16s$douche_frequency == 'No'] = 0
  
  # sexual_partners
  metadata_16s$sexual_partners_month[metadata_16s$sexual_partners_month == '3_to_5'] = 4
  metadata_16s$sexual_partners_month[metadata_16s$sexual_partners_month == '6_to_10'] = 8
  metadata_16s$sexual_partners_month[metadata_16s$sexual_partners_month == '11_to_20'] = 16
  
  # yogurt, milk, cheese, ice_cream
  metadata_16s[metadata_16s == '1_2_servings_per_day'] = 4
  metadata_16s[metadata_16s == '3_4_servings_per_day'] = 5
  metadata_16s[metadata_16s == '5_servings_per_day'] = 6
  metadata_16s[metadata_16s == 'Less_than_1_serving_per_week'] = 1
  metadata_16s[metadata_16s == '1_2_servings_per_week'] = 2
  metadata_16s[metadata_16s == '3_4_servings_per_week'] = 3
  
  # weight_change
  metadata_16s[metadata_16s == 'My_weight_has_NOT_changed_since_my_last_visit'] = 0
  metadata_16s[metadata_16s == 'I_have_LOST_weight_since_my_last_visit'] = -1
  metadata_16s[metadata_16s == 'I_have_GAINED_weight_since_my_last_visit'] = 1
  
  # weight_lost, weight_gained
  metadata_16s[metadata_16s == '0_5_lbs'] = 1
  metadata_16s[metadata_16s == '10_15_lbs'] = 3
  metadata_16s[metadata_16s == '5_10_lbs'] = 2
  metadata_16s[metadata_16s == 'More_than_15_lbs'] = 4
  
  # physical_activity_vigorous, physical_activity_moderate, physical_activity_light
  metadata_16s[metadata_16s == '0_times'] = 0
  metadata_16s[metadata_16s == '1_2_times'] = 1
  metadata_16s[metadata_16s == '3_4_times'] = 2
  metadata_16s[metadata_16s == '5_6_times'] = 3
  metadata_16s[metadata_16s == '7_times'] = 4
  
  # alcohol_frequency
  metadata_16s[metadata_16s == '3_5_times'] = 4
  metadata_16s[metadata_16s == '6_9_times'] = 7
  metadata_16s[metadata_16s == '10_19_times'] = 15
  metadata_16s[metadata_16s == '20_times'] = 20
  
  
  # uti_lifetime_number, 
  metadata_16s[metadata_16s == '2_to_4'] = 3
  metadata_16s[metadata_16s == '3_to_5'] = 4
  metadata_16s[metadata_16s == '6_to_10'] = 8
  metadata_16s[metadata_16s == '11_to_20'] = 15
  
  
  # douche_frequency, vaginal_sex_frequency, vaginal_penetration_frequency, received_oral_frequency, gave_oral_frequency, anal_sex_frequency, vaginal_lubrication_frequency
  metadata_16s[metadata_16s == '1_3_times_per_month'] = 2
  metadata_16s[metadata_16s == '2_6_times_per_week'] = 4
  metadata_16s[metadata_16s == 'Less_than_once_a_month'] = 1
  metadata_16s[metadata_16s == 'Never'] = 0
  metadata_16s[metadata_16s == 'Once_a_day'] = 5
  metadata_16s[metadata_16s == 'Once_a_week'] = 3
  
  # baby_time_in_hospital
  metadata_16s[metadata_16s == 'Less_than_24_hours'] = 1
  metadata_16s[metadata_16s == '24_to_48_hours'] = 2
  metadata_16s[metadata_16s == '3_days'] = 3
  metadata_16s[metadata_16s == '4_days'] = 4
  metadata_16s[metadata_16s == 'My_baby_is_still_in_the_hospital'] = 5
  metadata_16s[metadata_16s == 'My_baby_was_stillborn'] = NA
  
  # pregnancy_feelings
  metadata_16s[metadata_16s == 'It_was_one_of_the_happiest_times_of_my_life'] = 1
  metadata_16s[metadata_16s == 'It_was_generally_a_happy_time_with_few_exceptions'] = 2
  metadata_16s[metadata_16s == 'It_wasnt_any_more_or_less_difficult_than_before_I_was_pregnant'] = 3
  metadata_16s[metadata_16s == 'It_was_generally_a_lot_harder_than_before_I_was_pregnant'] = 4
  metadata_16s[metadata_16s == 'It_was_one_of_the_hardest_times_of_my_life'] = 5
  metadata_16s[metadata_16s == 'Other'] = NA
  
  # anal_sex_last_time
  metadata_16s[metadata_16s == 'More_than_1_year_ago_or_never'] = 1
  metadata_16s[metadata_16s == 'More_than_1_month_ago'] = 2
  metadata_16s[metadata_16s == '1_3_weeks_ago'] = 3
  metadata_16s[metadata_16s == '3_7_days_ago'] = 4
  metadata_16s[metadata_16s == '1_2_days_ago'] = 5
  
  # income
  metadata_16s[metadata_16s == '15_000_19_999'] = 2
  metadata_16s[metadata_16s == '20_000_39_999'] = 3
  metadata_16s[metadata_16s == '40_000_59_999'] = 4
  metadata_16s[metadata_16s == '60_000_79_999'] = 5
  metadata_16s[metadata_16s == '80_000_or_more'] = 6
  metadata_16s[metadata_16s == 'Under_15_000'] = 1
  
  # education
  metadata_16s[metadata_16s == '2_year_College_Degree'] = 3
  metadata_16s[metadata_16s == '4_year_College_Degree'] = 4
  metadata_16s[metadata_16s == 'Doctoral_or_Professional_Degree'] = 6
  metadata_16s[metadata_16s == 'High_School_GED'] = 2
  metadata_16s[metadata_16s == 'Less_than_High_School'] = 1
  metadata_16s[metadata_16s == 'Masters_Degree'] = 5
  metadata_16s[metadata_16s == 'Some_College'] = NA
  
  # prenatal_care_start
  metadata_16s[metadata_16s == '0_4_weeks_after_conception'] = 3
  metadata_16s[metadata_16s == '4_weeks_to_the_end_of_the_1st_trimester'] = 2
  metadata_16s[metadata_16s == '4_weeks_to_the_end_of_the_first_trimester'] = 2
  metadata_16s[metadata_16s == 'Before_conception'] = 4
  metadata_16s[metadata_16s == 'During_the_2nd_trimester'] = 1
  metadata_16s[metadata_16s == 'During_the_second_trimester'] = 1
  metadata_16s$SampleID_ori = NULL
  
  # gave_oral_last_time
  metadata_16s[metadata_16s == 'More_than_1_year_ago_or_never'] = 1
  metadata_16s[metadata_16s == 'More_than_1_month_ago'] = 2
  metadata_16s[metadata_16s == '1_3_weeks_ago'] = 3
  metadata_16s[metadata_16s == '3_7_days_ago'] = 4
  metadata_16s[metadata_16s == '1_2_days_ago'] = 5
  metadata_16s[metadata_16s == 'Within_the_last_24_hours'] = 6
  
  
  
}

# get all missing information when participant ID and visit time the same
{
  {
    metadata_16s$days_rel2birth[metadata_16s$KitType == 'Baby_Birth' & is.na(metadata_16s$days_rel2birth)] = 0
    
    metadata_instruction = read.csv('metadata_instruction copy 2.csv')
    metadata_16s = metadata_16s[,colnames(metadata_16s) %in% metadata_instruction$Metadata_subject_.Unless_otherwise_stated._listed_subjects_are_collected_from_mothers.]
    metadata_16s$ParticipantID = str_remove_all(metadata_16s$ParticipantID_Mombaby,'_.*')
    
    metadata_mom = metadata_16s[metadata_16s$Mombaby == 'Mom',]
    metadata_baby = metadata_16s[metadata_16s$Mombaby == 'Baby',]
    
    participant_list = unique(metadata_mom$ParticipantID)
    
    for (a in 1: length(participant_list)) {
      n = which(metadata_mom$ParticipantID == participant_list[a])
      metadata = metadata_mom[n,]
      for (b in 1: ncol(metadata)) {
        if (colnames(metadata)[b] %in% metadata_instruction$Metadata_subject_.Unless_otherwise_stated._listed_subjects_are_collected_from_mothers.[metadata_instruction$Same_in_all_visit == 'y'] & 
            length(unique(metadata[,b])[!is.na(unique(metadata[,b]))]) == 1) {
          metadata[,b] = unique(metadata[,b])[!is.na(unique(metadata[,b]))]
        }
      }
      metadata_mom[n,] = metadata
    }
    
    participant_list = unique(metadata_baby$ParticipantID)
    
    for (a in 1: length(participant_list)) {
      n = which(metadata_baby$ParticipantID == participant_list[a])
      metadata = metadata_baby[n,]
      for (b in 1: ncol(metadata)) {
        if (colnames(metadata)[b] %in% metadata_instruction$Metadata_subject_.Unless_otherwise_stated._listed_subjects_are_collected_from_mothers.[metadata_instruction$Same_in_all_visit == 'y'] & 
            length(unique(metadata[,b])[!is.na(unique(metadata[,b]))]) == 1) {
          metadata[,b] = unique(metadata[,b])[!is.na(unique(metadata[,b]))]
        }
      }
      metadata_baby[n,] = metadata
    }
    
  }
  
  metadata_16s = rbind(metadata_mom,metadata_baby)
  reads_table_16s = reads_table_16s[metadata_16s$SampleID]
  
  
  # get all the days_rel2birth
  metadata_mom = metadata_16s[metadata_16s$Mombaby == 'Mom',]
  metadata_baby = metadata_16s[metadata_16s$Mombaby == 'Baby',]
  
  participant_list = unique(metadata_mom$ParticipantID)
  
  for (a in 1: length(participant_list)) {
    n = which(metadata_mom$ParticipantID == participant_list[a])
    metadata = metadata_mom[n,]
    visit_time = unique(metadata$VisitNum)
    
    for (b in 1: length(visit_time)) {
      x = unique(metadata$days_rel2birth[metadata$VisitNum == visit_time[b]])
      x = x[!is.na(x)]
      if (length(x) !=0) {
        metadata$days_rel2birth[metadata$VisitNum == visit_time[b]] = x
      }
    }
    metadata_mom[n,] = metadata
  }
  
  participant_list = unique(metadata_baby$ParticipantID)
  
  for (a in 1: length(participant_list)) {
    n = which(metadata_baby$ParticipantID == participant_list[a])
    metadata = metadata_baby[n,]
    visit_time = unique(metadata$VisitNum)
    
    for (b in 1: length(visit_time)) {
      x = unique(metadata$days_rel2birth[metadata$VisitNum == visit_time[b]])
      x = x[!is.na(x)]
      if (length(x) !=0) {
        metadata$days_rel2birth[metadata$VisitNum == visit_time[b]] = x
      }
    }
    metadata_baby[n,] = metadata
  }
  
  metadata_16s = rbind(metadata_mom, metadata_baby)
  
  
  metadata_16s = metadata_16s[!is.na(metadata_16s$days_rel2birth),]
  reads_table_16s = reads_table_16s[,colnames(reads_table_16s) %in% metadata_16s$SampleID]
  reads_table_16s = reads_table_16s[metadata_16s$SampleID]
  
  # remove baby more than 2 day; only 22 samples 
  metadata_16s$days_rel2birth <- as.numeric(as.character(metadata_16s$days_rel2birth))
  keep = metadata_16s$days_rel2birth > 2 & !is.na(metadata_16s$days_rel2birth)
  sum(keep)
  keep = !keep
  reads_table_16s = reads_table_16s[,keep]
  metadata_16s = metadata_16s[keep,]
  
  # remove mom after delivery to make MB the same as MV and MR in sample collection time
  keep = metadata_16s$days_rel2birth >= 0 & !is.na(metadata_16s$days_rel2birth) &
    metadata_16s$Mombaby == 'Mom' 
  sum(keep)
  keep = !keep
  reads_table_16s = reads_table_16s[,keep]
  metadata_16s = metadata_16s[keep,]
  
  
  metadata_16s$Preterm[metadata_16s$ga_at_delivery < 259] = 'Yes'
  metadata_16s$Preterm[metadata_16s$ga_at_delivery >= 259] = 'No'
  
}

# remove outliers by grubbs.test and put ordinal data back to ordinal form
{
  metadata_num = as.data.frame(matrix(data = NA, nrow = nrow(metadata_16s), ncol =0))
  metadata_cha = as.data.frame(matrix(data = NA, nrow = nrow(metadata_16s), ncol =0))
  n =1
  m =1
  metadata_16s[metadata_16s == NaN] = NA
  metadata_16s[,1:ncol(metadata_16s)] <- lapply(metadata_16s[,1:ncol(metadata_16s)],as.character)
  
  c = c('0','1','2','3','4','5','6','7','8','9','.','-')
  
  for (a in 1:ncol(metadata_16s)) {
    b = metadata_16s[,a]
    b = b[!is.na(b)]
    b = strsplit(b,'*')
    b = unlist(b)
    b = unique(b)
    keep = b %in% c
    
    if (sum(!keep) == 0) {
      metadata_16s[,a] <- as.numeric(metadata_16s[,a])
      metadata_num <- cbind(metadata_num,metadata_16s[,a] )
      colnames(metadata_num)[n] = colnames(metadata_16s)[a]
      n=n+1
    } else {
      metadata_cha <- cbind(metadata_cha,metadata_16s[,a])
      colnames(metadata_cha)[m] = colnames(metadata_16s)[a]
      m=m+1
    }
  }
  
  for (a in 1:ncol(metadata_num)) {

    data = metadata_num[,a]
    
    if (length(unique(data[!is.na(data)])) < 2) {
      next
    }
    
    if (length(unique(data[!is.na(data)])) > 6) {
      # remove outliers
      {
        library(outliers)
          data_1 = grubbs.test(data)
          data_1 = data_1$statistic
          
          if (length(data_1) > 0) {
            data_2 = data
            data[data %in% data_1] = NA
          }

        metadata_num[,a] = data
      }
      next
    }
    
    data_1 = as.data.frame(table(data))
    data_1 = data_1[order(data_1$data),]
    
    data[!is.na(data) & data == data_1$data[1]] = 'A';data[!is.na(data) & data == data_1$data[2]] = 'B'
    data[!is.na(data) & data == data_1$data[3]] = 'C';data[!is.na(data) & data == data_1$data[4]] = 'D'
    data[!is.na(data) & data == data_1$data[5]] = 'E';data[!is.na(data) & data == data_1$data[6]] = 'F'

    metadata_num[,a] = data

  }

  metadata_16s = cbind(metadata_cha,metadata_num)
  
  rm(metadata,metadata_baby,metadata_mom,data, metadata_cha,metadata_num,data_1)
}

# modify taxa name
{
  row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Lachnospiraceae_G-9_bacterium_HMT_924')] = 'BVAB1'
  row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Sneathia_amnii_Not_Validly_Published')] = 'Sneathia_amnii'
  row.names(reads_table_16s)[which(row.names(reads_table_16s) == '1-68_')] = 'Tissierellaceae_1_68_'
  row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'WAL_1855D_')] = 'Tissierellaceae_WAL_1855D_'
  row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Rs-045_OTU_3')] = 'TM7_H1'
  row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'Rs-045_')] = 'TM7_H1_'
  row.names(reads_table_16s)[which(row.names(reads_table_16s) == 'ph2_')] = 'Tissierellaceae_ph2_'
  row.names(reads_table_16s)[which(row.names(reads_table_16s) == '258ds10_OTU_')] = 'Fibrobacteria_258ds10_OTU_'
  
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'Lachnospiraceae_G-9_bacterium_HMT_924')] = 'BVAB1'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'Lachnospiraceae_G-9_bacterium_HMT_924_BT')] = 'BVAB1_BT'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'Sneathia_amnii_Not_Validly_Published')] = 'Sneathia_amnii'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Sneathia_amnii')] = 'k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;g__Sneathia;s__amnii'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'Sneathia_amnii_Not_Validly_Published_BT')] = 'Sneathia_amnii_BT'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Sneathia_amnii_BT')] = 'k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;g__Sneathia;s__amnii'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == '1-68_OTU_12')] = 'Tissierellaceae_1_68_OTU_12'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_1_68_OTU_12')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'WAL_1855D_OTU_17')] = 'Tissierellaceae_WAL_1855D_OTU_17'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_WAL_1855D_OTU_17')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'WAL_1855D_OTU_17_BT')] = 'Tissierellaceae_WAL_1855D_OTU_17_BT'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_WAL_1855D_OTU_17_BT')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'Rs-045_OTU_3')] = 'TM7_H1'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'TM7_H1')] = 'k__Bacteria;p__TM7;c__TM7-3;o__;f__;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'ph2_OTU_4')] = 'Tissierellaceae_ph2_OTU_4'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_ph2_OTU_4')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == '1-68_OTU_17')] = 'Tissierellaceae_1_68_OTU_17'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_1_68_OTU_17')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'WAL_1855D_OTU_29')] = 'Tissierellaceae_WAL_1855D_OTU_29'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_WAL_1855D_OTU_29')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == '1-68_OTU_6')] = 'Tissierellaceae_1_68_OTU_6'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_1_68_OTU_6')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == '258ds10_OTU_4_BT')] = 'Fibrobacteria_258ds10_OTU_4_BT'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Fibrobacteria_258ds10_OTU_4_BT')] = 'k__Bacteria;p__Fibrobacteres;c__Fibrobacteria;o__;f__;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == '1-68_OTU_16')] = 'Tissierellaceae_1_68_OTU_16'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_1_68_OTU_16')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'WAL_1855D_OTU_19')] = 'Tissierellaceae_WAL_1855D_OTU_19'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_WAL_1855D_OTU_19')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'WAL_1855D_OTU_27')] = 'Tissierellaceae_WAL_1855D_OTU_27'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_WAL_1855D_OTU_27')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'WAL_1855D_OTU_8')] = 'Tissierellaceae_WAL_1855D_OTU_8'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_WAL_1855D_OTU_8')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == 'ph2_OTU_13')] = 'Tissierellaceae_ph2_OTU_13'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_ph2_OTU_13')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  taxonomy_16s$V3[which(taxonomy_16s$V3 == '1-68_OTU_18')] = 'Tissierellaceae_1_68_OTU_18'
  taxonomy_16s$V2[which(taxonomy_16s$V3 == 'Tissierellaceae_1_68_OTU_18')] = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Tissierellaceae;g__;s__'
  
}

# match mom and baby samples
{
  {
    # give flags
    setwd('/Users/binzhu/Desktop/Mom-baby/')
    
    # flag
    color1=c('#E43B2D','#0070FF','#4DAF4A','#FEF840','#F27E33','#F180BF')
    
    metadata_16s$Flag = NA
    metadata_16s$Flag = metadata_16s$SampleType
    metadata_16s$Flag[metadata_16s$SampleType == "BCKD" & metadata_16s$days_rel2birth == 2 & !is.na(metadata_16s$days_rel2birth)] = "BCKD_birth_2_days"
    metadata_16s$Flag[metadata_16s$SampleType == "BRCD" & metadata_16s$days_rel2birth == 2 & !is.na(metadata_16s$days_rel2birth)] = "BRCD_birth_2_days"
    metadata_16s$Flag[metadata_16s$SampleType == "BS1D" & metadata_16s$days_rel2birth == 2 & !is.na(metadata_16s$days_rel2birth)] = "BS1D_birth_2_days"
    
    metadata_16s$Flag[metadata_16s$SampleType == "BCKD" & metadata_16s$days_rel2birth == 1 & !is.na(metadata_16s$days_rel2birth)] = "BCKD_birth_1_days"
    metadata_16s$Flag[metadata_16s$SampleType == "BRCD" & metadata_16s$days_rel2birth == 1 & !is.na(metadata_16s$days_rel2birth)] = "BRCD_birth_1_days"
    metadata_16s$Flag[metadata_16s$SampleType == "BS1D" & metadata_16s$days_rel2birth == 1 & !is.na(metadata_16s$days_rel2birth)] = "BS1D_birth_1_days"
    
    metadata_16s$Flag[metadata_16s$SampleType == "BCKD" & metadata_16s$KitType == 'Baby_Birth' & metadata_16s$days_rel2birth == 0 & !is.na(metadata_16s$days_rel2birth)] = 'BCKD_birth_0_day'
    metadata_16s$Flag[metadata_16s$SampleType == "BRCD" & metadata_16s$KitType == 'Baby_Birth' & metadata_16s$days_rel2birth == 0 & !is.na(metadata_16s$days_rel2birth)] = 'BRCD_birth_0_day'
    metadata_16s$Flag[metadata_16s$SampleType == "BS1D" & metadata_16s$KitType == 'Baby_Birth' & metadata_16s$days_rel2birth == 0 & !is.na(metadata_16s$days_rel2birth)] = 'BS1D_birth_0_day'
    
    metadata_16s$Flag[metadata_16s$SampleType == "MV1D"] = 'MV'
    metadata_16s$Flag[metadata_16s$SampleType == "MRCD"] = 'MR'
    metadata_16s$Flag[metadata_16s$SampleType == "MCKD"] = 'MB'
    
    keep = !is.na(metadata_16s$Flag); sum(keep)
    reads_table_16s = reads_table_16s[,keep]
    metadata_16s = metadata_16s[keep,]
    
    participant_id_unique = unique(metadata_16s$ParticipantID)
    
    keep = colSums(reads_table_16s) >= 5000
    sum(keep)
    metadata_16s = metadata_16s[keep,]
    reads_table_16s = reads_table_16s[,keep]
    
    # remove mom samples that do not have paired baby
    baby_participant_list = unique(metadata_16s$ParticipantID[metadata_16s$Mombaby == 'Baby'])
    keep = metadata_16s$ParticipantID %in% baby_participant_list
    sum(keep)
    metadata_16s = metadata_16s[keep,]
    reads_table_16s = reads_table_16s[,keep]
    
    mom_participant_list = unique(metadata_16s$ParticipantID[metadata_16s$Mombaby == 'Mom'])
    keep = metadata_16s$ParticipantID %in% mom_participant_list
    sum(keep)
    metadata_16s = metadata_16s[keep,]
    reads_table_16s = reads_table_16s[,keep]
    
    
    
    
    # select maternal microbiomes at the last visit
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
    keep = rowSums(!is.na(output_table)) > 3 ; sum(keep)
    output_table = output_table[keep,]
    
    sample_list = output_table
    sample_list$participant = NULL
    sample_list = gather(output_table)
    sample_list = sample_list[sample_list$key != 'participant',]
    sample_list = sample_list[!is.na(sample_list$value),]
    
    reads_table_16s = reads_table_16s[,colnames(reads_table_16s) %in% sample_list$value]
    metadata_16s = metadata_16s[metadata_16s$SampleID %in% sample_list$value, ]
    
    reads_table_cytokines = reads_table_cytokines[,colnames(reads_table_cytokines) %in% sample_list$value]
    reads_table_lipid = reads_table_lipid[,colnames(reads_table_lipid) %in% sample_list$value]

    keep = sapply(1:ncol(metadata_16s), function(j) (length(unique(metadata_16s[,j])))); keep = keep >1; sum(keep)
    metadata_16s = metadata_16s[,keep]
  }

  # get taxa
  {
    # get main taxa list in 6 microbiomes
    type_all = unique(metadata_16s$SampleType)
    
    for (a in 1:length(type_all)) {
      keep = metadata_16s$SampleType == type_all[a]
      reads_table = reads_table_16s[,keep]
      
      # species and sample total reads thresholds and get taxonomy_16s
      keep <- colSums(reads_table) > 5000   # sample total reads >= 5000
      reads_table = reads_table[,keep]
      
      keep <- rowSums(reads_table) >= ncol(reads_table)   # species total reads >= sample number
      sum(keep)
      reads_table <- reads_table[keep,]
      reads_table = prepare_reads_table_2(reads_table, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8)
      
      if (type_all[a] == "BRCD") {
        reads_table_16s_BRCD = reads_table
      } else if (type_all[a] == "BS1D") {
        reads_table_16s_BS1D = reads_table
      } else if (type_all[a] == "BCKD") {
        reads_table_16s_BCKD = reads_table
      } else if (type_all[a] == "MCKD") {
        reads_table_16s_MCKD = reads_table
      } else if (type_all[a] == "MV1D") {
        reads_table_16s_MV1D = reads_table
      } else if (type_all[a] == "MRCD") {
        reads_table_16s_MRCD = reads_table
      } 
    }
    
    taxa_list = c(row.names(reads_table_16s_BRCD),row.names(reads_table_16s_BS1D),
                  row.names(reads_table_16s_BCKD),row.names(reads_table_16s_MCKD),
                  row.names(reads_table_16s_MV1D),row.names(reads_table_16s_MRCD))
    taxa_list = unique(taxa_list)
    
    # remove samples not passed threshold
    passed_sample_list = c(colnames(reads_table_16s_BRCD),colnames(reads_table_16s_BS1D),
                           colnames(reads_table_16s_BCKD),colnames(reads_table_16s_MCKD),
                           colnames(reads_table_16s_MV1D),colnames(reads_table_16s_MRCD))
    
    keep = metadata_16s$SampleID %in% passed_sample_list
    sum(keep)
    metadata_16s = metadata_16s[keep,]
    reads_table_16s = reads_table_16s[,keep]
    
    keep = metadata_16s$Mombaby == 'Mom' & metadata_16s$days_rel2birth >=0 & !is.na(metadata_16s$days_rel2birth); sum(keep)
    metadata_16s = metadata_16s[!keep,]
    reads_table_16s = reads_table_16s[,!keep]
    
    
    # remove samples that do not have paired samples
    baby_participant_list = unique(metadata_16s$ParticipantID[metadata_16s$Mombaby == 'Baby'])
    keep = metadata_16s$ParticipantID %in% baby_participant_list
    sum(keep)
    metadata_16s = metadata_16s[keep,]
    reads_table_16s = reads_table_16s[,keep]
    
    mom_participant_list = unique(metadata_16s$ParticipantID[metadata_16s$Mombaby == 'Mom'])
    keep = metadata_16s$ParticipantID %in% mom_participant_list
    sum(keep)
    metadata_16s = metadata_16s[keep,]
    reads_table_16s = reads_table_16s[,keep]
    
    keep = colnames(reads_table_16s_MCKD) %in% colnames(reads_table_16s)
    reads_table_16s_MCKD = reads_table_16s_MCKD[,keep]
    
    keep = colnames(reads_table_16s_MRCD) %in% colnames(reads_table_16s)
    reads_table_16s_MRCD = reads_table_16s_MRCD[,keep]
    
    keep = colnames(reads_table_16s_MV1D) %in% colnames(reads_table_16s)
    reads_table_16s_MV1D = reads_table_16s_MV1D[,keep]
    
    keep = colnames(reads_table_16s_BCKD) %in% colnames(reads_table_16s)
    reads_table_16s_BCKD = reads_table_16s_BCKD[,keep]
    
    keep = colnames(reads_table_16s_BRCD) %in% colnames(reads_table_16s)
    reads_table_16s_BRCD = reads_table_16s_BRCD[,keep]
    
    keep = colnames(reads_table_16s_BS1D) %in% colnames(reads_table_16s)
    reads_table_16s_BS1D = reads_table_16s_BS1D[,keep]
    
    # taxa threshold
    keep = row.names(reads_table_16s) %in% taxa_list
    sum(keep)
    reads_table_16s = reads_table_16s[keep,]
    
    # total reads threshold
    keep = colSums(reads_table_16s) >= 5000
    sum(keep)
    reads_table_16s = reads_table_16s[,keep]
    metadata_16s = metadata_16s[keep,]
    
    # get taxonomy_16s
    taxonomy_16s = taxonomy_16s[,c(3,2)]
    colnames(taxonomy_16s) = c('Taxa','taxonomy_16s')
    taxonomy_16s = taxonomy_16s[taxonomy_16s$Taxa %in% row.names(reads_table_16s),]
    keep = as.data.frame(rowSums(reads_table_16s)); keep$V2 = NA; keep = keep[order(keep$`rowSums(reads_table_16s)`, decreasing = T),]; keep = row.names(keep)
    reads_table_16s = reads_table_16s[match(keep, row.names(reads_table_16s)),]; row.names(reads_table_16s)
    taxonomy_16s = taxonomy_16s[match(keep,taxonomy_16s$Taxa),]
    rm(metadata_16s_cha,metadata_16s_num,metadata_16s_cha,metadata_16s_num, reads_table,taxonomy_16s_2, metadata,metadata_2)
  }
  
  # match again after removing some samples
  {
    # select maternal microbiomes at the last visit
    keep = !is.na(metadata_16s$Flag); sum(keep)
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
    
    output_table = output_table[output_table$participant !='P31350',]
    which(rowSums(!is.na(output_table))<4)
    
    sample_list = output_table
    sample_list$participant = NULL
    sample_list = gather(output_table)
    sample_list = sample_list[sample_list$key != 'participant',]
    sample_list = sample_list[!is.na(sample_list$value),]
    
    reads_table_16s = reads_table_16s[,colnames(reads_table_16s) %in% sample_list$value]
    metadata_16s = metadata_16s[metadata_16s$SampleID %in% sample_list$value, ]
    
    reads_table_cytokines = reads_table_cytokines[,colnames(reads_table_cytokines) %in% sample_list$value]
    reads_table_lipid = reads_table_lipid[,colnames(reads_table_lipid) %in% sample_list$value]
    
    keep = colnames(reads_table_16s_MCKD) %in% sample_list$value; sum(keep)
    reads_table_16s_MCKD = reads_table_16s_MCKD[,keep]
    
    keep = colnames(reads_table_16s_MRCD) %in% sample_list$value
    reads_table_16s_MRCD = reads_table_16s_MRCD[,keep]
    
    keep = colnames(reads_table_16s_MV1D) %in% sample_list$value
    reads_table_16s_MV1D = reads_table_16s_MV1D[,keep]
    
    keep = colnames(reads_table_16s_BCKD) %in% colnames(reads_table_16s)
    reads_table_16s_BCKD = reads_table_16s_BCKD[,keep]
    
    keep = colnames(reads_table_16s_BRCD) %in% colnames(reads_table_16s)
    reads_table_16s_BRCD = reads_table_16s_BRCD[,keep]
    
    keep = colnames(reads_table_16s_BS1D) %in% colnames(reads_table_16s)
    reads_table_16s_BS1D = reads_table_16s_BS1D[,keep]
  }
  
}

metadata_16s$`delivery_vaginal ` = NULL
colnames(metadata_16s) = str_remove_all(colnames(metadata_16s),' ')

sum(metadata_16s$Flag == "BCKD_birth_0_day" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BCKD_birth_1_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BCKD_birth_2_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BRCD_birth_0_day" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BRCD_birth_1_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BRCD_birth_2_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BS1D_birth_0_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BS1D_birth_1_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "BS1D_birth_2_days" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MB" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MR" & !is.na(metadata_16s$Flag))
sum(metadata_16s$Flag == "MV" & !is.na(metadata_16s$Flag))

rm(sample_list,metadata,metadata_2, reads_table, metadata_16s_ori)

# convert cytokine and lipids reads table to log10 values
{
  reads_table_lipid[reads_table_lipid ==0] = 0.000001
  reads_table_lipid = log10(reads_table_lipid)
#  row.names(reads_table_lipid) = paste0('Log10_',row.names(reads_table_lipid))
  
  reads_table_cytokines[reads_table_cytokines ==0] = 0.000001
  reads_table_cytokines = log10(reads_table_cytokines)
#  row.names(reads_table_cytokines) = paste0('Log10_',row.names(reads_table_cytokines))
}

rm(list = ls.str(mode = 'numeric'))
rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'logical'))





