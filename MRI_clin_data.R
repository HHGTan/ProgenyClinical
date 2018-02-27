#Load packages
library(plyr);library(dplyr)
library(tidyr)
library(reshape2);library(readxl)
library(survival);library(ggplot2)
#library(lme4); library(effects)


# Set functions
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# Set work directory
setwd("/Volumes/Promise_Pegasus/Harold/MRI_clinical_master/Data")
# Load data format
# Direct export from progeny --> Netherlands/ MRI_master|MRI_long|MRI_ecas
#raw_clin_wide <- read_excel("MRI_clin_20171222.xlsx",sheet=1)
#raw_alsfrs_long <- read_excel("MRI_ALSFRS_20171222.xlsx",sheet=1)
#raw_ecas_wide <- read_excel("MRI_ECAS_20171222.xlsx",sheet=1)

#raw_clin_wide <- read_excel("PAN_clin_20180112.xlsx",sheet=1)
#raw_alsfrs_long <- read_excel("PAN_ALSFRS_20180112.xlsx",sheet=1)
#raw_ecas_wide <- read_excel("PAN_ECAS_20180112.xlsx",sheet=1)

raw_clin_wide <- read_excel("PANMR_clin_20180206.xlsx",sheet=1)
raw_alsfrs_long <- read_excel("PANMR_ALSFRS_20180208.xlsx",sheet=1)
raw_ecas_wide <- read_excel("PANMR_ECAS_20180206.xlsx",sheet=1)



# ECAS and ALSFRS should be less than x numbers away of scan.
cutoff <- 90


##################
## Clean data.. ##
##################
# To be cleaned:
# A few missing fields for MRI subjecttype
# ALSFRS 1111 values
# 

raw_clin_wide$``Deelname MRI 3T``

#######################
## RAW CLINICAL DATA ##
#######################

# Superquickanddirty: Als hij geen mrisubjecttype heeft haalt hij het uit NL pat/cont. 
# Om beter FCO's er uit te halen moet ook FALS velden worden meegenomen
# tzt ook alle mimics uitsplitsen.
raw_clin_wide$All_Diagnosis <-sapply(ifelse(!is.na(raw_clin_wide$Diagnosis),raw_clin_wide$Diagnosis,
                                     ifelse(!is.na(raw_clin_wide$`MRI 3T: Subjecttype`), raw_clin_wide$`MRI 3T: Subjecttype`,
                                            ifelse(raw_clin_wide$`NL - Patient/Control` == "Controle", "Controle",NA))),simpleCap)


# Calculate ages at each event
raw_clin_wide$Age_onset <- as.numeric((as.Date(raw_clin_wide$`Date of onset`) - as.Date(raw_clin_wide$`Date of birth`))/365.25)
raw_clin_wide$Age_diagnosis <- as.numeric((as.Date(raw_clin_wide$`Date of diagnosis`) - as.Date(raw_clin_wide$`Date of birth`))/365.25)
raw_clin_wide$Diagnostic_delay <- as.numeric((as.Date(raw_clin_wide$`Date of diagnosis`) - as.Date(raw_clin_wide$`Date of onset`))/365.25*12)

raw_clin_wide$Dead <- as.numeric(!is.na(raw_clin_wide$`Date of Death`))
raw_clin_wide$Surv_temp <- ifelse(!is.na(raw_clin_wide$`Date of Death`), raw_clin_wide$`Date of Death`,raw_clin_wide$`AFLOOP: Datum check` )
raw_clin_wide$Surv_onset <- as.numeric(as.Date(raw_clin_wide$Surv_temp) - as.Date(raw_clin_wide$`Date of onset`))/365.25*12


# Explore survival
temp4 <- temp6 <- NULL
for (i in c("PLS", "ALS", "PMA", "Other, Namely:")) {
  df1 <- filter(raw_clin_wide, All_Diagnosis==i)
  temp1 <- coxph(Surv(Surv_onset, Dead==1) ~ 1, df1)
  temp2 <- survfit(temp1)
  temp3 <-  cbind.data.frame(survival = temp2$surv, time = temp2$time)
  temp3$Diagnosis <- i
  temp4 <- rbind(temp4,temp3)
  
  temp5 <- c(i, unlist(quantile(temp2, probs = c(0.25, 0.50, 0.75))$quantile))
  temp6 <- rbind(temp6,temp5)
}

ggplot(temp4, aes(time, survival, col=Diagnosis)) +geom_point(shape=3,size=2) + geom_path(size=0.5)+
  ggtitle("Kaplan-Meier Survival Curve")+ theme(plot.title=element_text(size=20, face="bold"))+
  theme(axis.title.x=element_text(size=18), axis.title.y=element_text(size=18))+
  scale_color_discrete(name="P") + theme(legend.title=element_text(size=16))


##############
## MRI DATA ##
##############

# Create new Long dataframe per MRI scan.
mri_wide <- raw_clin_wide[,c("ALS number","FollowUp1_datum_scan", "FollowUp2_datum_scan", "FollowUp3_datum_scan", 
               "FollowUp4_datum_scan", "FollowUp5_datum_scan", "FollowUp6_datum_scan")]

mri_long <- melt(mri_wide, id.vars = "ALS number", variable.name = "scan_fu", value.name = "scan_datum")
mri_long$scan_datum <- as.Date(mri_long$scan_datum,format="%Y-%m-%d")
mri_long$scan_fu <- gsub("FollowUp|_datum_scan", "", mri_long$scan_fu)


###############
## ECAS DATA ##
###############

# Reformat ECAS column names
raw_ecas_wide <- raw_ecas_wide  # As for now (01/12/2017, FU6 data is not correctly available in progeny)
colnames(raw_ecas_wide)[grep("FU|ALS number",invert = T,colnames(raw_ecas_wide))] <- paste0("FU1: ",colnames(raw_ecas_wide)[grep("FU|ALS number",invert = T,colnames(raw_ecas_wide))]  )
colnames(raw_ecas_wide) <- gsub("ECAS 1", "ECAS", colnames(raw_ecas_wide))
colnames(raw_ecas_wide) <- gsub("(Datum.*) (FU.*)","\\2: \\1",colnames(raw_ecas_wide))
colnames(raw_ecas_wide) <- gsub("-|_"," ", colnames(raw_ecas_wide))
colnames(raw_ecas_wide)[-1] <- sapply(colnames(raw_ecas_wide)[-1], simpleCap)
colnames(raw_ecas_wide) <- gsub("Interia","Inertia",colnames(raw_ecas_wide))

# Reformat from wide to long format using tidyr
ecas_long <-  raw_ecas_wide %>% 
  gather(v, value, grep("FU[1-9]",colnames(raw_ecas_wide)))%>% 
  separate(v, c("ECAS_FU", "col"), sep=": ") %>% 
  spread(col, value)

ecas_long$`Datum ECAS` <- as.Date(ecas_long$`Datum ECAS`,format = "%Y-%m-%d")
domains_ecas <- c("ECAS ALS NONSPECIFIC", "ECAS ALS SPECIFIC", "ECAS Total Score", 
                  "Executive", "Fluency Free + Fixed", "Language", "Memory", "Visuospatial")
ecas_long[domains_ecas] <- lapply(ecas_long[domains_ecas],as.numeric)

for( i in 1:nrow(ecas_long)) {
  if ( is.na(ecas_long[i,"Symptoms"]) & is.na(ecas_long[i,"Total Behaviour"]) ) {
    ecas_long[i,c("Behaviour Disinhibition", "Apathy Inertia", "Hyperorality Altered Food", 
                  "Loss Of Sympathy", "Perseverative Stereotyped") ] <- NA
  } else {
    for (j in c("Behaviour Disinhibition", "Apathy Inertia", "Hyperorality Altered Food",
                "Loss Of Sympathy", "Perseverative Stereotyped")){
      if(is.na(ecas_long[i,j])){
        ecas_long[i,j] <- "N"
      } 
      if(ecas_long[i,j] =="Yes"){
        ecas_long[i,j] <- "Y"
      }
    }
  }
}
ecas_bhv <- ecas_long[,c("ALS number","Behaviour Disinhibition", "Apathy Inertia", "Hyperorality Altered Food", 
                         "Loss Of Sympathy", "Perseverative Stereotyped", "Symptoms", "Total Behaviour")]




# Kies de bijpassende ECAS per subject, per scan.
mri_ecas_match <-ldply(unique(mri_long$`ALS number`[which(!is.na(mri_long$`ALS number`))]),function(i){
  #subj_mri_dates <- filter(mri_long, `ALS number`== i) #MRI datums
  #subj_ecas_dates <- filter(ecas_long, `ALS number`== i) #ECAS datums
  
  #workaround omdat dplyr niet werkt.
  subj_mri_dates <- mri_long[which(mri_long$`ALS number`==i),]  #MRI datums
  subj_ecas_dates <- ecas_long[which(ecas_long$`ALS number`==i),] #ECAS datums
  
  
  if (length(grep(i, ecas_long$`ALS number`)) ==0   ){
    subj_ecas_dates[1,] <- NA
  }
  
  #bekijk per ECAS het verschil in dagen met alle MRI scans 
  #outer(subj_ecas_dates$`Datum ECAS`,subj_mri_dates$scan_datum, "-")
  mat1<- t(sapply(subj_ecas_dates$`Datum ECAS`, function(j){
    return(j - subj_mri_dates$scan_datum)
  }))
  
  rownames(mat1) <- subj_ecas_dates$ECAS_FU
  
  #Als het verschil meer dan de cutoff is(in dagen), wordt de ECAS niet gebruikt.
  mat2 <- ifelse( abs(mat1) > cutoff,NA, abs(mat1))
  mat3 <- apply(mat2, 1, function(x) {if (all(is.na(x))) {NA}  else {which.min(x)} }) 

  if (length(which(!is.na(unique(mat3))))!=length(which(!is.na(mat3)))){
    dup_fu <- names(table(mat3))[which(table(mat3)>1)]
    for (m in dup_fu) {
      fu_scandif <- mat2[,as.integer(m)]
    rm_duplicate <- ifelse ( (min(fu_scandif,na.rm = T) != fu_scandif) & !is.na(fu_scandif),NA,1)
    mat3 <- rm_duplicate * mat3
    }
  }
  
  #werkt niet voor NAs
  x3<- subj_mri_dates 
  x3$ECAS_FU <- sapply(x3$scan_fu, function (k){
    if (length(which(k==mat3))==0) {
      out_FU <- NA
    } else {
      out_FU <- names(mat3)[which(k==mat3)]
    }
    return(out_FU)
  })
  
  
  if (length(which(!is.na(unique(mat3))))!=length(which(!is.na(mat3)))){
    x3$ECAS_FU <- rep(i,6)
  }
  
  return(x3)
})

data_ecas <- join_all(list(mri_ecas_match,ecas_long),by = c("ALS number", "ECAS_FU"), type = "left", match = "all")
data_ecas$ECAS_FU <- gsub("FU|: ", "", data_ecas$ECAS_FU)




#######
#Next step: ALSFRS

##### VOEG TOE: if empty, dan geen datum gebruiken. 
#If complete dan is het normaal dat total ontbreekt. (optioneel)



for(i in 1:length(raw_alsfrs_long$`ALS number`)){
  if (is.na(raw_alsfrs_long$`ALS number`[i])) {
    raw_alsfrs_long$`ALS number`[i] <- raw_alsfrs_long$`ALS number`[(i-1)]
  } 
}
ALSFRS_col <- grep("ALS-FRS-R.*",colnames(raw_alsfrs_long), value = T)
ALSFRS_subf <- grep("ALS-FRS-R.[0-9].*",colnames(raw_alsfrs_long), value = T)

raw_alsfrs_long$checkAF <- apply(raw_alsfrs_long[,ALSFRS_subf],1, function(i){
  fields <-  sum(!is.na(i))
  logical <- ifelse(fields == 12 | fields == 0 , 0 , 1)
  return(logical)
})
raw_alsfrs_long$checkAF2  <- ifelse(raw_alsfrs_long$checkAF == 1 & !is.na(raw_alsfrs_long$`Table ALS-FRS-R.Total ALS-FRS-R SCORE:`),1,0)


alsfrs_long <- as.data.frame(raw_alsfrs_long[,c("ALS number",ALSFRS_col)])
alsfrs_long$`Table ALS-FRS-R.Date of ALS-FRS-R` <- as.Date(alsfrs_long$`Table ALS-FRS-R.Date of ALS-FRS-R`, format = "%Y-%m-%d")
subscores_alfrs <- grep("ALS-FRS-R.[0-9T]", colnames(alsfrs_long),value = T)
alsfrs_long[subscores_alfrs] <- lapply(alsfrs_long[subscores_alfrs],as.numeric)


mri_alsfrs_match <-ldply(unique(mri_long$`ALS number`[which(!is.na(mri_long$`ALS number`))]),function(i){
  #subj_mri_dates <- filter(mri_long, `ALS number`== i) #MRI datums
  #subj_alsfrs_dates <- filter(alsfrs_long, `ALS number`== i) #ALSFRS datums

  #workaround
  subj_mri_dates <- mri_long[which(mri_long$`ALS number`==i),]  #MRI datums
  subj_alsfrs_dates <- alsfrs_long[which(alsfrs_long$`ALS number`==i),] #ALSFRS datums
  
  
    
  if (length(grep(i, alsfrs_long$`ALS number`)) ==0   ){
    subj_alsfrs_dates[1,] <- NA
  }
  
  #bekijk per ALSFRS het verschil in dagen met alle MRI scans 
  mat1<- t(sapply(subj_alsfrs_dates$`Table ALS-FRS-R.Date of ALS-FRS-R`, function(j){
    return(j - subj_mri_dates$scan_datum)
  }))
  
  rownames(mat1) <- as.character(subj_alsfrs_dates$`Table ALS-FRS-R.Date of ALS-FRS-R`)
  
  #Als het verschil meer dan cutoff dagen is, wordt de ALSFRS niet gebruikt.
  mat2 <- ifelse( abs(mat1) > cutoff,NA, abs(mat1))
  mat3 <- apply(mat2, 1, function(x) {if (all(is.na(x))) {NA}  else {which.min(x)} }) 
  
  if (length(which(!is.na(unique(mat3))))!=length(which(!is.na(mat3)))){
    dup_fu <- names(table(mat3))[which(table(mat3)>1)]
    for (m in dup_fu) {
      fu_scandif <- mat2[,as.integer(m)]
      rm_duplicate <- ifelse ( (min(fu_scandif,na.rm = T) != fu_scandif) & !is.na(fu_scandif),NA,1)
      mat3 <- rm_duplicate * mat3
    }
  }
  
  x3<- subj_mri_dates 
  x3$`Table ALS-FRS-R.Date of ALS-FRS-R` <- sapply(x3$scan_fu, function (k){
    if (length(which(k==mat3))==0) {
      out_FU <- NA
    } else {
      out_FU <- names(mat3)[which(k==mat3)]
    }
    return(out_FU)
  })
  
  
  if (length(which(!is.na(unique(mat3))))!=length(which(!is.na(mat3)))){
    x3$`Table ALS-FRS-R.Date of ALS-FRS-R` <- rep(i,6)
  }
  
  return(x3)
})
mri_alsfrs_match$`Table ALS-FRS-R.Date of ALS-FRS-R` <- as.character(mri_alsfrs_match$`Table ALS-FRS-R.Date of ALS-FRS-R`, format="%Y-%m-%d") # Yes it's weird, in the next bit we'll join_all on character and Date classes, yet somehow this seems to work.
data_alsfrs <- join_all(list(mri_alsfrs_match,alsfrs_long), by=c("ALS number","Table ALS-FRS-R.Date of ALS-FRS-R"), type="left", match="first")


###############
## ALS-FTD-Q ##
###############

AFQ_col <- grep("ALS-FTD-Q.*",colnames(raw_alsfrs_long), value = T)
AFQ1 <- raw_alsfrs_long[,c("ALS number",AFQ_col[1:2])]
AFQ2 <-raw_alsfrs_long[,c("ALS number",AFQ_col[3:4])]
colnames(AFQ2)<- colnames(AFQ1)
AFQall <-rbind(AFQ1,AFQ2)
AFQall <- filter(AFQall, !is.na(`Datum van invullen ALS-FTD-Q`) | !is.na(`Totaal score ALS-FTD-Q`))

match2scan <- function(df, id, measure){
  ids <- unique(df[,id])
  
}

mri_alsfrs_match <-ldply(unique(mri_long$`ALS number`[which(!is.na(mri_long$`ALS number`))]),function(i){
  #subj_mri_dates <- filter(mri_long, `ALS number`== i) #MRI datums
  #subj_alsfrs_dates <- filter(alsfrs_long, `ALS number`== i) #ALSFRS datums
  
  #workaround
  subj_mri_dates <- mri_long[which(mri_long$`ALS number`==i),]  #MRI datums
  subj_alsfrs_dates <- alsfrs_long[which(alsfrs_long$`ALS number`==i),] #ALSFRS datums
  
  
  
  if (length(grep(i, alsfrs_long$`ALS number`)) ==0   ){
    subj_alsfrs_dates[1,] <- NA
  }
  
  #bekijk per ALSFRS het verschil in dagen met alle MRI scans 
  mat1<- t(sapply(subj_alsfrs_dates$`Table ALS-FRS-R.Date of ALS-FRS-R`, function(j){
    return(j - subj_mri_dates$scan_datum)
  }))
  
  rownames(mat1) <- as.character(subj_alsfrs_dates$`Table ALS-FRS-R.Date of ALS-FRS-R`)
  
  #Als het verschil meer dan cutoff dagen is, wordt de ALSFRS niet gebruikt.
  mat2 <- ifelse( abs(mat1) > cutoff,NA, abs(mat1))
  mat3 <- apply(mat2, 1, function(x) {if (all(is.na(x))) {NA}  else {which.min(x)} }) 
  
  if (length(which(!is.na(unique(mat3))))!=length(which(!is.na(mat3)))){
    dup_fu <- names(table(mat3))[which(table(mat3)>1)]
    for (m in dup_fu) {
      fu_scandif <- mat2[,as.integer(m)]
      rm_duplicate <- ifelse ( (min(fu_scandif,na.rm = T) != fu_scandif) & !is.na(fu_scandif),NA,1)
      mat3 <- rm_duplicate * mat3
    }
  }
  
  x3<- subj_mri_dates 
  x3$`Table ALS-FRS-R.Date of ALS-FRS-R` <- sapply(x3$scan_fu, function (k){
    if (length(which(k==mat3))==0) {
      out_FU <- NA
    } else {
      out_FU <- names(mat3)[which(k==mat3)]
    }
    return(out_FU)
  })
  
  
  if (length(which(!is.na(unique(mat3))))!=length(which(!is.na(mat3)))){
    x3$`Table ALS-FRS-R.Date of ALS-FRS-R` <- rep(i,6)
  }
  
  return(x3)
})
mri_alsfrs_match$`Table ALS-FRS-R.Date of ALS-FRS-R` <- as.character(mri_alsfrs_match$`Table ALS-FRS-R.Date of ALS-FRS-R`, format="%Y-%m-%d") # Yes it's weird, in the next bit we'll join_all on character and Date classes, yet somehow this seems to work.
data_alsfrs <- join_all(list(mri_alsfrs_match,alsfrs_long), by=c("ALS number","Table ALS-FRS-R.Date of ALS-FRS-R"), type="left", match="first")


##############
## All data ##
##############

##Merge data together
data_long <- join_all(list(data_ecas,data_alsfrs), by=c("ALS number", "scan_fu"), type="left", match = "all")

data_full <- join_all(list(data_long,raw_clin_wide[,grep("FollowUp",colnames(raw_clin_wide),invert=T)]),by = "ALS number", type="left", match = "all")
data_full$Age_MRI <- as.numeric((data_full$scan_datum - as.Date(data_full$`Date of birth`))/365.25)
data_full$Disease_duration <- as.numeric(data_full$scan_datum - as.Date(data_full$`Date of onset`))/365.25 * 12
data_full$FTD <- ifelse(data_full$`ALS plus`=="FTD","Yes","No")

data_select <- data_full[,
c("ALS number", "scan_fu", "scan_datum", "ECAS_FU", "Datum ECAS", "Table ALS-FRS-R.Date of ALS-FRS-R",
  "Age_MRI", "Disease_duration", "Sex (gender at birth)__1",
  "All_Diagnosis", "Diagnose 4-12 en 0", "Other diagnosis", "Vorige diagnose", "Vorige diagnose anders",
  "Age_onset", "Age_diagnosis", "Diagnostic_delay", "Surv_onset", "Dead", 
  "Site of Onset", "Side spinal onset",
  "REG_Bulb_CMN", "REG_Cerv_CMN", "REG_Lumb_CMN", "REG_Thor_CMN", 
  "REG_ Bulb_PMN", "REG_Cerv_PMN", "REG_Lumb_PMN", "REG_Thor_PMN", 
  "ALS plus", "ALS plus anders", "FTD", "Type: Familiair of Sporadisch", 
  
  
  "Table ALS-FRS-R.Total ALS-FRS-R SCORE:", 
  "Table ALS-FRS-R.1. Speech", "Table ALS-FRS-R.2. Salivation", 
  "Table ALS-FRS-R.3. Swallowing", "Table ALS-FRS-R.4. Handwriting", 
  "Table ALS-FRS-R.5a. Cutting Food and Handling Utensils", "Table ALS-FRS-R.5b. Cutting Food and Handling Utensils", 
  "Table ALS-FRS-R.6. Dressing and Hygiene", "Table ALS-FRS-R.7. Turning in Bed and Adjusting Bed Clothes", 
  "Table ALS-FRS-R.8. Walking", "Table ALS-FRS-R.9. Climbing Stairs", 
  "Table ALS-FRS-R.10. Dyspnea", "Table ALS-FRS-R.11. Orthopnea", 
  "Table ALS-FRS-R.12. Respiratory Insufficiency", "Table ALS-FRS-R.Notes", 
  "Table ALS-FRS-R.Studie", "015 MRI 3T", "Deelname MRI 3T", "NL - Patient/Control", 
  
  "ECAS Total Score", "ECAS ALS SPECIFIC", "ECAS ALS NONSPECIFIC", 
  "Language", "Fluency Free + Fixed", "Executive", "Memory", "Visuospatial",
  "ALS number","Behaviour Disinhibition", "Apathy Inertia", "Hyperorality Altered Food", 
  "Loss Of Sympathy", "Perseverative Stereotyped", "Symptoms", "Total Behaviour","Psychosis"
)]







# TO DO:
# - Spreidingslijst (MRI_long)
# - Add VFI? FAB? (MRI_long)
# 



# Misschien hier stoppen, volgende gedeelte in persoonlijke datasets maken.
#


