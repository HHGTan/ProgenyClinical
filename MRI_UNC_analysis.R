# UNC analysis, building on MRI_clin_data.R
# This part is study specific.

# Load libraries
library(plyr);library(dplyr)
library(tidyr)
library(reshape2);library(readxl)
library(survival);library(ggplot2)
library(lme4); library(effects)
library(tableone)

# Load clinical data from MRI_clin_data.R
#data_select <- ....
# Change a few colnames
colnames(data_select)[colnames(data_select) == "ALS number"] <- "ALSnr"
colnames(data_select)[colnames(data_select) == "Sex (gender at birth)__1"] <- "Gender"

data_select$Gender <- as.factor(data_select$Gender)
data_select$All_Diagnosis <- as.factor(data_select$All_Diagnosis)


# Load genotype data
setwd("/Volumes/Promise_Pegasus/Harold/MRI_clinical_master/Genotypes")

# Instead of below code, read Ruben's excel file
unc_raw_data <- read_excel(path="UNC13A_RvE_100118.xlsx",sheet = 1)
#unc_data <- data.frame(ALSnr=unc_raw_data$ID,genotype=unc_raw_data$UNC)
#unc_data$riskvariant <- as.factor(ifelse(unc_data$genotype == "G/G", "1", "0"))
#unc_data$dose <- ifelse(unc_data$genotype=="G/G",2, ifelse(unc_data$genotype=="G/T",1,0))
#unc_data$UNC13A <- ifelse(is.na(unc_data$genotype),"NA", as.character(unc_data$genotype))
#unc_data$riskvariant <- relevel(factor(unc_data$riskvariant), ref="0")
#unc_data$genotype <- relevel(factor(unc_data$genotype), ref="T/T")
unc_raw_data <- data.frame(ALSnr=unc_raw_data$ID,genotype=unc_raw_data$UNC)
unc_raw_data$genotype <- gsub("T","A", gsub("G","C",unc_raw_data$genotype))


# OLD UNC data from text file
unc_data <- read.table("20180131_UNC13_v3.txt", header = T, sep = ",")
test_unc <- join_all(list(unc_data,unc_raw_data), by="ALSnr",type="full", match="all")

unc_data <- test_unc
unc_data$riskvariant <- as.factor(ifelse(unc_data$genotype == "C/C", "1", "0"))
unc_data$dose <- ifelse(unc_data$genotype=="C/C",2, ifelse(unc_data$genotype=="C/A",1,0))
unc_data$UNC13A <- ifelse(is.na(unc_data$genotype),"NA", as.character(unc_data$genotype))
unc_data$riskvariant <- relevel(factor(unc_data$riskvariant), ref="0")
unc_data$genotype <- relevel(factor(unc_data$genotype), ref="A/A")

c9_data <- read.table("C9status_20171222.txt",header = T,sep = "\t",na.strings = c("",NA))
colnames(c9_data)[1] <- "ALSnr"


# Load imaging data
setwd("/Volumes/Promise_Pegasus/Harold/MRI")

T1_1 <- read.table("T1_NL1_HR/NL1_T1_aparc_cross_ct.txt",header=T,sep="")
T1_2 <- read.table("T1_NL1_HR/NL1_T1_BA.thresh_cross_ct.txt",header=T,sep="")
T1_3 <- read.table("T1_NL1_HR/NL1_T1_subcort_cross_volume.txt",header=T,sep="")
T1_all <- join_all(list(T1_1,T1_2,T1_3), type="left", by="ALSnr", match="first")
colnames(T1_all)[1] <- "ALSnr_MRI"
T1_all$th_asymmetrie <- ((T1_all$lh_MeanThickness_thickness - T1_all$rh_MeanThickness_thickness)*2 /
                           ( T1_all$lh_MeanThickness_thickness + T1_all$rh_MeanThickness_thickness))

DTI_1 <- read.table("DTI_NL1/NL1_trac_cvs_stats.txt",header = T)
DTI_2 <- read.table("DTI_NL1/NL1_Ptx_stats.txt",header=T)
DTI_select1 <- DTI_1[,c("ALSnr", grep("(FA|RD)_Avg_Weight", colnames(DTI_1),value = T))]
DTI_select2 <- DTI_2[,c(grep("thickness", colnames(DTI_2),value = T,invert = T))]
DTI_all <- join_all(list(DTI_select1,DTI_select2), by="ALSnr", type="left",match="all")
colnames(DTI_all)[1] <- "ALSnr_MRI"



# Combine all datasets (Clinical, Genetic, Imaging)
data_clingen <- join_all(list(data_select, unc_data,c9_data), by="ALSnr", type="left", match="all")
data_clingen$ALSnr_MRI <- paste(data_clingen$ALSnr,data_clingen$scan_fu,sep = "_")

dmaster <- join_all(list(data_clingen,T1_all,DTI_all), by="ALSnr_MRI",type="left", match="all")


#Balint <- read.table("/Users/htan4/Documents/C9_MET_pipeline/Balints_PLSers.txt", header=T)
######################################################################################################
# UNC13A analysis

# Select study group

# Cross sectional (x = cross)
dALSx <- filter(dmaster, (All_Diagnosis=="ALS") & scan_fu==1 )
dALSx <- filter(dALSx, !is.na(lh_bankssts_thickness))
dALSx <- filter(dALSx, !is.na(genotype))
dALSx <- filter(dALSx, is.na(C9orf72.mutated)|C9orf72.mutated=="No")


# Baseline characteristiscs
factorVars <- c("riskvariant")
vars <- c("Gender", "Age_MRI", "Disease_duration", "Age_onset", "Age_diagnosis",  "Surv_onset",
          "Site_of_Onset", "Type_Familiair_of_Sporadisch", "FTD",
          "C9orf72.mutated",  "UNC13A", "Table_ALS_FRS_R.Total_ALS_FRS_R_SCORE", "ALSFRS_slope",
          "ECAS_Total_Score", "ECAS_totcut","ECAS_ALS_SPECIFIC", "ECAS_speccut","ECAS_ALS_NONSPECIFIC","ECAS_nonscut")
varsCont <- c("Age_MRI", "Age_onset", "Age_diagnosis", "Disease_duration", "Surv_onset",  "Table_ALS_FRS_R.Total_ALS_FRS_R_SCORE","ALSFRS_slope", 
              "ECAS_Total_Score", "ECAS_ALS_SPECIFIC", "ECAS_ALS_NONSPECIFIC")


dtab <- dALSx
colnames(dtab) <- gsub(":","", gsub("-| ","_",colnames(dtab)))
dtab$ALSFRS_slope <- -(48-dtab$Table_ALS_FRS_R.Total_ALS_FRS_R_SCORE)/dtab$Disease_duration
dtab$ECAS_totcut <- ifelse(dtab$ECAS_Total_Score < 96,"1","0")
dtab$ECAS_speccut <- ifelse(dtab$ECAS_ALS_SPECIFIC < 69,"1","0")
dtab$ECAS_nonscut <- ifelse(dtab$ECAS_ALS_NONSPECIFIC < 24,"1","0")

tableOne <- CreateTableOne(vars = vars, strata = c("riskvariant"), data = dtab, factorVars = factorVars)
print(tableOne, nonnormal= varsCont, minMax=T, cramVars = c("Sex", "C9orf72.mutated"), missings = TRUE)



# Survival
temp4 <- temp6 <- NULL
for (i in c("0", "1")) {
  df1 <- filter(dtab, riskvariant==i)
  temp1 <- coxph(Surv(Surv_onset, Dead==1) ~ Gender + Age_onset + Site_of_Onset, df1)
  temp2 <- survfit(temp1)
  temp3 <-  cbind.data.frame(survival = temp2$surv, time = temp2$time/12)
  if (i=="0"){
    temp3$UNC13A <- "AA + AC"
  } else {
    temp3$UNC13A <- "CC"
  }
  
  temp4 <- rbind(temp4,temp3)
  
  temp5 <- c(i, unlist(quantile(temp2, probs = c(0.25, 0.50, 0.75))$quantile))
  temp6 <- rbind(temp6,temp5)
}

ggplot(temp4, aes(time, survival, col=UNC13A)) +geom_point(shape=3,size=2) + geom_path(size=0.5)+
  #  ggtitle("Kaplan-Meier Survival Curve")+ theme(plot.title=element_text(size=20, face="bold"))+
  theme(axis.title.x=element_text(size=18), axis.title.y=element_text(size=18))+
  scale_color_discrete(name="UNC13A genotype") + theme(legend.title=element_text(size=16)) +xlim(0,10) +
  ylab("Survival (%)") + xlab("Duration after disease onset (yr)")



######################################################################################################
# tableone


dtest <- filter(d9, scan_fu==1)
dtest$UNC13A <- ifelse(is.na(dtest$genotype),"NA", as.character(dtest$genotype))
dALS <- filter(dtest$All_Diagnosis=="ALS")
dtab <- filter(dtest, All_Diagnosis =="ALS"| All_Diagnosis=="Controle")

colnames(dtab) <- gsub(":","", gsub("-| ","_",colnames(dtab)))
tab$FTD<- ifelse(dtab$ALS_plus=="FTD","Yes","No")


factorVars <- c("All_Diagnosis")
vars <- c("Gender", "Age_MRI", "Age_onset", "Age_diagnosis", "Disease_duration", "Surv_onset", 
          "Site of onset", "Type_Familiair_of_Sporadisch", 
          "C9orf72.mutated",  "UNC13A", "Table_ALS_FRS_R.Total_ALS_FRS_R_SCORE", 
          "ECAS_Total_Score", "ECAS_bin","ECAS_ALS_SPECIFIC", "ECAS_sbin","ECAS_ALS_NONSPECIFIC","ECAS_nbin")
varsCont <- c("Age_MRI", "Age_onset", "Age_diagnosis", "Disease_duration", "Surv_onset",  "Table_ALS_FRS_R.Total_ALS_FRS_R_SCORE", 
              "ECAS_Total_Score", "ECAS_ALS_SPECIFIC", "ECAS_ALS_NONSPECIFIC")
dtab$C9orf72.mutated <- factor(dtab$C9orf72.mutated, level=c("No","Yes"))
CreateTableOne(vars = vars, strata = "All_Diagnosis", data = dtab, factorVars = factorVars)
tableOne <- CreateTableOne(vars = vars, strata = c("All_Diagnosis"), data = dtab, 
                           factorVars = factorVars)
tableOne
print(tableOne, nonnormal= varsCont,minMax=T)
print(tableOne, nonnormal= varsCont, minMax=T, cramVars = c("Sex", "C9orf72.mutated"), missings = TRUE)



dtab2 <- filter(dtab,All_Diagnosis=="ALS")
dtab2 <- filter(dtab,All_Diagnosis=="ALS" & (is.na(C9orf72.mutated) | C9orf72.mutated=="No" )& !is.na(lh_bankssts_thickness))
factorVars <- c("riskvariant")
vars <- c("Gender", "Age_MRI", "Age_onset", "Age_diagnosis", "Disease_duration", "Surv_onset", 
          "Site_of_Onset", "Type_Familiair_of_Sporadisch", 
          "C9orf72.mutated",  "UNC13A", "Table_ALS_FRS_R.Total_ALS_FRS_R_SCORE", "ALSFRS_slope",
          "ECAS_Total_Score", "ECAS_bin","ECAS_ALS_SPECIFIC", "ECAS_sbin","ECAS_ALS_NONSPECIFIC","ECAS_nbin")
varsCont <- c("Age_MRI", "Age_onset", "Age_diagnosis", "Disease_duration", "Surv_onset",  "Table_ALS_FRS_R.Total_ALS_FRS_R_SCORE","ALSFRS_slope", 
              "ECAS_Total_Score", "ECAS_ALS_SPECIFIC", "ECAS_ALS_NONSPECIFIC")
#CreateTableOne(vars = vars, strata = "riskvariant", data = dtab2, factorVars = factorVars)
tableOne <- CreateTableOne(vars = vars, strata = c("riskvariant"), data = dtab2, 
                           factorVars = factorVars)
print(tableOne, nonnormal= varsCont, minMax=T, cramVars = c("Sex", "C9orf72.mutated"), missings = TRUE)


test1 <- ldply(unique(ecas_long$`ALS number`), function (i){
  #print(i)
  df1 <- as.data.frame(filter(ecas_long, `ALS number`==i & !is.na(`ECAS Total Score`)))
  df1$`ECAS Total Score` <- as.numeric(df1$`ECAS Total Score`)
  if (nrow(df1) <2){
    return(data.frame(ALSnr=i,ECAS_slope=NA))
  } else {
    ecas_first <- which(df1$`Datum ECAS` == min(df1$`Datum ECAS`,na.rm = T))
    ecas_last <- which(df1$`Datum ECAS` == max(df1$`Datum ECAS`,na.rm = T))
    if ((length(ecas_first) + length(ecas_first)) >2 ){
      return(data.frame(ALSnr=i,ECAS_slope=NA))
    } else {
      diff_days <- as.numeric(max(df1$`Datum ECAS`,na.rm = T) - min(df1$`Datum ECAS`,na.rm = T))
      ecas_diff <- df1[ecas_last,"ECAS Total Score"]-df1[ecas_first,"ECAS Total Score"]
      ecas_slope <- ecas_diff/(diff_days/365.25)
      return(data.frame(ALSnr=i,ECAS_slope=ecas_slope))
    }
  }
})
dtab2 <- join_all(list(dtab2,test1),type="left", by="ALSnr", match="all")



## testzone
# Ecas Behaviour

ecasPresBhv <- filter(ecas_long, !is.na(Symptoms) | !is.na(`Total Behaviour`))
ECAS_FIRST <- ldply(unique(ecasPresBhv$`ALS number`), function (i){
  df1 <- as.data.frame(filter(ecasPresBhv, `ALS number`==i ))
  df1$ECAS_first <- ifelse(is.na(df1$`Datum ECAS`),NA,ifelse(df1$`Datum ECAS`==min(df1$`Datum ECAS`,na.rm=T) & !is.na(df1$`Datum ECAS`),"1","0"))
  #  if()
  return(df1[c("ALS number", "ECAS_FU", "ECAS_first")])
})
ecasBhv <- join_all(list(ecasPresBhv, ECAS_FIRST), by=c("ALS number", "ECAS_FU"), type="left", match = "all")
unc_dataB<-unc_data
colnames(unc_dataB)[1] <- "ALS number"
ecasBhv <- join_all(list(ecasBhv, unc_dataB), by="ALS number",type="left", match="all")
ecasBhv <- filter(ecasBhv, !is.na(riskvariant))
ecasBhv1 <-filter(ecasBhv, ECAS_first==1)

ggplot(ecasBhv1,aes(riskvariant,fill=`Total Behaviour`)) + geom_bar()
prop.table(table(ecasBhv3$Symptoms, ecasBhv3$riskvariant),2)


ggplot(ecasBhv1,aes(riskvariant,fill=`Total Behaviour`)) + geom_bar()





######################################################################################################
#
d11 <- select(dALSx, ALSnr, Age_MRI, Disease_duration, Gender, riskvariant, genotype, dose, All_Diagnosis,eTIV,
              lh_bankssts_thickness:lh_MeanThickness_thickness,rh_bankssts_thickness:rh_MeanThickness_thickness,
              Left.Thalamus.Proper:Left.Pallidum, Left.Hippocampus, Left.Amygdala, Left.Accumbens.area,
              Right.Thalamus.Proper:Right.Accumbens.area,
              Left.Lateral.Ventricle, Left.Inf.Lat.Vent, Right.Lateral.Ventricle, Right.Inf.Lat.Vent,
              left_Hippocampal_tail:right_Whole_hippocampus, fmajor_RD_Avg_Weight:RD_rh_BA4p)
volROI <- c("Left.Thalamus.Proper", "Left.Caudate", "Left.Putamen", "Left.Pallidum", 
  "Left.Hippocampus", "Left.Amygdala", "Left.Accumbens.area", "Right.Thalamus.Proper", 
  "Right.Caudate", "Right.Putamen", "Right.Pallidum", "Right.Hippocampus", 
  "Right.Amygdala", "Right.Accumbens.area", "Left.Lateral.Ventricle", 
  "Left.Inf.Lat.Vent", "Right.Lateral.Ventricle", "Right.Inf.Lat.Vent", 
  "left_Hippocampal_tail", "left_subiculum", 
  "left_CA1", "left_hippocampal.fissure", "left_presubiculum", 
  "left_parasubiculum", "left_molecular_layer_HP", "left_GC.ML.DG", 
  "left_CA3", "left_CA4", "left_fimbria", "left_HATA", "left_Whole_hippocampus", 
  "right_Hippocampal_tail", "right_subiculum", "right_CA1", "right_hippocampal.fissure", 
  "right_presubiculum", "right_parasubiculum", "right_molecular_layer_HP", 
  "right_GC.ML.DG", "right_CA3", "right_CA4", "right_fimbria", 
  "right_HATA", "right_Whole_hippocampus")

a1 <- which(colnames(d11)=="lh_bankssts_thickness")
#a2 <- which(colnames(d11)=="Left.Cerebellum.Cortex")
#a2 <- which(colnames(d11)=="right_Whole_hippocampus")
a2 <- which(colnames(d11)=="RD_rh_BA4p")

#d11 <- filter(d11,!is.na(genotype) & !is.na(lh_bankssts_thickness))

l1 <- lapply(a1:a2, function(i){
  require(broom)
  if(colnames(d11)[i] %in% volROI){
    df1 <- select_(d11, "Age_MRI", "Gender", "riskvariant", "eTIV", brain_region=colnames(d11)[i])
    m1 <- lm(brain_region ~ Age_MRI + Gender + riskvariant + eTIV , data=df1)
  } else {
    df1 <- select_(d11, "Age_MRI", "Gender", "riskvariant", brain_region=colnames(d11)[i])
    m1 <- lm(brain_region ~ Age_MRI + Gender + riskvariant, data=df1)
    }
  
  
  s1 <- tidy(m1)
  out1 <- filter(s1, grepl("riskvariant", term))
  out1$brain_region <- rep(colnames(d11)[i], nrow(out1))
  return(out1)
})
l1_out1 <- do.call(rbind.data.frame, l1)
#arrange(l1_out1, p.value)
View(arrange(l1_out1, p.value))
#write.table(l1_out1,"output/ALS_dose.txt",col.names = T, row.names = F,quote=F)




l1 <- lapply(a1:a2, function(i){
  require(broom)
    df1 <- select_(d11, "Age_MRI", "Gender", "riskvariant", brain_region=colnames(d11)[i])
    m1 <- lm(brain_region ~ Age_MRI + Gender + riskvariant, data=df1)
  
  
  s1 <- tidy(m1)
  out1 <- filter(s1, grepl("riskvariant", term))
  out1$brain_region <- rep(colnames(d11)[i], nrow(out1))
  return(out1)
})
l1_out1 <- do.call(rbind.data.frame, l1)
#arrange(l1_out1, p.value)
View(arrange(l1_out1, p.value))






#interactions
l1 <- lapply(a1:a2, function(i){
  require(broom)
  df1 <- select_(d11, "Age_MRI", "Gender", "riskvariant", "All_Diagnosis", brain_region=colnames(d11)[i])
  m1 <- lm(brain_region ~ Age_MRI + Gender + riskvariant*All_Diagnosis, data=df1)
  s1 <- tidy(m1)
  out1 <- filter(s1, grepl("riskvariant", term))
  out1$brain_region <- rep(colnames(d11)[i], nrow(out1))
  return(out1)
})
l1_out2 <- do.call(rbind.data.frame, l1)
#arrange(l1_out2, p.value)
#View(arrange(l1_out2, p.value))
write.table(arrange(l1_out1, p.value),"output/AC_riskvariant.txt",col.names = T, row.names = F,quote=F)



# Longitudinal MRI analysis
dALSl <- filter(dmaster, (All_Diagnosis=="ALS"))
dALSl <- filter(dALSl, !is.na(genotype))
dALSl <- filter(dALSl, is.na(C9orf72.mutated)|C9orf72.mutated=="No")

# per region
region <- "lh_inferiortemporal_thickness"
formlong <- formula(paste0("`",region,"`","~ Age_MRI + Gender  + riskvariant+ (1|ALSnr) + ( 0+ Disease_duration| ALSnr)"))
m1 <- lmer(formlong , dALSl ) 

# for all regions
thROI <-  c("lh_bankssts_thickness", 
"lh_caudalanteriorcingulate_thickness", "lh_caudalmiddlefrontal_thickness", 
"lh_cuneus_thickness", "lh_entorhinal_thickness", "lh_fusiform_thickness", 
"lh_inferiorparietal_thickness", "lh_inferiortemporal_thickness", 
"lh_isthmuscingulate_thickness", "lh_lateraloccipital_thickness", 
"lh_lateralorbitofrontal_thickness", "lh_lingual_thickness", 
"lh_medialorbitofrontal_thickness", "lh_middletemporal_thickness", 
"lh_parahippocampal_thickness", "lh_paracentral_thickness", "lh_parsopercularis_thickness", 
"lh_parsorbitalis_thickness", "lh_parstriangularis_thickness", 
"lh_pericalcarine_thickness", "lh_postcentral_thickness", "lh_posteriorcingulate_thickness", 
"lh_precentral_thickness", "lh_precuneus_thickness", "lh_rostralanteriorcingulate_thickness", 
"lh_rostralmiddlefrontal_thickness", "lh_superiorfrontal_thickness", 
"lh_superiorparietal_thickness", "lh_superiortemporal_thickness", 
"lh_supramarginal_thickness", "lh_frontalpole_thickness", "lh_temporalpole_thickness", 
"lh_transversetemporal_thickness", "lh_insula_thickness", "lh_MeanThickness_thickness", 

 "rh_bankssts_thickness", "rh_caudalanteriorcingulate_thickness", 
"rh_caudalmiddlefrontal_thickness", "rh_cuneus_thickness", "rh_entorhinal_thickness", 
"rh_fusiform_thickness", "rh_inferiorparietal_thickness", "rh_inferiortemporal_thickness", 
"rh_isthmuscingulate_thickness", "rh_lateraloccipital_thickness", 
"rh_lateralorbitofrontal_thickness", "rh_lingual_thickness", 
"rh_medialorbitofrontal_thickness", "rh_middletemporal_thickness", 
"rh_parahippocampal_thickness", "rh_paracentral_thickness", "rh_parsopercularis_thickness", 
"rh_parsorbitalis_thickness", "rh_parstriangularis_thickness", 
"rh_pericalcarine_thickness", "rh_postcentral_thickness", "rh_posteriorcingulate_thickness", 
"rh_precentral_thickness", "rh_precuneus_thickness", "rh_rostralanteriorcingulate_thickness", 
"rh_rostralmiddlefrontal_thickness", "rh_superiorfrontal_thickness", 
"rh_superiorparietal_thickness", "rh_superiortemporal_thickness", 
"rh_supramarginal_thickness", "rh_frontalpole_thickness", "rh_temporalpole_thickness", 
"rh_transversetemporal_thickness", "rh_insula_thickness", "rh_MeanThickness_thickness")

library(lmerTest)
test2 <- ldply(thROI, function(i){
  formlong <- formula(paste0("`",i,"`","~ Age_MRI + Gender  + riskvariant+ (1|ALSnr) + ( 0+ Disease_duration| ALSnr)"))
  m1 <- lmer(formlong , dALSl ) 
  coefs <- summary(m1)$coefficients
  out <- c(i,coefs[grepl("riskvariant", rownames(coefs)), c("Estimate","Pr(>|t|)")])
  return(out)
})


# Get unbiased p-values through bootstrapping
mysumm1 <- function(.){b=getME(.,"beta")}
nSim=10000

t1 <- mysumm1(m1)
t2 <- data.frame(t1)
set.seed(123)
boo01 <- bootMer(m1, mysumm1, nsim=nSim)
bCI1 <- boot.ci(boo01, index=5, type="norm")
SE1 <- (bCI1$normal[3]-bCI1$normal[2])/(2*qnorm(0.975))
z1 <- as.vector(bCI1$t0)/SE1
rawP <- 2*pnorm(-abs(z1))



######################################################################################################

#ECAS analysis
ecas_data <- join_all(list(ecas_long, raw_clin_wide), by="ALS number", type="left", match="all")
ecas_data$Age_ecas <- as.numeric((ecas_data$`Datum ECAS` - as.Date(ecas_data$`Date of birth`))/365.25)
ecas_data$onset2ecas <- as.numeric((ecas_data$`Datum ECAS` - as.Date(ecas_data$`Date of onset`))/365.25*12)

unc_data2 <- unc_data
colnames(unc_data2)[1] <- "ALS number"
c9_data2 <- c9_data
colnames(c9_data2)[1] <- "ALS number"

d2 <- join_all(list(ecas_data,unc_data2,c9_data2), by="ALS number", type="left", match="all")
#d2 <- filter(d2, `ALS number` %in% d11$ALSnr)
d2$riskvariant <- as.factor(d2$riskvariant)
library(lme4);library(lmerTest)
colnames(d2)[grepl("ISCED",colnames(d2))] <- "ISCED"
d2$ISCED <- as.factor(d2$ISCED)
d2$Visuospatial <- as.numeric(d2$Visuospatial)
d2$Memory <- as.numeric(d2$Memory)
d2$Executive <- as.numeric(d2$Executive)
d2$`Fluency Free + Fixed` <- as.numeric(d2$`Fluency Free + Fixed`)
d2$Language <- as.numeric(d2$Language)
colnames(d2)[which(colnames(d2)=="Sex (gender at birth)__1")] <- "Gender"
d2$Gender <- as.factor(d2$Gender)
d2$UNC13A_Genotype <- as.factor(ifelse(d2$riskvariant=="0","AA+AC","CC"))

ECAS_FIRST <- ldply(unique(ecas_long$`ALS number`), function (i){
  df1 <- as.data.frame(filter(ecas_long, `ALS number`==i ))
  df1$ECAS_first <- ifelse(is.na(df1$`Datum ECAS`),NA,ifelse(df1$`Datum ECAS`==min(df1$`Datum ECAS`,na.rm=T) & !is.na(df1$`Datum ECAS`),"1","0"))
#  if()
  return(df1[c("ALS number", "ECAS_FU", "ECAS_first")])
})


# Of mogelijk juist mean slope gebruiken?
ECAS_SLOPE <- ldply(unique(ecas_long$`ALS number`), function (i){
  #print(i)
  df1 <- as.data.frame(filter(ecas_long, `ALS number`==i & !is.na(`ECAS Total Score`)))
  df1$`ECAS Total Score` <- as.numeric(df1$`ECAS Total Score`)
  if (nrow(df1) <2){
    return(data.frame(ALSnr=i,ECAS_slope=NA))
  } else {
    ecas_first <- which(df1$`Datum ECAS` == min(df1$`Datum ECAS`,na.rm = T))
    ecas_last <- which(df1$`Datum ECAS` == max(df1$`Datum ECAS`,na.rm = T))
    if ((length(ecas_first) + length(ecas_first)) >2 ){
      return(data.frame(ALSnr=i,ECAS_slope=NA))
    } else {
      diff_days <- as.numeric(max(df1$`Datum ECAS`,na.rm = T) - min(df1$`Datum ECAS`,na.rm = T))
      ecas_diff <- df1[ecas_last,"ECAS Total Score"]-df1[ecas_first,"ECAS Total Score"]
      ecas_slope <- ecas_diff/(diff_days/365.25)
      return(data.frame(ALSnr=i,ECAS_slope=ecas_slope))
    }
  }
})


# ECAS learning effect: number of attempts, number of alterations?






d2 <- join_all(list(d2,ECAS_FIRST), by=c("ALS number","ECAS_FU"),type="left", match="first")
#ECAS crosssectional
d3 <- filter(d2, All_Diagnosis=="ALS" & ECAS_first=="1")
for (i in c("ECAS Total Score", "ECAS ALS SPECIFIC", "ECAS ALS NONSPECIFIC", 
       "Executive", "Fluency Free + Fixed", "Language", "Memory", "Visuospatial")){
  form1 <- formula(paste0("`",i,"` ~  Age_ecas + Gender + ISCED + riskvariant"))
  cm1 <- lm(form1,d3)
  print(i)
  print(summary(cm1))
  plot(effect("riskvariant",cm1))
  readline(prompt="Press [enter] to continue")
}

is.num <- sapply(test100, is.numeric)
test100[is.num] <- lapply(test100[is.num], round, 2)
ecross <- ldply(c("ECAS Total Score", "ECAS ALS SPECIFIC", "ECAS ALS NONSPECIFIC", 
            "Executive", "Fluency Free + Fixed", "Language", "Memory", "Visuospatial"),
            function(i){
  form1 <- formula(paste0("`",i,"` ~  Age_ecas + Gender + ISCED + riskvariant"))
  cm1 <- lm(form1,d3)
  sum1 <- summary(cm1)$coefficients
  ef1 <- effect("riskvariant",cm1)
  ef2 <- data.frame(ef1)
  is.num <- sapply(ef2, is.numeric)
  ef2[is.num] <- lapply(ef2[is.num], round, 2)
  rv0 <- ef2[ef2$riskvariant=="0",]
  prv0 <- paste0(rv0$fit," (",rv0$lower,"-",rv0$upper,")")
  rv1 <- ef2[ef2$riskvariant=="1",]
  prv1 <- paste0(rv1$fit," (",rv1$lower,"-",rv1$upper,")")
  pval <- sum1[nrow(sum1),ncol(sum1)]
  return(data.frame(Domain=i,"AA/AC"=prv0, CC=prv1,pval=pval))
})

is.num <- sapply(test100, is.numeric)
test100[is.num] <- lapply(test100[is.num], round, 2)

form1 <- formula(paste0("`",i,"` ~  Age_ecas + ISCED + riskvariant"))
cm1 <- lm(form1,d3)
ef <- effect("riskvariant",cm1)
plot(ef)



#ECAS LONGITUDINAL
ECAS_cols <- c("ECAS Total Score", "ECAS ALS SPECIFIC", "ECAS ALS NONSPECIFIC", 
               "Executive", "Fluency Free + Fixed", "Language", "Memory", "Visuospatial")
d2[ECAS_cols] <- lapply(d2[ECAS_cols], as.numeric)
i=ECAS_cols[8]
formlong <- formula(paste0("`",i,"`","~ Age_ecas + Gender + ISCED + log(onset2ecas)*UNC13A_Genotype+ (1 + log(onset2ecas) | `ALS number`)"))
m1 <- lmer(formlong , d2 ) # isket? age?,  
ef <- effect("log(onset2ecas)*UNC13A_Genotype",m1,quantiles=c(0.25,0.5,0.75))
x <- as.data.frame(ef)

ggplot(x, aes(y=fit, x=log(onset2ecas), col=UNC13A_Genotype, fill=UNC13A_Genotype)) + geom_line() +
  geom_ribbon(colour=NA, alpha=0.1, aes(ymin=lower,ymax=upper)) + 
  geom_rug(data=ef$data,aes(y=NULL),sides="b") + ylab(i) + xlab("log(Time since onset (mth))") +
  theme(axis.title.x=element_text(size=18), axis.title.y=element_text(size=18))+
  scale_color_discrete(name="UNC13A genotype") +scale_fill_discrete(name="UNC13A genotype") + theme(legend.title=element_text(size=16))



d2$dose <- as.factor(d2$dose)
m1 <- lmer(Visuospatial ~ Age_ecas + Gender + ISCED+ log(onset2ecas)*dose + (1 + log(onset2ecas) | `ALS number`), d2 ) # isket? age?,  
ef <- effect("log(onset2ecas)*dose",m1)
x <- as.data.frame(ef)

ggplot(x, aes(y=fit, x=log(onset2ecas), col=dose, fill=dose)) + geom_line() +
  geom_ribbon(colour=NA, alpha=0.1, aes(ymin=lower,ymax=upper)) + 
  geom_rug(data=ef$data,aes(y=NULL),sides="b")





# NO log transform

m1 <- lmer(`ECAS Total Score` ~ Age_ecas + ISCED  + onset2ecas + riskvariant + ( onset2ecas | `ALS number`), d2 ) # isket? age?,  
ef <- effect("riskvariant",m1)
x <- as.data.frame(ef)

ggplot(x, aes(y=fit, x=riskvariant, col=riskvariant, fill=riskvariant)) + geom_line() +
  geom_ribbon(colour=NA, alpha=0.1, aes(ymin=lower,ymax=upper)) + 
  geom_rug(data=ef$data,aes(y=NULL),sides="b")

ggplot(d2,aes(y=`ECAS Total Score`, x=log(onset2ecas), col=riskvariant)) + geom_point()





#C9 longitudinaal
d2$C9 <- ifelse(d2$C9orf72.mutated =="Yes" & !is.na(d2$C9orf72.mutated),"Yes","No")
d2$C9 <- as.factor(d2$C9)
m1 <- lmer(`ECAS Total Score` ~ Age_ecas + ISCED  + log(onset2ecas)*C9 + ( log(onset2ecas) | `ALS number`), d2 ) # isket? age?,  
ef <- effect("log(onset2ecas)*C9",m1)
x <- as.data.frame(ef)

ggplot(x, aes(y=fit, x=log(onset2ecas), col=C9, fill=C9)) + geom_line() +
  geom_ribbon(colour=NA, alpha=0.1, aes(ymin=lower,ymax=upper)) + 
  geom_rug(data=ef$data,aes(y=NULL),sides="b")







table(d2$riskvariant,d2$`ALS plus`)
#Controls
d1 <- filter(ecas_data, All_Diagnosis=="Controle")
unc_data2 <- unc_data
colnames(unc_data2)[1] <- "ALS number"
d2 <- join_all(list(d1,unc_data2), by="ALS number", type="left", match="all")
d2$`ECAS Total Score` <- as.numeric(d2$`ECAS Total Score`)
d2$riskvariant <- as.factor(d2$riskvariant)
library(lme4);library(lmerTest)
colnames(d2)[grepl("ISCED",colnames(d2))] <- "ISCED"
d2$ISCED <- as.factor(d2$ISCED)

m1 <- lmer(`ECAS Total Score` ~ Age_ecas + ISCED + log(onset2ecas)*riskvariant + (1 + log(onset2ecas) | `ALS number`), d2 ) # isket? age?,  
ef <- effect("log(onset2ecas)*riskvariant",m1,quantiles=c(0.25,0.5,0.75))
x <- as.data.frame(ef)

ggplot(x, aes(y=fit, x=log(onset2ecas), col=riskvariant, fill=riskvariant)) + geom_line() +
  geom_ribbon(colour=NA, alpha=0.1, aes(ymin=lower,ymax=upper)) + 
  geom_rug(data=ef$data,aes(y=NULL),sides="b")







# TEST
ecas_data <- join_all(list(ecas_long, raw_clin_wide), by="ALS number", type="left", match="all")
ecas_data$Age_ecas <- as.numeric((ecas_data$`Datum ECAS` - as.Date(ecas_data$`Date of birth`))/365.25)
ecas_data$onset2ecas <- as.numeric((ecas_data$`Datum ECAS` - as.Date(ecas_data$`Date of onset`))/365.25*12)

d1 <- filter(ecas_data, All_Diagnosis=="ALS"|All_Diagnosis=="Controle")
dunc_data2 <- unc_data
colnames(unc_data2)[1] <- "ALS number"
d2 <- join_all(list(d1,unc_data2), by="ALS number", type="left", match="all")
d2$`ECAS Total Score` <- as.numeric(d2$`ECAS Total Score`)
d2$riskvariant <- as.factor(d2$riskvariant)
d2$dose <- as.factor(d2$dose)
d2$All_Diagnosis <- as.factor(d2$All_Diagnosis)
m1 <- lmer(`ECAS Total Score` ~ Age_ecas + onset2ecas*riskvariant*All_Diagnosis + (1 + onset2ecas | `ALS number`), d2 ) # isket? age?,  
ef <- effect("onset2ecas*riskvariant*All_Diagnosis",m1)
x <- as.data.frame(ef)

ggplot(x, aes(y=fit, x=onset2ecas, col=riskvariant, fill=riskvariant)) + geom_line() +
  geom_ribbon(colour=NA, alpha=0.1, aes(ymin=lower,ymax=upper)) + 
  geom_rug(data=ef$data,aes(y=NULL),sides="b")


# Useful functions
#write.table(arrange(l1_out1, p.value),"/Users/htan4/Documents/UNC13_test.txt",col.names = T, row.names = F, quote = F,sep=",")
#View(raw_clin_wide[raw_clin_wide$`ALS number` %in% c("ALS08456","ALS25181"),])




######################################################################################################

temp4 <- temp6 <- NULL
for (i in c("PLS", "ALS", "PMA")) {
  df1 <- filter(raw_clin_wide, All_Diagnosis==i)
  temp1 <- coxph(Surv(Surv_onset, Dead==1) ~ 1 , df1)
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
  scale_color_discrete(name="P") + theme(legend.title=element_text(size=16)) + xlim(0,120)



temp4 <- temp6 <- NULL
for (i in c("0", "1")) {
  df1 <- filter(d10, riskvariant==i)
  temp1 <- coxph(Surv(Surv_onset, Dead==1) ~ Gender + Age_onset + `Site of Onset`, df1)
  temp2 <- survfit(temp1)
  temp3 <-  cbind.data.frame(survival = temp2$surv, time = temp2$time/12)
  if (i=="0"){
    temp3$UNC13A <- "AA + AC"
  } else {
    temp3$UNC13A <- "CC"
  }
  
  temp4 <- rbind(temp4,temp3)
  
  temp5 <- c(i, unlist(quantile(temp2, probs = c(0.25, 0.50, 0.75))$quantile))
  temp6 <- rbind(temp6,temp5)
}

ggplot(temp4, aes(time, survival, col=UNC13A)) +geom_point(shape=3,size=2) + geom_path(size=0.5)+
#  ggtitle("Kaplan-Meier Survival Curve")+ theme(plot.title=element_text(size=20, face="bold"))+
  theme(axis.title.x=element_text(size=18), axis.title.y=element_text(size=18))+
  scale_color_discrete(name="UNC13A genotype") + theme(legend.title=element_text(size=16)) +xlim(0,10) +
  ylab("Survival (%)") + xlab("Duration after disease onset (yr)")








dALS <- filter(d9, All_Diagnosis=="ALS" & scan_fu==1 )
temp4 <- temp6 <- NULL
for (i in c("Yes", "No")) {
  df1 <- filter(dALS, C9orf72.mutated==i)
  temp1 <- coxph(Surv(Surv_onset, Dead==1) ~ Gender + Age_onset + `Site of Onset`, df1)
  temp2 <- survfit(temp1)
  temp3 <-  cbind.data.frame(survival = temp2$surv, time = temp2$time/12)
  temp3$riskvariant <- i
  temp4 <- rbind(temp4,temp3)
  
  temp5 <- c(i, unlist(quantile(temp2, probs = c(0.25, 0.50, 0.75))$quantile))
  temp6 <- rbind(temp6,temp5)
}

ggplot(temp4, aes(time, survival, col=riskvariant)) +geom_point(shape=3,size=2) + geom_path(size=0.5)+
  ggtitle("Kaplan-Meier Survival Curve")+ theme(plot.title=element_text(size=20, face="bold"))+
  theme(axis.title.x=element_text(size=18), axis.title.y=element_text(size=18))+
  scale_color_discrete(name="P") + theme(legend.title=element_text(size=16)) + xlim(0,10)









lx <- lapply(grep("_thickness",colnames(dtab2),value = T), function(i){
  form <- formula(paste(i, "~Age_MRI + Gender + ECAS_slope*riskvariant"))
  model <- lm(form,dtab2)
  summ <- summary(model)$coefficients
  pval <- summ[nrow(summ),ncol(summ)]
  return(pval)
  #eff <- effect("ECAS_slope*riskvariant")
})
View(data.frame(Region=grep("_thickness",colnames(dtab2),value = T),pval=unlist(lx)))




