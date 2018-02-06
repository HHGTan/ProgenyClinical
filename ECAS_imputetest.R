


ecas_data <- join_all(list(ecas_long, raw_clin_wide), by="ALS number", type="left", match="all")
colnames(ecas_data)[which(colnames(ecas_data)=="ALS number")] <- "ALSnr"
ecas_data2 <- join_all(list(ecas_data,unc_data,c9_data), by="ALSnr", type="left", match="all")

ecas_data2 <-  filter(ecas_data2, All_Diagnosis=="ALS" & (C9orf72.mutated=="No" | is.na(C9orf72.mutated)) &
                        (!is.na(`Datum ECAS`) | !is.na(`ECAS ALS NONSPECIFIC`) | !is.na(`ECAS ALS SPECIFIC`) | 
                                                           !is.na(`ECAS Total Score`) | !is.na(`Executive`) | !is.na(`Fluency Free + Fixed`) | 
                                                           !is.na(`Language`) | !is.na(`Memory`) | !is.na(`Visuospatial`) ) )

ecas_data2$Age_ECAS <-  as.numeric(ecas_data2$`Datum ECAS` - as.Date(ecas_data2$`Date of birth`))/365.25

ed.mis <- ecas_data2[,c("ALSnr", "ECAS_FU", "Age_ECAS", "ECAS ALS NONSPECIFIC", 
                        "ECAS ALS SPECIFIC", "ECAS Total Score", "Executive", "Fluency Free + Fixed", 
                        "Language", "Memory", "Visuospatial",
                        "Date of birth", "Sex (gender at birth)__1", 
                        "ISCED 1997", "All_Diagnosis", "Age_onset", "Age_diagnosis","riskvariant"
)]
colnames(ed.mis)[which(colnames(ed.mis)=="Sex (gender at birth)__1")] <- "Gender"
colnames(ed.mis)[which(colnames(ed.mis)=="ISCED 1997")] <- "ISCED"
colnames(ed.mis)[which(colnames(ed.mis)=="Fluency Free + Fixed")] <- "Fluency"

colnames(ed.mis)[which(colnames(ed.mis)=="ECAS ALS NONSPECIFIC")] <- "ECAS_nspec"
colnames(ed.mis)[which(colnames(ed.mis)=="ECAS ALS SPECIFIC")] <- "ECAS_spec"
colnames(ed.mis)[which(colnames(ed.mis)=="ECAS Total Score")] <- "ECAS_Total"


library(mice)
library(VIM)
mice_plot <- aggr(ed.mis, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(ed.mis), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))

#ed.mis2 <- subset(ed.mis, select = -c(ISCED, `ALS number`, ECAS_FU, `Datum ECAS`, `Date of birth`, Gender, All_Diagnosis))

ed.mis[sapply(ed.mis,is.character)] <- lapply(ed.mis[sapply(ed.mis,is.character)], as.factor)
ed.mis$ISCED <- as.factor(ifelse(ed.mis$ISCED=="1111",NA,as.character(ed.mis$ISCED)))


init = mice(ed.mis, maxit=0) 
meth = init$method
predM = init$predictorMatrix
predM[,c("ALSnr","All_Diagnosis")] <- 0
meth[c("ISCED","riskvariant")]="logreg" 

imputed_Data <- mice(ed.mis, m=10, maxit = 50, method = meth,predictorMatrix=predM, seed = 500)

fit <- with(data = imputed_Data, exp = lm(ECAS_Total ~ Age_onset + ECAS_spec + ECAS_nspec)) 
fit <- with(data = imputed_Data, exp = lmer(ECAS_Total ~ Gender + ISCED + Age_ECAS*riskvariant + (0 + Age_ECAS | ALSnr) )) 

#Age_ecas + Gender + ISCED + onset2ecas*riskvariant + (0 + Age_ECAS | ALSnr)
            
#combine results of all 5 models
combine <- pool(fit)
summary(combine)


completeData <- complete(imputed_Data,2)




#test1 <- filter(dmaster, (All_Diagnosis=="ALS"| All_Diagnosis=="Controle") & scan_fu==1)
#test2 <- test1[,c("ALSnr","All_Diagnosis")]

#test3  <- read_excel(path="/Volumes/Promise_Pegasus/Harold/Test_env/Copy\ of\ UNC13A_RvE_100118.xlsx",sheet=1)
#colnames(test3)[1] <- "ALSnr"
#colnames(test3)[5] <- "All_Diagnosis"

#test4 <- join_all(list(test3,test1), by="ALSnr",type="full",match="first")
#writexl::write_xlsx(test4[,c("ALSnr","All_Diagnosis")],path = "/Volumes/Promise_Pegasus/Harold/Test_env/IDs_UNCneeded.xlsx",col_names = T)

#length(unique(test4$ALSnr))


#test5 <- filter(ecas_long, (!is.na(`Datum ECAS`) | !is.na(`ECAS Total Score`)) ) 
#unique_ecasid <- unique(test5$`ALS number`)

#test6 <- filter(raw_clin_wide, (`ALS number` %in% unique_ecasid) & Diagnosis=="ALS")
#test8 <- filter(test6, (is.na(`Date of birth`) | is.na(`Date of onset`) | is.na(`ISCED 1997`)) & Diagnosis=="ALS" )
#test9 <- test8[,c("ALS number" ,"Date of onset", "ISCED 1997")]


#test8 <- filter(raw_clin_wide, (is.na(`Date of birth`) | is.na(`Date of onset`) | is.na(`ISCED 1997`)) & Diagnosis=="ALS" )
#test9 <- test8[,c("ALS number" ,"Date of birth", "Date of onset", "ISCED 1997")]

#writexl::write_xlsx(test9,path = "/Volumes/Promise_Pegasus/Harold/Test_env/ECASanalysis_missings.xlsx",col_names = T)




