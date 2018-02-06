#load packages
library(plyr); library(dplyr)

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

d11 <- read_excel("/Volumes/Promise_Pegasus/Harold/Test_env/ALSUNC_20180103.xlsx",sheet=1)
d11$Gender <- as.factor(d11$Gender)
#selecteer data zo dat de variabele waar je in geinteresseerd bent slechts 2 levels kan bevatten (bijv ALS vs CO)
#d2 <- filter(d1, Group!="Thoracic", C9!="long")
d2 <- filter(d11,!is.na(riskvariant)&!is.na(lh_bankssts_thickness))
d2$riskvariant <- factor(d2$riskvariant)
d2$riskvariant <- relevel(d2$riskvariant, ref="0")
#d3 <- d2
#d3$Group <- factor(d3$Group)
a1 <- which(colnames(d11)=="lh_bankssts_thickness")
a2<- which(colnames(d11)=="RD_rh_BA4p")

#maak nieuwe permuatie functie
#maak 3 datasets
#dataset 1 bevat uitkomst (dat wat voor de tilde staat), van alle verschillende modellen
#dataset 2 bevat covariates
#dataset 3 bevat variabele die je wilt shuffelen (in dit geval vaak case vs control)
dd1 <- select(d2, colnames(d2)[a1:a2])
#dd1 <- select(dd1, -matches("precentral_|FA_CC_BA4p"))
dd2 <- select(d2, Age_MRI, Gender)
dd2v <- select(d2, Age_MRI, Gender, eTIV)
dd3 <- select(d2, riskvariant)
df1 <- cbind.data.frame(dd3, dd1, dd2v) #verander deze regel ip niet (omdat perm1 op dit format rekent)

#voor de veryfast option heb je matrices nodig
mm1 <- data.matrix(dd1)
mm2 <- data.matrix(dd2v)
mm3 <- data.matrix(dd3)
mat <- cbind(rep(1, nrow(mm2)), mm3, mm2) #dit is design matrix (=model matrix), verander ip niet omdat perm1 op dit format rekent

#maak formules voor lineaire modellen (alleen nodig als veryfast=FALSE is)
form1 <- lapply(1:length(dd1), function(j){
  outcome1 <- colnames(dd1)[j]
  if (outcome1 %in% volROI){
    form0 <- as.formula(paste0(outcome1, "~", paste(colnames(dd2v), collapse="+"), "+", colnames(dd3)))
  } else {
    form0 <- as.formula(paste0(outcome1, "~", paste(colnames(dd2), collapse="+"), "+", colnames(dd3)))
  }
  return(form0)
})

#functie om p-waardes uit veryfast option te berekenen
#pval1 is bedoeld voor modellen die met lm.fit gefit zijn
#het is gebaseerd op de summary functie van lm. Functie is ip nog experimenteel.
pval1 <- function(linear_model){
  rss1 <- sum(linear_model$residuals^2) #RSS
  resvar1 <- rss1/linear_model$df.residual #residual variance
  
  p <- linear_model$rank
  Qr <- linear_model$qr
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar1)
  est <- linear_model$coefficients[Qr$pivot[p1]]
  tval <- est/se
  pval <- 2*pt(-abs(tval), linear_model$df.residual)
  return(pval)
}

#lineaire modellen kunnen nog sneller gefit worden met .lm.fit.
#daarvoor heb ik onderstaand pval2 geschreven. Hiervoor geldt hetzelfde als voor pval1 (en is eigenlijk
#zelfs nog wat meer experimenteel), maar ook nog sneller :-)
pval2 <- function(linear_model){
  rss1 <- sum(linear_model$residuals^2) #RSS
  rdf1 <- length(linear_model$residuals)-linear_model$rank
  resvar1 <- rss1/rdf1 #residual variance
  
  p <- linear_model$rank
  Qr <- linear_model
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar1)
  est <- linear_model$coefficients[Qr$pivot[p1]]
  tval <- est/se
  pval <- 2*pt(-abs(tval), rdf1)
  return(pval)
}

#permuteren
#nu zo gemaakt dat je bij factoren met >2 levels gewoon apart permuteerd en
#dus onderstaande functie 2 moet runnen na gebruik van relevel
#later ook optie "speed" toegevoegd.
#speed heeft 3 opties (om zeer snel lineaire modellen te kunnen fitten):
# normal      normal (default), is meest robust en traagst
# fast        gebruikt lm.fit (wat experimenteel), is ongeveer 14x sneller dan normal
# veryfast    gebruikt .lm.fit (meest experimenteel), is ongeveer 24x sneller dan normal
# verifieer resultaten van fast en veryfast altijd even met normal (omdat het nog experimenteel is)
#ip zou uit normal, fast en veryfast hetzelfde moeten komen
library(parallel)
nperm <- 10000
set.seed(123)
system.time(perm1 <- mclapply(0:nperm, mc.cores=12, function(i, speed="normal"){
  if(speed=="veryfast"){
    if(i==0){
      out1 <- lapply(1:ncol(mm1), function(j){
        #LET OP! de veryfast optie maakt gebruik van .lm.fit
        missing1 <- which(is.na(mm1[,j]))
        if(length(missing1)>0){
          mat <- mat[-missing1,]
          lm_fit1 <- .lm.fit(mat, mm1[-missing1,j])
        }else{
          lm_fit1 <- .lm.fit(mat, mm1[,j])
        }
        pval <- pval2(lm_fit1)
        pval <- as.vector(pval[2])
        return(pval)
      })
    }else{
      mat[,2] <- sample(mat[,2])
      out0 <- lapply(1:ncol(mm1), function(j){
        #LET OP! de veryfast optie maakt gebruik van .lm.fit
        missing1 <- which(is.na(mm1[,j]))
        if(length(missing1)>0){
          mat <- mat[-missing1,]
          lm_fit1 <- .lm.fit(mat, mm1[-missing1,j])
        }else{
          lm_fit1 <- .lm.fit(mat, mm1[,j])
        }
        pval <- pval2(lm_fit1)
        pval <- as.vector(pval[2])
        return(pval)
      })
      out1 <- min(unlist(out0))
    }
    out2 <- out1
    return(out2)
  }else if(speed=="fast"){
    if(i==0){
      out1 <- lapply(1:ncol(mm1), function(j){
        #LET OP! de veryfast optie maakt gebruik van lm.fit
        missing1 <- which(is.na(mm1[,j]))
        if(length(missing1)>0){
          mat <- mat[-missing1,]
          lm_fit1 <- lm.fit(mat, mm1[-missing1,j])
        }else{
          lm_fit1 <- lm.fit(mat, mm1[,j])
        }
        pval <- pval1(lm_fit1)
        pval <- as.vector(pval[2])
        return(pval)
      })
    }else{
      mat[,2] <- sample(mat[,2])
      out0 <- lapply(1:ncol(mm1), function(j){
        #LET OP! de veryfast optie maakt gebruik van lm.fit
        missing1 <- which(is.na(mm1[,j]))
        if(length(missing1)>0){
          mat <- mat[-missing1,]
          lm_fit1 <- lm.fit(mat, mm1[-missing1,j])
        }else{
          lm_fit1 <- lm.fit(mat, mm1[,j])
        }
        pval <- pval1(lm_fit1)
        pval <- as.vector(pval[2])
        return(pval)
      })
      out1 <- min(unlist(out0))
    }
    out2 <- out1
    return(out2)
  }else{
    if(i==0){
      out1 <- lapply(1:length(form1), function(j){
        lm1 <- lm(form1[[j]], data=df1)
        lm2 <- summary(lm1)$coefficients
        #rows1 <- grep(colnames(dd3), rownames(lm2))
        #lm3 <- lm2[rows1,ncol(lm2)]
        lm3 <- lm2[nrow(lm2), ncol(lm2)]
        return(lm3)
      })
    }else{
      df1[,1] <- sample(df1[,1])
      out0 <- lapply(1:length(form1), function(j){
        lm1 <- lm(form1[[j]], data=df1)
        lm2 <- summary(lm1)$coefficients
        #rows1 <- grep(colnames(dd3), rownames(lm2))
        #lm3 <- lm2[rows1,ncol(lm2)]
        lm3 <- lm2[nrow(lm2), ncol(lm2)]
        return(lm3)
      })
      out1 <- min(unlist(out0))
    }
  }

  out2 <- out1
  return(out2)
}))

origp <- unlist(perm1[[1]]) #orig p-values
permp <- unlist(perm1[-1]) #permutated p-values
#empP <- lapply(1:length(origp), function(i){
#  out1 <- (length(which(permp<origp[i]))+1)/(length(permp)+1)
#  return(out1)
#})
#empP
empP <- lapply(1:length(origp), function(i){
  out1 <- (length(which(permp<origp[i])))/(length(permp))
  return(out1)
})
unlist(empP)
hist(unlist(perm1[-1]))


