#Evaluation of gene expression-based predictors of lymph node metastasis in bladder cancer

library(nestedcv)
library(caret)
library(glmnet)
library(ranger)
library(pROC)

#TCGA data filtered for low expression genes in Count Data
TCGA_Filtering_NoExp <-
  read.delim(file = "TCGA_TMM.txt",
             sep = "\t",
             row.names = 1)

#########
#Read and prepare main data (TCGA)
TCGA_GEX <-
  read.delim(
    file = "TCGA_BLCA_TMM406_filtered_symbols_metadata.txt",
    sep = "\t",
    row.names = 1,
    skip = 4
  )
TCGA_header <-
  read.delim(
    file = "TCGA_BLCA_TMM406_filtered_symbols_metadata.txt",
    sep = "\t",
    nrows = 1,
    row.names = 1
  )
colnames(TCGA_GEX) <- colnames(TCGA_header)
rm(TCGA_header)
TCGA_GEX[1:5, 1:5]
table(rownames(TCGA_Filtering_NoExp) %in% rownames(TCGA_GEX))
TCGA_GEX <- TCGA_GEX[rownames(TCGA_Filtering_NoExp),]

Metadata <-
  read.delim(
    file = "TCGA_BLCA_TMM406_filtered_symbols_metadata.txt",
    sep = "\t",
    row.names = 1,
    nrows = 4
  )
Metadata <- as.data.frame(t(Metadata))
Metadata$pN_status_YN <- Metadata$pN_status

Metadata$pN_status_YN[Metadata$pN_status_YN == 0] <- "NO"
Metadata$pN_status_YN[Metadata$pN_status_YN == 1] <- "YES"

#Change gene symbol dashes to underscore for ranger compatibility
Original_rownames <- rownames(TCGA_GEX)
RF_rownames <- gsub("-", "_", Original_rownames)
rownames(TCGA_GEX) <- RF_rownames

#Create a scaled version (used when applying Lund2017 models to the TCGA cohort)
TCGA_GEX_var <- apply(TCGA_GEX, 1, var)

TCGA_GEX_scaled <-
  t(apply(log2(TCGA_GEX[c(TCGA_GEX_var != 0),] + 1), 1, scale))

colnames(TCGA_GEX_scaled) <- colnames(TCGA_GEX)
x_scaled <- t(TCGA_GEX_scaled)

#Read LN predictor signatures
Signatures <-
  read.delim(file = "Signatures.txt", sep = "\t", header = F)
Signatures_Original <- Signatures
Signatures[, 2] <-
  gsub("-", "_", Signatures[, 2]) #Change gene symbol dashes to underscore for ranger compatibility
Signatures[, 2][!(Signatures[, 2] %in% rownames(TCGA_GEX))] #Check if any signature genes are missing

Signatures_Filtered <-
  Signatures[c(Signatures[, 2] %in% rownames(TCGA_GEX)),]

Signatures_Filtered_split <-
  split(Signatures_Filtered[, 2], Signatures_Filtered[, 1])
#Signatures <- split(Signatures[, 2], Signatures[, 1])



#Read and prepare main data (s400/Lund2017)
s400_GEX <-
  read.delim(
    file = "Lund2017_pN_uncentered_raw_metadata.txt.txt",
    sep = "\t",
    row.names = 1,
    skip = 14
  )
s400_header <-
  read.delim(
    file = "Lund2017_pN_uncentered_raw_metadata.txt.txt",
    sep = "\t",
    nrows = 1,
    row.names = 1
  )
s400_Metadata <-
  read.delim(
    file = "Lund2017_pN_uncentered_raw_metadata.txt.txt",
    sep = "\t",
    row.names = 1,
    nrows = 14
  )

colnames(s400_GEX) <- colnames(s400_Metadata)
table(rownames(TCGA_Filtering_NoExp) %in% rownames(s400_GEX))
s400_GEX <- s400_GEX[rownames(TCGA_Filtering_NoExp),]
s400_Original_rownames <- rownames(s400_GEX)
s400_RF_rownames <-
  gsub("-", "_", s400_Original_rownames) #Change gene symbol dashes to underscore for ranger compatibility
rownames(s400_GEX) <- s400_RF_rownames

s400_Metadata <- as.data.frame(t(s400_Metadata))

s400_Metadata$pN_YN <- s400_Metadata$pN
s400_Metadata$pN_YN[s400_Metadata$pN_YN == 0] <- "NO"
s400_Metadata$pN_YN[s400_Metadata$pN_YN == 1] <- "YES"


#Create a scaled version (used when applying TCGA models to the Lund2017 cohort)
s400_GEX_scaled <- t(apply(s400_GEX, 1, scale))
colnames(s400_GEX_scaled) <- colnames(s400_GEX)
s400x_scaled <- t(s400_GEX_scaled)




#Define the TCGA training Data (x) and Node status (y) // Lund2017 training Data (s400x) and Node status (s400y) variables for downstream analysis
x <- as.matrix(t(TCGA_GEX))
y <- factor(Metadata$pN_status_YN)

s400x <- as.matrix(t(s400_GEX))
s400y <- factor(s400_Metadata$pN_YN)



##########Build predictors
#Note, seed is not set. Repeated runs will be highly similar but match exactly to reported results.

### glmnet - TCGA
res.rtx_TCGA_all <- lapply(1:14, function(SigN) {
  nestcv.glmnet(
    y = y,
    x = log2(x[, Signatures_Filtered_split[[SigN]]] + 1) ,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )
})
names(res.rtx_TCGA_all) <- names(Signatures_Filtered_split)


##TCGA 10 times 5x cross-validation
# Create fixed folds
# folds <- repeatfolds(y, repeats = 10, n_outer_folds = 5)
#save(folds,file="folds.Rdata")
load("folds.Rdata") #load the saved fixed fold

#Set number of utilized cores (cv.cores =12) appropriate to user
repcv_fixed <- list()
for (SigN in c(1:14)) {
  repcv_fixed[[SigN]] <-
    nestcv.glmnet(
      y = y,
      x = log2(x[, Signatures_Filtered_split[[SigN]]] + 1),
      family = "binomial",
      n_outer_folds = 5,
      n_inner_folds = 5,
      outer_train_predict = TRUE,
      cv.cores = 12,
      alphaSet = seq(0, 1, 0.05)
    ) |> repeatcv(10, repeat_folds = folds, rep.cores = 1)
}
names(repcv_fixed) <- names(Signatures_Filtered_split)

#Calculate mean and 95% CI
repcv_fixed_summary <-
  do.call("rbind", lapply(1:14, function(i) {
    cbind(as.data.frame(broom::tidy(t.test(
      repcv_fixed[[i]]$result[, 1]
    )))[c(1, 5, 6)],
    as.data.frame(broom::tidy(t.test(
      repcv_fixed[[i]]$result[, 3]
    )))[c(1, 5, 6)])
  }))

rownames(repcv_fixed_summary) <- names(Signatures_Filtered_split)
repcv_fixed_summary <-
  cbind(do.call("c", lapply(1:14, function(i) {
    signif(pROC::auc(res.rtx_TCGA_all[[i]]$roc), 3)
  })), repcv_fixed_summary)
colnames(repcv_fixed_summary) <-
  c(
    "NestedAUC",
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )
#repcv_fixed

write.table(
  round(repcv_fixed_summary, 3),
  file = "TCGA_10xCV_repcv_fixed_summary.txt",
  sep = "\t",
  quote = F
)


#Train Random Forest (ranger) models with 5x5 cross-validation approach
#TCGA RF 5x CV
{
  AllDataCVx5x5 <- lapply(1:14, function(SigN) {
    DataTrain <- log2(x[, Signatures_Filtered_split[[SigN]]] + 1)
    tgrid <- expand.grid(
      .mtry = ceiling(sqrt(ncol(DataTrain))),
      .splitrule = "gini",
      .min.node.size = c(1)
    )
    
    
    model_caret <- train(
      x = DataTrain,
      y = y,
      method = "ranger",
      trControl = trainControl(
        method = "repeatedcv",
        number = 5,
        repeats = 5,
        verboseIter = T,
        classProbs = T,
        savePred = T
      ),
      tuneGrid = tgrid,
      num.trees = 1000,
      importance = 'impurity'
    )
    return(model_caret)
  })
  
  
  names(AllDataCVx5x5) <- names(Signatures_Filtered_split)
  
  #Summarize all CV results for mean and 95% CI (AUC and Balanced Accuracy)
  bob <- do.call("rbind", lapply(names(AllDataCVx5x5), function(x) {
    cbind(AllDataCVx5x5[[x]]$resample, x)
  }))
  rownames(bob) <- paste0(bob$x, "_", bob$Resample)
  
  
  bab <- lapply(names(AllDataCVx5x5), function(x) {
    bab1 <-
      split(AllDataCVx5x5[[x]]$pred, AllDataCVx5x5[[x]]$pred$Resample)
    X4 <- do.call("rbind", lapply(bab1, function(x2) {
      x3 <-
        caret::confusionMatrix(
          data = factor(x2$pred, levels = c("NO", "YES")),
          reference = factor(x2$obs, levels = c("NO", "YES")),
          mode = "everything"
        )
      return(c(x3$overall, x3$byClass, x))
    }))
    return(X4)
  })
  bab2 <- do.call("rbind", bab)
  rownames(bab2) <- paste0(bab2[, 19], "_", rownames(bab2))
  bab2 <- bab2[,-19]
  class(bab2) <- "numeric"
  bab2 <- as.data.frame(bab2)
  bab2 <- bab2[rownames(bob),]
  bab2$x <- bob$x
  
  TCGA_RF_cv <- bab2
}

#Get average and 95% CI for Balanced Accuracy
RF_ttest_TCGA <-
  round(do.call("rbind", lapply(split(TCGA_RF_cv, TCGA_RF_cv$x), function(x) {
    return(as.data.frame(broom::tidy(t.test(
      x$'Balanced Accuracy'
    )))[c(1, 5, 6)])
  })), 3)
colnames(RF_ttest_TCGA) <-
  c("5x5_RF_CV_BalancedAccuracy", "CI_lower", "CI_upper")

# write.table(RF_ttest_TCGA,
#             file = "TCGA_RF_BalancedAccuracy.txt",
#             sep = "\t",
#             quote = F)


#Get average and 95% CI for AUC
TCGA_RF_CVauc <- do.call("rbind", lapply(1:14, function(k) {
  print(k)
  tmp <-
    split(AllDataCVx5x5[[k]]$pred, AllDataCVx5x5[[k]]$pred$Resample)
  as.data.frame(broom::tidy(t.test(do.call(
    "c", lapply(1:length(tmp), function(l) {
      #print(paste(k,l))
      pROC::roc(tmp[[l]]$obs, tmp[[l]]$YES)$auc
    })
  ))))[c(1, 5, 6)]
}))
rownames(TCGA_RF_CVauc) <- names(Signatures_Filtered_split)
colnames(TCGA_RF_CVauc) <-
  c("5x5_RF_CV_AUC", "CI_lower", "CI_upper")

# write.table(
#   round(TCGA_RF_CVauc, 3),
#   file = "TCGA_RF_CV_AUC.txt",
#   sep = "\t",
#   quote = F
# )

write.table(round(cbind(TCGA_RF_CVauc, RF_ttest_TCGA), 3),
            file = "TCGA_5x5_CV_RF_summary.txt",
            sep = "\t",
            quote = F)



################################################################################
################################################################################
################################################################################

###Main Function s400
res.rtx_s400_all <- lapply(1:14, function(SigN) {
  nestcv.glmnet(
    y = s400y,
    x = s400x[, Signatures_Filtered_split[[SigN]]] ,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 0.95, 0.05)
  )
})


##s400  10 times 5x cross-validation
# CV folds s400 Fixed
#folds_s400 <- repeatfolds(s400y, repeats = 10, n_outer_folds = 5)
#save(folds_s400,file="folds_s400.Rdata")
load("folds_s400.Rdata") #load the saved fixed fold

repcv_s400_fixed <- list()

#Set number of utilized cores (cv.cores =12) appropriate to user
for (SigN in c(1:14)) {
  repcv_s400_fixed[[SigN]] <-
    nestcv.glmnet(
      y = s400y,
      x = s400x[, Signatures_Filtered_split[[SigN]]],
      family = "binomial",
      n_outer_folds = 5,
      n_inner_folds = 5,
      outer_train_predict = TRUE,
      cv.cores = 12,
      alphaSet = seq(0, 1, 0.05)
    ) |> repeatcv(10, repeat_folds = folds_s400, rep.cores = 1)
}
names(repcv_s400_fixed) <- names(Signatures_Filtered_split)
repcv_s400_fixed_backup <- repcv_s400_fixed

#note!! #some signatures performed so poorly that regularization shrink all coefficients to zero. If rerunning, check which models that fail and exclude these!
#In the publication run, signatures 5,7,9 ("Lu_2023","Luo_2021","Mitra_2014") did not yield viable models.
repcv_s400_fixed_summary <-
  do.call("rbind", lapply(c(1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 14), function(i) {
    cbind(as.data.frame(broom::tidy(t.test(
      repcv_s400_fixed[[i]]$result[, 1]
    )))[c(1, 5, 6)],
    as.data.frame(broom::tidy(t.test(
      repcv_s400_fixed[[i]]$result[, 3]
    )))[c(1, 5, 6)])
  }))

#5,7,9

rownames(repcv_s400_fixed_summary) <-
  names(Signatures_Filtered_split)[c(1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 14)]
names(Signatures_Filtered_split)[c(5, 7, 9)]
repcv_s400_fixed_summary <-
  cbind(do.call("c", lapply(c(1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 14), function(i) {
    signif(pROC::auc(res.rtx_s400_all[[i]]$roc), 3)
  })), repcv_s400_fixed_summary)
colnames(repcv_s400_fixed_summary) <-
  c(
    "NestedAUC",
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )
repcv_s400_fixed

write.table(
  round(repcv_s400_fixed_summary, 3),
  file = "Lund2017_10xCV_repcv_fixed_summary.txt",
  sep = "\t",
  quote = F
)


#Train Random Forest (ranger) models with 5x5 cross-validation approach
#s400 RF 5x5 CV
{
  AllDataCVx5x5_s400 <- lapply(1:14, function(SigN) {
    DataTrain <- s400x[, Signatures_Filtered_split[[SigN]]]
    tgrid <- expand.grid(
      .mtry = ceiling(sqrt(ncol(DataTrain))),
      .splitrule = "gini",
      .min.node.size = c(1)
    )
    
    
    model_caret <- train(
      x = DataTrain,
      y = s400y,
      method = "ranger",
      trControl = trainControl(
        method = "repeatedcv",
        number = 5,
        repeats = 5,
        verboseIter = T,
        classProbs = T,
        savePred = T
      ),
      tuneGrid = tgrid,
      num.trees = 1000,
      importance = 'impurity'
    )
    return(model_caret)
  })
  
  
  names(AllDataCVx5x5_s400) <- names(Signatures_Filtered_split)
  
  #Summarize all CV results
  bob <-
    do.call("rbind", lapply(names(AllDataCVx5x5_s400), function(x) {
      cbind(AllDataCVx5x5_s400[[x]]$resample, x)
    }))
  rownames(bob) <- paste0(bob$x, "_", bob$Resample)
  
  
  bab <- lapply(names(AllDataCVx5x5_s400), function(x) {
    bab1 <-
      split(AllDataCVx5x5_s400[[x]]$pred,
            AllDataCVx5x5_s400[[x]]$pred$Resample)
    X4 <- do.call("rbind", lapply(bab1, function(x2) {
      x3 <-
        caret::confusionMatrix(
          data = factor(x2$pred, levels = c("NO", "YES")),
          reference = factor(x2$obs, levels = c("NO", "YES")),
          mode = "everything"
        )
      return(c(x3$overall, x3$byClass, x))
    }))
    return(X4)
  })
  bab2 <- do.call("rbind", bab)
  rownames(bab2) <- paste0(bab2[, 19], "_", rownames(bab2))
  bab2 <- bab2[,-19]
  class(bab2) <- "numeric"
  bab2 <- as.data.frame(bab2)
  bab2 <- bab2[rownames(bob),]
  bab2$x <- bob$x
  
  s400_RF_cv <- bab2
}

#Get average and 95% CI for Balanced Accuracy
RF_ttest_s400 <-
  round(do.call("rbind", lapply(split(s400_RF_cv, s400_RF_cv$x), function(x) {
    return(as.data.frame(broom::tidy(t.test(
      x$'Balanced Accuracy'
    )))[c(1, 5, 6)])
  })), 3)
colnames(RF_ttest_s400) <-
  c("5x5_RF_CV_BalancedAccuracy", "CI_lower", "CI_upper")

# write.table(RF_ttest_s400,
#             file = "Lund2017_RF_CV_BalancedAccuracy.txt",
#             sep = "\t",
#             quote = F)


#Get average and 95% CI for AUC
s400_RF_CVauc <- do.call("rbind", lapply(1:14, function(k) {
  print(k)
  tmp <-
    split(AllDataCVx5x5_s400[[k]]$pred,
          AllDataCVx5x5_s400[[k]]$pred$Resample)
  as.data.frame(broom::tidy(t.test(do.call(
    "c", lapply(1:length(tmp), function(l) {
      #print(paste(k,l))
      pROC::roc(tmp[[l]]$obs, tmp[[l]]$YES)$auc
    })
  ))))[c(1, 5, 6)]
}))
rownames(s400_RF_CVauc) <- names(Signatures_Filtered_split)
colnames(s400_RF_CVauc) <-
  c("5x5_RF_CV_AUC", "CI_lower", "CI_upper")

# write.table(
#   round(s400_RF_CVauc, 3),
#   file = "Lund2017_RF_CV_AUC.txt",
#   sep = "\t",
#   quote = F
# )

write.table(round(cbind(s400_RF_CVauc, RF_ttest_s400), 3),
            file = "Lund2017_5x5_CV_RF_summary.txt",
            sep = "\t",
            quote = F)
#############################################

#some signatures performed so poorly that regularization shrink all coefficients to zero. If rerunning code, check which models that fail and exclude these!
#TCGA
which(do.call("c",lapply(1:14,function(i){is.null(nrow(res.rtx_TCGA_all[[i]]$final_coef))})))
names(Signatures_Filtered_split)[which(do.call("c",lapply(1:14,function(i){is.null(nrow(res.rtx_TCGA_all[[i]]$final_coef))})))]
#s400 (Lund2017)
which(do.call("c",lapply(1:14,function(i){is.null(nrow(res.rtx_s400_all[[i]]$final_coef))})))
names(Signatures_Filtered_split)[which(do.call("c",lapply(1:14,function(i){is.null(nrow(res.rtx_s400_all[[i]]$final_coef))})))]


#Plot Final Model AUC, and TCGA-models applied to Lund2017 data
pdf("AUC_FinalModelInTCGA_vs_Lund2017.pdf",
    width = 30,
    height = 15)
par(mfrow = c(4, 8),
    mar = c(4, 4, 4, 4) + .1,
    pty = "s")
for (i in c(1:14)) {
  print(i)
  #Exclude the models where all variables were rejected
    if (i %in% c(9)) { 
    plot.new()
    plot.new()
    next
  } 
  boxplot(
    scale(colSums(
      t(x_scaled)[rownames(res.rtx_TCGA_all[[i]]$final_coef)[-1], , drop = FALSE] *
        res.rtx_TCGA_all[[i]]$final_coef[-c(1), 1]
    )) ~ y,
    main = paste(names(Signatures_Filtered_split)[i], "\n", "Training: TCGA"),
    ylab = "Calculated Score: Sum of Exp*Coef",
    xlab = "LN positive"
  )
  legend("bottomright",
         legend = paste0("AUC = ", round(signif(
           pROC::auc(y, colSums(
             t(x_scaled)[rownames(res.rtx_TCGA_all[[i]]$final_coef)[-1], , drop = FALSE] *
               res.rtx_TCGA_all[[i]]$final_coef[-c(1), 1]
           ))
         ), 3)),
         bty = 'n')
  boxplot(
    scale(colSums(
      t(s400x_scaled)[rownames(res.rtx_TCGA_all[[i]]$final_coef)[-1], , drop =
                        FALSE] * res.rtx_TCGA_all[[i]]$final_coef[-c(1), 1]
    )) ~ s400y,
    main = paste(names(Signatures_Filtered_split)[i], "\n", "Test: Lund2017"),
    ylab = "Calculated Score: Sum of Exp*Coef",
    xlab = "LN positive"
  )
  legend("bottomright",
         legend = paste0("AUC = ", round(signif(
           pROC::auc(s400y, colSums(
             t(s400x_scaled)[rownames(res.rtx_TCGA_all[[i]]$final_coef)[-1], , drop =
                               FALSE] * res.rtx_TCGA_all[[i]]$final_coef[-c(1), 1]
           ))
         ), 3)),
         bty = 'n')
  
}
dev.off()







##################################
#Train models on the TCGA dataset, using either internal glmnet feature selection, or with preselection of top 100 features though t.test, wilcox.test, or Random Forest (ranger) feature importance. Feature selection is nested (i.e., only performed in the current CV train fold).

#no filter
res.rtx_TCGA_PC_noFilter <-
  nestcv.glmnet(
    y,
    log2(x + 1) ,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )

summary(res.rtx_TCGA_PC_noFilter)

repcv_res.rtx_TCGA_PC_noFilter <-
  nestcv.glmnet(
    y = y,
    x = log2(x + 1),
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    family = "binomial",
    alphaSet = seq(0, 1, 0.05),
    cv.cores = 12
  ) |> repeatcv(10, repeat_folds = folds, rep.cores = 1)

repcv_res.rtx_TCGA_PC_noFilter_summary <-
  cbind(as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_TCGA_PC_noFilter$result[, 1])
  ))[c(1, 5, 6)], as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_TCGA_PC_noFilter$result[, 3])
  ))[c(1, 5, 6)])

names(repcv_res.rtx_TCGA_PC_noFilter_summary) <-
  c(
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )

#ttest
res.rtx_TCGA_PC_ttest <-
  nestcv.glmnet(
    y,
    log2(x + 1) ,
    filterFUN = ttest_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )
summary(res.rtx_TCGA_PC_ttest)
repcv_res.rtx_TCGA_PC_ttest <-
  nestcv.glmnet(
    y = y,
    x = log2(x + 1),
    filterFUN = ttest_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    alphaSet = seq(0, 1, 0.05),
    cv.cores = 12
  ) |> repeatcv(10, repeat_folds = folds, rep.cores = 1)
repcv_res.rtx_TCGA_PC_ttest_summary <-
  cbind(as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_TCGA_PC_ttest$result[, 1])
  ))[c(1, 5, 6)], as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_TCGA_PC_ttest$result[, 3])
  ))[c(1, 5, 6)])

names(repcv_res.rtx_TCGA_PC_ttest_summary) <-
  c(
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )

#wilcox.test
res.rtx_TCGA_PC_wilcox <-
  nestcv.glmnet(
    y,
    log2(x + 1) ,
    filterFUN = wilcoxon_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )

summary(res.rtx_TCGA_PC_wilcox)

repcv_res.rtx_TCGA_PC_wilcox <-
  nestcv.glmnet(
    y = y,
    x = log2(x + 1),
    filterFUN = wilcoxon_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    alphaSet = seq(0, 1, 0.05),
    cv.cores = 12
  ) |> repeatcv(10, repeat_folds = folds, rep.cores = 1)

repcv_res.rtx_TCGA_PC_wilcox_summary <-
  cbind(as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_TCGA_PC_wilcox$result[, 1])
  ))[c(1, 5, 6)], as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_TCGA_PC_wilcox$result[, 3])
  ))[c(1, 5, 6)])

names(repcv_res.rtx_TCGA_PC_wilcox_summary) <-
  c(
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )

#ranger_filter
res.rtx_TCGA_PC_ranger <-
  nestcv.glmnet(
    y,
    log2(x + 1) ,
    filterFUN = ranger_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )

summary(res.rtx_TCGA_PC_ranger)

repcv_res.rtx_TCGA_PC_ranger <-
  nestcv.glmnet(
    y = y,
    x = log2(x + 1),
    filterFUN = ranger_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    alphaSet = seq(0, 1, 0.05),
    cv.cores = 12
  ) |> repeatcv(10, repeat_folds = folds, rep.cores = 1)

repcv_res.rtx_TCGA_PC_ranger_summary <-
  cbind(as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_TCGA_PC_ranger$result[, 1])
  ))[c(1, 5, 6)], as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_TCGA_PC_ranger$result[, 3])
  ))[c(1, 5, 6)])

names(repcv_res.rtx_TCGA_PC_ranger_summary) <-
  c(
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )



##################################################
#Train models on the Lund2017 dataset, using either internal glmnet feature selection, or with preselection of top 100 features though t.test, wilcox.test, or Random Forest (ranger) feature importance. Feature selection is nested (i.e., only performed in the current CV train fold).

#Lund2017
#no filter
res.rtx_s400_PC_noFilter <-
  nestcv.glmnet(
    y = s400y,
    x = s400x ,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )

summary(res.rtx_s400_PC_noFilter)

repcv_res.rtx_s400_PC_noFilter <-
  nestcv.glmnet(
    y = s400y,
    x = s400x,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    family = "binomial",
    alphaSet = seq(0, 1, 0.05),
    cv.cores = 12
  ) |> repeatcv(10, repeat_folds = folds_s400, rep.cores = 1)

repcv_res.rtx_s400_PC_noFilter_summary <-
  cbind(as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_s400_PC_noFilter$result[, 1])
  ))[c(1, 5, 6)], as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_s400_PC_noFilter$result[, 3])
  ))[c(1, 5, 6)])

names(repcv_res.rtx_s400_PC_noFilter_summary) <-
  c(
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )

#ttest
res.rtx_s400_PC_ttest <-
  nestcv.glmnet(
    y = s400y,
    x = s400x ,
    filterFUN = ttest_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )

summary(res.rtx_s400_PC_ttest)

repcv_res.rtx_s400_PC_ttest <-
  nestcv.glmnet(
    y = s400y,
    x = s400x,
    filterFUN = ttest_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    alphaSet = seq(0, 1, 0.05),
    cv.cores = 12
  ) |> repeatcv(10, repeat_folds = folds_s400, rep.cores = 1)

repcv_res.rtx_s400_PC_ttest_summary <-
  cbind(as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_s400_PC_ttest$result[, 1])
  ))[c(1, 5, 6)], as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_s400_PC_ttest$result[, 3])
  ))[c(1, 5, 6)])

names(repcv_res.rtx_s400_PC_ttest_summary) <-
  c(
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )

#wilcox.test
res.rtx_s400_PC_wilcox <-
  nestcv.glmnet(
    y = s400y,
    x = s400x ,
    filterFUN = wilcoxon_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )

summary(res.rtx_s400_PC_wilcox)

repcv_res.rtx_s400_PC_wilcox <-
  nestcv.glmnet(
    y = s400y,
    x = s400x,
    filterFUN = wilcoxon_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100, p_cutoff = NULL),
    family = "binomial",
    alphaSet = seq(0, 1, 0.05),
    cv.cores = 12
  ) |> repeatcv(10, repeat_folds = folds_s400, rep.cores = 1)

repcv_res.rtx_s400_PC_wilcox_summary <-
  cbind(as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_s400_PC_wilcox$result[, 1])
  ))[c(1, 5, 6)], as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_s400_PC_wilcox$result[, 3])
  ))[c(1, 5, 6)])

names(repcv_res.rtx_s400_PC_wilcox_summary) <-
  c(
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )

#ranger_filter
res.rtx_s400_PC_ranger <-
  nestcv.glmnet(
    y = s400y,
    x = s400x ,
    filterFUN = ranger_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100),
    family = "binomial",
    cv.cores = 12,
    alphaSet = seq(0, 1, 0.05)
  )

summary(res.rtx_s400_PC_ranger)

repcv_res.rtx_s400_PC_ranger <-
  nestcv.glmnet(
    y = s400y,
    x = s400x,
    filterFUN = ranger_filter,
    n_outer_folds = 5,
    n_inner_folds = 5,
    outer_train_predict = TRUE,
    filter_options = list(nfilter = 100),
    family = "binomial",
    alphaSet = seq(0, 1, 0.05),
    cv.cores = 12
  ) |> repeatcv(10, repeat_folds = folds_s400, rep.cores = 1)

repcv_res.rtx_s400_PC_ranger_summary <-
  cbind(as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_s400_PC_ranger$result[, 1])
  ))[c(1, 5, 6)], as.data.frame(broom::tidy(
    t.test(repcv_res.rtx_s400_PC_ranger$result[, 3])
  ))[c(1, 5, 6)])

names(repcv_res.rtx_s400_PC_ranger_summary) <-
  c(
    "10xNestedAUC",
    "CI_lower",
    "CI_upper",
    "10xBalancedAccuracy",
    "CI_lower",
    "CI_upper"
  )




#Summarize 10x repeated 5x CV results for TCGA and s400 (Lund2017) datasets. This is an assessment of internal performance in training data without information leakage.
FilteredFeaturesCVresults <- round(
  rbind(
    repcv_res.rtx_TCGA_PC_noFilter_summary,
    repcv_res.rtx_TCGA_PC_ttest_summary,
    repcv_res.rtx_TCGA_PC_wilcox_summary,
    repcv_res.rtx_TCGA_PC_ranger_summary,
    
    repcv_res.rtx_s400_PC_noFilter_summary,
    repcv_res.rtx_s400_PC_ttest_summary,
    repcv_res.rtx_s400_PC_wilcox_summary,
    repcv_res.rtx_s400_PC_ranger_summary
  ),
  3
)

rownames(FilteredFeaturesCVresults) <-
  c(
    "TCGA_GLMNET_filter",
    "TCGA_t.test100",
    "TCGA_wilcoxon100",
    "TCGA_ranger100",
    "Lund2017_GLMNET_filter",
    "Lund2017_t.test100",
    "Lund2017_wilcoxon100",
    "Lund2017_ranger100"
  )

write.table(
  FilteredFeaturesCVresults,
  file = "FilteredFeaturesCVresults.txt",
  sep = "\t",
  quote = F
)




########Testing signatures
#############################################################
#############################################################
#res.rtx_TCGA_PC_noFilter
#res.rtx_TCGA_PC_ttest
#res.rtx_TCGA_PC_wilcox
#res.rtx_TCGA_PC_ranger

#Plot each model as fitted on training data (note that this fit and AUC-value is not a reflection of prediction performance in the training set (which is assessed in FilteredFeaturesCVresults above) , but only illustrates how the selected variables separate the classes in the training set). Application to the test dataset (TCGA-> Lund2017, and Lund2017->TCGA) IS a prediction performance evaluation, which  shows poor results across all models.
#Top row: panel pairs of: TCGA (training fit) :: Lund2017 (test)
#Bottom row: panel pairs of: Lund2017 (training fit) :: TCGA (test)
pdf(
  "InternalFeatureSelection_InTCGA_vs_Lund2017.pdf",
  width = 30,
  height = 15
)
par(mfrow = c(4, 8),
    mar = c(4, 4, 4, 4) + .1,
    pty = "s")

boxplot(
  scale(colSums(
    t(x_scaled)[rownames(res.rtx_TCGA_PC_noFilter$final_coef)[-1], , drop =
                  FALSE] * res.rtx_TCGA_PC_noFilter$final_coef[-c(1), 1]
  )) ~ y,
  main = paste("GLMNET filter", "\n", "Training: TCGA"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(y, colSums(
           t(x_scaled)[rownames(res.rtx_TCGA_PC_noFilter$final_coef)[-1], , drop =
                         FALSE] * res.rtx_TCGA_PC_noFilter$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')
boxplot(
  scale(colSums(
    t(s400x_scaled)[rownames(res.rtx_TCGA_PC_noFilter$final_coef)[-1], , drop =
                      FALSE] * res.rtx_TCGA_PC_noFilter$final_coef[-c(1), 1]
  )) ~ s400y,
  main = paste("GLMNET filter", "\n", "Test: Lund2017"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(s400y, colSums(
           t(s400x_scaled)[rownames(res.rtx_TCGA_PC_noFilter$final_coef)[-1], , drop =
                             FALSE] * res.rtx_TCGA_PC_noFilter$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')

boxplot(
  scale(colSums(
    t(x_scaled)[rownames(res.rtx_TCGA_PC_ttest$final_coef)[-1], , drop = FALSE] *
      res.rtx_TCGA_PC_ttest$final_coef[-c(1), 1]
  )) ~ y,
  main = paste("t.test 100", "\n", "Training: TCGA"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(y, colSums(
           t(x_scaled)[rownames(res.rtx_TCGA_PC_ttest$final_coef)[-1], , drop = FALSE] *
             res.rtx_TCGA_PC_ttest$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')
boxplot(
  scale(colSums(
    t(s400x_scaled)[rownames(res.rtx_TCGA_PC_ttest$final_coef)[-1], , drop =
                      FALSE] * res.rtx_TCGA_PC_ttest$final_coef[-c(1), 1]
  )) ~ s400y,
  main = paste("t.test 100", "\n", "Test: Lund2017"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(s400y, colSums(
           t(s400x_scaled)[rownames(res.rtx_TCGA_PC_ttest$final_coef)[-1], , drop =
                             FALSE] * res.rtx_TCGA_PC_ttest$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')

boxplot(
  scale(colSums(
    t(x_scaled)[rownames(res.rtx_TCGA_PC_wilcox$final_coef)[-1], , drop = FALSE] *
      res.rtx_TCGA_PC_wilcox$final_coef[-c(1), 1]
  )) ~ y,
  main = paste("wilcoxon 100", "\n", "Training: TCGA"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(y, colSums(
           t(x_scaled)[rownames(res.rtx_TCGA_PC_wilcox$final_coef)[-1], , drop = FALSE] *
             res.rtx_TCGA_PC_wilcox$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')
boxplot(
  scale(colSums(
    t(s400x_scaled)[rownames(res.rtx_TCGA_PC_wilcox$final_coef)[-1], , drop =
                      FALSE] * res.rtx_TCGA_PC_wilcox$final_coef[-c(1), 1]
  )) ~ s400y,
  main = paste("wilcoxon 100", "\n", "Test: Lund2017"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(s400y, colSums(
           t(s400x_scaled)[rownames(res.rtx_TCGA_PC_wilcox$final_coef)[-1], , drop =
                             FALSE] * res.rtx_TCGA_PC_wilcox$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')

boxplot(
  scale(colSums(
    t(x_scaled)[rownames(res.rtx_TCGA_PC_ranger$final_coef)[-1], , drop = FALSE] *
      res.rtx_TCGA_PC_ranger$final_coef[-c(1), 1]
  )) ~ y,
  main = paste("ranger 100", "\n", "Training: TCGA"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(y, colSums(
           t(x_scaled)[rownames(res.rtx_TCGA_PC_ranger$final_coef)[-1], , drop = FALSE] *
             res.rtx_TCGA_PC_ranger$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')
boxplot(
  scale(colSums(
    t(s400x_scaled)[rownames(res.rtx_TCGA_PC_ranger$final_coef)[-1], , drop =
                      FALSE] * res.rtx_TCGA_PC_ranger$final_coef[-c(1), 1]
  )) ~ s400y,
  main = paste("ranger 100", "\n", "Test: Lund2017"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(s400y, colSums(
           t(s400x_scaled)[rownames(res.rtx_TCGA_PC_ranger$final_coef)[-1], , drop =
                             FALSE] * res.rtx_TCGA_PC_ranger$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')

boxplot(
  scale(colSums(
    t(s400x_scaled)[rownames(res.rtx_s400_PC_noFilter$final_coef)[-1], , drop =
                      FALSE] * res.rtx_s400_PC_noFilter$final_coef[-c(1), 1]
  )) ~ s400y,
  main = paste("GLMNET filter", "\n", "Training: Lund2017"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(s400y, colSums(
           t(s400x_scaled)[rownames(res.rtx_s400_PC_noFilter$final_coef)[-1], , drop =
                             FALSE] * res.rtx_s400_PC_noFilter$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')
boxplot(
  scale(colSums(
    t(x_scaled)[rownames(res.rtx_s400_PC_noFilter$final_coef)[-1], , drop =
                  FALSE] * res.rtx_s400_PC_noFilter$final_coef[-c(1), 1]
  )) ~ y,
  main = paste("GLMNET filter", "\n", "Test: TCGA"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(y, colSums(
           t(x_scaled)[rownames(res.rtx_s400_PC_noFilter$final_coef)[-1], , drop =
                         FALSE] * res.rtx_s400_PC_noFilter$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')


boxplot(
  scale(colSums(
    t(s400x_scaled)[rownames(res.rtx_s400_PC_ttest$final_coef)[-1], , drop =
                      FALSE] * res.rtx_s400_PC_ttest$final_coef[-c(1), 1]
  )) ~ s400y,
  main = paste("t.test 100", "\n", "Training: Lund2017"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(s400y, colSums(
           t(s400x_scaled)[rownames(res.rtx_s400_PC_ttest$final_coef)[-1], , drop =
                             FALSE] * res.rtx_s400_PC_ttest$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')
boxplot(
  scale(colSums(
    t(x_scaled)[rownames(res.rtx_s400_PC_ttest$final_coef)[-1], , drop = FALSE] *
      res.rtx_s400_PC_ttest$final_coef[-c(1), 1]
  )) ~ y,
  main = paste("t.test 100", "\n", "Test: TCGA"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(y, colSums(
           t(x_scaled)[rownames(res.rtx_s400_PC_ttest$final_coef)[-1], , drop = FALSE] *
             res.rtx_s400_PC_ttest$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')

boxplot(
  scale(colSums(
    t(s400x_scaled)[rownames(res.rtx_s400_PC_wilcox$final_coef)[-1], , drop =
                      FALSE] * res.rtx_s400_PC_wilcox$final_coef[-c(1), 1]
  )) ~ s400y,
  main = paste("wilcoxon 100", "\n", "Training: Lund2017"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(s400y, colSums(
           t(s400x_scaled)[rownames(res.rtx_s400_PC_wilcox$final_coef)[-1], , drop =
                             FALSE] * res.rtx_s400_PC_wilcox$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')
boxplot(
  scale(colSums(
    t(x_scaled)[rownames(res.rtx_s400_PC_wilcox$final_coef)[-1], , drop = FALSE] *
      res.rtx_s400_PC_wilcox$final_coef[-c(1), 1]
  )) ~ y,
  main = paste("wilcoxon 100", "\n", "Test: TCGA"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(y, colSums(
           t(x_scaled)[rownames(res.rtx_s400_PC_wilcox$final_coef)[-1], , drop = FALSE] *
             res.rtx_s400_PC_wilcox$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')

boxplot(
  scale(colSums(
    t(s400x_scaled)[rownames(res.rtx_s400_PC_ranger$final_coef)[-1], , drop =
                      FALSE] * res.rtx_s400_PC_ranger$final_coef[-c(1), 1]
  )) ~ s400y,
  main = paste("ranger 100", "\n", "Training: Lund2017"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(s400y, colSums(
           t(s400x_scaled)[rownames(res.rtx_s400_PC_ranger$final_coef)[-1], , drop =
                             FALSE] * res.rtx_s400_PC_ranger$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')
boxplot(
  scale(colSums(
    t(x_scaled)[rownames(res.rtx_s400_PC_ranger$final_coef)[-1], , drop = FALSE] *
      res.rtx_s400_PC_ranger$final_coef[-c(1), 1]
  )) ~ y,
  main = paste("ranger 100", "\n", "Test: TCGA"),
  ylab = "Calculated Score: Sum of Exp*Coef",
  xlab = "LN positive"
)
legend("bottomright",
       legend = paste0("AUC = ", round(signif(
         pROC::auc(y, colSums(
           t(x_scaled)[rownames(res.rtx_s400_PC_ranger$final_coef)[-1], , drop = FALSE] *
             res.rtx_s400_PC_ranger$final_coef[-c(1), 1]
         ))
       ), 3)),
       bty = 'n')


dev.off()


#Save the identified genes and coefficients
write.table(
  res.rtx_TCGA_PC_noFilter$final_coef[-1, 1, drop = FALSE],
  file = "TCGA_GLMNET_filter_genes.txt",
  sep = "\t",
  quote = F
)
write.table(
  res.rtx_TCGA_PC_ttest$final_coef[-1, 1, drop = FALSE],
  file = "TCGA_t.test_filter_genes.txt",
  sep = "\t",
  quote = F
)
write.table(
  res.rtx_TCGA_PC_wilcox$final_coef[-1, 1, drop = FALSE],
  file = "TCGA_wilcoxon_filter_genes.txt",
  sep = "\t",
  quote = F
)
write.table(
  res.rtx_TCGA_PC_ranger$final_coef[-1, 1, drop = FALSE],
  file = "TCGA_ranger_filter_genes.txt",
  sep = "\t",
  quote = F
)

write.table(
  res.rtx_s400_PC_noFilter$final_coef[-1, 1, drop = FALSE],
  file = "Lund2017_GLMNET_filter_genes.txt",
  sep = "\t",
  quote = F
)
write.table(
  res.rtx_s400_PC_ttest$final_coef[-1, 1, drop = FALSE],
  file = "Lund2017_t.test_filter_genes.txt",
  sep = "\t",
  quote = F
)
write.table(
  res.rtx_s400_PC_wilcox$final_coef[-1, 1, drop = FALSE],
  file = "Lund2017_wilcoxon_filter_genes.txt",
  sep = "\t",
  quote = F
)
write.table(
  res.rtx_s400_PC_ranger$final_coef[-1, 1, drop = FALSE],
  file = "Lund2017_ranger_filter_genes.txt",
  sep = "\t",
  quote = F
)
