# Johansson et al. 2025 Node-predictor Evaluation
R code and data for evaluating RNA-based predictors of lymph node status in muscle-invasive bladder cancer (MIBC). Includes validation of 12 published predictors and development of new models using gene expression data from two cohorts, revealing limitations in predicting nodal spread with bulk tumor RNA profiles.

## Background
The presence of cancer in pelvic lymph nodes removed during radical surgery for muscleinvasive
bladder cancer (MIBC) is a key determinant of patient outcome. It would be beneficial to
predict node status preoperatively to tailor the use of neoadjuvant chemotherapy and extent of lymph
node dissection. Of 12 published node status predictors based on tumor RNA expression signatures,
none have been successfully validated in subsequent reports.

## Objective 
We aimed to validate all published node status predictors and evaluate new prediction
models in MIBC.

## Methods
Gene expression data and node status from two MIBC cohorts were used to test 12 published
node-predictive signatures. The overlap in differential expression was examined across the two datasets,
and new prediction models were tested in cross-validation and by application to the independent cohort.
Results: Published node status predictors performed either no better, or only slightly better than chance
in the two independent validation datasets (maximum AUC 0.59 and 0.65, and maximum balanced
accuracy 0.54 and 0.57). Among very few genes and signatures differentially expressed in the same
direction in both data sets we identified upregulation of interferon-response signatures in node negative
cases. Transcriptomic predictors trained in one dataset performed poorly when applied to the
independent dataset (AUC 0.60-0.62).

## Conclusions
In this systematic evaluation, neither the 12 published signatures nor our own models
reached an adequate performance for clinical node status prediction in independent data. This indicates
that the biological determinants of nodal spread are poorly captured by bulk tumor RNA expression
profiles.