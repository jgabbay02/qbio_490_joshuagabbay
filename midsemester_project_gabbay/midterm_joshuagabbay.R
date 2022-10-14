
dir.create("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs")

knitr::opts_knit$set(root.dir = normalizePath("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs")) 

if (!require(survival)){ 
  install.packages("survival")
}
library(survival)

if (!require(survminer)){
  install.packages("survminer")
}
library(survminer)

if (!require(ggplot2)){ 
  install.packages("ggplot2")
}
library(ggplot2)

if (!require(TCGAbiolinks)){ 
  install.packages("TCGAbiolinks")
}
library(TCGAbiolinks)

if (!require(BiocManager)){
  install.packages("BiocManager")
}
library(BiocManager)
BiocManager::install("maftools")
library(maftools)

#GDC on clinical data
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml" )
GDCdownload(clinical_query)
clinical <- GDCprepare_clinic(clinical_query, clinical.info = "patient")

#GDC to prepare clinical radiation data
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")

#Creating a her2 status column in clinical, then cleaning out all NA values and all "borderline" values
clinical$her2_status <- ifelse(clinical$her2_immunohistochemistry_level_result == "1+" | clinical$her2_immunohistochemistry_level_result == "0", "negative", ifelse(clinical$her2_immunohistochemistry_level_result == "2+", "borderline", "positive"))
her2_na_mask <- ifelse(clinical$her2_immunohistochemistry_level_result == "", F, T)
her2_cleaned_clinical <- clinical[her2_na_mask, ]
her2_borderline_mask <- ifelse(her2_cleaned_clinical$her2_status == "borderline", F, T)
her2_cleaned_clinical <- her2_cleaned_clinical[her2_borderline_mask, ]
her2_mask <- ifelse(her2_cleaned_clinical$her2_status == "negative", F, T)

#use her2_mask to subset for only negative patients
her2_negative_clinical <- her2_cleaned_clinical[!her2_mask, ]

#use her2_mask to subset for only positive patients
her2_positive_clinical <- her2_cleaned_clinical[her2_mask, ]

#making a survival time column for survival plots
her2_cleaned_clinical$survival_time <- ifelse(is.na(her2_cleaned_clinical$days_to_death), her2_cleaned_clinical$survival_time <- her2_cleaned_clinical$days_to_last_followup, her2_cleaned_clinical$survival_time <- her2_cleaned_clinical$days_to_death)

#making a death event (T/F) column for survival plots
her2_cleaned_clinical$death_event <- ifelse(her2_cleaned_clinical$vital_status == "Alive", her2_cleaned_clinical$death_event <- FALSE, her2_cleaned_clinical$death_event <- TRUE)

#initializing a survival object
surv_object_her2 <- Surv(time = her2_cleaned_clinical$survival_time, event = her2_cleaned_clinical$death_event)

#creating a fit object
her2_fit <- survfit(surv_object_her2 ~ her2_cleaned_clinical$her2_status, data = her2_cleaned_clinical)

#formats and creates KM plot
survplot_her2 = ggsurvplot(her2_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/her2_KM_plot.jpg")
KM_plot_her2 = survplot_her2$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_her2
dev.off()

#creating a boxplot using radiation dosage and anatomic treatment site
radiation_mask <- ifelse(clinical.rad$radiation_dosage == "", F, T)
boxplot_clinical <- clinical.rad[radiation_mask, ]
site_mask <- ifelse(clinical.rad$anatomic_treatment_site == "", F, T)
boxplot_clinical <- boxplot_clinical[site_mask, ]

jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/boxplot_clinical_radiation_dosage_v_site.jpg")
boxplot(boxplot_clinical$radiation_dosage~boxplot_clinical$anatomic_treatment_site)
dev.off()

#MAF

#changing colnames of clinical so it is readable by MAF
colnames(her2_cleaned_clinical)[ colnames(her2_cleaned_clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

#GDC on MAF data
maf_query <- GDCquery(project = "TCGA-BRCA", data.category = "Simple Nucleotide Variation", access = "open", data.type = "Masked Somatic Mutation", workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query)
maf <- GDCprepare(maf_query) 
maf_object <- read.maf(maf = maf, clinicalData = her2_cleaned_clinical, isTCGA = TRUE)

#creating an oncoplot with clinical data of her2 status
jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/her2_status_onco_top10genes.jpg")
oncoplot(maf = maf_object, top = 10, clinicalFeatures = "her2_status")
dev.off()

#storing barcodes to subset MAF data
positive_barcodes <- her2_cleaned_clinical[her2_cleaned_clinical$her2_status == "positive", 1]
positive_maf <- subsetMaf(maf = maf_object, tsb = positive_barcodes)

negative_barcodes <- her2_cleaned_clinical[her2_cleaned_clinical$her2_status == "negative", 1]
negative_maf <- subsetMaf(maf = maf_object, tsb = negative_barcodes)

#creating (and saving) our coOncoplot for HER2 positive and HER2 negative patients
jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/coOnco_HER2status.jpg")
coOncoplot(m1 = positive_maf, 
           m2 = negative_maf, 
           m1Name = 'HER2 Positive Genetic Mutations', 
           m2Name = 'HER2 Negative Genetic Mutations')
dev.off()

#creating a co-lollipop plot for the PIK3CA gene based on HER2 status
jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/pik3ca_her2status_colollipop.jpg")
lollipopPlot2(m1 = positive_maf, m2 = negative_maf, m1_name = 'HER2 Positive Patient HER2 Mutation Layout', m2_name = 'HER2 Negative Patient HER2 Mutation Layout', gene = "PIK3CA")
dev.off()

#creating and cleaning survival time and death event columns in MAF data
maf_object@clinical.data$survival_time <- ifelse(is.na(maf_object@clinical.data$days_to_death), maf_object@clinical.data$survival_time <- maf_object@clinical.data$days_to_last_followup, maf_object@clinical.data$survival_time <- maf_object@clinical.data$days_to_death)
maf_object@clinical.data$death_event <- ifelse(maf_object@clinical.data$vital_status == "Alive", maf_object@clinical.data$death_event <- FALSE, maf_object@clinical.data$death_event <- TRUE)

#graphing a MAF survival plot based on PIK3CA mutations
jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/pik3ca_mafsurvival.jpg")
mafSurvival(maf = maf_object, genes = "PIK3CA", time = "survival_time", Status = "death_event", isTCGA = TRUE)
dev.off()

#Saving 4 CSV files
write.csv(clinical, "/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/midterm_clinical.csv")
write.csv(clinical.rad, "/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/midterm_clinical_radiation.csv")
write.csv(her2_cleaned_clinical, "/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/midterm_her2_cleaned_clinical.csv")
write.csv(boxplot_clinical, "/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/midsemester_project_gabbay/outputs/midterm_boxplot_clinical_radiation.csv")
























