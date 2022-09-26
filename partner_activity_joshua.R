
setwd("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/analysis_data")

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


clinical <- read.csv("brca_clinical_data.csv")

clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")

clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")

clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")

sum(is.na(clinical.rad$anatomic_treatment_site))

surgery <- clinical[ ,c(1,34)] # creates a data frame with the patient barcodes and breast carcinoma surgery data
surgery <- surgery[surgery$breast_carcinoma_surgical_procedure_name != "", ] # takes out all empty values
surgery$breast_carcinoma_surgical_procedure_name <- as.character(surgery$breast_carcinoma_surgical_procedure_name) 
rad.region <- clinical.rad[ ,c(1,11)] # creates a data frame witht the patient barcodes and anatomic treatment site data
rad.region <- rad.region[rad.region$anatomic_treatment_site != "", ] # takes out all empty values

region.v.surgery <- rad.region # create new data frame
region.v.surgery$breast_carcinoma_surgical_procedure_name <- NA # add breast carcinoma column to this dataframe

for (i in 1:length(region.v.surgery$anatomic_treatment_site)){ # goes through each row and matches patient barcodes up between the two data sources
  surgery.entry <- surgery[which(surgery$bcr_patient_barcode == region.v.surgery$bcr_patient_barcode[i]),]
  if(nrow(surgery.entry) != 0){
    region.v.surgery[i,3] <- surgery.entry[1,2]
  }
}

region.v.surgery <- region.v.surgery[!is.na(region.v.surgery$breast_carcinoma_surgical_procedure_name),] #removes all na values
region.v.surgery$anatomic_treatment_site <- as.character(region.v.surgery$anatomic_treatment_site) 

rad.site <- c(rep("Primary Tumor Field", 4), rep("Regional site", 4), rep("Local Recurrence", 4), 
          rep("Distant Recurrence", 4), rep("Distant site", 4)) # for x value of ggplot, different categories with 4 columns each
surgery.type <- rep(c("Lumpectomy", "Simple Mastectomy", "Modified Radical Mastectomy", "Other"), 5) # 4 different types of columns

test <- data.frame(rad.site, surgery.type, count = NA) # creates a test dataframe for plot

for (i in 1: length(test$rad.site)){ # goes through each row and counts how many times each breast carcinoma surgery appears
  all.entries <- region.v.surgery[which((region.v.surgery$anatomic_treatment_site == test[i,1]) 
                                        & (region.v.surgery$breast_carcinoma_surgical_procedure_name == test[i,2])), ]
  test[i,3] <- nrow(all.entries)
}

jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/analysis_data/region_km.jpg")
ggplot(test, aes(fill=surgery.type, y=test[ ,3], x=rad.site)) + geom_bar(position ="dodge", stat="identity"); # plots the two categorical variables
dev.off() # saves plot as jpeg

# KM PLOT

clinical$survival_time <- ifelse(!is.na(clinical$days_to_death), clinical$days_to_death, clinical$days_to_last_followup) # creates a survival time column

clinical$death_event <- ifelse(as.character(clinical$vital_status) == "Alive", clinical$death_event <- FALSE, clinical$death_event <- TRUE) #creates a death event column

surv_object_surgery <- Surv(time = clinical$survival_time, event = clinical$death_event) # creates a surv object for surgery

surgery_fit <- survfit(surv_object_surgery ~ clinical$breast_carcinoma_surgical_procedure_name, data = clinical) # creates a fit object for surgery

survplot_surgery = ggsurvplot(surgery_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

KM_plot_surgery = survplot_surgery$plot + theme_bw() +  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5))

jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/analysis_data/surgery_km.jpg")
KM_plot_surgery # plots the KM plot for surgery data
dev.off() # saves as jpeg

clinical.rad$survival_time <- NA # creates an empty survival time column in clinical.rad
clinical.rad$death_event <- NA # creates an empty death event column in clinical.rad

for (i in 1: length(clinical.rad$bcr_patient_barcode)){ # fills in the survival time column and death event column for patients in clinical.rad
  exists <- clinical[which(clinical$bcr_patient_barcode == clinical.rad$bcr_patient_barcode[i]),]
  if (nrow(exists) != 0){
    clinical.rad[i, 21] <- exists[1, 116]
    clinical.rad[i, 22] <- exists[1, 117]
  }
}

surv_object_region <- Surv(time = clinical.rad$survival_time, event = clinical.rad$death_event) # creates a surv object for region

region_fit <- survfit(surv_object_region ~ clinical.rad$anatomic_treatment_site, data = clinical.rad) # creates a fit object for region

survplot_region = ggsurvplot(region_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")


KM_plot_region = survplot_region$plot + theme_bw() + theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

jpeg("/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/analysis_data/region_km.jpg")
KM_plot_region # plots the KM plot for surgery data
dev.off() # saves as jpeg

write.csv(clinical.rad, "/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/analysis_data/brca_clinical_radiation_data.csv", row.names = FALSE)
write.csv(clinical.drug, "/Users/joshuagabbay/Desktop/QBIO/qbio_490_joshuagabbay/analysis_data/brca_clinical_drug_data.csv", row.names =  FALSE)


