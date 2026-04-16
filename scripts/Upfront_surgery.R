rm(list=ls())
library(data.table)
library(survival) 
library(survminer)
library(gridtext)
library(ggplot2)
library(scales)
library(officer)
library(tidyverse)
library(gtsummary)
library(flextable)
library(vautils)
library(coxphf)
library(pheatmap)
library(rcompanion)

########################
## reporting function ##
########################
report_surv_extra <- function(KM_curve, cox_fit, times = c(24, 36, 48)) {
  
  # --- 1. Auto-Detect Reserved Keywords ---
  reserved_keywords <- c("ACT", "MDT", "Treatment", "RX", "Chemo")
  
  term_names <- if (inherits(cox_fit, "coxphf")) {
    all.vars(cox_fit$formula)[-1] 
  } else {
    attr(cox_fit$terms, "term.labels")
  }
  
  is_special_mode <- any(sapply(reserved_keywords, function(k) grepl(k, term_names[1], ignore.case = TRUE)))
  
  # --- 2. Cox Model Reporting ---
  if (inherits(cox_fit, "coxphf")) {
    hr_label <- "Firth HR"
    HR      <- exp(cox_fit$coefficients[1])
    lowerCI <- cox_fit$ci.lower[1]
    upperCI <- cox_fit$ci.upper[1]
    pval    <- cox_fit$prob[1]
  } else {
    hr_label <- "HR"
    csum    <- summary(cox_fit)
    HR      <- csum$conf.int[1, "exp(coef)"]
    lowerCI <- csum$conf.int[1, "lower .95"]
    upperCI <- csum$conf.int[1, "upper .95"]
    pval    <- csum$coefficients[1, "Pr(>|z|)"]
  }
  
  # --- 3. Conditional Inversion ---
  if (is_special_mode && HR > 1) {
    hr_label <- paste(hr_label, "(inverted)")
    old_low <- lowerCI
    HR      <- 1 / HR
    lowerCI <- 1 / upperCI
    upperCI <- 1 / old_low
  }
  
  label_text <- sprintf("%s = %.2f (%.2f - %.2f); Cox p = %.3g",
                        hr_label, HR, lowerCI, upperCI, pval)
  cat(label_text, "\n\n")
  
  # --- 4. Median Survival and Events (Updated with 95% CI) ---
  km_table <- summary(KM_curve)$table
  cat("Group Statistics (Events & Median Survival [95% CI]):\n")
  
  if (is.null(KM_curve$strata)) {
    n_total <- KM_curve$n
    n_ev    <- sum(KM_curve$n.event)
    pct_ev  <- (n_ev / n_total) * 100
    med     <- km_table["median"]
    m_low   <- km_table["0.95LCL"]
    m_up    <- km_table["0.95UCL"]
    
    med_text <- if(is.na(med)) "NA" else sprintf("%.2f (%.2f - %.2f)", med, m_low, m_up)
    
    cat(sprintf("All: Events = %.1f%% (%d/%d); Median = %s\n", 
                pct_ev, as.integer(n_ev), as.integer(n_total), med_text))
  } else {
    for (i in 1:nrow(km_table)) {
      grp_name  <- sub(".*=", "", rownames(km_table)[i])
      n_total   <- KM_curve$n[i]
      n_ev      <- km_table[i, grep("event", colnames(km_table))] 
      pct_ev    <- (n_ev / n_total) * 100
      med       <- km_table[i, "median"]
      m_low     <- km_table[i, "0.95LCL"]
      m_up      <- km_table[i, "0.95UCL"]
      
      med_text <- if(is.na(med)) "NA" else sprintf("%.2f (%.2f - %.2f)", med, m_low, m_up)
      
      cat(sprintf("%s: Events = %.1f%% (%d/%d); Median = %s\n", 
                  grp_name, pct_ev, as.integer(n_ev), as.integer(n_total), med_text))
    }
  }
  cat("\n")
  
  # --- 5. Landmark Survival ---
  s_all <- summary(KM_curve, times = times, extend = FALSE)
  actual_strata <- if(is.null(s_all$strata)) rep("All", length(s_all$time)) else sub(".*=", "", as.character(s_all$strata))
  expected_strata <- if(is.null(KM_curve$strata)) "All" else sub(".*=", "", names(KM_curve$strata))
  
  actual_df <- data.frame(strata = actual_strata, time = s_all$time, surv = s_all$surv,
                          lower = s_all$lower, upper = s_all$upper, stringsAsFactors = FALSE)
  expected_df <- expand.grid(strata = expected_strata, time = times, stringsAsFactors = FALSE)
  res_df <- merge(expected_df, actual_df, by = c("strata", "time"), all.x = TRUE)
  res_df <- res_df[order(res_df$strata, res_df$time), ]
  
  lines_out <- apply(res_df, 1, function(row) {
    if (is.na(row["surv"])) {
      sprintf("%s @%s mo: NA (NA-NA)", trimws(row["strata"]), trimws(row["time"]))
    } else {
      sprintf("%s @%s mo: %.3f (%.3f-%.3f)", trimws(row["strata"]), trimws(row["time"]), 
              as.numeric(row["surv"]), as.numeric(row["lower"]), as.numeric(row["upper"]))
    }
  })
  cat("Landmark Survival:\n")
  cat(paste(lines_out, collapse = "\n"), "\n")
  
  invisible(label_text)
}
#########################
## load and clean data ##
#########################

all_data <- flexread("Galaxy Liver Mets Data_20260131 ESMO GI ver2.3.xlsx - FANTASTIC plus.tsv")
dup <- duplicated(names(all_data)) |
  duplicated(names(all_data), fromLast = TRUE)
unique(names(all_data)[dup])
all_data <- all_data[, !dup, with = FALSE]

all_data=all_data[study_name=='Galaxy']

if ("RFS.ESMOGI" %in% names(all_data)) {
  setnames(all_data, "RFS.ESMOGI", "RFS.ESMO.GI")
}
all_data[NumLiverMetsGroup=='竕･2',NumLiverMetsGroup:='≥2']
all_data[SizeLiverMetsmmGroup=='竕･50',SizeLiverMetsmmGroup:='≥50']
all_data[Hepatectomy=='major',Hepatectomy:='Major']
all_data[Hepatectomy=='minor',Hepatectomy:='Minor']
all_data[,Postop.Complication:=as.character(Postop.Complication)]
all_data[Postop.Complication=='0',Postop.Complication:='No']
all_data[Postop.Complication=='1',Postop.Complication:='Yes']
all_data[Postop.Complication=='no',Postop.Complication:='No']
all_data[Postop.Complication=='yes',Postop.Complication:='Yes']
all_data[Oxaliplatin.History=='no',Oxaliplatin.History:='No']
all_data[Oxaliplatin.History=='yes',Oxaliplatin.History:='Yes']

# table(all_data$ctDNA.MRD)
# table(all_data$LiverMets)
cols_date <- c("Enrollment.Date","Surgery.Date", "Date.of.mets.surgery", 
               "Date.last.FU","Date.of.death","Date.of.surgery.for.recurrence","Date.of.recurrence",
               "DFS2.event.date","PFS2.event.date","start.of.ACT","Date.of.end.of.ACT")
all_data[, (cols_date) := lapply(.SD, as.Date, format = "%m/%d/%Y"),
   .SDcols = cols_date]

all_data[DFS.ESMO.GI=='29,07',.(as.numeric(Date.last.FU-Surgery.Date)/30,Date.last.FU-Surgery.Date)]
all_data[DFS.ESMO.GI=='29,07',DFS.ESMO.GI:=29.07]
all_data$DFS.ESMO.GI=as.numeric(all_data$DFS.ESMO.GI)

all_data[PFS2=='(2.80)',PFS2:='2.80']
all_data[,PFS2:=as.numeric(PFS2)]



########################
## check invalid data ##
########################

# if(F){
table(all_data[,.(MDT,is.na(Date.of.surgery.for.recurrence))])
all_data[MDT==FALSE & !is.na(Date.of.surgery.for.recurrence),Date.of.surgery.for.recurrence:=NA] # fix based on Kataoka e-mail title "Update analysis of CLM" at 8:40 pm Jan 31
all_data[MDT==FALSE & !is.na(Date.of.surgery.for.recurrence)]
all_data[(MDT==FALSE | is.na(MDT) )& !is.na(DFS2)]

all_data[DFS.Event==FALSE & (!is.na(DFS2) | !is.na(PFS2))]

all_data[ (DFS2.event==FALSE | is.na(DFS2.event) ) & !is.na(DFS2.event.date)]
all_data[ (PFS2.event==FALSE | is.na(PFS2.event) ) & !is.na(PFS2.event.date)]


## ACT 
all_data[start.of.ACT>Date.of.end.of.ACT,]
all_data[start.of.ACT<Date.of.mets.surgery,.(PtsID,start.of.ACT,Date.of.mets.surgery)]
table(all_data[,.(is.na(start.of.ACT),is.na(Date.of.end.of.ACT))])

## different Surgery.Date and Date.of.mets.surgery
all_data[Surgery.Date!=Date.of.mets.surgery,.(PtsID,Surgery.Date,Date.of.mets.surgery)]

## shorter DFS
all_data[DFS.months>DFS.ESMO.GI & DFS.Event==TRUE,
         .(PtsID,DFS.Event,old_DFS=round(DFS.months*30),new_DFS=round(DFS.ESMO.GI*30),DFS_2_surgery=Date.of.recurrence-Date.of.mets.surgery,followup_2_surgery=Date.last.FU-Date.of.mets.surgery)]
all_data[DFS.months>DFS.ESMO.GI & DFS.Event==FALSE,
         .(PtsID,DFS.Event,old_DFS=round(DFS.months*30),new_DFS=round(DFS.ESMO.GI*30),DFS_2_surgery=Date.of.recurrence-Date.of.mets.surgery,followup_2_surgery=Date.last.FU-Date.of.mets.surgery)]

## accuracy of current DFS 
all_data[(!is.na(DFS.ESMO.GI) & DFS.Event==TRUE & (abs(round(as.numeric(Date.of.recurrence-Date.of.mets.surgery))-round(DFS.ESMO.GI*30))>1)) |
           (is.na(Date.of.recurrence) & (abs(round(as.numeric(Date.of.death-Date.of.mets.surgery))-round(DFS.ESMO.GI*30))>1)),
         .(PtsID,DFS.Event,DFS=round(DFS.ESMO.GI*30),OS.Event,DFS_2_surgery=Date.of.recurrence-Date.of.mets.surgery,followup_2_surgery=Date.last.FU-Date.of.mets.surgery,OS_2_surgery=Date.of.death-Date.of.mets.surgery)]
all_data[!is.na(DFS.ESMO.GI) & DFS.Event==FALSE & (abs(round(as.numeric(Date.last.FU-Date.of.mets.surgery))-round(DFS.ESMO.GI*30))>1),
         .(PtsID,DFS.Event,DFS=round(DFS.ESMO.GI*30),DFS_2_surgery=Date.of.recurrence-Date.of.mets.surgery,followup_2_surgery=Date.last.FU-Date.of.mets.surgery)]


## accuracy of current DFS2 ** need to add NA case for recurrence
all_data[(!is.na(DFS2) & DFS2.event==TRUE & (abs(round(as.numeric(DFS2.event.date-Date.of.surgery.for.recurrence))-round(DFS2*30))>1)) |
           (is.na(DFS2.event.date) & (abs(round(as.numeric(Date.of.death-Date.of.surgery.for.recurrence))-round(DFS2*30))>1)),
         .(PtsID,DFS2.event,DFS2=round(DFS2*30),OS.Event,DFS2_2_MDT=DFS2.event.date-Date.of.surgery.for.recurrence,followup_2_MDT=Date.last.FU-Date.of.surgery.for.recurrence,OS_2_surgery=Date.of.death-Date.of.surgery.for.recurrence)]
all_data[!is.na(DFS2) & DFS2.event==FALSE & (abs(round(as.numeric(Date.last.FU-Date.of.surgery.for.recurrence))-round(DFS2*30))>1),
         .(PtsID,DFS2.event,DFS2=DFS2*30,DFS2_2_MDT=DFS2.event.date-Date.of.surgery.for.recurrence,followup_2_MDT=Date.last.FU-Date.of.surgery.for.recurrence)]


## accuracy of current PFS2
all_data[(!is.na(PFS2) & PFS2.event==TRUE & (abs(round(as.numeric(PFS2.event.date-Date.of.recurrence))-round(PFS2*30))>1)) |
           (is.na(PFS2.event.date) & (abs(round(as.numeric(Date.of.death-Date.of.recurrence))-round(PFS2*30))>1)),
         .(PtsID,PFS2.event,PFS2.event.date,Date.of.recurrence,PFS2=PFS2*30,PFS2_2_recur=PFS2.event.date-Date.of.recurrence,followup_2_recur=Date.last.FU-Date.of.recurrence)]
all_data[!is.na(PFS2) & PFS2.event==FALSE & (abs(round(as.numeric(Date.last.FU-Date.of.recurrence))-round(PFS2*30))>1),
         .(PtsID,PFS2.event,PFS2=PFS2*30,PFS2_2_recur=PFS2.event.date-Date.of.recurrence,followup_2_recur=Date.last.FU-Date.of.recurrence)]

## accuracy of OS
all_data[!is.na(OS.ESMO.GI) & OS.Event==TRUE & (abs(round(as.numeric(Date.of.death-Date.of.mets.surgery))-round(OS.ESMO.GI*30))>1),
         .(PtsID,OS.Event,OS=round(OS.ESMO.GI*30),OS_2_surgery=Date.of.death-Date.of.mets.surgery,followup_2_surgery=Date.last.FU-Date.of.mets.surgery)]

all_data[!is.na(OS.ESMO.GI) & OS.Event==FALSE & (abs(round(as.numeric(Date.last.FU-Date.of.mets.surgery))-round(OS.ESMO.GI*30))>1),
         .(PtsID,OS.Event,OS=round(OS.ESMO.GI*30),OS_2_surgery=Date.of.death-Date.of.mets.surgery,followup_2_surgery=Date.last.FU-Date.of.mets.surgery)]

## OS before DFS
all_data[DFS.ESMO.GI > OS.ESMO.GI,.(PtsID,DFS.Event,OS.Event,DFS.ESMO.GI,OS.ESMO.GI,
                                    Date.of.mets.surgery,Date.of.recurrence,Date.of.death,
                                    Date.last.FU,death_calculate=(Date.of.death-Date.of.mets.surgery),
                                    relapse_calculate=(Date.of.recurrence-Date.of.mets.surgery)
                                    )]

all_data[Date.of.recurrence > Date.of.death,.(PtsID,DFS.Event,OS.Event,DFS.ESMO.GI,OS.ESMO.GI,
                                    Date.of.mets.surgery,Date.of.recurrence,Date.of.death,
                                    Date.last.FU,death_calculate=(Date.of.death-Date.of.mets.surgery),
                                    relapse_calculate=(Date.of.recurrence-Date.of.mets.surgery)
)]
## old vs new OS
all_data[OS.ESMO.GI < OS.MRD.months & OS.Event==TRUE & (abs(round(as.numeric(Date.of.death-Date.of.mets.surgery))-round(OS.ESMO.GI*30))>1),
         .(PtsID,OS.Event,OS=round(OS.ESMO.GI*30),OS_2_surgery=Date.of.death-Date.of.mets.surgery,followup_2_surgery=Date.last.FU-Date.of.mets.surgery)]

## check number of death in old data
table(all_data[study_name=='Galaxy',OS.Event])
old_data=fread('~/Downloads/Galaxy Data_20240603 Complete Dataset.csv')
table(old_data[PtsID %in% all_data[study_name=='Galaxy',PtsID],OS.Event])

# }
#############################
## add additional variable ##
#############################
summary(all_data$ctDNA.MRD.Time) #69/7 =9.8

all_data[as.numeric(DFS.ESMO.GI)*30==1847.1,.(Surgery.Date,Date.of.mets.surgery,Date.last.FU,Date.last.FU-Date.of.mets.surgery)]

all_data[,DFS_0Dlandmark:=((as.numeric(DFS.ESMO.GI)*30)-0)/360*12]
all_data[,RFS_0Dlandmark:=((as.numeric(RFS.ESMO.GI)*30)-0)/360*12]
all_data[,OS_0Dlandmark:=((as.numeric(OS.ESMO.GI)*30)-0)/360*12]                                                
all_data[,DFS_70Dlandmark:=((as.numeric(DFS.ESMO.GI)*30)-70)/360*12]
all_data[,RFS_70Dlandmark:=((as.numeric(RFS.ESMO.GI)*30)-70)/360*12]
all_data[,OS_70Dlandmark:=((as.numeric(OS.ESMO.GI)*30)-70)/360*12]
all_data[,DFS_48Dlandmark:=((as.numeric(DFS.ESMO.GI)*30)-48)/360*12]
all_data[,RFS_48Dlandmark:=((as.numeric(RFS.ESMO.GI)*30)-48)/360*12]
all_data[,OS_48Dlandmark:=((as.numeric(OS.ESMO.GI)*30)-48)/360*12]
all_data[,DFS_98Dlandmark:=((as.numeric(DFS.ESMO.GI)*30)-98)/360*12]
all_data[,RFS_98Dlandmark:=((as.numeric(RFS.ESMO.GI)*30)-98)/360*12]
all_data[,OS_98Dlandmark:=((as.numeric(OS.ESMO.GI)*30)-98)/360*12]
all_data[,DFS_9Mlandmark:=(as.numeric(DFS.ESMO.GI)*30)/360*12-9]
all_data[,RFS_9Mlandmark:=(as.numeric(RFS.ESMO.GI)*30)/360*12-9]
all_data[,OS_9Mlandmark:=(as.numeric(OS.ESMO.GI)*30)/360*12-9]
all_data[,DFS_12.5Mlandmark:=(as.numeric(DFS.ESMO.GI)*30)/360*12-12.5]
all_data[,RFS_12.5Mlandmark:=(as.numeric(RFS.ESMO.GI)*30)/360*12-12.5]
all_data[,OS_12.5Mlandmark:=(as.numeric(OS.ESMO.GI)*30)/360*12-12.5]
all_data[,DFS2:=(as.numeric(DFS2)*30)/360*12]
all_data[,PFS2:=(as.numeric(PFS2)*30)/360*12]

all_data[!is.na(Date.of.recurrence),Surgery2Relapse_mrel:=as.numeric(Date.of.recurrence-Date.of.mets.surgery)/360*12]
all_data[!is.na(Date.of.surgery.for.recurrence),Surgery2MDT_mrel:=as.numeric(Date.of.surgery.for.recurrence-Date.of.mets.surgery)/360*12]
all_data[!is.na(Date.of.surgery.for.recurrence),Relapse2MDT_mrel:=as.numeric(Date.of.surgery.for.recurrence-Date.of.recurrence)/360*12]
all_data[!is.na(Date.of.end.of.ACT),Surgery2endACT_mrel:=as.numeric(Date.of.end.of.ACT-Date.of.mets.surgery)/360*12]
all_data[!is.na(start.of.ACT),Surgery2startACT_mrel:=as.numeric(start.of.ACT-Date.of.mets.surgery)/360*12]

all_data[!is.na(Date.of.recurrence) & OS.ESMO.GI>Surgery2Relapse_mrel,OS2.ESMO.GI:=OS.ESMO.GI-Surgery2Relapse_mrel]
all_data[!is.na(OS2.ESMO.GI),OS2.ESMO.GI_4Mlandmark:=OS2.ESMO.GI-4]
all_data[!is.na(PFS2),PFS2_4Mlandmark:=PFS2-4]
all_data[!is.na(OS2.ESMO.GI),OS2.ESMO.GI_28.1Mlandmark:=OS2.ESMO.GI-28.1]
all_data[!is.na(PFS2),PFS2_28.8Mlandmark:=PFS2-28.1]



###############################
## molecular positivity rate ##
###############################
a <- names(all_data)
a[grepl("ctdna", a, ignore.case = TRUE)]

all_data[,.(ctDNA.Baseline,ctDNA.MRD)] %>% summary()

table(all_data[,.(ctDNA.Baseline)])
prop.table(table(all_data[,.(ctDNA.Baseline)]))

Kataoka_baseline_positive=all_data[ctDNA.Baseline=='POSITIVE',.(PtsID)]
table(all_data[,.(ctDNA.MRD)])

load('~/Google Drive/My Drive/Galaxy_data/Galaxy_clinical2025Jul15_prep2025Aug22.RData')
dt.pat[cirID %in% all_data$PtsID,.(p_Baseline)] %>% table()
GALAXY_baseline_positive=dt.pat[cirID %in% all_data$PtsID & p_Baseline==TRUE,.(cirID)] 
setdiff(Kataoka_baseline_positive$PtsID,GALAXY_baseline_positive$cirID)



all_sig_test=dt.events[cirID %in% all_data$PtsID & evType=='Test/ctDNA/Signatera' & evValue %in% c('POSITIVE','NEGATIVE'),.(cirID,evDate,evValue,MTM=evValue2,p_dabsRef)] 
all_sig_test=all_data[,.(cirID=PtsID,Date.of.mets.surgery)][all_sig_test,on='cirID']
unique(all_sig_test[,.(cirID,Date.of.mets.surgery,p_dabsRef)])[,.(Date.of.mets.surgery==p_dabsRef)] %>% table() # we wil use Date.of.mets.surgery as ref for entire study
all_sig_test[,ref2test:=as.numeric(evDate-Date.of.mets.surgery)]

all_sig_test[ref2test<=0 ] %>% dim()

pre_surgican_data=all_sig_test[ref2test<=0]
unique(pre_surgican_data[,.(cirID,evValue)])[,.(evValue)] %>% table()
# NEGATIVE POSITIVE 
# 3      185 

unique(pre_surgican_data[,.(cirID,evValue)])[,.(evValue)] %>% table() %>% prop.table()

pre_surgican_data=unique(pre_surgican_data[,.(PtsID=cirID,pre_surgical_ctDNA=evValue,MTM=as.numeric(MTM))])
all_data=pre_surgican_data[all_data,on='PtsID']
all_data[,cohort_def:='Upfront Surgery']

all_data[,followup_drel:=(Date.last.FU-Date.of.mets.surgery)/30]

write.table(all_data,'UpfrontSurgery.tsv',quote=F,row.names=F,sep='\t')



############
## Table1 ##
############

circ_data_subset <- all_data %>%
  select(
    Age,
    Gender,
    ECOG,
    PrimSite,
    Mets.Type,
    Hepatectomy,
    NumLiverMetsGroup,
    SizeLiverMetsmmGroup,
    pT,
    pN,
    ACT,
    Oxaliplatin.History,
    Postop.Complication,
    BRAF.V600E,
    RAS,
    MSI,
    # RFS.Event,
    # OS.months,
    ctDNA.MRD,
    ctDNA.presurgical=pre_surgical_ctDNA,
    OS.Event,
    DFS.Event) %>%
  mutate(
    Age = as.numeric(Age),
    Gender = factor(Gender, levels = c("Male", "Female")),
    ECOG = factor(ECOG, levels = c(0, 1)),
    PrimSite = factor(PrimSite, levels = c("Right-sided colon", "Left-sided colon")),
    Mets.Type = factor(Mets.Type, levels = c("Synchronous", "Metachronous")),
    Hepatectomy = factor(Hepatectomy, levels = c("Minor", "Major")),
    NumLiverMetsGroup = factor(NumLiverMetsGroup, levels = c("1", "≥2")),
    SizeLiverMetsmmGroup = factor(SizeLiverMetsmmGroup, levels = c("<50", "≥50")),
    pT = factor(pT, levels = c("T1-T2", "T3-T4")),
    pN = factor(pN, levels = c("N0", "N1-N2")),
    ACT = factor(ACT, levels = c("TRUE", "FALSE"), labels = c("Adjuvant Chemotherapy", "Observation")),
    Oxaliplatin.History = factor(Oxaliplatin.History, levels = c("No", "Yes")),
    Postop.Complication = factor(Postop.Complication, levels = c("No", "Yes")),
    BRAF.V600E = factor(BRAF.V600E, levels = c("WT", "MUT"), labels = c("BRAF wt", "BRAF V600E")),
    RAS = factor(RAS, levels = c("WT", "MUT"), labels = c("RAS wt", "RAS mut")),
    MSI = factor(MSI, levels = c("MSS", "MSI-High"))
    # RFS.Event = factor(RFS.Event, levels = c("TRUE", "FALSE"), labels = c("Recurrence", "No Recurrence")),
    # OS.months = as.numeric(OS.months)
    )
table1 <- circ_data_subset %>%
  tbl_summary(
    by=ctDNA.MRD,
    type = list(
      Oxaliplatin.History ~ "categorical",
      Postop.Complication ~ "categorical"
    ),    statistic = list(
      all_continuous() ~ "{median} ({min} - {max})",
      all_categorical() ~ "{n} ({p}%)")) |>
  add_p(
    test = list(
      all_continuous()  ~ "wilcox.test",
      all_categorical() ~ "fisher.test"
    )) |>
  # add_q() |>
  add_overall()

table1

table_flextable=as_flex_table(table1)
doc <- read_docx() %>%
  body_add_flextable(table_flextable)

print(doc,target='TableS3A_UpfrontSurgery_cohort_description.docx')



###################
## MRD OS/DFS KM ##
###################

## DFS
circ_data <- all_data[all_data$ctDNA.MRD!="",]
circ_data <- circ_data[circ_data$DFS_70Dlandmark>=0,]
circ_datadf <- as.data.frame(circ_data)

survfit(Surv(time = circ_data$DFS_70Dlandmark, event = circ_data$DFS.Event)~ctDNA.MRD, data = circ_data)
surv_object <-Surv(time = circ_data$DFS_70Dlandmark, event = circ_data$DFS.Event)
KM_curve <- survfit(surv_object ~ ctDNA.MRD, data = circ_data,conf.int=0.95,conf.type="log-log") 
g=ggsurvplot(KM_curve, data = circ_data, pval = FALSE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("#46AC46","#004677"), title="DFS - ctDNA MRD window", ylab= "Disease-Free Survival", xlab="Time from Landmark Time point (Months)", legend.labs=c("ctDNA Negative", "ctDNA Positive"), legend.title="")
g
pdf('Figure1A_UpfrontSurgery_MRD_DFS.pdf')
print(g, newpage = FALSE)
dev.off()
summary(KM_curve, times= c(24, 30, 36))
circ_data$ctDNA.MRD <- factor(circ_data$ctDNA.MRD, levels=c("NEGATIVE","POSITIVE"))
cox_fit <- coxph(surv_object ~ ctDNA.MRD, data=circ_data) 
report_surv_extra(KM_curve, cox_fit)


## OS
circ_data <- all_data[all_data$ctDNA.MRD!="",]
circ_data <- circ_data[circ_data$OS_70Dlandmark>=0,]

survfit(Surv(time = circ_data$OS_70Dlandmark, event = circ_data$OS.Event)~ctDNA.MRD, data = circ_data)
surv_object <-Surv(time = circ_data$OS_70Dlandmark, event = circ_data$OS.Event)
KM_curve <- survfit(surv_object ~ ctDNA.MRD, data = circ_data,conf.int=0.95,conf.type="log-log") 
g=ggsurvplot(KM_curve, data = circ_data, pval = FALSE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("#46AC46","#004677"), title="OS - ctDNA MRD window", ylab= "Overall Survival", xlab="Time from Landmark Time point (Months)", legend.labs=c("ctDNA Negative", "ctDNA Positive"), legend.title="")
g
pdf('Figure1B_UpfrontSurgery_MRD_OS.pdf')
print(g, newpage = FALSE)
dev.off()
summary(KM_curve, times= c(24, 30, 36))
circ_data$ctDNA.MRD <- factor(circ_data$ctDNA.MRD, levels=c("NEGATIVE","POSITIVE"))
cox_fit <- coxph(surv_object ~ ctDNA.MRD, data=circ_data) 
report_surv_extra(KM_curve, cox_fit)


###############################################
## Surveillance non-time-dependent OS/DFS KM ##
###############################################

## DFS
circ_data <- all_data[all_data$ctDNA.Surveillance!="",]
circ_data <- circ_data[circ_data$DFS_70Dlandmark>=0,]

survfit(Surv(time = circ_data$DFS_70Dlandmark, event = circ_data$DFS.Event)~ctDNA.Surveillance, data = circ_data)
surv_object <-Surv(time = circ_data$DFS_70Dlandmark, event = circ_data$DFS.Event)
KM_curve <- survfit(surv_object ~ ctDNA.Surveillance, data = circ_data,conf.int=0.95,conf.type="log-log") 
g=ggsurvplot(KM_curve, data = circ_data, pval = FALSE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("#46AC46","#004677"), title="DFS - ctDNA surveillance window", ylab= "Disease-Free Survival", xlab="Time from Landmark Time point (Months)", legend.labs=c("ctDNA Negative", "ctDNA Positive"), legend.title="")
g
pdf('FigureS2A_UpfrontSurgery_Surv_DFS.pdf')
print(g, newpage = FALSE)
dev.off()
summary(KM_curve, times= c(24, 30, 36))
circ_data$ctDNA.Surveillance <- factor(circ_data$ctDNA.Surveillance, levels=c("NEGATIVE","POSITIVE"))
cox_fit <- coxph(surv_object ~ ctDNA.Surveillance, data=circ_data) 
report_surv_extra(KM_curve, cox_fit)



## OS
circ_data <- all_data[all_data$ctDNA.Surveillance!="",]
circ_data <- circ_data[circ_data$OS_70Dlandmark>=0,]
circ_datadf <- as.data.frame(circ_data)

survfit(Surv(time = circ_data$OS_70Dlandmark, event = circ_data$OS.Event)~ctDNA.Surveillance, data = circ_data)
surv_object <-Surv(time = circ_data$OS_70Dlandmark, event = circ_data$OS.Event)
KM_curve <- survfit(surv_object ~ ctDNA.Surveillance, data = circ_data,conf.int=0.95,conf.type="log-log") 
g=ggsurvplot(KM_curve, data = circ_data, pval = FALSE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("#46AC46","#004677"), title="OS - ctDNA surveillance window", ylab= "Overall Survival", xlab="Time from Landmark Time point (Months)", legend.labs=c("ctDNA Negative", "ctDNA Positive"), legend.title="")
g
pdf('FigureS2B_UpfrontSurgery_Surv_OS.pdf')
print(g, newpage = FALSE)
dev.off()
summary(KM_curve, times= c(24, 30, 36))
circ_data$ctDNA.Surveillance <- factor(circ_data$ctDNA.Surveillance, levels=c("NEGATIVE","POSITIVE"))
cox_fit <- coxph(surv_object ~ ctDNA.Surveillance, data=circ_data) 
report_surv_extra(KM_curve, cox_fit)


############################################
## Surveillance time-dependent OS/DFS cox ##
############################################

circ_data <- all_data
dim(circ_data)
dim(all_sig_test)
length(unique(all_sig_test$cirID))
all_sig_test=all_sig_test[order(cirID,evDate)]
all_sig_test[,PtsID:=cirID]
merge_all_test=all_sig_test[circ_data,on='PtsID']
dim(merge_all_test)
surveillance_test=merge_all_test
dim(surveillance_test)
# surveillance_test[,.(Date.of.recurrence,Date.of.mets.surgery,ref2test,Date.of.end.of.ACT)]
surveillance_test=surveillance_test[is.na(Date.of.recurrence) | (Date.of.recurrence-Date.of.mets.surgery) >=ref2test]

surveillance_test=surveillance_test[(ACT==FALSE & ref2test>70) | (ACT==TRUE & (ref2test-(Date.of.end.of.ACT-Date.of.mets.surgery)>=28))]
dim(surveillance_test)
surveillance_test


# 1. Sort and Remove Duplicates
setkey(surveillance_test, cirID, evDate)

# Identify duplicates
dup_indices <- duplicated(surveillance_test, by = c("cirID", "evDate"))
if (any(dup_indices)) {
  message("Duplicate records found and removed for the following IDs:")
  print(surveillance_test[dup_indices, .(cirID, evDate)])
  surveillance_test <- unique(surveillance_test, by = c("cirID", "evDate"))
}

# 2. Define Outcome Dates
# We need to establish the terminal date for DFS and OS for each patient
surveillance_test[, `:=`(
  # DFS Outcome: Recurrence date if exists, otherwise last FU
  dfs_end_date = as.Date(ifelse(!is.na(Date.of.recurrence), 
                                Date.of.recurrence, 
                                Date.last.FU), origin = "1970-01-01"),
  # OS Outcome: Death date if exists, otherwise last FU
  os_end_date = as.Date(ifelse(!is.na(Date.of.death), 
                               Date.of.death, 
                               Date.last.FU), origin = "1970-01-01")
)]

# 3. Create Time-Dependent Intervals
# We normalize time to "Days since Surgery"
surveillance_test[, t_start := as.numeric(evDate - Date.of.mets.surgery)]

# Calculate stop times (next evDate or the terminal outcome date)
surveillance_test[, next_evDate := shift(evDate, type = "lead"), by = cirID]

# 4. Generate tdDFS and tdOS segments
# For the last test of each patient, the 'stop' time is the outcome date
surveillance_test[, `:=`(
  tdDFS_stop = as.numeric(fifelse(!is.na(next_evDate), next_evDate, dfs_end_date) - Date.of.mets.surgery),
  tdOS_stop  = as.numeric(fifelse(!is.na(next_evDate), next_evDate, os_end_date) - Date.of.mets.surgery)
)]

# 5. Filter "Collisions" (Zero-length intervals)
# Remove intervals where the test date is the same as the outcome date (duration 0)
surveillance_test <- surveillance_test[tdDFS_stop > t_start | tdOS_stop > t_start]

# 6. Set Surveillance Labels and Events
surveillance_test[, `:=`(
  # 3. Label for the interval is the ctDNA result at the start of that window
  surveillance.interval.label = evValue,
  
  # 4. tdDFS.event logic: 
  # Set to 1 only if it's the last interval for a patient who had an event
  tdDFS.event = as.integer(DFS.Event == TRUE & is.na(next_evDate)),
  
  # 4. tdOS.event logic:
  # Set to 1 only if it's the last interval for a patient who had an event
  tdOS.event = as.integer(OS.Event == TRUE & is.na(next_evDate))
)]

# Cleanup temporary columns
surveillance_test[, c("next_evDate", "dfs_end_date", "os_end_date") := NULL]

# View result
head(surveillance_test[, .(cirID, t_start, tdDFS_stop, surveillance.interval.label, tdDFS.event)])

# Set levels so NEGATIVE is the reference group
surveillance_test[, surveillance.interval.label := factor(
  surveillance.interval.label, 
  levels = c("NEGATIVE", "POSITIVE")
)]

## DFS

# 1. Create the Time-Dependent Surv object
# Note: Surv(time, time2, event)
surv_obj_dfs <- Surv(
  time = surveillance_test$t_start, 
  time2 = surveillance_test$tdDFS_stop, 
  event = surveillance_test$tdDFS.event
)

# 2. Fit the Cox Model
cox_fit_dfs <- coxph(surv_obj_dfs ~ surveillance.interval.label, data = surveillance_test)

# 3. Create the KM Curve
# survfit handles the start-stop format automatically
KM_curve_dfs <- survfit(surv_obj_dfs ~ surveillance.interval.label, data = surveillance_test)

# 4. Report results using your function
report_surv_extra(KM_curve_dfs, cox_fit_dfs)


## OS

# 1. Create the Time-Dependent Surv object for OS
surv_obj_os <- Surv(
  time = surveillance_test$t_start, 
  time2 = surveillance_test$tdOS_stop, 
  event = surveillance_test$tdOS.event
)

# 2. Fit the Cox Model
cox_fit_os <- coxph(surv_obj_os ~ surveillance.interval.label, data = surveillance_test)

# 3. Create the KM Curve
KM_curve_os <- survfit(surv_obj_os ~ surveillance.interval.label, data = surveillance_test)

# 4. Report results
report_surv_extra(KM_curve_os, cox_fit_os)


#############
## MRD MVA ##
#############
## option
circ_data <- all_data[all_data$ctDNA.MRD!="",]


circ_data <- circ_data[circ_data$DFS_70Dlandmark>=0,]
# circ_data <- circ_data[circ_data$OS_70Dlandmark>=0,]

circ_data$ctDNA.MRD <- factor(circ_data$ctDNA.MRD, levels=c("NEGATIVE","POSITIVE"), labels = c("Negative", "Positive"))
circ_data$Gender <- factor(circ_data$Gender, levels = c("Female", "Male"))
circ_data$Age.Group <- factor(circ_data$Age.Group, levels = c("1", "2"), labels = c("<70", ">70"))
circ_data$PrimSite <- factor(circ_data$PrimSite, levels = c("Right-sided colon", "Left-sided colon"), labels = c("Right-sided", "Left-sided"))
circ_data$NumLiverMetsGroup <- factor(circ_data$NumLiverMetsGroup, levels = c("1", "≥2"))
circ_data$SizeLiverMetsmmGroup <- factor(circ_data$SizeLiverMetsmmGroup, levels = c("<50", "≥50"))
circ_data$ECOG <- factor(circ_data$ECOG, levels = c("0", "1"))
circ_data$pT <- factor(circ_data$pT, levels = c("T1-T2", "T3-T4"))
circ_data$pN <- factor(circ_data$pN, levels = c("N0", "N1-N2"))
circ_data$ACT <- factor(circ_data$ACT, levels = c("FALSE", "TRUE"), labels = c("Observation", "Chemotherapy"))
circ_data$Postop.Complication <- factor(circ_data$Postop.Complication, levels = c("No", "Yes"))
circ_data$Oxaliplatin.History <- factor(circ_data$Oxaliplatin.History, levels = c("No", "Yes"))
circ_data$Mets.Type <- factor(circ_data$Mets.Type, levels = c("Metachronous", "Synchronous"))
circ_data$MSI <- factor(circ_data$MSI, levels = c("MSS", "MSI-High"), labels = c("MSS", "MSI-High"))
circ_data$BRAF.V600E <- factor(circ_data$BRAF.V600E, levels = c("WT", "MUT"), labels = c("Wild-Type", "V600E"))
circ_data$RAS <- factor(circ_data$RAS, levels = c("WT", "MUT"), labels = c("Wild-Type", "Mutant"))


surv_object <- Surv(time = circ_data$DFS_70Dlandmark, event = circ_data$DFS.Event)
# surv_object <- Surv(time = circ_data$OS_70Dlandmark, event = circ_data$OS.Event)
# 
## option
cox_fit <- coxph(surv_object ~ ctDNA.MRD + Age.Group + PrimSite + pT + pN + RAS + Mets.Type + Postop.Complication + SizeLiverMetsmmGroup + NumLiverMetsGroup  + Oxaliplatin.History + ACT, data=circ_data)
# summary(cox_fit)

## option
g=ggforest(cox_fit, data = circ_data, main = "Multivariate Regression Model for DFS", refLabel = "Reference Group")
# g=ggforest(cox_fit, data = circ_data, main = "Multivariate Regression Model for OS", refLabel = "Reference Group")
g
pdf('Figure1C_UpfrontSurgery_MVA_MRD_DFS.pdf')
# pdf('Figure1D_UpfrontSurgery_MVA_MRD_OS.pdf')
print(g, newpage = FALSE)
dev.off()





#########################
## ACT KM by MRD group ##
#########################

## DFS

## option
circ_data <- all_data[ctDNA.MRD=="POSITIVE",]
# circ_data <- all_data[ctDNA.MRD=="NEGATIVE",]


## option
circ_data <- circ_data[circ_data$DFS_70Dlandmark>=0,]
survfit(Surv(time = circ_data$DFS_70Dlandmark, event = circ_data$DFS.Event)~ACT, data = circ_data)
surv_object <-Surv(time = circ_data$DFS_70Dlandmark, event = circ_data$DFS.Event)



KM_curve <- survfit(surv_object ~ ACT, data = circ_data,conf.int=0.95,conf.type="log-log") 

## option
g=ggsurvplot(KM_curve, data = circ_data, pval = FALSE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("#004677","#46AC46"), title="DFS - ctDNA MRD Positive ACT vs Observation", ylab= "Disease-Free Survival", xlab="Time from Landmark Time point (Months)", legend.labs=c("Observation", "ACT"), legend.title="")
# g=ggsurvplot(KM_curve, data = circ_data, pval = FALSE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("#004677","#46AC46"), title="DFS - ctDNA MRD Negative ACT vs Observation", ylab= "Disease-Free Survival", xlab="Time from Landmark Time point (Months)", legend.labs=c("Observation", "ACT"), legend.title="")

g
pdf('Figure3A_UpfrontSurgery_ACT_ctDNAPos_DFS.pdf')
# pdf('Figure3B_UpfrontSurgery_ACT_ctDNANeg_DFS.pdf')
print(g, newpage = FALSE)
dev.off()

# summary(KM_curve)
circ_data$ACT <- factor(circ_data$ACT, levels=c("TRUE","FALSE"))

## option
# cox_fit <- coxph(surv_object ~ ACT, data=circ_data)
#1# full model
cox_fit <- coxph(surv_object ~ ACT + Age.Group + Gender + RAS +
                   SizeLiverMetsmmGroup + NumLiverMetsGroup + Oxaliplatin.History +
                   Mets.Type, data=circ_data)


# ggforest(cox_fit)
report_surv_extra(KM_curve, cox_fit)




## OS

## option
# circ_data <- all_data[all_data$ctDNA.MRD=="POSITIVE",]
circ_data <- all_data[all_data$ctDNA.MRD=="NEGATIVE",]


## option
circ_data <- circ_data[circ_data$OS_70Dlandmark>=0,]
survfit(Surv(time = circ_data$OS_70Dlandmark, event = circ_data$OS.Event)~ACT, data = circ_data)
surv_object <-Surv(time = circ_data$OS_70Dlandmark, event = circ_data$OS.Event)


KM_curve <- survfit(surv_object ~ ACT, data = circ_data,conf.int=0.95,conf.type="log-log") 

## option
# g=ggsurvplot(KM_curve, data = circ_data, pval = FALSE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("#004677","#46AC46"), title="OS - ctDNA MRD Positive ACT vs Observation", ylab= "Overall Survival", xlab="Time from Landmark Time point (Months)", legend.labs=c("Observation", "ACT"), legend.title="")
g=ggsurvplot(KM_curve, data = circ_data, pval = FALSE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("#004677","#46AC46"), title="OS - ctDNA MRD Negative ACT vs Observation", ylab= "Overall Survival", xlab="Time from Landmark Time point (Months)", legend.labs=c("Observation", "ACT"), legend.title="")

g
# pdf('Figure3C_UpfrontSurgery_ACT_ctDNAPos_OS.pdf')
pdf('Figure3D_UpfrontSurgery_ACT_ctDNANeg_OS.pdf')
print(g, newpage = FALSE)
dev.off()

# summary(KM_curve)
circ_data$ACT <- factor(circ_data$ACT, levels=c("TRUE","FALSE"))


## option
# cox_fit <- coxph(surv_object ~ ACT, data=circ_data)
#1# full model
cox_fit <- coxph(surv_object ~ ACT + Age.Group + Gender + RAS +
                   SizeLiverMetsmmGroup + NumLiverMetsGroup + Oxaliplatin.History +
                   Mets.Type, data=circ_data)


report_surv_extra(KM_curve, cox_fit)



####################
## recur location ##
####################

circ_data <- all_data[RFS.Event == "TRUE"]
circ_data=circ_data[Rec.Site!='']
circ_data[Rec.Site == 'liver', Rec.Site := 'Liver']
circ_data[Rec.Site == 'lung', Rec.Site := 'Lung']
circ_data[Rec.Site == 'Peritoneum & Others', Rec.Site := 'Peritoneum']
circ_data$ctDNA.MRD=factor(circ_data$ctDNA.MRD,levels=c('POSITIVE','NEGATIVE'))

table_df <- circ_data[!is.na(ctDNA.MRD) & !is.na(Rec.Site), 
                      .(Freq = .N), 
                      by = .(ctDNA.MRD, Rec.Site)]

table_df[, Rec.Site := factor(Rec.Site, 
                              levels = c("Liver","Lung", "Peritoneum", "Local/LN","Brain"))]

table_df[, Total := sum(Freq), by = ctDNA.MRD]
table_df[, Percentage := Freq / Total]
table_df[, Label := ifelse(Freq > 0, paste0(Freq, " (", scales::percent(Percentage, accuracy = 1), ")"), NA)]

contingency_table <- xtabs(Freq ~ ctDNA.MRD + Rec.Site, data = table_df)
chi_square_test <- chisq.test(contingency_table)
p_value <- chi_square_test$p.value

g=ggplot(table_df, aes(x = ctDNA.MRD, y = Percentage, fill = Rec.Site)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Label), 
            position = position_stack(vjust = 0.5), 
            na.rm = TRUE, 
            color = "black", size = 7) +
  scale_fill_viridis_d(option = "plasma", direction = -1, name = "Recurrence Site") +
  theme_minimal() +
  labs(title = "Patients with Radiological Recurrence", 
       x = "ctDNA at the MRD Window", 
       y = "Patients (%)", 
       caption = paste("Chi-squared test p-value: ", format.pval(p_value, digits=3))) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

g
pdf('FigureS4A_UpfrontSurgery_location.pdf')
print(g, newpage = FALSE)
dev.off()


