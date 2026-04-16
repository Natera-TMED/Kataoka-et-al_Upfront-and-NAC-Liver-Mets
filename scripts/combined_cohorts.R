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
source("/Users/cfungtammasan/gitlab/translationaloncology/tmed_stat.R")



non_NAC=fread('UpfrontSurgery.tsv')
non_NAC=non_NAC[,.(PtsID,
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
           Margins,
           ctDNA.MRD,
           ctDNA.presurgical=pre_surgical_ctDNA,
           ctDNA.Surveillance,
           OS.Event,
           DFS.Event,
           RFS.Event,
           cohort_def,
           `Follow up month`=followup_drel
)]
dim(non_NAC)
NAC=fread('NAC.tsv')
NAC=NAC[,.(PtsID,
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
                   Margins,
                   ctDNA.MRD,
                   ctDNA.presurgical=pre_surgical_ctDNA,
                   ctDNA.Surveillance,
                   OS.Event,
                   DFS.Event,
                   RFS.Event,
                   cohort_def,
                   `Follow up month`=followup_drel
           
)]
dim(NAC)

combined_data <- rbind(non_NAC, NAC)
combined_data$PtsID=NULL
combined_data$cohort_def=factor(combined_data$cohort_def,levels=c('Upfront Surgery','NAC'))
combined_data[ctDNA.Surveillance=='',ctDNA.Surveillance:=NA]
median(combined_data$`Follow up month`)

table1 <- combined_data %>%
  tbl_summary(
    by=cohort_def,
    type = list(
      Oxaliplatin.History ~ "categorical",
      Postop.Complication ~ "categorical",
      ACT ~ "categorical"
    ),    statistic = list(
      all_continuous() ~ "{median} ({min} - {max})",
      all_categorical() ~ "{n} ({p}%)")) |>
  add_p(
    test = list(
      all_continuous()  ~ "wilcox.test",
      all_categorical() ~ "fisher.test"
    )) |>
  gtsummary::add_q(method = "bonferroni") |>   
  add_overall()

table1

table_flextable=as_flex_table(table1)
doc <- read_docx() %>%
  body_add_flextable(table_flextable)

# print(doc,target='TableS1_both_cohort_description.docx')

combined_data <- rbind(non_NAC, NAC)

combined_data[(!is.na(DFS2) & DFS2.event==TRUE & (abs(round(as.numeric(DFS2.event.date-Date.of.surgery.for.recurrence))-round(DFS2*30))>1)) |
           (is.na(DFS2.event.date) & (abs(round(as.numeric(Date.of.death-Date.of.surgery.for.recurrence))-round(DFS2*30))>1)),
         .(PtsID,DFS2.event,DFS2=round(DFS2*30),OS.Event,DFS2_2_MDT=DFS2.event.date-Date.of.surgery.for.recurrence,followup_2_MDT=Date.last.FU-Date.of.surgery.for.recurrence,OS_2_surgery=Date.of.death-Date.of.surgery.for.recurrence,cohort_def)]
combined_data[!is.na(DFS2) & DFS2.event==FALSE & (abs(round(as.numeric(Date.last.FU-Date.of.surgery.for.recurrence))-round(DFS2*30))>1),
         .(PtsID,DFS2.event,DFS2=DFS2*30,DFS2_2_MDT=DFS2.event.date-Date.of.surgery.for.recurrence,followup_2_MDT=Date.last.FU-Date.of.surgery.for.recurrence,cohort_def)]





load('~/Google Drive/My Drive/Galaxy_data/Galaxy_clinical2025Jul15_prep2025Aug22.RData')



all_sig_test=dt.events[cirID %in% combined_data$PtsID & evType=='Test/ctDNA/Signatera' & evValue %in% c('POSITIVE','NEGATIVE'),.(PtsID=cirID,evDate,evValue)] 
all_sig_test=combined_data[,.(PtsID,Date.of.surgery.for.recurrence=as.Date(Date.of.surgery.for.recurrence),cohort_def,DFS2.event,DFS2,Date.last.FU=as.Date(Date.last.FU))][all_sig_test,on='PtsID']
all_sig_test[,MDT2test:=as.numeric(evDate-Date.of.surgery.for.recurrence)]

all_sig_test[MDT2test>0 ,.(ctDNA=evValue,cohort_def)]  %>% table()
all_sig_test[MDT2test>14 , PtsID] %>% unique() %>% length()
all_sig_test[MDT2test>14 & MDT2test<(12*7) ,.(ctDNA=evValue,cohort_def)]  %>% table()
all_sig_test[MDT2test>14 & MDT2test<(12*7), PtsID] %>% unique() %>% length()
all_sig_test[MDT2test>0 , PtsID] %>% unique() %>% length()

all_sig_test=all_sig_test[MDT2test>14]
# all_sig_test=all_sig_test[MDT2test>14 & MDT2test<(12*7)]

all_sig_test[,postMDT_ctDNA_cat:=paste0(evValue,collapse=''),by=PtsID]

all_sig_test[,post_MDT_ctDNA:=fcase(postMDT_ctDNA_cat %like% 'POSITIVE','POSITIVE',
                                    postMDT_ctDNA_cat %like% 'NEGATIVE','NEGATIVE',
                                    default=NA
                                    )]
unique(all_sig_test[,.(PtsID,cohort_def,post_MDT_ctDNA)])
all_sig_test[order(PtsID,evDate)][evValue=='POSITIVE']
all_sig_test[order(PtsID,evDate)][evValue=='NEGATIVE',.(PtsID,evDate,(Date.last.FU-Date.of.surgery.for.recurrence)/30)]
all_sig_test[order(PtsID,evDate)][evValue=='POSITIVE']


postMDT_DFS2_data=unique(all_sig_test[,.(PtsID,DFS2.event,DFS2,post_MDT_ctDNA)])

             
surv_object <-Surv(time = postMDT_DFS2_data$DFS2, event = postMDT_DFS2_data$DFS2.event)

KM_curve <- survfit(surv_object ~ post_MDT_ctDNA, data = postMDT_DFS2_data,conf.int=0.95,conf.type="log-log")

g=ggsurvplot(KM_curve, data = postMDT_DFS2_data, pval = TRUE, conf.int = FALSE, risk.table = TRUE, break.time.by=6, palette=c("blue","red"), title="DFS2", ylab= "Disease-Free Survival after recur", xlab="Time from MDT (Months)", legend.title="")

cox_fit <- coxphf(surv_object ~ post_MDT_ctDNA, data=postMDT_DFS2_data)
g
pdf('Figure5_DFS2.pdf')
print(g, newpage = FALSE)
dev.off()

report_surv_extra(KM_curve, cox_fit)
