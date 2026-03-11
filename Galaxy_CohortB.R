library(survival)
library(survminer)
library(ggplot2)
library(ggfortify)
library(gtsummary)
library(tidycmprsk)
library(lubridate)
library(dplyr)
library(coxphf)
library(readxl)


rm(list=ls()) # clears all environments plots etc so starts the script clean
setwd("~/Documents/Clinical Data Mining/CRC/Rectal/CohortB/")
getwd()

#DFS & OS (switch between DFS.Month or OS.Month and updated DFS.Event and OS.Event)
{mrd_outcome <- read_excel("DFS_OS_KM.xlsx")
mrd_outcome <- mrd_outcome %>%
  filter(Eligible == TRUE)

mrd_outcome


km_mrd <- with(mrd_outcome[which(mrd_outcome$PostNAT!= "NA"),], Surv(mrd_outcome$OS.Month, mrd_outcome$OS.Event))
km_mrd_fit <- survfit(Surv(OS.Month, OS.Event) ~ PostNAT, data = mrd_outcome)
surv_pvalue(km_mrd_fit) #p=0.47
#Plot version 2

ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = TRUE,
           #fun = "cumhaz",
           risk.table = TRUE,
           size = 1,
           # linetype = "strata",
           #linetypes= c("solid", "solid"),
           censor = TRUE,
           palette = c("#253494",
                       "#238443"
           ),
           xlab = "Months from Landmark",
           break.time.by = 3,
           #xlim = c(0, 42),
           ylab = "OS Probability", 
           legend = "none",
           legend.title = "Post-NAT",
           legend.labs = c("Negative", "Positive"),
           ylim = c(0, 1))
#+ ggtitle("Rectal Baseline")
# surv.median.line = c("hv"))
summary(km_mrd_fit)$table





#color schemes
#253494 is blue
#238443 is green
#"#8BC400" is natera green
#"#008BCE" is natera blue

#Calculate %PFS/OS at 12 mo, 24 mo
summary(km_mrd_fit, times= c(24,36))
#Calculate Median DFS and median OS for both groups
survfit(Surv(time = mrd_outcome$DFS.Month, event = mrd_outcome$DFS.Event)~PostNAT, data = mrd_outcome)


#exp(coef) is the HR
#lower .95 is lower bound CI
#upper .95 is upper bound CI


# coxfit hazard ratio for just two factors if there are multiple lines in KM curve
# Filter the data to include only two groups
mrd_subset <- mrd_outcome %>% dplyr::filter(PostNAT %in% c("0", "1")) %>%
  droplevels()   # <-- NEW: drop levels 2–5 so they don't appear as "reference"

# Relevel the grouping variable to set the reference level
mrd_subset$PostNAT <- factor(mrd_subset$PostNAT)
mrd_subset$PostNAT <- relevel(factor(mrd_subset$PostNAT), ref = "0")

# Create the survival object
surv_object <- Surv(time = mrd_subset$OS.Month, event = mrd_subset$OS.Event)

# Fit the Cox model
cox_fit <- coxphf(surv_object ~ PostNAT, data = mrd_subset)

#view the HR
ggforest(cox_fit,data = mrd_outcome)

# View the summary for Hazard Ratio, CI, and p-value
summary(cox_fit)
}

#pCR (switch between DFS.Month or OS.Month and updated DFS.Event and OS.Event)
{
mrd_outcome <- readxl::read_excel("DFS_OS_KM.xlsx") %>%
  mutate(
    Eligible = case_when(
      as.character(Eligible) %in% c("TRUE","True","1","Yes","YES") ~ TRUE,
      as.character(Eligible) %in% c("FALSE","False","0","No","NO") ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  filter(Eligible == TRUE) %>%
  # ↓↓↓ MINIMAL ADD: drop rows with missing or non-positive DFS time. update to DFS and OS depending on endpoint
  filter(!is.na(OS.Month) & OS.Month > 0)

mrd_outcome

# fix NA test and avoid re-referencing mrd_outcome$ inside with()
km_mrd <- with(mrd_outcome[!is.na(mrd_outcome$pCR), ], Surv(OS.Month, OS.Event))
km_mrd_fit <- survfit(Surv(OS.Month, OS.Event) ~ pCR, data = mrd_outcome)
surv_pvalue(km_mrd_fit) # p=0.47

ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = TRUE,
           risk.table = TRUE,
           size = 1,
           censor = TRUE,
           palette = c("#253494", "#238443"),
           xlab = "Months from Surgery",
           break.time.by = 3,
          # xlim = c(0, 42),
           ylab = "OS Probability",
           legend = "none",
           legend.title = "Path Response",
           legend.labs = c("pCR", "non-pCR"),
           ylim = c(0, 1))
summary(km_mrd_fit)$table


#color schemes
#253494 is blue
#238443 is green
#"#8BC400" is natera green
#"#008BCE" is natera blue

#Calculate %PFS/OS at 12 mo, 24 mo
summary(km_mrd_fit, times= c(24,36))
#Calculate Median DFS and median OS for both groups
survfit(Surv(time = mrd_outcome$DFS.Month, event = mrd_outcome$DFS.Event)~pCR, data = mrd_outcome)


#exp(coef) is the HR
#lower .95 is lower bound CI
#upper .95 is upper bound CI


# coxfit hazard ratio for just two factors if there are multiple lines in KM curve
# Filter the data to include only two groups
mrd_subset <- mrd_outcome %>% dplyr::filter(pCR %in% c("0", "1")) %>%
  droplevels()   # <-- NEW: drop levels 2–5 so they don't appear as "reference"

# Relevel the grouping variable to set the reference level
mrd_subset$pCR <- factor(mrd_subset$pCR)
mrd_subset$pCR <- relevel(factor(mrd_subset$pCR), ref = "0")

# Create the survival object
surv_object <- Surv(time = mrd_subset$DFS.Month, event = mrd_subset$DFS.Event)

# Fit the Cox model
cox_fit <- coxphf(surv_object ~ pCR, data = mrd_subset)

#view the HR
ggforest(cox_fit,data = mrd_outcome)

# View the summary for Hazard Ratio, CI, and p-value
summary(cox_fit)
}

#pCR & NAC (switch between DFS.Month or OS.Month and updated DFS.Event and OS.Event)
{
  mrd_outcome <- read_excel("pCR_NAC_KM.xlsx")
  mrd_outcome <- mrd_outcome %>%
    filter(Eligible == TRUE)
  
  mrd_outcome
  
  
  km_mrd <- with(mrd_outcome[which(mrd_outcome$PostNAT.pCR!= "NA"),], Surv(mrd_outcome$OS.Month, mrd_outcome$OS.Event))
  km_mrd_fit <- survfit(Surv(OS.Month, OS.Event) ~ PostNAT.pCR, data = mrd_outcome)
  surv_pvalue(km_mrd_fit) #p=0.47
  #Plot version 2
  
  ggsurvplot(km_mrd_fit,
             conf.int = FALSE,
             pval = TRUE,
             #fun = "cumhaz",
             risk.table = TRUE,
             size = 1,
             # linetype = "strata",
             #linetypes= c("solid", "solid"),
             censor = TRUE,
             palette = c("#253494",
                         "#238443","orange"
             ),
             xlab = "Months from Surgery",
             break.time.by = 3,
             #xlim = c(0, 39),
             ylab = "DFS Probability", 
             legend = "none",
             legend.title = "PostNAT & Path Response",
             legend.labs = c("ctDNA- & pCR", "ctDNA- & nonpCR", "ctDNA+ & nonpCR"),
             ylim = c(0, 1))
  #+ ggtitle("Rectal Baseline")
  # surv.median.line = c("hv"))
  summary(km_mrd_fit)$table
  
  
  

  
  #color schemes
  #253494 is blue
  #238443 is green
  #"#8BC400" is natera green
  #"#008BCE" is natera blue
  
  #Calculate %PFS/OS at 12 mo, 24 mo
  summary(km_mrd_fit, times= c(24,36))
  #Calculate Median PFS and median OS for both groups
  survfit(Surv(time = mrd_outcome$OS.Month, event = mrd_outcome$OS.Event)~PostNAT.pCR, data = mrd_outcome)
  
  
  #exp(coef) is the HR
  #lower .95 is lower bound CI
  #upper .95 is upper bound CI
  
  
  # coxfit hazard ratio for just two factors if there are multiple lines in KM curve
  # Filter the data to include only two groups
  mrd_subset <- mrd_outcome %>% dplyr::filter(PostNAT.pCR %in% c("1", "2")) %>%
    droplevels()   # <-- NEW: drop levels 2–5 so they don't appear as "reference"
  
  # Relevel the grouping variable to set the reference level
  mrd_subset$PostNAT.pCR <- factor(mrd_subset$PostNAT.pCR)
  mrd_subset$PostNAT.pCR <- relevel(factor(mrd_subset$PostNAT.pCR), ref = "1")
  
  # Create the survival object
  surv_object <- Surv(time = mrd_subset$DFS.Month, event = mrd_subset$DFS.Event)
  
  # Fit the Cox model
  cox_fit <- coxphf(surv_object ~ PostNAT.pCR, data = mrd_subset)
  
  #view the HR
  ggforest(cox_fit,data = mrd_outcome)
  
  # View the summary for Hazard Ratio, CI, and p-value
  summary(cox_fit)
}

#DFS_4W (MRD)
{
mrd_outcome <- read_excel("DFS_4W_KM.xlsx")
mrd_outcome <- mrd_outcome %>%
  filter(Eligible == TRUE)

mrd_outcome


km_mrd <- with(mrd_outcome[which(mrd_outcome$w4.ctDNA!= "NA"),], Surv(mrd_outcome$DFS.Month, mrd_outcome$DFS.Event))
km_mrd_fit <- survfit(Surv(DFS.Month, DFS.Event) ~ w4.ctDNA, data = mrd_outcome)
surv_pvalue(km_mrd_fit) #p=0.47
#Plot version 2

ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = TRUE,
           #fun = "cumhaz",
           risk.table = TRUE,
           size = 1,
           # linetype = "strata",
           #linetypes= c("solid", "solid"),
           censor = TRUE,
           palette = c("#253494",
                       "#238443"
           ),
           xlab = "Months from Landmark: MRD Timepoint",
           break.time.by = 3,
          # xlim = c(0, 42),
           ylab = "DFS Probability", 
           legend = "none",
           legend.title = "MRD window",
           legend.labs = c("Negative", "Positive"),
           ylim = c(0, 1))
#+ ggtitle("Rectal Baseline")
# surv.median.line = c("hv"))
summary(km_mrd_fit)$table





#color schemes
#253494 is blue
#238443 is green
#"#8BC400" is natera green
#"#008BCE" is natera blue

#Calculate %PFS/OS at 12 mo, 24 mo
summary(km_mrd_fit, times= c(12,24,36))
#Calculate Median PFS and median OS for both groups
survfit(Surv(time = mrd_outcome$DFS.Month, event = mrd_outcome$DFS.Event)~w4.ctDNA, data = mrd_outcome)


#exp(coef) is the HR
#lower .95 is lower bound CI
#upper .95 is upper bound CI


# coxfit hazard ratio for just two factors if there are multiple lines in KM curve
# Filter the data to include only two groups
mrd_subset <- mrd_outcome %>% dplyr::filter(w4.ctDNA %in% c("0", "1")) %>%
  droplevels()   # <-- NEW: drop levels 2–5 so they don't appear as "reference"

# Relevel the grouping variable to set the reference level
mrd_subset$w4.ctDNA <- factor(mrd_subset$w4.ctDNA)
mrd_subset$w4.ctDNA <- relevel(factor(mrd_subset$w4.ctDNA), ref = "0")

# Create the survival object
surv_object <- Surv(time = mrd_subset$DFS.Month, event = mrd_subset$DFS.Event)

# Fit the Cox model
cox_fit <- coxph(surv_object ~ w4.ctDNA, data = mrd_subset)

#view the HR
ggforest(cox_fit,data = mrd_outcome)

# View the summary for Hazard Ratio, CI, and p-value
summary(cox_fit)
}

#DFS_4w (MRD) & pCR
{
mrd_outcome <- read_excel("DFS_4W_KM.xlsx")
mrd_outcome <- mrd_outcome %>%
  filter(Eligible == TRUE)

mrd_outcome


km_mrd <- with(mrd_outcome[which(mrd_outcome$w4.ctDNA.pCR!= "NA"),], Surv(mrd_outcome$DFS.Month, mrd_outcome$DFS.Event))
km_mrd_fit <- survfit(Surv(DFS.Month, DFS.Event) ~ w4.ctDNA.pCR, data = mrd_outcome)
surv_pvalue(km_mrd_fit) #p=0.47
#Plot version 2

ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = TRUE,
           #fun = "cumhaz",
           risk.table = TRUE,
           size = 1,
           # linetype = "strata",
           #linetypes= c("solid", "solid"),
           censor = TRUE,
           palette = c("#253494",
                       "#238443", "orange"
           ),
           xlab = "Months from Landmark: MRD Timepoint",
           break.time.by = 3,
          # xlim = c(0, 42),
           ylab = "DFS Probability", 
           legend = "none",
           legend.title = "MRD window",
           legend.labs = c("ctDNA (-) & pCR", "ctDNA (-) & non-pCR", "ctDNA (+) & non-pCR"),
           ylim = c(0, 1))
#+ ggtitle("Rectal Baseline")
# surv.median.line = c("hv"))
summary(km_mrd_fit)$table





#color schemes
#253494 is blue
#238443 is green
#"#8BC400" is natera green
#"#008BCE" is natera blue

#Calculate %PFS/OS at 12 mo, 24 mo
summary(km_mrd_fit, times= c(12,24,36))
#Calculate Median PFS and median OS for both groups
survfit(Surv(time = mrd_outcome$DFS.Month, event = mrd_outcome$DFS.Event)~w4.ctDNA.pCR, data = mrd_outcome)


#exp(coef) is the HR
#lower .95 is lower bound CI
#upper .95 is upper bound CI


# coxfit hazard ratio for just two factors if there are multiple lines in KM curve
# Filter the data to include only two groups
mrd_subset <- mrd_outcome %>% dplyr::filter(w4.ctDNA.pCR %in% c("1", "2")) %>%
  droplevels()   # <-- NEW: drop levels 2–5 so they don't appear as "reference"

# Relevel the grouping variable to set the reference level
mrd_subset$w4.ctDNA.pCR <- factor(mrd_subset$w4.ctDNA.pCR)
mrd_subset$w4.ctDNA.pCR <- relevel(factor(mrd_subset$w4.ctDNA.pCR), ref = "1")

# Create the survival object
surv_object <- Surv(time = mrd_subset$DFS.Month, event = mrd_subset$DFS.Event)

# Fit the Cox model
cox_fit <- coxph(surv_object ~ w4.ctDNA.pCR, data = mrd_subset)

#view the HR
ggforest(cox_fit,data = mrd_outcome)

# View the summary for Hazard Ratio, CI, and p-value
summary(cox_fit)
}

#DFS_postNAC
{
  mrd_outcome <- read_excel("postNAC_DFS_KM.xlsx")
  mrd_outcome <- mrd_outcome %>%
    filter(Eligible == TRUE)
  
  mrd_outcome
  
  
  km_mrd <- with(mrd_outcome[which(mrd_outcome$postNAC.ctDNA!= "NA"),], Surv(mrd_outcome$DFS.Month, mrd_outcome$DFS.Event))
  km_mrd_fit <- survfit(Surv(DFS.Month, DFS.Event) ~ postNAC.ctDNA, data = mrd_outcome)
  surv_pvalue(km_mrd_fit) #p=0.47
  #Plot version 2
  
  ggsurvplot(km_mrd_fit,
             conf.int = FALSE,
             pval = TRUE,
             #fun = "cumhaz",
             risk.table = TRUE,
             size = 1,
             # linetype = "strata",
             #linetypes= c("solid", "solid"),
             censor = TRUE,
             palette = c("#253494",
                         "#238443"
             ),
             xlab = "Months from Landmark",
             break.time.by = 3,
             #xlim = c(0, 39),
             ylab = "DFS Probability", 
             legend = "none",
             legend.title = "Post-NAC",
             legend.labs = c("Negative", "Positive"),
             ylim = c(0, 1))
  #+ ggtitle("Rectal Baseline")
  # surv.median.line = c("hv"))
  summary(km_mrd_fit)$table
  
  
  
  
  
  #color schemes
  #253494 is blue
  #238443 is green
  #"#8BC400" is natera green
  #"#008BCE" is natera blue
  
  #Calculate %PFS/OS at 12 mo, 24 mo
  summary(km_mrd_fit, times= c(24,36))
  #Calculate Median PFS and median OS for both groups
  survfit(Surv(time = mrd_outcome$DFS.Month, event = mrd_outcome$DFS.Event)~postNAC.ctDNA, data = mrd_outcome)
  
  
  #exp(coef) is the HR
  #lower .95 is lower bound CI
  #upper .95 is upper bound CI
  
  
  # coxfit hazard ratio for just two factors if there are multiple lines in KM curve
  # Filter the data to include only two groups
  mrd_subset <- mrd_outcome %>% dplyr::filter(postNAC.ctDNA %in% c("0", "1")) %>%
    droplevels()   # <-- NEW: drop levels 2–5 so they don't appear as "reference"
  
  # Relevel the grouping variable to set the reference level
  mrd_subset$postNAC.ctDNA <- factor(mrd_subset$postNAC.ctDNA)
  mrd_subset$postNAC.ctDNA <- relevel(factor(mrd_subset$postNAC.ctDNA), ref = "0")
  
  # Create the survival object
  surv_object <- Surv(time = mrd_subset$DFS.Month, event = mrd_subset$DFS.Event)
  
  # Fit the Cox model
  cox_fit <- coxph(surv_object ~ postNAC.ctDNA, data = mrd_subset)
  
  #view the HR
  ggforest(cox_fit,data = mrd_outcome)
  
  # View the summary for Hazard Ratio, CI, and p-value
  summary(cox_fit)
}

#DFS_postCRT
{
mrd_outcome <- read_excel("postCRT_DFS_KM.xlsx")
mrd_outcome <- mrd_outcome %>%
  filter(Eligible == TRUE)

mrd_outcome


km_mrd <- with(mrd_outcome[which(mrd_outcome$PostCRT.ctDNA!= "NA"),], Surv(mrd_outcome$DFS.Month, mrd_outcome$DFS.Event))
km_mrd_fit <- survfit(Surv(DFS.Month, DFS.Event) ~ PostCRT.ctDNA, data = mrd_outcome)
surv_pvalue(km_mrd_fit) #p=0.47
#Plot version 2

ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = TRUE,
           #fun = "cumhaz",
           risk.table = TRUE,
           size = 1,
           # linetype = "strata",
           #linetypes= c("solid", "solid"),
           censor = TRUE,
           palette = c("#253494",
                       "#238443"
           ),
           xlab = "Months from Landmark",
           break.time.by = 3,
           #xlim = c(0, 39),
           ylab = "DFS Probability", 
           legend = "none",
           legend.title = "Post-CRT",
           legend.labs = c("Negative", "Positive"),
           ylim = c(0, 1))
#+ ggtitle("Rectal Baseline")
# surv.median.line = c("hv"))
summary(km_mrd_fit)$table





#color schemes
#253494 is blue
#238443 is green
#"#8BC400" is natera green
#"#008BCE" is natera blue

#Calculate %PFS/OS at 12 mo, 24 mo
summary(km_mrd_fit, times= c(24,36))
#Calculate Median PFS and median OS for both groups
survfit(Surv(time = mrd_outcome$DFS.Month, event = mrd_outcome$DFS.Event)~PostCRT.ctDNA, data = mrd_outcome)


#exp(coef) is the HR
#lower .95 is lower bound CI
#upper .95 is upper bound CI


# coxfit hazard ratio for just two factors if there are multiple lines in KM curve
# Filter the data to include only two groups
mrd_subset <- mrd_outcome %>% dplyr::filter(PostCRT.ctDNA %in% c("0", "1")) %>%
  droplevels()   # <-- NEW: drop levels 2–5 so they don't appear as "reference"

# Relevel the grouping variable to set the reference level
mrd_subset$PostCRT.ctDNA <- factor(mrd_subset$PostCRT.ctDNA)
mrd_subset$PostCRT.ctDNA <- relevel(factor(mrd_subset$PostCRT.ctDNA), ref = "0")

# Create the survival object
surv_object <- Surv(time = mrd_subset$DFS.Month, event = mrd_subset$DFS.Event)

# Fit the Cox model
cox_fit <- coxph(surv_object ~ PostCRT.ctDNA, data = mrd_subset)

#view the HR
ggforest(cox_fit,data = mrd_outcome)

# View the summary for Hazard Ratio, CI, and p-value
summary(cox_fit)
}

#DFS_Surveillance
{
mrd_outcome <- read_excel("Surveillance.xlsx")
mrd_outcome <- mrd_outcome %>%
  filter(Eligible == TRUE)

mrd_outcome


km_mrd <- with(mrd_outcome[which(mrd_outcome$Surveillance.Status!= "NA"),], Surv(mrd_outcome$DFS.Month.Landmark.PostDefTx, mrd_outcome$DFS.Event))
km_mrd_fit <- survfit(Surv(DFS.Month.Landmark.PostDefTx, DFS.Event) ~ Surveillance.Status, data = mrd_outcome)
surv_pvalue(km_mrd_fit) #p=0.47
#Plot version 2

ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = TRUE,
           #fun = "cumhaz",
           risk.table = TRUE,
           size = 1,
           # linetype = "strata",
           #linetypes= c("solid", "solid"),
           censor = TRUE,
           palette = c("#253494",
                       "#238443"
           ),
           xlab = "Months Post Definitive Therapy",
           break.time.by = 3,
           #xlim = c(0, 39),
           ylab = "DFS Probability", 
           legend = "none",
           legend.title = "Surveillance",
           legend.labs = c("Negative", "Positive"),
           ylim = c(0, 1))
#+ ggtitle("Rectal Baseline")
# surv.median.line = c("hv"))
summary(km_mrd_fit)$table





#color schemes
#253494 is blue
#238443 is green
#"#8BC400" is natera green
#"#008BCE" is natera blue

#Calculate %PFS/OS at 12 mo, 24 mo
summary(km_mrd_fit, times= c(12,24,36))
#Calculate Median PFS and median OS for both groups
survfit(Surv(time = mrd_outcome$DFS.Month.Landmark.PostDefTx, event = mrd_outcome$DFS.Event)~Surveillance.Status, data = mrd_outcome)


#exp(coef) is the HR
#lower .95 is lower bound CI
#upper .95 is upper bound CI


# coxfit hazard ratio for just two factors if there are multiple lines in KM curve
# Filter the data to include only two groups
mrd_subset <- mrd_outcome %>% dplyr::filter(Surveillance.Status %in% c("0", "1")) %>%
  droplevels()   # <-- NEW: drop levels 2–5 so they don't appear as "reference"

# Relevel the grouping variable to set the reference level
mrd_subset$Surveillance.Status <- factor(mrd_subset$Surveillance.Status)
mrd_subset$Surveillance.Status <- relevel(factor(mrd_subset$Surveillance.Status), ref = "0")

# Create the survival object
surv_object <- Surv(time = mrd_subset$DFS.Month.Landmark.PostDefTx, event = mrd_subset$DFS.Event)

# Fit the Cox model
cox_fit <- coxph(surv_object ~ Surveillance.Status, data = mrd_subset)

#view the HR
ggforest(cox_fit,data = mrd_outcome)

# View the summary for Hazard Ratio, CI, and p-value
summary(cox_fit)
}

#Time Varying covariate for surveillance HR calculation
{
df <- read.csv("ctdna_timevarying_intervals_with_draw_date_123.csv")

# Order intervals by patient and time
df <- df %>%
  arrange(patient_id, t_start, t_stop)

# Ensure dfs_event = 1 only on the last interval for each patient
df <- df %>%
  group_by(patient_id) %>%
  mutate(dfs_event = ifelse(row_number() == n() & any(dfs_event == 1), 1, 0)) %>%
  ungroup()

# Fit Cox model with time-varying covariate
cox_model <- coxph(
  Surv(t_start, t_stop, dfs_event) ~ ctdna_result + cluster(patient_id),
  data = df
)

# Display results
print(summary(cox_model))

# Optionally display baseline hazard directly (without saving)
base_haz <- basehaz(cox_model, centered = FALSE)
head(base_haz)
}

#Multivariate
{
rm(list=ls()) # clears all environments plots etc so starts the script clean

setwd("~/Documents/Clinical Data Mining/CRC/Rectal/CohortB/")
multi_data <- read_excel("MVA.xlsx")
multi_data <- subset(multi_data, !is.na(ctDNA))
multi_data <- subset(multi_data, Eligible == TRUE)
multi_data <- subset(multi_data, ctDNA %in% c(0, 1))

multi_datadf <- as.data.frame(multi_data)

multi_datadf$Gender <- factor(multi_datadf$Gender, levels = c("0", "1"), labels = c("Female", "Male"))
multi_datadf$Age <- factor(multi_datadf$Age, levels = c("0", "1"), labels = c("<70", "≥70"))
multi_datadf$Tumor <- factor(multi_datadf$`Location`, levels = c("0", "1"), labels = c("Lower Rectum", "Upper Rectum"))
multi_datadf$Neoadjuvant <- factor(multi_datadf$Therapy.Type, levels = c("CRT", "NAC"), labels = c("CRT", "NAC"))
multi_datadf$T.Stage <- factor(multi_datadf$T.Stage, levels = c("0", "1"), labels = c("≤T3", "T4"))
multi_datadf$N.Stage <- factor(multi_datadf$N.Stage, levels = c("0", "1"), labels = c("Negative", "Positive"))
multi_datadf$Adjuvant <- factor(multi_datadf$ACT, levels = c("0", "1"), labels = c("observation", "ACT"))
multi_datadf$RAS <- factor(multi_datadf$RAS, levels = c("0", "1"), labels = c("WT", "MUT"))
multi_datadf$BRAF <- factor(multi_datadf$BRAF, levels = c("0", "1"), labels = c("WT", "MUT"))
multi_datadf$ctDNA <- factor(multi_datadf$ctDNA, levels = c("0", "1"), labels = c("Negative", "Positive"))
multi_datadf$ECOG <- factor(multi_datadf$ECOG, levels = c("0", "1"), labels = c("0", "1"))


surv_object<-Surv(time = multi_datadf$DFS.Month, event = multi_datadf$DFS.Event) 
cox_fit <- coxph(surv_object ~ctDNA + Gender + Age + ECOG+ Tumor + Neoadjuvant + Adjuvant +T.Stage + N.Stage  + RAS, data=multi_datadf) 
# change coxph to coxphf if need to apply firth correction for some factors 
ggforest(cox_fit, data = multi_datadf, main = "Multivariate Regression Model for DFS", refLabel = "Reference Group")
test.ph <- cox.zph(cox_fit)
}
