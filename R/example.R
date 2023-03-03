require(haven)
source("R/library.R")

# prepare data ------------------------------------------------------------

# data2021_ori = read_sas("~/K/Projects/Oversights/Data/FACREP/PYF/shr_data_2021_unq_hosp.sas7bdat")

data2021_ori = read_sas("K:/Projects/Oversights/Data/FACREP/PYF/shr_data_2021_unq_hosp.sas7bdat")
data = data2021_ori
data2021 = data2021_ori
set.seed(22110122)
missing_ind = ceiling(runif(400, 0, 912559))
data2021[missing_ind[1:100], "T_start"] = NA
data2021[missing_ind[101:200], "t_end"] = NA
data2021[missing_ind[201:300], "Pfemale"] = NA
data2021[missing_ind[301:400], "provfs"] = NA

input_var = c("Pdiab", "Pdismiss", "Pfemale", "bmi_cat", "Age_period",
              "ashd1", "othcardiac1",
              "carfail", "noambul",
              "pulmon", "notrans", "cancer", "diabetes", "pvasc", "cva", "smoke", "alcoh", "drug", "inci_one",
              "inci_miss",
              paste0("ccs_", c(1:91)), "pcom_miss", "PYF_period_ESRD",
              "p_day_ma",
              "NH_type_lt90", "NH_type_ge90")
result = fe.data.prep(data2021, T1.char = "T_start", T2.char = "t_end", event.prefix = "hosp_beg_", event.num = 37,
             event.char = NULL, nevents.char = NULL, Zti.char = input_var, Ztd.char = NULL,
             ttd.char = c("covid_start", "covid_start10days", "covid_post_start", "covid_persistent_start"),
             prov.char = "provfs")

