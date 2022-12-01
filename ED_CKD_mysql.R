##### This code is part of the phenotype analysis. 
##### It is one of three sources of diagnoses; i.e. GP records, the other two being and self-reported diagnoses.
##### The two phenotypes that we will be looking into here, are earache, rash, joint ache/stiffness and ANCA associated vasculitis (ANCA).
##### the diagnosis codes have been saved in .csv files that are also available on this repository.
##### the output of this code is a file that will subsequently be transformed in STATA for easier collation with diagnoses from other sources (code: GP_diagnoses.do).

library(RMySQL)
library(dplyr)
library(tcltk)
library(data.table)
lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)  

pwd <- .rs.askForPassword("Database Password:")


abs_path <- "/slade/home/tj358/pilot"

setwd(abs_path)

file_list <- c("ed_readcodes.txt", 
               "ckd_readcodes.txt", 
               "diabetes_readcodes.txt",
               "cvd_readcodes.txt",
               "hypertension_readcodes.txt",
               "chol_readcodes.txt",
               "ihd_readcodes.txt")

for (file in file_list) {
  
  con=dbConnect(MySQL(),dbname="UKB_GP_RECORDS",password=pwd, user='tj358',host="slade.ex.ac.uk") 
  
  path_to_file <- file.path(abs_path, file)
  
  read_codes <- read.delim(path_to_file, sep = "\t", quote = "", header=F)
  
  read_list <- paste(paste0("'",as.list(read_codes$V1),"'"),collapse = ",") 
  
  query <-"SELECT * FROM UKB_GP_RECORDS.gp_clinical_230K_171019 WHERE READ_3 IN (%s) OR READ_2 IN (%s);"
  
  format_query <- sprintf(query,read_list, read_list)
  
  df = dbSendQuery(con, format_query) %>% fetch(rsa, n=-1)
  
  write.table(df, file = gsub("_readcodes.txt", "_gpcases.tsv", path_to_file), sep= "\t", quote=FALSE, row.names=FALSE) 
  
  lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)
  
}


lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)  

con=dbConnect(MySQL(),dbname="UKB_GP_RECORDS",password=pwd, user='tj358',host="slade.ex.ac.uk") 


query <- "select distinct(n_eid) from gp_clinical_230K_171019;"
query2 <- "select distinct(n_eid) from gp_registrations_230K_171019;"

df = dbSendQuery(con, query) %>% fetch(rsa, n=-1)
df2 = dbSendQuery(con, query2) %>% fetch(rsa, n=-1)

lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)

df_add <- as.data.frame(df2[!df2$n_eid %in% df$n_eid,])
df_add$n_eid <- df_add$`df2[!df2$n_eid %in% df$n_eid, ]`
df_add$`df2[!df2$n_eid %in% df$n_eid, ]` <- NULL
df_all <- rbind(df, df_add)

write.csv(df_all, "~/pilot/primarycaredata_available.csv", row.names = F)

