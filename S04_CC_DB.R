#____________________________________________01/06/2023___________________________________________
# Codes to generate Table 2 of lipids MS1: Average numbers of carbons and double bonds per acyl chain across 
# lipids in each subclass in Day 1 and Day 19 B.tryoni males 
# Calculating the mean and standard error for subclass abundance across samples by line by age by time

library("tidyverse")
library("ggrepel")
library("ggpubr")
library("readxl")
library("gtools")
library("broom")
library("ggpubr")
library("Rmisc")

rm(list=ls())

#_________________Outlier function to be used, when required______________________________ 
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

#__________________Create a list of ambiguous lipid species from Shirleen's csv file______
ALS <- read_csv("data/Lipids to filter ploty.csv") %>% 
  dplyr::select(3) %>% 
  unique() %>% 
  c()
#______________ Load samples and remove outliers from batch02 ______________________________
# correct names of incorrectly labeled sample names
'%!in%' <- function(x,y)!('%in%'(x,y)) # a function to create opposite of %in% for filtering out

B01 <- read_csv("data/CD_results_Batch01_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))%>% 
  dplyr::filter(samples != "ct_n09")

B02 <- read_csv("data/CD_results_Batch02_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(samples != "S06_95") %>% 
  dplyr::mutate(samples = recode(samples,  syd_n87 = "syd_o87",syd_n104 = "syd_o104"))  

B03 <- read_csv("data/CD_results_Batch03_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>%
  dplyr::filter(samples != "syd_n134") %>% # remove outlier
  dplyr::mutate(samples = recode(samples, ct_n161 = "ct_o161", S06_171 = "s06_171",
                                 syd_160 = "syd_o160", syd_n138 = "syd_o138", syd_n169 = "syd_o169")) %>% 
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))

B04 <- read_csv("data/CD_results_Batch04_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(samples != "s06_217") %>% # remove outlier
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))

#_________________ calculate total area of all peaks in each sample and calculate their mean ___________________
F01 <- function(x){
  x %>% 
    dplyr::group_by(samples) %>% 
    dplyr::mutate(S.Area = sum(Area)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(samples, S.Area) %>% 
    unique() %>% 
    dplyr::mutate(M.S.Area = mean(S.Area)) %>% 
    dplyr::select(M.S.Area) %>% 
    slice(1)
}

MSA.B01 <- F01(x = B01)
MSA.B02 <- F01(x = B02)
MSA.B03 <- F01(x = B03)
MSA.B04 <- F01(x = B04)

NP.B01 <- B01 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B01))) %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

NP.B02 <- B02 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B02))) %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100)  

NP.B03 <- B03 %>% 
  dplyr::mutate(N.Area = Area)  %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

NP.B04 <- B04 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B04)))  %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

DF1 <- NP.B01 %>% 
  dplyr::bind_rows(NP.B02) %>% 
  dplyr::bind_rows(NP.B03) %>%
  dplyr::bind_rows(NP.B04) 
#______________________________________________________________________________________________________

#__________________________ Create main data frame _________________________________
Variables <- read_excel("data/Variables.xlsx") %>% 
  dplyr::select(ID, Weight, Time2) %>% 
  dplyr::rename(samples = ID) 

DF <- DF1 %>% 
  dplyr::left_join(., Variables, by = "samples") %>% 
  dplyr::filter(Weight < 7.5) %>% 
  dplyr::filter(Weight > 1) %>% 
  dplyr::select(Name, Weight, samples, N.Area, count, class, species, batch, DB_new, class_new) %>%
  dplyr::mutate(Bond = DB_new) %>% 
  dplyr::mutate(SubClass = class_new) %>% 
  dplyr::select(-DB_new, -class_new) %>% 
  tidyr::separate(samples, into = c("Line", "Cage"), remove = FALSE) %>% 
  dplyr::mutate(Line = replace(Line, Line == "S06", "s06")) %>% 
  dplyr::mutate(Class = class) %>% 
  dplyr::mutate(Batch = batch) %>%
  dplyr::select(-batch) %>% 
  dplyr::rename(Samples = samples) %>% 
  dplyr::mutate(Age = if_else(Batch %in% c("B01", "B02"), "19 Day", "1 Day")) %>% 
  dplyr::mutate(Time = str_extract(Cage, "[a-z]+")) %>% 
  dplyr::mutate(Time = replace_na(Time, "o")) %>% 
  dplyr::mutate(Time = if_else(Time %in% c("o"), "Old", "New")) %>% 
  dplyr::mutate(Time = replace(Time, Line == "s06", "Older")) %>% 
  dplyr::select(-count, -class) %>% 
  tidyr::separate(col = species, into = c("Carbon", NA), ":", remove = FALSE) %>% 
  mutate(Carbon = as.double(Carbon))

#___________________

# Import standards data

df.old <- read_csv("Table_ConcAllStandards.csv") %>% # this is data from 1st run
  dplyr::select(Name, Class, area, Conc, Batch) %>%
  dplyr::filter(Class != "Cer" & Class != "FA" & Class != "SM" & Class != "PCp") %>% 
  dplyr::group_by(Name, Class, Conc, Batch) %>% 
  dplyr::summarise(Area = mean(area)) %>% # calculates total area for any class that had more than one standard 
  dplyr::ungroup() %>%
  unique()  %>% 
  dplyr::select(Conc, Name, Class, Area)

# Add manually calculate exp 2 CL areas to DF
DF.CL <- tibble(Conc = c(0.01, 0.10, 1.00, 10.00), 
                Name = rep("CL 72:4",4), Class = rep("CL", 4),
                Area = c(1023.129827,32570.02475,579994.8638,45993575.15))

DF.Stds <- rbind(df.old, DF.CL)

Regression <- DF.Stds %>% split(~Class) %>%
  map(~lm(log10(Area) ~ log10(Conc), data = .x)$coefficient) ### Obtain all the regression coefficient

Regression_coef <- tibble(Class2 = names(Regression) ,
                          bind_rows(lapply(Regression, as.data.frame.list))) %>% ### 
  dplyr::rename(Intercept =  "X.Intercept.",  Slope = "log10.Conc.") %>% 
  mutate(Class2 = if_else(Class2 =="PEp", "PE p", Class2))

#____________create main DF_____________________

DF.main <- DF %>%
  ungroup() %>% 
  mutate(Class2 = if_else(grepl(" e", SubClass), Class, SubClass)) %>% ## All ether lipids use ester lipid standards
  left_join(Regression_coef, by = "Class2") %>%
  mutate(Conc = 10^((log10(N.Area)-Intercept)/Slope)) %>%
  mutate(LogAreaCheck = Intercept + Slope*(log10(Conc)), AreaCheck = 10^LogAreaCheck) %>%  ### Just to double check
  dplyr::select(-N.Area, -Class2, -Intercept, -Slope, -LogAreaCheck,-AreaCheck) 

#____________________create dataframe for mean carbon and double bond calculation_____________________

Complete_DF_new <- DF.main %>% ungroup() %>% complete(nesting(SubClass, Bond, Carbon, Name), nesting(Samples, Age, Line, Time),
                                                 fill = list(Conc = 0)) %>%
  select(Samples, Age, Line, Time, SubClass, Name, Conc, Bond, Carbon)

### DATA is a data frame containing the Conc (but must be separated by SubClass)
### RemoveUnresolved, if TRUE remove from calculation the non-disambiguated lipids
Compute_Average <- function(DATA, RemoveUnresolved = FALSE){
  DATA$Name <- gsub(".*\\((.+)\\).*", "\\1", DATA$Name)
  NameSplit <- str_split(DATA$Name,"_") ### Split the lipid species... creates a list of vectors
  SCNum <- NameSplit %>% map(~length(.)) %>% unlist() ### determine the number of side chains
  DATA$SCNum <- SCNum
  MaxSCNum <- max(SCNum) ### Should give the actual number of side chains in the subclass
  DATA$Disambiguated <- ifelse(DATA$SCNum == MaxSCNum, TRUE, FALSE) ### If SC number is less than maximum side chain number, it must be Unresolved
  DATA$SCNum <- MaxSCNum ### Change all SCNum to the actual number which is MaxSCNum (for later processing purposes)
  if (RemoveUnresolved){
    return(DATA %>% filter(Disambiguated == TRUE) %>% ### remove Unresolved lipid species
             dplyr::group_by(Samples, SubClass) %>% 
             dplyr::mutate(TotalAbundance = sum(Conc)) %>% 
             dplyr::mutate(AveAbundance = if_else(TotalAbundance == 0, 0,Conc/TotalAbundance)) %>% 
             dplyr::ungroup() %>% 
             dplyr::group_by(Age, Line, Time, SubClass, Name, SCNum, Disambiguated, Bond, Carbon) %>% 
             dplyr::summarise(Average=mean(AveAbundance))
    )
  }
  else {
    return(DATA %>% 
             dplyr::group_by(Samples, SubClass) %>% 
             dplyr::mutate(TotalAbundance = sum(Conc)) %>% 
             dplyr::mutate(AveAbundance = if_else(TotalAbundance == 0, 0,Conc/TotalAbundance)) %>% 
             dplyr::ungroup() %>% 
             dplyr::group_by(Age, Line, Time, SubClass, Name, SCNum, Disambiguated, Bond, Carbon) %>% 
             dplyr::summarise(Average=mean(AveAbundance))
    )
  }
}

Complete_List <- Complete_DF_new %>% split(~SubClass) %>% # must split into list by SubClass to allow using the function
  map(~Compute_Average(., RemoveUnresolved = FALSE))

# Creates a list after removing the unresolved lipid species
Complete_removeUnresolved <- Complete_DF_new %>% split(~SubClass) %>% map(~Compute_Average(., RemoveUnresolved = TRUE))

### NameSplit is a vector. Lipid species that are non-disambiguated/unresolved will only have one value
### This will temporarily duplicate the non-disambiguated/unresolved side chain for ease of processing purposes
### The values will be corrected in the later script

FillSC <- function(NameSplit, MaxSCNum){
  if (length(NameSplit) < MaxSCNum & length(NameSplit) == 1){
    Duplicate <- gsub("(.+)e+", "\\1", NameSplit)
    NewName <- c(NameSplit, rep(Duplicate, MaxSCNum - 1))
  } else {
    NewName <- NameSplit
  }
  return(NewName)
}


### DATA is the average abundance of lipid species within one subclass
### Averages were computed by age, line and time
###
SideChains <- function(DATA){
  NameSplit <- str_split(DATA$Name, "_") %>% 
    map(~FillSC(.,unique(DATA$SCNum))) ### Duplicate non-disambiguated side chains
  for (i in 1:unique(DATA$SCNum)){
    DATA[,paste("SC",i,sep="_")] <- NameSplit %>% map(i) %>% unlist() ### Fill in columns of side chain
  }
  NewDat <- DATA %>% 
    pivot_longer(cols = starts_with("SC_"), names_to = "SideChainNum", values_to = "SCType") %>%  ### Pivot longer SC columns into one
    dplyr::mutate(Type = if_else(grepl("e", SCType), "Ether",
                                 if_else(grepl("p", SCType), "Vinyl", "Ester"))) %>% 
    dplyr::mutate(CarbonNum = gsub("(\\d+):(\\d+)[ep]*$", "\\1", SCType) %>% as.numeric()) %>%
    dplyr::mutate(DoubleBond = gsub("(\\d+):(\\d+)[ep]*$", "\\2", SCType) %>% as.numeric()) %>%
    dplyr::mutate(CarbonNum = if_else(Disambiguated == FALSE, CarbonNum/SCNum, CarbonNum)) %>%
    dplyr::mutate(DoubleBond = if_else(Disambiguated == FALSE, DoubleBond/SCNum, DoubleBond)) 
  return(NewDat)
}

### Standard error

wtd.sterror <- function(Vect, Weight){
  NormWt <- Weight/sum(Weight) ### normalised weights
  Mean <- Hmisc::wtd.mean(Vect, NormWt)
  Numerator <- sum(NormWt*(Vect - Mean)^2)
  Denominator <- 1 - sum(NormWt^2)
  Variance <- Numerator/Denominator
  #print(Variance)
  SquareSE <- Variance*sum(NormWt^2) ###
  #print(SquareSE)
  return(sqrt(SquareSE))
}

wtd.sterror2 <- function(Vect, Weight){
  NormWt <- Weight/sum(Weight)
  Mean <- Hmisc::wtd.mean(Vect, NormWt)
  Variance <- Hmisc::wtd.var(Vect, NormWt, normwt = TRUE)
  #print(Variance)
  SquareSE <- Variance*sum(NormWt^2)
  #print(SquareSE)
  return(sqrt(SquareSE))
}

Test1 <- Complete_List %>% map(~SideChains(.))
### Use Hmisc weighted mean and normalised weighted variance
### normalised weights such that sum of weight = 1, 
# necessary since we can't know for sure the true abundance within a given population

Test <- Complete_List %>% map(~SideChains(.)) %>% dplyr::bind_rows() %>% 
  dplyr::group_by(Age,Line,Time,SubClass) %>% 
  dplyr::summarise(MEAN_Carbon = Hmisc::wtd.mean(CarbonNum, Average), 
                   SE_Carbon = wtd.sterror(CarbonNum, Average),  ###
                   MEAN_DB = Hmisc::wtd.mean(DoubleBond, Average), 
                   SE_DB = wtd.sterror(DoubleBond, Average)) %>% 
  ungroup()

Test2 <- Complete_removeUnresolved %>% map(~SideChains(.)) %>% dplyr::bind_rows() %>% 
  dplyr::group_by(Age,Line,Time,SubClass) %>% 
  dplyr::summarise(MEAN_Carbon_rm = Hmisc::wtd.mean(CarbonNum, Average), 
                   SE_Carbon_rm = wtd.sterror(CarbonNum, Average), 
                   MEAN_DB_rm = Hmisc::wtd.mean(DoubleBond, Average), 
                   SE_DB_rm = wtd.sterror(DoubleBond, Average)) %>% 
  ungroup()


Final <- full_join(Test,Test2) %>% 
  filter(Line=="s06")

EtherVinyl <- Complete_List %>% map(~SideChains(.)) %>% dplyr::bind_rows() %>% filter(Type != "Ester") %>% 
  dplyr::group_by(Age,Line,Time,SubClass) %>% 
  dplyr::summarise(MEAN_Carbon = Hmisc::wtd.mean(CarbonNum, Average), 
                   SE_Carbon = wtd.sterror(CarbonNum, Average),  ###
                   MEAN_DB = Hmisc::wtd.mean(DoubleBond, Average), 
                   SE_DB = wtd.sterror(DoubleBond, Average))%>% 
  ungroup()

EtherVinyl2 <- Complete_removeUnresolved %>% map(~SideChains(.)) %>% dplyr::bind_rows() %>% filter(Type != "Ester") %>%
  dplyr::group_by(Age,Line,Time,SubClass) %>% 
  dplyr::summarise(MEAN_Carbon_rm = Hmisc::wtd.mean(CarbonNum, Average), 
                   SE_Carbon_rm = wtd.sterror(CarbonNum, Average), 
                   MEAN_DB_rm = Hmisc::wtd.mean(DoubleBond, Average), 
                   SE_DB_rm = wtd.sterror(DoubleBond, Average))%>% 
  ungroup()


FinalEtherVinyl <- full_join(EtherVinyl,EtherVinyl2)%>% 
  filter(Line=="s06")

Ester <- Complete_List %>% map(~SideChains(.)) %>% dplyr::bind_rows() %>% filter(Type == "Ester") %>% 
  dplyr::group_by(Age,Line,Time,SubClass) %>% 
  dplyr::summarise(MEAN_Carbon = Hmisc::wtd.mean(CarbonNum, Average), 
                   SE_Carbon = wtd.sterror(CarbonNum, Average),  ###
                   MEAN_DB = Hmisc::wtd.mean(DoubleBond, Average), 
                   SE_DB = wtd.sterror(DoubleBond, Average))%>% 
  ungroup()

Ester2 <- Complete_removeUnresolved %>% map(~SideChains(.)) %>% dplyr::bind_rows() %>% filter(Type == "Ester") %>%
  dplyr::group_by(Age,Line,Time,SubClass) %>% 
  dplyr::summarise(MEAN_Carbon_rm = Hmisc::wtd.mean(CarbonNum, Average), 
                   SE_Carbon_rm = wtd.sterror(CarbonNum, Average), 
                   MEAN_DB_rm = Hmisc::wtd.mean(DoubleBond, Average), 
                   SE_DB_rm = wtd.sterror(DoubleBond, Average))%>% 
  ungroup()

FinalEster <- full_join(Ester,Ester2)%>% 
  filter(Line=="s06")

write.csv(Final, "WeightedMeanVariance_Carbon&DoubleBond.csv", row.names = FALSE)
write.csv(FinalEtherVinyl, "WeightedMeanVariance_EtherVinyl_Carbon&DoubleBond.csv", row.names = FALSE)
write.csv(FinalEster, "WeightedMeanVariance_Ester_Carbon&DoubleBond.csv", row.names = FALSE)

#________________________________________END_________________________________________

