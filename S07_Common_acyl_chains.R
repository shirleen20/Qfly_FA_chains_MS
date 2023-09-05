#_______________________________02/06/2023_______________________________

# Analysis to identify common acyl chains present in lipid subclasses and their % abundance in Day1 and Day19 B.tryoni males
# Use output from here to generate S9, S10, S11 for lipid MS1 and main text figure 4 which is the heat map we generate in S09.R

library("tidyverse")
library("readxl")
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
  mutate(LogAreaCheck = Intercept + Slope*(log10(Conc)), AreaCheck = 10^LogAreaCheck) %>% ### Just to double check
  select(Name, Weight, Samples, Line, Cage, Conc, species, Carbon, Bond, SubClass, Class, Batch, Age, Time)

#__________________________________________________________

# start analysis

#First remove subclasses PSe and DGe from DF because these two subclasses have all lipid species that do not have acyl chain identified
DF <- DF.main %>% 
  dplyr::filter(SubClass != "PS e" & SubClass != "DG e")

#grep - finds things, gives index 1, 3, 5 not showing 1 is apple, 3 is pear uless you specify value =T) check help ?gsub or grep
#gsub is a substitution

### Function to create for each subclass data frame, all side chains separated into columns
# X = data frame to process
# Y = Subclass level as string

F1 <- function(X, Y){
  Subclass <- X %>% 
    filter(SubClass == Y)
  length1 <- length(Subclass)  # The starting number of columns
  Names <- colnames(Subclass) # Save column names
  Subclass_SC <- gsub(".*\\((.+)\\).*","\\1",Subclass$Name) %>%  ### keep only string within "( )"
    str_remove_all("[ep]") %>%           
    strsplit("_") ### split by "_"
  SCNumber <- Subclass_SC %>% map(length) %>% unlist ### find the number of side chains in each lipid for each individual (want to remove those with no side chains identified)
  RemoveSpecies <- which(SCNumber < max(SCNumber)) ### vector of indices of those without side chains identified. They will have length less than the maximum number of side chains
  Subclass_SC <- do.call(rbind, Subclass_SC) %>% as.data.frame() ### 
  Subclass <- cbind(Subclass, Subclass_SC) ### Join the new columns to the old data frame
  length2 <- length(Subclass) ### New number of columns
  colnames(Subclass) <- c(Names, paste("SC", 1:(length2-length1), sep = ""))  ### Rename new colunms with SC1, SC2 ... etc length2-length1 is the total number of new columns
  if (length(RemoveSpecies) == 0){
    return(Subclass)
  }
  else{
    return(Subclass[-RemoveSpecies,]) ### Output new data frame after removing rows with unindentified side chains
  }
}


# X is the data frame to process
# Subclass is the name of the subclass variable as a string (default: "SubClass")
F2 <- function(X, Subclass = "SubClass"){
  Subclass_names <- unique(X[[Subclass]]) ### All the unique Subclasses
  Out_list <- setNames(replicate(length(Subclass_names), data.frame()), Subclass_names) ### create a list of empty data frames with names Subclass_names with length of the number of subclasses 
  # if using vector(mode = "list", length = length(Subclass_names)) ... required renaming objects in Out_list, see next line
  # names(Out_list) <- Subclass_names ## rename each object within the list with the subclass names
  for (i in Subclass_names){ ### For each separate subclass
    Out_list[[i]] <- F1(X, i) ### in each object of that ith subclass, run F1 function on data frame X with subclass i
  }
  return(Out_list)
}

LIST1 <- F2(X = DF)

F3 <- function(X, Samples = "Samples", SubClass = "SubClass", Age = "Age", Time = "Time"){
  LIST <- F2(X, Subclass = SubClass)
  Subclass_names <- unique(X[[SubClass]])
  Sample_names <- unique(X[[Samples]])
  Ages <- X[[Age]]
  Times <- X[[Time]]
  OriginalLength <- length(X)
  Outdf <- data.frame(Samples=c(), Age = c(), Line = c(), Time = c(), SubClass = c(), SideChain = c(), N_SC = c(), P1_SC = c(), Abundance_SC = c(), Proportion = c())
  SampleNoLipid <- data.frame(SubClass = c(), Samples = c(), Age = c(), Time = c(), Status = c())
  #print(SampleNoLipid)
  for (i in Subclass_names){
    NumberSC <- length(LIST[[i]]) - OriginalLength #NumberSC is length of the list i.e the nbr of columns
    for (j in Sample_names){
      df <- LIST[[i]] %>% filter(Samples == j)
      
      if (nrow(df) == 0){
        if (length(SampleNoLipid) == 0){
          
          SampleNoLipid <- tibble(SubClass_names = i, Samples = j,
                                  Age = Ages[X[["Samples"]] == j] %>%  unique(),
                                  Time = Times[X[["Samples"]] == j] %>%  unique(), 
                                  Status = "No lipids")
        }
        else {
          SampleNoLipid <- rbind(SampleNoLipid,c(i,j, Ages[X[["Samples"]] == j] %>% unique(),
                                                 Times[X[["Samples"]] == j] %>%  unique(), "No lipids"))
          #print(unique(X[["Age"]][X[["Samples"]]==j & X[["SubClass"]] == i]))
          
          #print(SampleNoLipid)
        }
        next
      }
      else {
        All_SC <- df[,-1:-OriginalLength] %>% as.matrix() %>% t() %>% as.vector()
        Samples2 = j
        Age2 = unique(df$Age)
        Line2 = unique(df$Line)
        Time2 = unique(df$Time)
        SubClass2 = i
        #rint(c(i,j))
        IntermediateDF <- data.frame(S = Samples2, A = Age2, L = Line2, T = Time2, SubClass = SubClass2, table(All_SC), P1 = 0)
        #print(IntermediateDF)
        names(IntermediateDF) <- c("Samples","Age","Line","Time", "SubClass","SideChain","N_SC", "P1_SC")
        IntermediateDF$P1_SC <- IntermediateDF$N_SC/sum(IntermediateDF$N_SC)
        #print(IntermediateDF)
        Weights <- rep(df$Conc, each = NumberSC)
        #print(Weights)
        #print(All_SC)
        #return(data.frame(SideChain = All_SC, wt = Weights))
        Abundance <- data.frame(SideChain = All_SC, wt = Weights) %>% dplyr::count(SideChain,wt=wt)
        #print(Abundance)
        names(Abundance) <- c("SideChain", "Abundance_SC")
        #Abundance_SC <- (data.frame(var = All_SC, wt = Weights) %>% count(var,wt=wt))$n
        #print(data.frame(var = All_SC, wt = Weights) %>% count(var,wt=wt))
        #IntermediateDF$Abundance_SC <- Abundance_SC
        IntermediateDF <- IntermediateDF %>% full_join(Abundance, by = "SideChain")
        IntermediateDF$Proportion <- (IntermediateDF$Abundance_SC/(sum(df$Conc)*NumberSC))*100
        Outdf <- rbind(Outdf, IntermediateDF)
      }
    }
  }
  return(list(Outdf, SampleNoLipid))
}

SCTable <- F3(DF)

### Include rows of zeroes where samples do not have a reading for SubClass by SideChain
Complete_SCTable <- SCTable[[1]] %>% complete(nesting(SubClass, SideChain), nesting(Samples,Age, Line, Time),
                                              fill = list(N_SC = 0, P1_SC = 0, Abundance_SC = 0, Proportion = 0)) %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, SideChain, N_SC, P1_SC, Abundance_SC, Proportion)

# add another column to dataframe SCTable called type that has lipid categories
phospholipids <- c("CL","PC", "PE", "PG", "PI", "PS","LPC")
etherlipidsPL <- c("PE p", "PE e", "PC e", "PS e")
etherlipidsNL <- c("TG e", "DG e")
neutrallipids <- c("TG", "DG")

DF1 <- Complete_SCTable %>% 
  dplyr::mutate(Type = ifelse(SubClass %in% phospholipids, "phospholipids", 
                              ifelse(SubClass %in% etherlipidsPL, "etherlipidsPL",
                                     ifelse(SubClass %in% etherlipidsNL, "etherlipidsNL",
                                            "neutrallipids"))))


#________________________________________________________

# calculate the mean proportions of common acyl chains present in each lipid subclass to generate Table S9 of lipid MS1 and the heatmap in script07.R
# (Output for Table S9 I have saved in script07.R after further pocessing the data to filter out data for s06 only)
Abundance_Subclass1 <- Complete_SCTable %>%  
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, SideChain, Proportion) %>% #the proportion here is mean computed earlier for a sample
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass", "SideChain"), measurevar = "Proportion") %>% 
  dplyr::select(Line, Time, Age, SubClass, SideChain, Proportion, se) %>% # here proportion refers to mean percentage
  dplyr::mutate(Proportion = format(round(as.numeric(.$Proportion), 4), nsmall = 4)) %>% 
  #dplyr::mutate(se = format(round(as.numeric(.$se), 4), nsmall = 4)) %>% 
  #unite("mean.se", Proportion:se, sep = "±") %>% 
  dplyr::select(Line, Time, Age, SubClass, SideChain, Proportion) 

sum(as.numeric((Abundance_Subclass1 %>% dplyr::filter(Line == "s06" & SubClass == "TG e" & Age == "1 Day"))$Proportion)) 

# check to see if each Line by Subclass by Age comes to 100%
# Abundance_Subclass1 %>% dplyr::filter(Line == "s06" & SubClass == "TG e" & Age == "19 Day")
# we see that TGe doesn't come to 100% because 2 samples don't have TGe present: check using this code
# View(Complete_SCTable %>% ungroup() %>% group_by(Samples, Age, Line, Time, SubClass) %>% dplyr::summarise(SUM=sum(Proportion)))

# Prefilter_AbundanceSubclass is used to generate Table S9: Common acyl chains present in lipid subclasses and their percent abundance in Day 1 and Day 19 B. tryoni males

Prefilter_AbundanceSubclass <- Abundance_Subclass1 %>% # we have called this prefilter because later we might apply 5% or 10% on column Proportions
  dplyr::mutate(Proportion = as.numeric(Proportion)) %>% 
  pivot_wider(names_from = c( Age, SubClass), values_from = Proportion) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% # replace NAs with 0
  dplyr::filter(Line == "s06") %>% 
  write_excel_csv("Final_Scripts/S06commonacylchains.csv")   


#_________________________________________________

# Calculate mean total percentage for each of the most abundant side chains using percentages for total lipids
# Column Total Lipids
FinalSC_stats1 <- DF1 %>% 
  dplyr::group_by(Samples, SideChain) %>% 
  dplyr::mutate(Numerator = sum(Abundance_SC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Samples) %>%
  dplyr::mutate(Denominator = sum(Abundance_SC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SideChain, Numerator, Denominator) %>% 
  dplyr::mutate(T.P = (Numerator/Denominator)*100) %>% # T.P is total percentage should come to 100 for each individual sample to check: sum((FinalSC_stats1 %>% dplyr::filter(Samples == "cbr_n04"))$T.P)
  unique() %>% 
  dplyr::select(-Numerator, -Denominator)

# Calculate mean total percentages for side chains in each of the four lipid categories
FinalSC_stats2 <- DF1 %>% 
  dplyr::group_by(Samples, SideChain, Type) %>% 
  dplyr::mutate(Numerator = sum(Abundance_SC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(Samples, Type) %>%
  dplyr::mutate(Denominator = sum(Abundance_SC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Line, Time, Age, Type, SideChain, Numerator, Denominator) %>% 
  dplyr::mutate(TP.type = (Numerator/Denominator)*100) %>% #T.P_type is total percentage should come to 100 for each individual category in a sample to check:  sum((FinalSC_stats2 %>% filter(Type == "neutrallipids" & Samples == "cbr_n04"))$TP.type)
  unique() %>% 
  dplyr::select(-Numerator, -Denominator) %>% 
  mutate_at(vars(TP.type), ~replace(., is.nan(.), 0)) # I added this line of codes to replace NaNs with 0

FinalSC2_summary <- FinalSC_stats2 %>% 
  dplyr::group_by(Age, Line, Time, SideChain, Type) %>% 
  dplyr::mutate(meanTP.type = mean(TP.type), seTP.type = sd(TP.type)/sqrt(length(TP.type))) %>% 
  dplyr::ungroup() %>% 
  #mutate(meanTP.type = format(round(.$meanTP.type, 5), nsmall = 5)) %>% 
  #mutate(seTP.type = format(round(.$seTP.type, 5), nsmall = 5)) %>% 
  dplyr::select(-Samples, -TP.type) %>% 
  unique()

# make dataframe  
Combine1 <- FinalSC_stats1 %>% 
  dplyr::group_by(Age, Line, Time, SideChain) %>% 
  dplyr::mutate(mean = mean(T.P), se = sd(T.P)/sqrt(length(T.P))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-Samples, -T.P) %>% 
  unique() %>% 
  dplyr::mutate(Type = "TotalLipids") %>% 
  dplyr::select(Age, Line, Time, Type, SideChain, mean, se)

Combine2 <- FinalSC2_summary %>% 
  dplyr::rename('mean' = meanTP.type, 'se' = seTP.type)

Combined <- rbind(Combine1, Combine2)

Combined %>% dplyr::filter(Age == "1 Day" & Line == "cbr" & Time == "New") %>% dplyr::select(-se) %>% 
  pivot_wider(names_from = "Type", values_from = "mean")

df_se <- Combined %>% dplyr::mutate(LineC = paste(Line,Age,Time,Type, sep = "_")) %>%  dplyr::select(-mean)

df_mean2 <- Combined %>% dplyr::mutate(LineC = paste(Line,Age,Time,Type, sep = "_")) %>% 
  dplyr::select(-se, -Type, -Age, -Line, -Time) %>% 
  pivot_wider(names_from = "LineC", values_from = "mean") %>% 
  dplyr::filter_at(vars(-SideChain), any_vars(.>2)) %>% 
  pivot_longer(col = -SideChain, names_to = "LineC", values_to = "mean") %>% 
  left_join(df_se, by = c("SideChain", "LineC")) %>% 
  dplyr::select(Age, Line, Time, Type, SideChain, mean, se)

Filtered_type <- df_mean2 %>%
  dplyr::filter(!is.na(mean) & !is.na(se)) %>% 
  dplyr::mutate(mean = format(round(as.numeric(.$mean), 4), nsmall = 4)) %>% 
  dplyr::mutate(se = format(round(as.numeric(.$se), 4), nsmall = 4)) %>% 
  unite("mean.se", mean:se, sep = "±") %>% 
  dplyr::select(Age, Line, Time, SideChain, Type, mean.se) %>% 
  unique() %>% 
  pivot_wider(names_from = SideChain, values_from = mean.se) # this retains data retains columns for those sidechains that were present after 2% filtering in any of the three types of lipids

# 1. Column neutral lipids (Table S10 in lipids MS1)

Filtered_type_NeutralLipids <- Filtered_type %>% 
  filter(Type == "neutrallipids") %>% 
  dplyr::filter(Line == "s06") %>% 
  write_excel_csv("S06NLchains.csv")

# 2. Column Phospholipids (Table S11 in lipids MS1)

Filtered_type_PhosphoLipids <- Filtered_type %>% 
  filter(Type == "phospholipids") %>%
  dplyr::filter(Line == "s06") %>% 
  write_excel_csv("S06PLchains.csv")   

# 3. Column ether lipids 

Filtered_type_EtherLipids <- Filtered_type %>% 
  filter(Type == "etherlipids") 


#________________________________________END_________________________________________
