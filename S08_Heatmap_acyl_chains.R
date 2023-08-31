#______________________________________________02/06/2023______________________________________________

# Generate Panel A and Panel B heat map for Figure 4 in lipid MS 1 using data obtained from Script07.R 
# on the mean proportions of common acyl chains present in each lipid subclass using dataframe saved as 
# Prefilter_AbundanceSubclass.csv

library(tidyverse)
library(ggplot2)
library(readxl)
library(viridis) # for color Brewer palette

rm(list=ls())

#___Create Panel A of Figure 4: Heatmap showing most common acyl chains in lipid subclasses based on <5% inclusion threshold in Day 1 and Day 19 B. tryoni males____________________

Abundance_Subclass <- as.data.frame(read_csv("MS1Standards/S06commonacylchains.csv"))

S06Abundance <- Abundance_Subclass %>% filter(Line == "s06") 

S06Abundance1 <- S06Abundance %>%
  mutate(SideChain = gsub(":","\\.", SideChain) %>% as.numeric() ) %>%  # Substitute ":" with "." and treat as numeric so that we can treat chains as number and arrange them in descending order, those smaller than 10 were appearing later when we treated them as character for e.g. 8:0 appeared after 10:0
  arrange(desc(SideChain)) %>%  # arrange in descending order
  mutate(SideChain = format(SideChain, nsmall = 1) %>% gsub(" ","",.) %>% gsub("\\.",":",.) ) %>%
  mutate(SideChain = as.factor(SideChain)) %>%
  mutate(SideChain = factor(SideChain, levels = unique(SideChain) %>%
                              gsub(":","\\.", .) %>%
                              as.numeric() %>%
                              sort(decreasing = TRUE) %>%
                              format(nsmall =1) %>%
                              gsub(" ","",.) %>%
                              gsub("\\.",":",.))) %>%
  #mutate(CountAboveFilter = rowSums(.[,-c(1:3)] > 0.2) != 0) %>% # apply 0.2% filtering
  mutate(CountAboveFilter = rowSums(.[,-c(1:3)] > 5) != 0) %>% # apply 5% filtering
  filter(CountAboveFilter == TRUE) %>%
  select(-CountAboveFilter) %>% 
  pivot_longer(cols = !Line:SideChain, names_to = "AcylChain", values_to="PercentageAbundance")

# __________Save dataframe S06Abundance1 to generate supplementary table S9 in lipid MS1 to present raw data used to construct heatmap________

Heatmap_rawdf <- S06Abundance1 %>%
  mutate(SideChain = paste(SideChain,"#", sep="")) #%>% 
# write_excel_csv("MS1Standards/Heatmap_rawdf.csv")

p <- S06Abundance1 %>%
  dplyr::filter(AcylChain %in% c("1 Day_DG","19 Day_DG","1 Day_TG","19 Day_TG","1 Day_CL","19 Day_CL","1 Day_PC", "19 Day_PC",
                      "1 Day_PE", "19 Day_PE","1 Day_PG", "19 Day_PG","1 Day_PI", "19 Day_PI","1 Day_PS", "19 Day_PS",
                      "1 Day_LPC", "19 Day_LPC")) %>% 
  #dplyr::mutate(PercentageAbundance = log10(PercentageAbundance)) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("1 Day_DG","19 Day_DG","1 Day_TG","19 Day_TG","1 Day_CL","19 Day_CL","1 Day_PC", "19 Day_PC",
                                                                    "1 Day_PE", "19 Day_PE","1 Day_PG", "19 Day_PG","1 Day_PI", "19 Day_PI","1 Day_PS", "19 Day_PS",
                                                                    "1 Day_LPC", "19 Day_LPC","1 Day_TG e", "19 Day_TG e", "1 Day_PC e", "19 Day_PC e","1 Day_PE e",
                                                                    "19 Day_PE e","1 Day_PE p", "19 Day_PE p"))) %>% 
  #ggplot(aes(x = AcylChain, y = SideChain, fill = log10(PercentageAbundance))) + #log transformed percentageabundance
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust=0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ # JO likes this white background
  #ggtitle("Heatmap showing most common acyl chains in lipid subclasses based on <5% inclusion threshold" ) +
  labs(x= "Lipid subclasses", y = "Acyl chains", fill = "Percentage abundance") +
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (8))) + 
  theme(axis.title = element_text(size = 8))+
  theme(legend.position="none") #remove legend title

# save plots below to create Panel A of Figure 4: Heatmap showing most common acyl chains in lipid subclasses based on <5% inclusion threshold 
# in Day1 and Day19 B. tryoni males

ggsave(plot = p, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/Common acyl chains in lipid subclasses.jpg")                                                              

# ggsave(plot = p, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/log_transfomed.jpg")

# Note we have used the untransformed percentage abundance scale as heat map (figure 4) in lipid MS1. 

#__________________________Create Panel B of Figure 4___________________________________

# 1st generate data for there lipids

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

df.old <- read_csv("MS1Standards/Table_ConcAllStandards.csv") %>% # this is data from 1st run
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
    #str_remove_all("[ep]") %>%           
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
        #print(c(i,j))
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

#___________________________________________________________________________________________________

# Calculate the proportion of ester and ether linked SC in the 4 EL subclasses

Complete_SCTable1 <- Complete_SCTable %>% 
  dplyr::mutate(Type = ifelse(grepl("p",SideChain), "vinylether", 
                              ifelse(grepl("e", SideChain), "ether", "ester"))) %>% 
  dplyr::group_by(Samples, Age, Line, Time, SubClass, Type) %>% 
  dplyr::mutate(Sum = sum(Proportion)) %>% 
  dplyr::mutate(Proportion2 = ifelse(Sum == 0, 0, Proportion*100/Sum)) #if its 0 output o, iof not 0 than do the calculations 

# Calculate the mean proportions of common acyl chains present in each lipid subclass

Abundance_Subclass1 <- Complete_SCTable1 %>%  
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, SideChain, Proportion2, Type) %>% #the proportion here is mean computed earlier for a sample
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass", "SideChain", "Type"), measurevar = "Proportion2") %>% 
  dplyr::select(Line, Time, Age, SubClass, SideChain, Type, Proportion2, se) %>% # here proportion refers to mean percentage
  dplyr::mutate(Proportion2 = format(round(as.numeric(.$Proportion2), 4), nsmall = 4)) %>% 
  dplyr::select(Line, Time, Age, SubClass, Type, SideChain, Proportion2) 

# check to see if each Line by Subclass by Age comes to 100%
sum(as.numeric((Abundance_Subclass1 %>% dplyr::filter(Line == "s06" & SubClass == "TG e" & Age == "1 Day" & Type == "ester"))$Proportion2)) 


Prefilter_AbundanceSubclass <- Abundance_Subclass1 %>% # we have called this prefilter because later we might apply 5% or 10% on column Proprotions
  dplyr::mutate(Proportion2 = as.numeric(Proportion2)) %>% 
  dplyr::mutate(CombinedSubClasses = ifelse(grepl("e|p",SideChain), ### Check if e or p present in SideChain
                                            paste(Age, SubClass, Type, sep = "_"), # If true, concatenate Age_SubClass_Type
                                            paste(Age, SubClass, sep = "_"))) %>% # else, Age_SubClass
  # dplyr::mutate(Combinedsubclasses = ifelse(length(unique(Type))==1, 
  #                                           paste(Age, SubClass, sep = "_"),
  #                                           paste(Age, SubClass, Type, sep = "_"))) %>% # we looking at variable type and when we grouby sc, line, time, age, than within unique if there's only one type of linkage eg ester only 
  dplyr::select(-Age, -SubClass, -Type) %>% 
  mutate(SideChain = gsub("[ep]","", SideChain)) %>% 
  pivot_wider(names_from = c(CombinedSubClasses), values_from = Proportion2) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%  #replace NAs with 0 
  dplyr::filter(Line == "s06") %>% 
  write_excel_csv("MS1Standards/S06commonacylchainsEL.csv") # save data to generate Panel B of heatmap 

#______________________________
# plot Panle B of Figure 4 heat map

rm(list=ls())

Abundance_Subclass <- as.data.frame(read_csv("MS1Standards/S06commonacylchainsEL.csv"))

S06Abundance <- Abundance_Subclass %>% ungroup() %>% select(-Time) 

S06Abundance1 <- S06Abundance %>%
  arrange(desc(SideChain)) %>%  # arrange in descending order
  mutate(SideChain = as.factor(SideChain)) %>%
  mutate(SideChain = factor(SideChain, levels = unique(SideChain))) %>%
  mutate(CountAboveFilter = rowSums(.[,-c(1:5)] > 5) != 0) %>% # apply 5% filtering
  filter(CountAboveFilter == TRUE) %>%
  select(-CountAboveFilter) %>% 
  pivot_longer(cols = !Line:SideChain, names_to = "AcylChain", values_to="PercentageAbundance") #so combined Sidechain as been put into Acylchain

# Save raw data used to generate Panel B of heat map
# Heatmap_rawdf <- S06Abundance1 %>%
#  mutate(SideChain = paste(SideChain,"#", sep="")) %>% 
#  pivot_wider(names_from = c(SideChain), values_from = PercentageAbundance) %>% 
#  write_excel_csv("MS1Standards/Heatmap_rawdfEL.csv")

#S06Abundance1 %>% select(AcylChain) %>% unique()

p1 <- S06Abundance1 %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG e", "19 Day_TG e", "1 Day_TG e_ether", "19 Day_TG e_ether",
                                 "1 Day_PC e", "19 Day_PC e", "1 Day_PC e_ether", "19 Day_PC e_ether", 
                                 "1 Day_PE e", "19 Day_PE e", "1 Day_PE e_ether", "19 Day_PE e_ether",
                                 "1 Day_PE p", "19 Day_PE p", "1 Day_PE p_vinylether", "19 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("1 Day_TG e", "19 Day_TG e", "1 Day_TG e_ether", "19 Day_TG e_ether",
                                                                    "1 Day_PC e", "19 Day_PC e", "1 Day_PC e_ether", "19 Day_PC e_ether", 
                                                                    "1 Day_PE e", "19 Day_PE e", "1 Day_PE e_ether", "19 Day_PE e_ether",
                                                                    "1 Day_PE p", "19 Day_PE p", "1 Day_PE p_vinylether", "19 Day_PE p_vinylether"))) %>% 
  #ggplot(aes(x = AcylChain, y = SideChain, fill = log10(PercentageAbundance))) + #log transformed percentage abundance
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  #ggtitle("Heatmap showing most common acyl chains in lipid subclasses based on <5% inclusion threshold") +
  labs(x= "Lipid subclasses", y = "Acyl/alkyl/alkenyl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (8))) + theme(axis.title = element_text(size = 8))

ggsave(plot = p1, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "MS1Standards/Common acyl, alkyl, alkenyl chains in lipid subclasses.jpg")                                                              

#ggsave(plot = p, width = 6.8, height = 3.2, units = "in", dpi = 300,filename = "log_transfomed.jpg")

#__________________________________END_______________________________________________
