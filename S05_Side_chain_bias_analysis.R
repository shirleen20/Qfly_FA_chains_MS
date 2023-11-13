#___________________________________________________02/062023_____________________________________________

# This consolidated script separates the ether and ester bond in ether lipids for side chain analysis 
# to generate Tables 3, S13-S15
# we divide them into short chains that  were less than 16C, medium chains that were 16-18C long and long 
# containing >18C. 
# double bond proportions in each subclass by dividing into three categories of saturated fatty acids s0 chains 
# with 0 DB, monounsaturated with 1 DB and polyunsaturated more than one DB

# Inclusion criteria in bold italic:
#  1)	Excess:
#  a.	Ratio of observed to expected is greater than 1.5 in at least one day OR magnitude difference between observed and expected is greater than 10 (%)
#  b.	At least one day must be statistically significant
#  c.	Observed percentage > 5(%) in both days
#  d.	Observed greater than expected in both days

# 2)	Deficit:
#  a.	Ratio of observed to expected is less than 0.67 in at least one day OR magnitude difference between observed and expected is greater than 10 (%)
#  b.	At least one day must be statistically significant
#  c.	Expected percentage > 5 (%) in both days
#  d.	Expected greater than observed in both days

# NB: We use outputs from here to generate figure 2 and 3 

library("tidyverse")
library("readxl")
library("Rmisc")
library("combinat") # combn and permn
library("gtools")   

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

#__________________________ Create data frame _________________________________
Variables <- read_excel("data/Variables.xlsx") %>% 
  dplyr::select(ID, Weight, Time2) %>% 
  dplyr::rename(samples = ID) 

DF.1 <- DF1 %>% 
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
  dplyr::mutate(Carbon = as.double(Carbon))

# _________________________________________________________

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

DF.main <- DF.1 %>%
  ungroup() %>% 
  mutate(Class2 = if_else(grepl(" e", SubClass), Class, SubClass)) %>% ## All ether lipids use ester lipid standards
  left_join(Regression_coef, by = "Class2") %>%
  mutate(Conc = 10^((log10(N.Area)-Intercept)/Slope)) %>%
  mutate(LogAreaCheck = Intercept + Slope*(log10(Conc)), AreaCheck = 10^LogAreaCheck) %>%  ### Just to double check
  dplyr::select(-N.Area, -Class2, -Intercept, -Slope, -LogAreaCheck,-AreaCheck) %>% 
  dplyr::relocate(Conc,.after = Cage)

# start analysis

#First remove subclasses PSe and DGe from DF because these two subclasses have all lipid species that do not have acyl chain identified
DF <- DF.main %>% 
  dplyr::filter(SubClass != "PS e" & SubClass != "DG e")

#grep - finds things, gives index 1, 3, 5 not showing 1 is apple, 3 is pear unless you specify value =T) check help ?g sub or grep
#g sub is a substitution

######## NEW code for F1 is not compatible with old codes using F3 and F4 functions (if F3/F4 functions to be resurrected, must rearrange output of F1 to place Merged and SC#_class before SC#)
### Function to create for each subclass data frame, all side chains separated into columns
# X = data frame to process
# Y = Subclass level as string
# type takes either "SCBond", "SC" or "Bond"... SCBond outputs combined SideChain and double bond info, SC is Side Chain only, Bond is dounble bond only
# keep "All", "Ethervinyl", "Ester"
F1 <- function(X, Y, type="SCBond", keep="All"){
  Subclass <- X %>% # Subclass is the dataframe
    dplyr::filter(SubClass == Y)
  length1 <- length(Subclass)  # The starting number of columns
  Names <- colnames(Subclass) # Save column names
  
  Subclass_SC <- gsub(".*\\((.+)\\).*","\\1",Subclass$Name) %>%  ### keep only string within "( )"
    #str_remove_all("[ep]") %>% 
    
    strsplit("_") ### split by "_"
  SCNumber <- Subclass_SC %>% map(length) %>% unlist ### find the number of side chains in each lipid for each individual (want to remove those with no side chains identified)
  
  RemoveSpecies <- which(SCNumber < max(SCNumber)) ### vector of indices of those without side chains identified. They will have length less than the maximum number of side chains
  Subclass_SC <- do.call(rbind, Subclass_SC) %>% as.tibble()### THis will have columns of side chains (e.g. if TGe, there will be 3 columns)
  #print(Subclass_SC)
  if (keep =="All"){
    Subclass_SC <- sapply(Subclass_SC, function(x){str_remove_all(x,"[ep]")}) %>% ### remove all e's and p's from every column of Subclass_SC
      as.tibble()### change from a matrix to a dataframe
  } else if (keep == "Ethervinyl"){
    EtherVinyl <- sapply(Subclass_SC, function(x) { sum(grepl("[ep]",x)) == length(x)}) ### which column has all e or p
    Subclass_SC <- sapply(Subclass_SC[,EtherVinyl],function(x){str_remove_all(x,"[ep]")}) %>% as.tibble()# keep column with e or p
    
  } else {
    Ester <- sapply(Subclass_SC, function(x) { !sum(grepl("[ep]",x)) == length(x)}) ### which column(s) does(do) not have e or p
    Subclass_SC <- sapply(Subclass_SC[,Ester],function(x){str_remove_all(x,"[ep]")}) %>% as.tibble()# keep column(s) with ester bonds (no e or p)
  }
  #print("After")
  #print(Subclass_SC)
  Subclass <- cbind(Subclass, Subclass_SC)%>% as.tibble() ### Join the new columns to the old data frame
  length2 <- length(Subclass) ### New number of columns
  #print(Subclass)
  colnames(Subclass) <- c(Names, paste("SC", 1:(length2-length1), sep = ""))  ### Rename new colunms with SC1, SC2 ... etc length2-length1 is the total number of new columns
  for (i in 1:(length2-length1)){ # length1 is columns from Name to Time, length 2 is after the new columns created based on # of SC present
    SClength <- str_split(pull(Subclass[,length1+i]),":") %>% map(1) %>% unlist() %>% as.numeric() ### Find the SC length for the ith side chain
    SCsaturation <- str_split(pull(Subclass[,length1+i]),":") %>% map(2) %>% unlist() %>% as.numeric()
    SClength2 <- ifelse(SClength < 16, "S", ifelse(SClength < 19, "M", "L"))
    SCsat2 <- ifelse(SCsaturation == 0, "0", ifelse(SCsaturation == 1, "1", "X"))
    if (type == "SCBond"){
      Subclass[paste0("SC",i,"_class")] <- paste0(SClength2, SCsat2)
    }
    else if(type == "SC"){
      Subclass[paste0("SC",i,"_class")] <- SClength2
    }
    else if(type == "Bond"){
      Subclass[paste0("SC",i,"_class")] <- SCsat2
    }
    else {
      stop("Incorrect 'type' specified, should either be 'SCBond', 'SC' or 'Bond'")
    }
  }
  Subclass <- Subclass %>% unite("Merged", "SC1_class":paste0("SC",length2-length1,"_class"), remove = FALSE)
  if (length(RemoveSpecies) == 0){
    return(Subclass)
  }
  else{
    return(Subclass[-RemoveSpecies,]) ### Output new data frame after removing rows with unindentified side chains
  }
}

F2 <- function(X, Subclass = "SubClass", type = "SCBond"){
  Subclass_names <- unique(X[[Subclass]]) ### All the unique Subclasses
  Out_list <- setNames(replicate(length(Subclass_names), data.frame()), Subclass_names) ### create a list of empty data frames with names Subclass_names with length of the number of subclasses 
  # if using vector(mode = "list", length = length(Subclass_names)) ... required renaming objects in Out_list, see next line
  # names(Out_list) <- Subclass_names ## rename each object within the list with the subclass names
  for (i in Subclass_names){ ### For each separate subclass
    Out_list[[i]] <- F1(X, i, type) ### in each object of that ith subclass, run F1 function on data frame X with subclass i
  }
  return(Out_list)
}

LIST1 <- F2(X = DF, type = "SCBond")
LIST2 <- F2(X = DF, type = "SC")
LIST3 <- F2(X = DF, type = "Bond")

DF_EtherVinyl <- DF %>% filter(SubClass %in% c( "PC e", "PE e", "PE p","TG e")) ### keep only these subclasses to look at ether, vinyl or ester side chain abundances
LIST1_Ethervinyl <- F2(X = DF_EtherVinyl, type = "SCBond")
LIST2_Ethervinyl <- F2(X = DF_EtherVinyl, type = "SC")
LIST3_Ethervinyl <- F2(X = DF_EtherVinyl, type = "Bond")
 
LIST1_Ester <- F2(X = DF_EtherVinyl, type = "SCBond")
LIST2_Ester <- F2(X = DF_EtherVinyl, type = "SC")
LIST3_Ester <- F2(X = DF_EtherVinyl, type = "Bond")

#View(LIST1$TG %>% select(Samples, SC1_class, SubClass, N.Area) %>%  complete(SC1_class, SubClass, Samples))

##### Summarise the Acylclass (short, medium or long, with level of un-saturation info coded as 0, 1 or more than 1 'X')
### Short = <16, Medium = 16-18, Long = >18
### 0 no double bond
SummariseLipid <- function(DF){
  NumberAcylChain <- sum(grepl("_class", colnames(DF))) ### Number of side chains
  for (i in 1:NumberAcylChain){
    New_dat <- DF %>% dplyr::select(Samples, Line, Age, Time, paste0("SC",i,"_class"), SubClass, Conc)  %>% 
      dplyr::rename(AcylClass=paste0("SC",i,"_class")) %>% ### rename SC#i_class to AcylClass
      dplyr::group_by(Samples, Line, Age, Time, SubClass, AcylClass) %>% 
      dplyr::summarise(SUMi=sum(Conc)) %>%  ungroup() %>%  ### SUMi is the sum Conc for each AcylClass by Sample/SubClass
      complete(nesting(Samples, Line, Age, Time), AcylClass, SubClass, fill = list(SUMi = 0)) ### fill in zeroes for samples with missing AcylClass
    names(New_dat) <- ifelse(names(New_dat) == "SUMi", paste0("SUM",i), names(New_dat) ) ### rename SUMi with unique ith value... i.e. SUM1, SUM2, etc
    if (i == 1){
      Final_dat <- New_dat
    }
    else {
      Final_dat <- full_join(Final_dat, New_dat, by = c("Samples", "Line", "Age", "Time", "AcylClass", "SubClass"))
    }
  }
  Final_dat <- Final_dat %>% replace(is.na(.), 0) %>% dplyr::mutate(SUM = rowSums(.[,(length(.)-(NumberAcylChain-1)):length(.)]) ) %>%  ## SUM across SUM1, SUM2, SUM3 etc...
    dplyr::group_by(Samples, Line, Age, Time, SubClass) %>%  dplyr::mutate(TotalSum = sum(SUM)) %>% ### TotalSum is the sum across acylclass
    dplyr::mutate(Proportion=SUM/TotalSum) %>% ungroup() ### Proportion of each acylclass for each sample
  return(list(Final_dat=Final_dat,
              summary=Final_dat %>% 
                dplyr::mutate(LineTimeAge = paste(Line, Time, Age, sep = "_")) %>% 
                dplyr::group_by(Line, Age, Time, LineTimeAge, SubClass, AcylClass) %>% 
                dplyr::summarise(MEAN=mean(Proportion), SE=sd(Proportion)/sqrt(n()), n=n()) %>% ungroup()))
}

### See example
# SummariseLipid(LIST2$CL %>% filter(Line == "s06"))

SummaryProportion <- LIST1 %>% map(~SummariseLipid(.)) ###### List containing data frames for each lipid subclass the proportion of acyl chain length and level of un-saturation
### SummaryProportion$DG$summary %>% select(LineTimeAge, n) %>% unique()
SummaryProportion2 <- LIST2 %>% map(~SummariseLipid(.))
SummaryProportion3 <- LIST3 %>% map(~SummariseLipid(.))

SummaryProportionEthervinyl <- LIST1_Ethervinyl %>% map(~SummariseLipid(.))
SummaryProportion2Ethervinyl <- LIST2_Ethervinyl %>% map(~SummariseLipid(.))
SummaryProportion3Ethervinyl <- LIST3_Ethervinyl %>% map(~SummariseLipid(.))

SummaryProportionEster <- LIST1_Ester  %>% map(~SummariseLipid(.))
SummaryProportion2Ester <- LIST2_Ester %>% map(~SummariseLipid(.))
SummaryProportion3Ester <- LIST3_Ester %>% map(~SummariseLipid(.))

### Will loop through this when performing the analysis
AcylChain_type <- c("S0","S1","SX","M0","M1","MX","L0","L1","LX")
AcylChain_type2 <- c("S","M","L")
AcylChain_type3 <- c("0","1","X")

AcylChain_typeEtherVinyl <- c("S0","S1","SX","M0","M1","MX","L0","L1","LX")
AcylChain_type2EtherVinyl <- c("S","M","L")
AcylChain_type3EtherVinyl <- c("0","1","X")

AcylChain_typeEster <- c("S0","S1","SX","M0","M1","MX","L0","L1","LX")
AcylChain_type2Ester <- c("S","M","L")
AcylChain_type3Ester <- c("0","1","X")

#______________ Create list containing all combination of strings when sampling n number of times_____________
# AcylChain_type --> Classification of acyl chains (Small, medium, long... 0, 1 or more than 1 double bond (X)) - or any combination of strings (as vector)
# NSc Number of side chains (or number of items)
CombinePerm <- function(AcylChain_type, NSc){
  AcylChain_typeOrder <- paste0(1:length(AcylChain_type),"-", AcylChain_type)
  ListContainingPerm <- sub(".*\\-(.+)", "\\1" , combinations(length(AcylChain_typeOrder),NSc,AcylChain_typeOrder, repeats.allowed = TRUE)) %>% # create all possible combintation when choosing 'NSc' times
    split(row(.)) %>% # create a list out of the matrix/array
    map(~permn(.)) %>% ### Within each unique combination, create all possible permutations
    map(~as.list(.) %>% map(~paste0(.,collapse="_"))) %>% ### Permutations are presented as lists (within the combination list) - paste together strings within permutations
    map(~unlist(.)) %>% ### unlist to create a vector of permutations (within list of combinations)
    map(~unique(.))
  names(ListContainingPerm) <- ListContainingPerm %>% map(~.[[1]]) %>%  unlist()### rename list of combinations with first element of each vector within the list
  return(ListContainingPerm)
}

#____________________________JOINED EXPECTATION (no conditional grouping)_______________________________
###### Wrap together as a function to work across all SubClass (using for loop) in LIST1
##### Also looping through 'AcylChain_type'
##### Output is a list of Subclasses each containing list of AcylChain_type each containing a data frame containing the observed and expected frequency of other acyl chain combinations 
##### REQUIREMENTS
### DFlist is the List of Data frames (one for each subclass) containing the Short/Medium/Long; 0/1/X classification and Conc values
### SumProp is the output of "SummariseLipid" containing summary of the observed proportions of individual acylchains
### AcylChain_type
### FILTER (Line filter)

ObsExpProportion_XGroup <- function(DFlist, SumProp, AcylChain_type, FILTER = ""){
  FinalList <- list()
  Leftover_SC <- list(CombinePerm(AcylChain_type, 1), CombinePerm(AcylChain_type, 2),CombinePerm(AcylChain_type, 3),CombinePerm(AcylChain_type, 4)) ### List of combination of acyl chain types
  for (i in names(DFlist)){
    #print(i)
    SCNum <- sum(grepl("_class", colnames(DFlist[[i]]))) ### Number of side chains
    # if (SCNum == 1){
    #   FinalList[[i]] <- list()
    #   next
    # }
    ############################ Computing the Expected proportion of combination of acylclass
    RelevantSummary <- SumProp[[i]][["summary"]] %>% dplyr::select(LineTimeAge, AcylClass, MEAN) ### Obtain the summary of individual acyl chain type observed mean frequencies, for ith subclass
    Lgroups <- length(unique(RelevantSummary$LineTimeAge)) ### Number of combination of variable groups (combinations of Line TIme and age) in the analysis
    IndividualExpected <- do.call(rbind, names(Leftover_SC[[SCNum]]) %>% strsplit("_")) %>%  as.data.frame()
    names(IndividualExpected) <- paste0("Class", 1:(SCNum))
    IndividualExpected <- do.call("rbind", replicate(Lgroups, IndividualExpected, simplify = FALSE))
    
    #print(c(i,Lgroups,length(names(Leftover_SC[[SCNum-1]]))))
    Expected <- tibble(LineTimeAge = rep(unique(RelevantSummary$LineTimeAge), each = length(names(Leftover_SC[[SCNum]]))),
                       Category=rep(names(Leftover_SC[[SCNum]]), Lgroups ) )
    Expected <- cbind(Expected, IndividualExpected)
    #print(RelevantSummary)
    #print(Expected)
    for (k in 1:(SCNum)){
      names(Expected)[length(Expected)-SCNum+1] <- "AcylClass"
      Expected <- left_join(Expected, RelevantSummary, by = c("LineTimeAge","AcylClass")) #%>% filter(is.na(MEAN) == FALSE) #############################
      names(Expected)[length(Expected)-SCNum] <- paste0("AcylClass",k)
      names(Expected)[length(Expected)] <- paste0("MEAN",k)
      #print(Expected)
      #print("SPACE")
    }
    Expected <- Expected %>% rowwise() %>% 
      dplyr::mutate(Denom = prod(factorial(table(c_across("AcylClass1":paste0("AcylClass",SCNum)))))) %>% 
      dplyr::mutate(EXPECTED=prod(c_across("MEAN1":paste0("MEAN",SCNum)))*factorial(SCNum)/Denom) %>% ungroup() %>%  ### row product (requires rowwise to vectorise prod)
      dplyr::select(LineTimeAge, Category, EXPECTED)
    ######################### Computing the observed proportion of combinations of acylclass
    AcylChain_list <- list()
    #for (j in AcylChain_type){
    #print(c(i,j))
    NEW <- DFlist[[i]] 
    #%>% filter(grepl(j, Merged)) %>% 
    #  dplyr::mutate(Leftover = sub(paste0(j,"_|_",j),"",Merged), AcylGroup=j) ### Replace only first occurrence of "Medium1_" or "_Medium1"
    NewCol <- rep(0, length(NEW$Merged))
    for (ii in names(Leftover_SC[[SCNum]])){
      NewCol <- ifelse(NEW$Merged %in% Leftover_SC[[SCNum]][[ii]], ii, NewCol)
    }
    NEW$Category <- NewCol
    #if (j %in% unique(SumProp[[i]][["summary"]]$AcylClass)){ ### If the acyl chain type is present in subclass i, compute this else nothing 
    CurrentSummary <- NEW %>% dplyr::select(Samples, Line, Age, Time, Category, SubClass, Conc) %>% 
      group_by(Samples, Line, Age, Time, SubClass, Category) %>% dplyr::summarise(SUM=sum(Conc)) %>% ungroup() %>% 
      complete(nesting(Samples, Line, Age, Time), Category, SubClass, fill = list(SUM = 0)) %>% 
      group_by(Samples, Line, Age, Time, SubClass) %>% dplyr::mutate(TotalSUM = sum(SUM)) %>% dplyr::mutate(Proportion=SUM/TotalSUM) %>% ungroup()
    FinalCurrSummary <- CurrentSummary %>% group_by(Line, Age, Time, SubClass, Category) %>% dplyr::summarise(ObservedMean=mean(Proportion), SE = sd(Proportion)/sqrt(length(Proportion)), n=n()) %>% 
      ungroup() %>% dplyr::mutate(LineTimeAge = paste(Line, Time, Age, sep = "_"))
    # print(c(CurrentSummary$Category))
    # print(c(FinalCurrSummary$Category))
    AcylChain_list <- full_join(FinalCurrSummary, Expected, by = c("LineTimeAge", "Category")) %>% 
      dplyr::mutate(Line = str_split(LineTimeAge, pattern = "_") %>% map(1) %>% unlist()) %>% 
      dplyr::mutate(Time = str_split(LineTimeAge, pattern = "_") %>% map(2) %>% unlist()) %>% 
      dplyr::mutate(Age = str_split(LineTimeAge, pattern = "_") %>% map(3) %>% unlist()) %>% 
      dplyr::mutate(SubClass = unique(SubClass)[1]) %>% 
      relocate(LineTimeAge, .after = SubClass) %>% 
      relocate(EXPECTED, .after = Category) %>% 
      replace_na(list(ObservedMean = 0, SE = 0, n = 0 )) %>% 
      dplyr::mutate(Zval = abs((ObservedMean-EXPECTED)/SE)) %>% 
      dplyr::mutate(pval = 1-pnorm(Zval)) %>% 
      group_by(Line, Time, Age) %>% ################## 
    dplyr::mutate(padj = p.adjust(pval, method = "bonferroni")) %>% 
      ungroup() %>% ###########################
    #dplyr::mutate(Direction = if_else(padj <= 0.05, if_else(ObservedMean > EXPECTED, "Greater", "Lesser"), "NS")) %>% 
    dplyr::mutate(Ratio = if_else(ObservedMean == 0 & EXPECTED == 0 , 1, (ObservedMean/EXPECTED)) ) %>% ########################
    dplyr::mutate(RatioSE = if_else(SE == 0 & EXPECTED == 0, 0, SE/EXPECTED)) %>%  ###################################
    #dplyr::mutate(Ratio = ObservedMean/EXPECTED ) %>% ########################
    #dplyr::mutate(RatioSE = SE/EXPECTED) %>%  ###################################
    dplyr::mutate(Significance = if_else(padj < 0.001, "***",
                                         if_else(padj < 0.01, "**",
                                                 if_else(padj < 0.05, "*",
                                                         if_else(padj <= 0.1, "~",""))))) ##################
    if (FILTER != ""){
      AcylChain_list <- AcylChain_list %>% dplyr::filter(Line == FILTER)
    } #######################################
    # }
    # else {
    #   AcylChain_list[[j]] <- ""
    # }
    #}
    FinalList[[i]] <- AcylChain_list
  }
  return(FinalList)
}

PROPORTION_LIST1_XGroup <- ObsExpProportion_XGroup(LIST1, SummaryProportion, AcylChain_type, FILTER = "s06")
PROPORTION_LIST2_XGroup <- ObsExpProportion_XGroup(LIST2, SummaryProportion2, AcylChain_type2, FILTER = "s06")
PROPORTION_LIST3_XGroup <- ObsExpProportion_XGroup(LIST3, SummaryProportion3, AcylChain_type3, FILTER = "s06") 
### Make sure the Expected values add up to one
PROPORTION_LIST2_XGroup$`CL` %>% group_by(LineTimeAge) %>% dplyr::select(LineTimeAge, Category, EXPECTED) %>% dplyr::mutate(SUM = sum(EXPECTED))

PROPORTION_LIST1_Ethervinyl_XGroup <- ObsExpProportion_XGroup(LIST1_Ethervinyl, SummaryProportionEthervinyl, AcylChain_typeEtherVinyl, FILTER = "s06")
#PROPORTION_LIST1_Ethervinyl_XGroup$`TG e` %>% dplyr::group_by(LineTimeAge) %>% dplyr::select(LineTimeAge, Category, EXPECTED) %>% dplyr::mutate(SUM = sum(EXPECTED))
PROPORTION_LIST2_Ethervinyl_XGroup <- ObsExpProportion_XGroup(LIST2_Ethervinyl, SummaryProportion2Ethervinyl, AcylChain_type2EtherVinyl, FILTER = "s06")
PROPORTION_LIST3_Ethervinyl_XGroup <- ObsExpProportion_XGroup(LIST3_Ethervinyl, SummaryProportion3Ethervinyl, AcylChain_type3EtherVinyl, FILTER = "s06") 

PROPORTION_LIST1_Ester_XGroup <- ObsExpProportion_XGroup(LIST1_Ester, SummaryProportionEster, AcylChain_typeEster, FILTER = "s06")
PROPORTION_LIST2_Ester_XGroup <- ObsExpProportion_XGroup(LIST2_Ester, SummaryProportion2Ester, AcylChain_type2Ester, FILTER = "s06")
#PROPORTION_LIST2_Ester_XGroup$`TG e` %>% dplyr::group_by(LineTimeAge) %>% dplyr::select(LineTimeAge, Category, EXPECTED) %>% dplyr::mutate(SUM = sum(EXPECTED))
PROPORTION_LIST3_Ester_XGroup <- ObsExpProportion_XGroup(LIST3_Ester, SummaryProportion3Ester, AcylChain_type3Ester, FILTER = "s06") 

TablingFunction <- function(LIST, AcylChaintype="AcylBond"){
  AcylChain_type <- switch(AcylChaintype,
                           AcylBond = c("S0"="1","S1"="2","SX"="3","M0"="4","M1"="5","MX"="6","L0"="7","L1"="8","LX"="9"),
                           AcylLength = c("S" = "1", "M" = "2", "L" = "3"),
                           Bond = c("0" = "1", "1" = "2", "X" = "3"))
  
  FINAL <- LIST %>% 
    bind_rows() %>%   
    dplyr::mutate(Class = case_when(SubClass %in% c('DG', 'TG') ~ "NL",
                                    SubClass %in% c("CL", "PE", "PC","PG","PI","PS","LPC") ~ "PL",
                                    SubClass %in% c("DG e","TG e") ~ "ELNL",
                                    SubClass %in% c("PE e","PE p","PC e","PS e") ~ "ELPL")) %>%
    relocate(Class, .before = SubClass) %>% 
    dplyr::mutate(Cat2 = str_replace_all(Category, AcylChain_type)) %>% 
    dplyr::arrange(factor(Class, levels = c("NL", "PL", "ELNL","ELPL")), SubClass, Cat2)   
  
  FINAL <- LIST %>% 
    bind_rows() %>%   
    dplyr::mutate(Class = case_when(SubClass %in% c('DG', 'TG') ~ "NL",
                                    SubClass %in% c("CL", "PE", "PC","PG","PI","PS","LPC") ~ "PL",
                                    SubClass %in% c("DG e","TG e") ~ "ELNL",
                                    SubClass %in% c("PE e","PE p","PC e","PS e") ~ "ELPL")) %>% 
    relocate(Class, .before = SubClass) %>% 
    dplyr::mutate(Cat2 = str_replace_all(Category, AcylChain_type)) %>% 
    dplyr::arrange(factor(Class, levels = c("NL", "PL", "ELNL","ELPL")), SubClass, Cat2)   
  
  
  FINAL <- FINAL %>% dplyr::select(Age, Class, SubClass, Category, Ratio, RatioSE, Significance, ObservedMean, EXPECTED) %>% 
    dplyr::mutate(ObservedMean = format(round(as.numeric(.$ObservedMean)*100, 5), nsmall = 5) %>% str_replace_all(" ","")) %>% ### remove all spaces after formatting to % and 1 decimal place
    dplyr::mutate(EXPECTED = format(round(as.numeric(.$EXPECTED)*100, 5), nsmall = 5) %>% str_replace_all(" ","")) %>%
    dplyr::mutate(Ratio = format(round(as.numeric(.$Ratio), 5), nsmall = 5)) %>% 
    dplyr::mutate(RatioSE = format(round(as.numeric(.$RatioSE), 5), nsmall = 5)) %>% 
    dplyr::mutate(Significance = replace_na(Significance, "")) %>% 
    unite("Ratio", Ratio:RatioSE, sep = " ±") %>%
    unite("Ratio", `Ratio`:Significance, sep = "") %>%
    dplyr::mutate(`Ratio` = if_else(grepl("NaN", Ratio), NA_character_, Ratio)) %>% 
    dplyr::mutate(`Ratio` = if_else(grepl("NA", Ratio), NA_character_, Ratio)) %>%
    pivot_wider(names_from = Age, values_from = c(`Ratio`, ObservedMean, EXPECTED) ) %>% 
    dplyr::filter(!is.na(`Ratio_1 Day`)|!is.na(`Ratio_19 Day`)) %>%
    relocate(`ObservedMean_1 Day`, .after = `Ratio_1 Day`) %>% 
    relocate(`EXPECTED_1 Day`, .after = `ObservedMean_1 Day`)
  return(FINAL)
}

AcylBond_Proportion <- TablingFunction(PROPORTION_LIST1_XGroup, "AcylBond") %>% #mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>% 
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19) 


AcylLength_Proportion <- TablingFunction(PROPORTION_LIST2_XGroup, "AcylLength") %>% dplyr::mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>%
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19) 


Bond_Proportion <- TablingFunction(PROPORTION_LIST3_XGroup, "Bond") %>% #dplyr::mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>%
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19) 

#####___repeat above for ELs separated into ether and ester chains______
AcylBond_ProportionEthervinyl <- TablingFunction(PROPORTION_LIST1_Ethervinyl_XGroup, "AcylBond") %>% #mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>% 
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19) 


AcylLength_ProportionEthervinyl <- TablingFunction(PROPORTION_LIST2_Ethervinyl_XGroup, "AcylLength") %>% dplyr::mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>%
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19) 


Bond_ProportionEthervinyl <- TablingFunction(PROPORTION_LIST3_Ethervinyl_XGroup, "Bond") %>% #dplyr::mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>%
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19) 

AcylBond_ProportionEster <- TablingFunction(PROPORTION_LIST1_Ester_XGroup, "AcylBond") %>% #mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>% 
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19) 


AcylLength_ProportionEster <- TablingFunction(PROPORTION_LIST2_Ester_XGroup, "AcylLength") %>% dplyr::mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>%
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19) 


Bond_ProportionEster <- TablingFunction(PROPORTION_LIST3_Ester_XGroup, "Bond") %>% #dplyr::mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::filter(SubClass != "LPC") %>%
  dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) %>% 
  dplyr::mutate(Diff1Day = as.numeric(format(round(as.numeric(.$Diff1Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Diff_19Day = as.numeric(format(round(as.numeric(.$Diff_19Day), 5), nsmall = 5))) %>% 
  dplyr::mutate(Day1 = if_else(as.numeric(`ObservedMean_1 Day`) == 0 & as.numeric(`EXPECTED_1 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_1 Day`),
                               "-",
                               paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) ) %>% 
  dplyr::mutate(Day19 = if_else(as.numeric(`ObservedMean_19 Day`) == 0 & as.numeric(`EXPECTED_19 Day`) == 0 & grepl("1.00.*±.*0.00", `Ratio_19 Day`),
                                "-",
                                paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) ) %>%  
  dplyr::select(Class, SubClass, Category, Day1, Day19)

#######___________


TablingFunctionMain <- function(LIST, AcylChaintype = "AcylBond", Greater = 1.5, Lesser = "", Expected = 5, Observed = 5,
                                Difference = 10){
  AcylChain_type <- switch(AcylChaintype,
                           AcylBond = c("S0"="1","S1"="2","SX"="3","M0"="4","M1"="5","MX"="6","L0"="7","L1"="8","LX"="9"),
                           AcylLength = c("S" = "1", "M" = "2", "L" = "3"),
                           Bond = c("0" = "1", "1" = "2", "X" = "3"))
  FINAL <- LIST %>% 
    bind_rows() %>%   
    dplyr::mutate(Class = case_when(SubClass %in% c('DG', 'TG') ~ "NL",
                                    SubClass %in% c("CL", "PE", "PC","PG","PI","PS","LPC") ~ "PL",
                                    SubClass %in% c("DG e","TG e") ~ "ELNL",
                                    SubClass %in% c("PE e","PE p","PC e","PS e") ~ "ELPL")) %>% 
    relocate(Class, .before = SubClass) %>% 
    dplyr::mutate(Cat2 = str_replace_all(Category, AcylChain_type)) %>% 
    dplyr::arrange(factor(Class, levels = c("NL", "PL", "ELNL","ELPL")), SubClass, Cat2)
  FINAL <- FINAL %>% 
    dplyr::select(Age, Class, SubClass, Category, Ratio, RatioSE, Significance, ObservedMean, EXPECTED) %>% 
    dplyr::mutate(ObservedMean = format(round(as.numeric(.$ObservedMean)*100, 1), nsmall = 1) %>% str_replace_all(" ","")) %>% ### remove all spaces after formatting to % and 1 decimal place
    dplyr::mutate(EXPECTED = format(round(as.numeric(.$EXPECTED)*100, 1), nsmall = 1) %>% str_replace_all(" ","")) %>%
    dplyr::mutate(Ratio = format(round(as.numeric(.$Ratio), 1), nsmall = 1)) %>% 
    dplyr::mutate(RatioSE = format(round(as.numeric(.$RatioSE), 1), nsmall = 1)) %>% 
    dplyr::mutate(Significance = replace_na(Significance, "")) %>% 
    unite("Ratio", Ratio:RatioSE, sep = " ±") %>%
    unite("Ratio", `Ratio`:Significance, sep = "") %>%
    dplyr::mutate(`Ratio` = if_else(grepl("NaN", Ratio), NA_character_, Ratio)) %>% 
    dplyr::mutate(`Ratio` = if_else(grepl("NA", Ratio), NA_character_, Ratio)) %>%
    pivot_wider(names_from = Age, values_from = c(`Ratio`, ObservedMean, EXPECTED) ) %>% 
    dplyr::filter(!is.na(`Ratio_1 Day`)|!is.na(`Ratio_19 Day`)) %>%
    relocate(`ObservedMean_1 Day`, .after = `Ratio_1 Day`) %>% 
    relocate(`EXPECTED_1 Day`, .after = `ObservedMean_1 Day`) %>% 
    dplyr::mutate(Day1Ratio = str_split(`Ratio_1 Day`, "±") %>% map(1) %>% unlist() %>% as.numeric()) %>% 
    dplyr::mutate(Day19Ratio = str_split(`Ratio_19 Day`, "±") %>% map(1) %>% unlist() %>% as.numeric()) %>% 
    dplyr::filter(as.numeric(`EXPECTED_1 Day`) > Expected & as.numeric(`EXPECTED_19 Day`) > Expected) %>% 
    dplyr::filter(as.numeric(`ObservedMean_1 Day`) > Observed & as.numeric(`ObservedMean_19 Day`) > Observed) %>% 
    dplyr::mutate(Diff1Day = as.numeric(`ObservedMean_1 Day`)-as.numeric(`EXPECTED_1 Day`)) %>% 
    dplyr::mutate(Diff_19Day = as.numeric(`ObservedMean_19 Day`)-as.numeric(`EXPECTED_19 Day`)) 
  if (Greater != Inf) {
    FINAL <- FINAL %>% dplyr::filter(as.numeric(Day1Ratio) >= Greater |as.numeric(Day19Ratio) >= Greater|
                                       as.numeric(Day1Ratio) < Lesser |as.numeric(Day19Ratio) < Lesser |
                                       as.numeric(Diff1Day) > Difference | as.numeric(Diff_19Day) > Difference) %>% 
      dplyr::filter(grepl("\\*|\\~", `Ratio_1 Day`) | grepl("\\*|\\~", `Ratio_19 Day`))  %>% 
      dplyr::filter(as.numeric(`ObservedMean_1 Day`)> as.numeric(`EXPECTED_1 Day`) & 
                      as.numeric(`ObservedMean_19 Day`) > as.numeric(`EXPECTED_19 Day`)) %>% 
      dplyr::select(-Day1Ratio, -Day19Ratio)
  }
  else {
    FINAL <- FINAL %>% dplyr::filter(as.numeric(Day1Ratio) < Lesser |as.numeric(Day19Ratio) < Lesser|
                                       as.numeric(Diff1Day) < -Difference | as.numeric(Diff_19Day) < -Difference ) %>% 
      dplyr::filter(grepl("\\*|\\~", `Ratio_1 Day`) | grepl("\\*|\\~", `Ratio_19 Day`))  %>% 
      dplyr::filter(as.numeric(`ObservedMean_1 Day`) < as.numeric(`EXPECTED_1 Day`) & 
                      as.numeric(`ObservedMean_19 Day`) < as.numeric(`EXPECTED_19 Day`)) %>%
      dplyr::select(-Day1Ratio, -Day19Ratio)
  }
  return(FINAL)
}


### Ratio greater than 1.5, in either or both 1 day, 19 day... must be statistically significant
## Observed frequency must be greater than 3% in at least on of the days (no filter on Expected frequency)
AcylBond_ProportionMain <- TablingFunctionMain(PROPORTION_LIST1_XGroup, AcylChaintype = "AcylBond", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")

### Ratio Less than 0.67 in at least one of 1 day or 19 day, and must be statistically significant, 
### Expected proportion must be greater than 3% in at least one of the days
AcylBond_ProportionMain2 <- TablingFunctionMain(PROPORTION_LIST1_XGroup, AcylChaintype = "AcylBond", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10)%>% 
  dplyr::mutate(Change = "Deficit")

#__________
#######for ether lipid separated into ether and ester chains
### Ratio greater than 1.5, in either or both 1 day, 19 day... must be statistically significant
## Observed frequency must be greater than 3% in at least on of the days (no filter on Expected frequency)
AcylBond_ProportionMainEthervinyl <- TablingFunctionMain(PROPORTION_LIST1_Ethervinyl_XGroup, AcylChaintype = "AcylBond", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")

AcylBond_ProportionMainEster <- TablingFunctionMain(PROPORTION_LIST1_Ester_XGroup, AcylChaintype = "AcylBond", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")


### Ratio Less than 0.67 in at least one of 1 day or 19 day, and must be statistically significant, 
### Expected proportion must be greater than 3% in at least one of the days
AcylBond_ProportionMain2Ethervinyl <- TablingFunctionMain(PROPORTION_LIST1_Ethervinyl_XGroup, AcylChaintype = "AcylBond", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10)%>% 
  dplyr::mutate(Change = "Deficit")

AcylBond_ProportionMain2Ester <- TablingFunctionMain(PROPORTION_LIST1_Ester_XGroup, AcylChaintype = "AcylBond", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10)%>% 
  dplyr::mutate(Change = "Deficit")

#________

# ### Merged

Merged0 <- rbind(AcylBond_ProportionMain,AcylBond_ProportionMain2) %>% #mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 


Merged0 <- Merged0 %>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)

#######for ether lipid separated into ether and ester chains
Merged0Ethervinyl <- rbind(AcylBond_ProportionMainEthervinyl,AcylBond_ProportionMain2Ethervinyl) %>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 

Merged0Ethervinyl<- Merged0Ethervinyl %>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)

Merged0Ester <- rbind(AcylBond_ProportionMainEster,AcylBond_ProportionMain2Ester)%>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 

Merged0Ester <- Merged0Ester %>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)

### Ratio greater than 1.5, in either or both 1 day, 19 day... must be statistically significant
## Observed frequency must be greater than 3% in at least on of the days (no filter on Expected frequency)
AcylLength_ProportionMain <- TablingFunctionMain(PROPORTION_LIST2_XGroup, AcylChaintype = "AcylLength", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")

### Ratio Less than 0.67 in at least one of 1 day or 19 day, and must be statistically significant, 
### Expected proportion must be greater than 3% in at least one of the days
AcylLength_ProportionMain2 <- TablingFunctionMain(PROPORTION_LIST2_XGroup, AcylChaintype = "AcylLength", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10)%>% 
  dplyr::mutate(Change = "Deficit")

#___________
### Ratio greater than 1.5, in either or both 1 day, 19 day... must be statistically significant
## Observed frequency must be greater than 3% in at least on of the days (no filter on Expected frequency)
AcylLength_ProportionMainEthervinyl <- TablingFunctionMain(PROPORTION_LIST2_Ethervinyl_XGroup, AcylChaintype = "AcylLength", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")

AcylLength_ProportionMainEster <- TablingFunctionMain(PROPORTION_LIST2_Ester_XGroup, AcylChaintype = "AcylLength", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")

### Ratio Less than 0.67 in at least one of 1 day or 19 day, and must be statistically significant, 
### Expected proportion must be greater than 3% in at least one of the days

AcylLength_ProportionMain2Ethervinyl <- TablingFunctionMain(PROPORTION_LIST2_Ethervinyl_XGroup, AcylChaintype = "AcylLength", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10)%>% 
  dplyr::mutate(Change = "Deficit")

AcylLength_ProportionMain2Ester <- TablingFunctionMain(PROPORTION_LIST2_Ester_XGroup, AcylChaintype = "AcylLength", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10)%>% 
  dplyr::mutate(Change = "Deficit")
#___________

### Merged

Merged1 <- rbind(AcylLength_ProportionMain,AcylLength_ProportionMain2) %>% #mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 

Merged1 <- Merged1 %>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)

#######for ether lipid separated into ether and ester chains
Merged1Ethervinyl <- rbind(AcylLength_ProportionMainEthervinyl,AcylLength_ProportionMain2Ethervinyl) %>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 

Merged1Ethervinyl <- Merged1Ethervinyl%>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)
  
Merged1Ester <- rbind(AcylLength_ProportionMainEster,AcylLength_ProportionMain2Ester) %>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 

Merged1Ester <- Merged1Ester%>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)

### Ratio greater than 1.5, in either or both 1 day, 19 day... must be statistically significant
## Observed frequency must be greater than 3% in at least on of the days (no filter on Expected frequency)

Bond_ProportionMain <- TablingFunctionMain(PROPORTION_LIST3_XGroup, AcylChaintype = "Bond", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")

### Ratio Less than 0.67 in at least one of 1 day or 19 day, and must be statistically significant, 
### Expected proportion must be greater than 3% in at least one of the days

Bond_ProportionMain2 <- TablingFunctionMain(PROPORTION_LIST3_XGroup, AcylChaintype = "Bond", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10) %>% 
  dplyr::mutate(Change = "Deficit")

#___________
### Ratio greater than 1.5, in either or both 1 day, 19 day... must be statistically significant
## Observed frequency must be greater than 3% in at least on of the days (no filter on Expected frequency)
Bond_ProportionMainEthervinyl <- TablingFunctionMain(PROPORTION_LIST3_Ethervinyl_XGroup, AcylChaintype = "Bond", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")

Bond_ProportionMainEster <- TablingFunctionMain(PROPORTION_LIST3_Ester_XGroup, AcylChaintype = "Bond", Greater = 1.5, Lesser = 0, Expected = -1, Observed = 5, Difference = 10) %>% 
  dplyr::mutate(Change = "Excess")

### Ratio Less than 0.67 in at least one of 1 day or 19 day, and must be statistically significant, 
### Expected proportion must be greater than 3% in at least one of the days

Bond_ProportionMain2Ethervinyl <- TablingFunctionMain(PROPORTION_LIST3_Ethervinyl_XGroup, AcylChaintype = "Bond", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10) %>% 
  dplyr::mutate(Change = "Deficit")

Bond_ProportionMain2Ester <- TablingFunctionMain(PROPORTION_LIST3_Ester_XGroup, AcylChaintype = "Bond", Greater = Inf, Lesser = 2/3, Expected = 5, Observed = -1, Difference = 10) %>% 
  dplyr::mutate(Change = "Deficit")

#___________

Merged2 <- rbind(Bond_ProportionMain, Bond_ProportionMain2) %>% #mutate(Category = str_replace_all(Category, "_","")) %>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 

Merged2 <- Merged2 %>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)

#######for ether lipid separated into ether and ester chains
Merged2Ethervinyl <- rbind(Bond_ProportionMainEthervinyl, Bond_ProportionMain2Ethervinyl)%>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 

Merged2Ethervinyl <- Merged2Ethervinyl %>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)


Merged2Ester <- rbind(Bond_ProportionMainEster, Bond_ProportionMain2Ester)%>% 
  dplyr::mutate(Day1 = paste0(`Ratio_1 Day`, " (", Diff1Day, ")")) %>% 
  dplyr::mutate(Day19 = paste0(`Ratio_19 Day`, " (", Diff_19Day, ")")) %>%  
  dplyr::select(Class, SubClass, Category, Change, Day1, Day19) 

Merged2Ester <- Merged2Ester%>% group_by(Class, SubClass, Change) %>% 
  dplyr::mutate(Category = paste(Category, collapse = "&&"),
                Day1 = paste(Day1, collapse = "&&"),
                Day19 = paste(Day19, collapse = "&&")) %>%
  ungroup() %>% 
  unique() %>% 
  pivot_wider(names_from = Change, values_from = c(Category, Day1, Day19)) %>% 
  relocate(Day1_Excess, .after=Category_Excess) %>% 
  relocate(Day19_Excess, .after=Day1_Excess)

#_________________________________________

# Save files for all lipids

write.csv(AcylBond_Proportion, "SideChainBias/Acyl_length_DoubleBond_Proportion.csv", row.names = FALSE)
write.csv(AcylLength_Proportion, "SideChainBias/Acyl_length_Proportion.csv", row.names = FALSE)
write.csv(Bond_Proportion, "SideChainBias/DoubleBond_Proportion.csv", row.names = FALSE)

write.csv(Merged0, "SideChainBias/AcylLengthBond_Proportion_Main_Greater1.5_Lesser0.67.csv", row.names = FALSE)
write.csv(Merged1, "SideChainBias/Acyl_length_Proportion_Main_Greater1.5_Lesser0.67.csv", row.names = FALSE)
write.csv(Merged2, "SideChainBias/DoubleBond_Proportion_Main_Greater1.5_Lesser0.67.csv", row.names = FALSE)

# Save files for ether lipid separated into ether and ester chains
write.csv(Merged0Ethervinyl, "SideChainBias/AcylLengthBond_Proportion_MainEthervinyl_Greater1.5_Lesser0.67.csv", row.names = FALSE)
write.csv(Merged1Ethervinyl, "SideChainBias/Acyl_length_Proportion_MainEthervinyl_Greater1.5_Lesser0.67.csv", row.names = FALSE)
write.csv(Merged2Ethervinyl, "SideChainBias/DoubleBond_Proportion_MainEthervinyl_Greater1.5_Lesser0.67.csv", row.names = FALSE)

write.csv(Merged0Ester, "SideChainBias/AcylLengthBond_Proportion_MainEster_Greater1.5_Lesser0.67.csv", row.names = FALSE)
write.csv(Merged1Ester, "SideChainBias/Acyl_length_Proportion_MainEster_Greater1.5_Lesser0.67.csv", row.names = FALSE)
write.csv(Merged2Ester, "SideChainBias/DoubleBond_Proportion_MainEster_Greater1.5_Lesser0.67.csv", row.names = FALSE)

#_______________________________________________________________________

# %s of individual chains in each subclass that were short vs medium vs long, and saturated, vs mono-unsaturated vs polyunsaturated
## Use objects 'SummaryProportion' which contain the output of SummariseLipid
### Individual proportions are calculated by line, age (1 or 19 days), time (new/old) and subclass. 
#### We use this data to generate Figure 2 of acyl chain: Patterns of acyl chain length and double bond distribution in various lipid subclasses in Day1 & Day19 B.tryoni males

FullSummarySCBond <- SummaryProportion %>% map(~.$summary) %>% bind_rows() 
FullSummarySC <- SummaryProportion2 %>% map(~.$summary) %>% bind_rows()
FullSummaryBond <- SummaryProportion3 %>% map(~.$summary) %>% bind_rows()

write.csv(FullSummarySCBond, "SideChainBias/IndividualSCBondcombined.csv", row.names = FALSE)
write.csv(FullSummarySC, "SideChainBias/IndividualSConly.csv", row.names = FALSE)
write.csv(FullSummaryBond, "SideChainBias/IndividualBondOnly.csv", row.names = FALSE)

#### We use this data to generate Figure 3 of acyl chain: Patterns of acyl chain length and double bond distribution in various lipid subclasses in Day1 & Day19 B.tryoni males

FullSummarySCBondEthervinyl <- SummaryProportionEthervinyl %>% map(~.$summary) %>% bind_rows() 
FullSummarySCEthervinyl <- SummaryProportion2Ethervinyl %>% map(~.$summary) %>% bind_rows()
FullSummaryBondEthervinyl <- SummaryProportion3Ethervinyl %>% map(~.$summary) %>% bind_rows()

write.csv(FullSummarySCBondEthervinyl, "SideChainBias/IndividualSCBondcombinedEthervinyl.csv", row.names = FALSE)
write.csv(FullSummarySCEthervinyl, "SideChainBias/IndividualSConlyEthervinyl.csv", row.names = FALSE)
write.csv(FullSummaryBondEthervinyl, "SideChainBias/IndividualBondOnlyEthervinyl.csv", row.names = FALSE)

FullSummarySCBondEster <- SummaryProportionEster %>% map(~.$summary) %>% bind_rows() 
FullSummarySCEster <- SummaryProportion2Ester  %>% map(~.$summary) %>% bind_rows()
FullSummaryBondEster <- SummaryProportion3Ester %>% map(~.$summary) %>% bind_rows()

write.csv(FullSummarySCBondEster, "SideChainBias/IndividualSCBondcombinedEster.csv", row.names = FALSE)
write.csv(FullSummarySCEster, "SideChainBias/IndividualSConlyEster.csv", row.names = FALSE)
write.csv(FullSummaryBondEster, "SideChainBias/IndividualBondOnlyEster.csv", row.names = FALSE)

#___________________________________END_______________________________________

