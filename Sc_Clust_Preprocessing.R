require(tidyverse)
require(stringr)
require(readxl)


#obtain Geng data, all probes and genes
df <- read_excel(path = "Supplemental Data Set 2.xlsx",sheet = "All Genes", skip = 1)


## keep only expression data
a<-df%>%select(-`Probe ID`, -Description)



#clean column names
require(janitor)


a<-a%>%clean_names()
colnames(a)

#get geneids to join to mapping file
a$atgeneid<-str_extract(string = a$agi_id,pattern = "At(M|C|[:digit:]{1})(g{1})([:digit:]{5})")


# read parvula v2.1 to v2.2 list
b<-read.delim(file = "SpV2.1-V2.2_rename.list", header = F)


# extract parvula gene ids
b$testV1<-str_extract(string = b$V1,pattern = "Sp(M|C|[:digit:]{1})(g{1})([:digit:]{5})")

b$testV2<-str_extract(string = b$V2,pattern = "Sp(M|C|[:digit:]{1})(g{1})([:digit:]{5})")

b$same<-b$testV1==b$testV2

#Check for ids that changed between versions

c<-b%>%filter(same==FALSE)


# read At to Sp file

AT_Sp <- read_excel("OrthologousPairing_CdsLengthFilAfterRcBlast.xlsx")


# require(devtools)
# install_github("KTMD-plant/gaRdenbox")


AT_Sp<-AT_Sp%>%select(AtSp_cdsFil_v1...1)


#create consistent id's
AT_Sp$atgeneid<-stringr::str_extract(string = AT_Sp$AtSp_cdsFil_v1...1,
                                         pattern =  
                                           "AT(M|C|[:digit:]{1})(G{1})([:digit:]{5})")


AT_Sp$spgeneid<-stringr::str_extract(string = AT_Sp$AtSp_cdsFil_v1...1,
                                         pattern =  
                                           "Sp(M|C|[:digit:]{1})(g{1})([:digit:]{5})")

#update to parvula ids to v2.2

AT_Sp$spgeneid<-str_replace(string = AT_Sp$spgeneid, pattern = "Sp1g13020",
                            replacement = "Sp1g13030")


AT_Sp$spgeneid<-str_replace(string = AT_Sp$spgeneid, pattern = "Sp4g09080",
                            replacement = "Sp4g09085")



AT_Sp$spgeneid<-str_replace(string = AT_Sp$spgeneid, pattern = "Sp7g05590",
                            replacement = "Sp7g05595")

#remove raw column
AT_Sp<-AT_Sp%>%select(-AtSp_cdsFil_v1...1)


#create columns to join to mapping file
AT_Sp$atgeneidsmall<-str_replace(string = AT_Sp$atgeneid,pattern = "AT",replacement = "At")
AT_Sp$atgeneidsmall<-str_replace(string = AT_Sp$atgeneidsmall,pattern = "G",replacement = "g")
AT_Sp$atgeneid<-AT_Sp$atgeneidsmall


#create mapping file

#The mapping file is a TAB delimited file in which the first row shows the names of the species and the first column shows the IDs of the orthologue groups (OGs). 
#Each OG includes zero, one, or many orthologous genes in each species' column split by commas.



mapfile<-data.frame("agi"=a$agi_id)

# mapfile$agi[1:length(a$agi_id)]<-a$agi_id
mapfile$atgeneid<-stringr::str_extract(string = mapfile$agi,
                                       pattern =  
                                         "At(M|C|[:digit:]{1})(g{1})([:digit:]{5})")

mapfile<-full_join(mapfile, AT_Sp)%>% select(agi, spgeneid)%>%unique()
mapfile$orthogroup<-rownames(mapfile)

mapfile<-mapfile%>%relocate(orthogroup)

write_delim(x = mapfile,file = "mapping.txt", delim = "\t")

# create geng input file
dfgeng<-a%>%select(-atgeneid)


write_delim(x = dfgeng,file = "pseudoGeng_preprocessed.tsv",delim = "\t")





#create parvula input data

#obtain counts; note: for clust FPKM etc would be better; do we have it somewhere?
parvula_gc <- read_delim("parvula_gc.csv",   #obtained from DE analysis; ids are already v2.2
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)

#obtain treatment data to give columns meaningfull names
coldat <- read_delim("RNAseq2.csv",delim = ",")


#combine treatment information
coldat$tandt<-paste(coldat$Treatment, coldat$Timepoint, coldat$Seq_Well, sep="")                      


colnames(parvula_gc)<-coldat$tandt


#change geneid to match new species format
parvula_gc$spgeneid<-str_replace(string = parvula_gc$NANANA,
                                 pattern = "Tp",
                                 replacement = "Sp")

#remove remaining column and move spgeneid to start
parvula_gc<-parvula_gc%>%select(-NANANA)%>%relocate(spgeneid)


#write file
write_delim(x = parvula_gc,file = "parvula_gc_prep.tsv",delim = "\t")
