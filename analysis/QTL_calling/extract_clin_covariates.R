library(data.table)
library(stringr)
library(ggplot2)
library(optparse)

#######################################################
### Author Gracie 2021
### 
### This script takes in the REDCAP export and 
### list of individuals, reformats genotype IDs
### to match redcap and extracts clinical covariates
#####################################################


###############
### READ IN ARGS
################

##accept commandline args
option_list = list(
make_option(c("-I", "--inds"), action="store", default=NULL, type='character',
              help="donor file, one line per donor"),
make_option(c("-M", "--meta"), action="store", default=NULL, type='character',
              help="metadata file"),
make_option(c("-O", "--out_file"), action="store", default=NULL, type='character',
              help="outfile"))

opt = parse_args(OptionParser(option_list=option_list))

donor_f=opt$inds
meta_f=opt$meta
outfile=opt$out_file

##read in list of individuals
inds_toget=as.data.frame(fread(donor_f,header=FALSE))$V1


##############
## MAIN ######
##############

#reformat to match REDCAP
# microbiome inds labeled MIB_0## in redcap but MIB_## in our genotype files
# ICC labeled _0### in redcap but ICC_X/S/M_#### in genotype files
t1=ifelse((nchar(inds_toget) == 6) & (substr(inds_toget,1,3)=="MIB"), str_replace(inds_toget,"MIB_","MIB_0"),inds_toget)
t2=ifelse((substr(t1,1,3)=="ICC"), substr(t1,6,length(t1)), t1)
t3=ifelse((substr(t2,1,1)=="_") & (nchar(t2) == 4), str_replace(t2,"_","_0"), t2)
mapping=as.data.frame(cbind(t3,inds_toget))
colnames(mapping)=c("participant_id","geno_ID")


## get metadata from READCAP, and filter
## Donor info is in multiple entries, get those entries and merge them
df=as.data.frame(fread(meta_f))
df_f=df[df$participant_id %in% t3,]
#split data types
df_enroll=df_f[df_f$redcap_event_name=="enrollment_arm_1",]
df_other=df_f[(df_f$redcap_event_name=="nonlongitudinal_vi_arm_1") ,]

##keep only relevant columns
df_other2=df_other[c('participant_id','age_visit')]
df_enroll2=df_enroll[c('participant_id','gender','ethnicity')]

#merge dataframes to get all data by individual
mdf=merge(df_enroll2,df_other2,by=c("participant_id"))

mdf=merge(mapping,mdf,by=c("participant_id"))
mdf=unique(mdf)
mdf[-1]

write.table(mdf[-1],outfile,quote=FALSE,row.names=FALSE)
