



#2015-09-17 this is the short-path calculations step - but now we do it parallelized with the myfunctions entries
rm(list=ls())
library(igraph)

source("fun_path_length_20150916.R")

load("Important R-images and cel files/2015-08-17 HGNC to Ensembl.rdata")
input<-read.table("/home/people/lasfol/Analysis/2014-12-19 Olink project/2015-02-05 from Eric input-data/2015-09-03 all loci all genes - ensembl and string-only.txt",sep="\t",header=T,stringsAsFactors=FALSE)
load("/home/people/lasfol/Analysis/2014-12-19 Olink project/2015-02-05 from Eric input-data/2015-04-27 trans-eQTL data.rdata")


load("/home/projects/allelic_imbalance/data/dataBulk/Miscellanous/2015-01-09 inweb/2015-08-10_string_v10/2015-08-12 string-ver10-cutoff-400.rdata")


load("/home/projects/allelic_imbalance/data/dataBulk/Miscellanous/2015-01-09 inweb/2015-01-13 inweb 5.5/2015-08-10 inweb-5-5.rdata")
inweb<-inweb[inweb[,"Confidence"]>0.2,]




iterations<-1000
checks<-c("basic length","eQTL-P cutoff", "eQTL-P cutoff unweighted-ends")


#making the string-nulls
l<-as.list(1:iterations)
names(l)<-as.character(1:iterations)
out<-qsubLapply(l,fun_path_length_20150916,transeQTL=transeQTL,master_network=string,checks=checks,simulate=TRUE,ncpus=50,jobName="rewire",verbose=3,waitTime=0.3)
pathLengths <- do.call(rbind,out)
save(pathLengths,file=paste("/home/people/lasfol/Analysis/2014-12-19 Olink project/2015-09-03 revisiting paths/2015-09-18_string_null.rdata",sep=""))

#making the string-paths
l<-as.list(1:nrow(input))
names(l)<-as.character(1:nrow(input))
out<-qsubLapply(l,fun_path_length_20150916,input=input,transeQTL=transeQTL,master_network=string,checks=checks,simulate=FALSE,hgnc_to_ensembl=hgnc_to_ensembl,ncpus=50,jobName="path",verbose=3,waitTime=0.3)
input_string <- do.call(rbind,out)
save(input_string,file=paste("/home/people/lasfol/Analysis/2014-12-19 Olink project/2015-09-03 revisiting paths/2015-09-18_string_path.rdata",sep=""))



#making the inweb-nulls
l<-as.list(1:iterations)
names(l)<-as.character(1:iterations)
out<-qsubLapply(l,fun_path_length_20150916,transeQTL=transeQTL,master_network=inweb,checks=checks,simulate=TRUE,ncpus=50,jobName="rewire",verbose=3,waitTime=0.3)
pathLengths <- do.call(rbind,out)
save(pathLengths,file=paste("/home/people/lasfol/Analysis/2014-12-19 Olink project/2015-09-03 revisiting paths/2015-09-18_inweb_null.rdata",sep=""))

#making the inweb-paths
l<-as.list(1:nrow(input))
names(l)<-as.character(1:nrow(input))
out<-qsubLapply(l,fun_path_length_20150916,input=input,transeQTL=transeQTL,master_network=inweb,checks=checks,simulate=FALSE,hgnc_to_ensembl=hgnc_to_ensembl,ncpus=50,jobName="path",verbose=3,waitTime=0.3)
input_inweb <- do.call(rbind,out)
save(input_inweb,file=paste("/home/people/lasfol/Analysis/2014-12-19 Olink project/2015-09-03 revisiting paths/2015-09-18_inweb_path.rdata",sep=""))




