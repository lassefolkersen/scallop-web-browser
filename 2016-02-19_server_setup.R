# 2016-02-19 server setup
#First copy off the usual shiny server AIM



#then get data
cd ~/data
mkdir 2015-12-01_olink_gwas_data
mkdir 2016-02-19_splits
cd 2015-12-01_olink_gwas_data
scp lasfol@computerome.cbs.dtu.dk:/home/people/lasfol/dataBulk/Miscellanous/2015-12-01_olink_gwas_data/* .




#for DNA-aspect, split in smaller sized chunks the data
# go to R	
rm(list=ls())
load("Important R-images and cel files/2015-09-28 myfunctions.rdata")

# in_dir<-"~/data/2015-12-01_olink_gwas_data"
# out_dir<-"~/data/2016-02-19_splits"

out_dir<-"~/2016-02-19_splits"
in_dir<-"~/dataBulk/Miscellanous/2015-12-01_olink_gwas_data"
files<-as.list(list.files(in_dir))
names(files)<-list.files(in_dir)



splitter<-function(file,in_dir,out_dir){
	window<-1000000
	#window size 1000000 will produce 229002 files, 17G -- maybe ok?
	pheno<-sub("\\..+$","",sub("^.+pheno","",file))
	chr<-sub("_pheno.+$","",sub("^olinkGWAS_chr","",file))
	print(paste("file",file,"processing: this is chr",chr,"pheno",pheno))
	data<-read.table(paste(in_dir,file,sep="/"),header=T)
	start<-floor(min(data[,"BP"])/window)*window
	end<-floor(max(data[,"BP"])/window)*window
	for(startPoints in seq(start,end,by=window)){
		print(startPoints)
		endPoints <- startPoints + window
		dh<-data[data[,"BP"] > startPoints & data[,"BP"] < endPoints,]
		dh<-dh[,c("SNP","CHR","BP","P")]
		fileName<-paste(out_dir,"/split_pheno_",pheno,"_chr",chr,"_",as.integer(startPoints/window),"_",as.integer(endPoints/window),".txt.gz",sep="")
		gz1 <- gzfile(fileName, "w")
		write.table(dh, gz1,sep="\t",col.names=T,row.names=F,quote=F)
		close(gz1)
	}
}


out<-qsubLapply(files, splitter, in_dir=in_dir,out_dir=out_dir,ncpus=50,jobName="split",verbose=3,waitTime=0.5)

#

cd /home/ubuntu/data/2016-02-19_splits
scp lasfol@computerome.cbs.dtu.dk:/home/people/lasfol/2016-02-19_splits/* .











#for protein-aspect, trim into only significant hits
#also prune, so that it only the top SNP in any given 1MB interval
# go to R	
rm(list=ls())
load("Important R-images and cel files/2015-09-28 myfunctions.rdata")


out_dir<-"~/2016-02-20_significant_bits"
in_dir<-"~/dataBulk/Miscellanous/2015-12-01_olink_gwas_data"
cutoff <- 1e-4
phenotypes<-as.list(sort(unique(sub("\\..+$","",sub("^.+pheno","",list.files(in_dir))))))
names(phenotypes)<-as.character(phenotypes)


trimmer<-function(pheno,in_dir,out_dir,cutoff){
	files<-list.files(in_dir,pattern=paste("pheno",pheno,"\\.txt",sep=""))
	print(paste("Pheno",pheno,"had",length(files),"files"))
	cols<-c("SNP","CHR","BP","P")
	data<-as.data.frame(matrix(nrow=0,ncol=length(cols),dimnames=list(NULL,cols)))
	for(file in files){
		print(file)
		d<-read.table(paste(in_dir,file,sep="/"),header=T)
		d<-d[d[,"P"]<cutoff,cols]
		d[,"BP_round"]<-round(d[,"BP"]/1000000)
		d<-d[order(d[,"BP_round"],d[,"P"]),]
		d<-d[!duplicated(d[,"BP_round"]),cols]
		data<-rbind(data,d)
	}
	fileName<-paste(out_dir,"/trimmed_pheno_",pheno,".txt.gz",sep="")
	gz1 <- gzfile(fileName, "w")
	write.table(data, gz1,sep="\t",col.names=T,row.names=F,quote=F)
	close(gz1)
}
out<-qsubLapply(phenotypes, trimmer, in_dir=in_dir,out_dir=out_dir,cutoff=cutoff, ncpus=84,jobName="trim",verbose=3,waitTime=0.5)


cd data/2016-02-20_significant_bits
scp lasfol@computerome.cbs.dtu.dk:/home/people/lasfol/2016-02-20_significant_bits/* .













#create symlinks for data download
mkdir /srv/shiny-server/www/
cd /srv/shiny-server/www/
	
for x in $(ls /home/ubuntu/data/2015-12-01_olink_gwas_data/)
do 
echo "/home/ubuntu/data/2015-12-01_olink_gwas_data/"
sudo ln -s /home/ubuntu/data/2015-12-01_olink_gwas_data/$x
done


cd /srv/shiny-server/www
sudo ln -s /srv/shiny-server/olink-improve/Olink_panel_IMPROVE_May_28th.txt

	

	
	
	
	
	
	
#getting names and positions of proteins in order
rm(list=ls())
library(biomaRt)
data<-read.table("Olink_panel_IMPROVE_May_28th.txt",sep="\t",header=T,stringsAsFactors=F)

data[,"short_name"]
data[data[,"gene"]%in%"BNP","gene"]<-"NPPC"
data[data[,"gene"]%in%"CTSL1","gene"]<-"CTSL"
data[data[,"gene"]%in%"IL8","gene"]<-"CXCL8"

rownames(data)<-data[,"gene"]
mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
attributes<-c("hgnc_symbol","start_position","end_position","chromosome_name")
trait_pos <- getBM(attributes=attributes, filters="hgnc_symbol",values=rownames(data),mart=mart)
trait_pos<-trait_pos[nchar(trait_pos[,"chromosome_name"])%in%1:2,]
rownames(trait_pos)<-trait_pos[,"hgnc_symbol"]
colnames(trait_pos)[4]<-"trait_chr"
colnames(trait_pos)[2]<-"trait_pos"


data<-cbind(data,trait_pos[rownames(data),])

save(data,file="2016-02-22_protein_pos_data.rdata")











#this is the shiny-server.conf file
ubuntu@ip-172-31-31-1:/srv/shiny-server$ sudo vi index.html
ubuntu@ip-172-31-31-1:/srv/shiny-server$ cd /etc/shiny-server/
	ubuntu@ip-172-31-31-1:/etc/shiny-server$ ls
shiny-server.conf
ubuntu@ip-172-31-31-1:/etc/shiny-server$ sudo vi shiny-server.conf
# Instruct Shiny Server to run applications as the user "shiny"
run_as shiny;

# Define a server that listens on port 3838
server {
	listen 80;
	
	# Define a location at the base URL
	location / {
		run_as ubuntu;
		
		# Host the directory of Shiny Apps stored in this directory
		site_dir /srv/shiny-server;
		
		# Log all Shiny output to files in this directory
		log_dir /var/log/shiny-server;
		
		# When a user visits the base URL rather than a particular application,
		# an index of the applications available in this directory will be shown.
		directory_index off;
	}
}




~
	
  
  
#2017-04-07  
#re-doing list of genes, because a ' left out the two last ones
#getting names and positions of proteins in order
rm(list=ls())
library(biomaRt)
data<-read.table("Olink_panel_IMPROVE_May_28th.txt",sep="\t",header=T,stringsAsFactors=F,quote="")

data[,"short_name"]
data[data[,"gene"]%in%"BNP","gene"]<-"NPPC"
data[data[,"gene"]%in%"CTSL1","gene"]<-"CTSL"
data[data[,"gene"]%in%"IL8","gene"]<-"CXCL8"

rownames(data)<-data[,"gene"]
mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
attributes<-c("hgnc_symbol","start_position","end_position","chromosome_name")
trait_pos <- getBM(attributes=attributes, filters="hgnc_symbol",values=rownames(data),mart=mart)
trait_pos<-trait_pos[nchar(trait_pos[,"chromosome_name"])%in%1:2,]
rownames(trait_pos)<-trait_pos[,"hgnc_symbol"]
colnames(trait_pos)[4]<-"trait_chr"
colnames(trait_pos)[2]<-"trait_pos"


data<-cbind(data,trait_pos[rownames(data),])

save(data,file="2017-04-07_protein_pos_data.rdata")


#comparing
rm(list=ls())
load("2017-04-07_protein_pos_data.rdata")
new<-data
load("2016-02-22_protein_pos_data.rdata")
old<-data

cor.test(new[rownames(old),"trait_pos"],old[rownames(old),"trait_pos"])
#good

na_in_new<-rownames(new)[is.na(new[,"trait_pos"])]
new[na_in_new,"trait_pos"]<-old[na_in_new,"trait_pos"]
new[na_in_new,"trait_chr"]<-old[na_in_new,"trait_chr"]
new[na_in_new,"end_position"]<-old[na_in_new,"end_position"]

save(data,file="2017-04-07_protein_pos_data.rdata")

