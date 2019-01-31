

# install.packages("XML")


rm(list=ls())
library(XML)
ICD10_file<-"C:/Users/FOLK/Documents/Work/Analysis/2018-02-13_network_ICD-10/icdClaML2016ens.xml"

ICD10_raw_xml<-xmlParse(ICD10_file)
ICD10 <- xmlToList(ICD10_raw_xml)





out_list<-list()
i_to_skip <- 9
for(i in 1:length(ICD10)){
  # print(names(ICD10[[i]]))
  if(i%in%i_to_skip)next #because
  if(!"Rubric"%in%names(ICD10[[i]])){next} #because then it's not an entry
  if(!"Label"%in%names(ICD10[[i]][["Rubric"]])){next} #but this never happens
  if(!"text"%in%names(ICD10[[i]][["Rubric"]][["Label"]])){next} #for some types of meta info. I think.
  if(!".attrs"%in%names(ICD10[[i]])){stop()} #this would be odd
  if(!"code"%in%names(ICD10[[i]][[".attrs"]])){stop()} #this would be odd
  
  text <- ICD10[[i]][["Rubric"]][["Label"]][["text"]]
  code <- ICD10[[i]][[".attrs"]][["code"]]
  
  #if superclass
  if(!"SuperClass"%in%names(ICD10[[i]])){
    #then it is a chapter - i.e. Roman Numeral
    superclass <- "top"
  }else{
    superclass <- ICD10[[i]][["SuperClass"]]
  }
  out_list[[as.character(i)]]<-data.frame(row.names=NULL,code,text,superclass)
}




#insert additional tops
out_list[["top"]] <- data.frame(code="top",text="Heading to Hospital",superclass="Feeling fine",stringsAsFactors = F)
# out_list[["Planning a baby"]] <- data.frame(code="Planning a baby",text="Planing a baby",superclass="Feeling fine",stringsAsFactors = F)

node<-do.call(rbind,out_list)
for(i in 1:ncol(node)){
  node[,i]<-as.character(node[,i])
}

#renaming chapter entries to something short
library(openxlsx)
newNames<-read.xlsx("diseaseNetwork/2018-02-26_shorter_names.xlsx",rowNames=T)
for(code in rownames(newNames)){
  node[node[,"code"]%in%code,"code"] <- newNames[code,"Easier.names"]
  node[node[,"superclass"]%in%code,"superclass"] <- newNames[code,"Easier.names"]
}




#making a edge list
edge_list <- list()
for(i in 1:nrow(node)){
  superclass <- as.character(node[i,"superclass"])
  if(!superclass %in% node[,"code"] & superclass!="Feeling fine"){
  }else{
    edge_list[[as.character(i)]]<-data.frame(a=superclass,b=node[i,"code"])
  }
}

edge<-do.call(rbind,edge_list)
for(i in 1:ncol(edge)){
  edge[,i]<-as.character(edge[,i])
}




#rename top to sick
edge[edge[,1]%in%"top",1]<-"Heading to Hospital"
edge[edge[,2]%in%"top",2]<-"Heading to Hospital"
node[node[,"code"]%in%"top","code"]<-"Heading to Hospital"
node[node[,"superclass"]%in%"top","superclass"] <-"Heading to Hospital"


a <- "Heading to Hospital"
a1<-union(grep(a,edge[,1]),grep(a,edge[,2]))
edge[a1,]
a2<-union(grep(a,node[,"code"]),grep(a,node[,"superclass"]))
node[a2,]


#making a node list
node<-node[node[,"code"] %in% c(edge[,1],edge[,2]),]

node<-rbind(
  node,
  data.frame(row.names="Feeling fine",code="Feeling fine",text="Feeling fine",superclass=NA,stringsAsFactors = F)
)


write.table(edge[,1:2],file="diseaseNetwork/2018-02-21_ICD10_edgelist.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(node,file="diseaseNetwork/2018-02-21_ICD10_nodelist.txt",sep="\t",col.names=T,row.names=F,quote=F)


#check there's complete match
all(node[,"code"] %in% c(edge[,1],edge[,2]))
all(c(edge[,1],edge[,2])%in%c(node[,"superclass"],node[,"code"]))

node[apply(is.na(node),1,sum)>0,]
edge[apply(is.na(edge),1,sum)>0,]





## Download and install the package
# install.packages("igraph")
rm(list=ls())
library(igraph)
edges_raw <- read.table("diseaseNetwork/2018-02-21_ICD10_edgelist.txt",sep="\t",stringsAsFactors = F,header=T)
nodes_raw <- read.table("diseaseNetwork/2018-02-21_ICD10_nodelist.txt",sep="\t",stringsAsFactors = F,comment.char="",quote="",row.names=1,header=T)



e<-graph_from_edgelist(as.matrix(edges_raw), directed = TRUE)
e<-set_vertex_attr(e,"niceName",value=nodes_raw[vertex_attr(e)[["name"]],"text"])

#add root distance
focus_node<-"Feeling fine"
i <- which(V(e)$name%in%focus_node)
a<-t(distances(e,v=i))[,1]
V(e)$distance <- a


save(e,file="diseaseNetwork/2018-02-21_igraph_object.rdata")




















#2018-03-16 adding sample count and date stamp to traits file
rm(list=ls())
#from here https://github.com/lassefolkersen/impute-me/blob/master/AllDiseases/2017-02-21_semi_curated_version_gwas_central.rdata
persnp_file<-"C:/Users/FOLK/Documents/Work/Bioinformatics/2015-08-17_impute-me/AllDiseases/2017-02-21_semi_curated_version_gwas_central.rdata"
load(persnp_file)


#from here https://github.com/lassefolkersen/impute-me/blob/master/AllDiseases/2017-02-21_trait_overoverview.rdata
traits_file<-"C:/Users/FOLK/Documents/Work/Bioinformatics/2015-08-17_impute-me/AllDiseases/2017-02-21_trait_overoverview.rdata"
load(traits_file)



for(trait in rownames(traits)){
  d<-data[data[,"study_id"]%in%trait,]
  
  sampleSize <- max(d[,"sampleSize"])
  
  DATE <- unique(d[,"DATE"])
  
  #stop if problem
  if(length(sampleSize)!=1 | length(DATE)!=1) stop("Found non-unique cases")
  
  
  #insert in trait
  traits[trait,"sampleSize"] <- sampleSize
  traits[trait,"DATE"] <- DATE
  
  
}



library(openxlsx)
write.xlsx("somename-this-is-a-good-overview-of-names.xlsx")











#2018-08-02 checking what could be added
rm(list=ls())
library(openxlsx)
load("AllDiseases/2018-05-28_trait_overoverview.rdata")
l<-read.xlsx("diseaseNetwork/2018-09-18_link_file.xlsx")
t<-traits
t<-t[t[,"most_recent"],]
t<-t[t[,"sampleSize"] > 5000,]
t<-t[!rownames(t) %in% l[,"gwas_code"],]
t<-t[t[,"disease"],]

rownames(t)

l[!l[,"gwas_code"] %in% rownames(t),2]
# grep("aneurysm",rownames(t),value=T)


