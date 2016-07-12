

fun_path_length_20150916<-function(i, master_network, input=NULL, checks = c("basic length","eQTL-P cutoff"),transeQTL=NULL, simulate=FALSE, hgnc_to_ensembl=NULL, niter = 10000){
  #This function searches for shortest paths in network-data, according to a set of 
  #pre-defined search-modes. It is written to be ready for multi-core processing using
  #any lapply-like function (hence the 'i' argument). the
    #i								a single integer - just corresponds to the line in input, see this (for use with qsubLappy)
  #input	 					a data frame  - for determining what path to investigate. Must contain columns indicating 'cis_genes' and 'trans_genes', which can be either genesymbol or ensembl
  #master_network		the network to use, in edge-list format
  #checks						a character vector specifying what checks to do
  #unWeightEnds			logical - if the ends (i.e. cis_gene and trans_gene) should have score 1 no matter what.
  #hgnc_to_ensembl	annotation information. Only used if input genes is of type ensembl
  #simulate					logical - indicating if the network should be randomized.
  
  
  
  
  
  library(igraph)
  
  if(!class(i)%in%c("integer","numeric"))stop("i must be class numeric")
  if(length(i)!=1)stop("i must be length 1")
  set.seed(i)
  
  if(class(simulate)!="logical")stop("simulate must be class logical")
  if(length(simulate)!=1)stop("simulate must be length 1")
  
  if(class(master_network)!="data.frame")stop("master_network must be class data.frame")
  if(length(grep("^Node1",colnames(master_network)))!=1)stop("master_network must have ONE column starting with Node1")
  if(length(grep("^Node2",colnames(master_network)))!=1)stop("master_network must have ONE column starting with Node2")
  node1name<-grep("^Node1",colnames(master_network),value=T)
  node2name<-grep("^Node2",colnames(master_network),value=T)
  if(any(master_network[,node2name]==master_network[,node1name]))stop(paste("master_network had",sum(master_network[,node2name]==master_network[,node1name]),"same-to-same edges"))
  if(any(master_network[,node2name]%in%"" | master_network[,node1name]%in%""))stop(paste("master_network had empty entries"))
  
  
  if(simulate){
    if(!is.null(input))stop("When simulating network, don't give an input argument")
  }else{
    if(class(input)!="data.frame")stop("input must be class data.frame")
    if(nrow(input)<i)stop("input now must be above the value of i")
    
    if(length(grep("cis_genes$",colnames(input)))!=1)stop("input must have ONE column ending with cis_genes")
    if(length(grep("trans_genes$",colnames(input)))!=1)stop("input must have ONE column ending with trans_genes")
    cis_gene_col<-grep("cis_genes$",colnames(input),value=T)
    trans_gene_col<-grep("trans_genes$",colnames(input),value=T)
    requiredCols<-c(cis_gene_col, trans_gene_col, "snp")
    # requiredCols<-c(cis_gene_col, trans_gene_col, "snp")
    if(!all(requiredCols%in%colnames(input)))stop(paste("These columns must be present in input:",paste(requiredCols,collapse=", ")))
  }
  
  
  
  if(class(checks)!="character")stop("checks must be class character")
  if(length(checks)<1)stop("checks must be of length 1 or more")
  implementedChecks<-c("basic length","eQTL-P raw","eQTL-P cutoff","eQTL-P cutoff unweighted-ends","eQTL-P cutoff mean-P","Edge-confidence","Edge-confidence cutoff","Degree centrality")
  if(!all(checks%in%implementedChecks))stop("Not all checks were recognised. Only these are allowed:",paste(implementedChecks,collapse=", "))
  
  if(any(c("eQTL-P cutoff","eQTL-P raw","eQTL-P cutoff unweighted-ends","eQTL-P cutoff mean-P")%in%checks)){
    if(is.null(transeQTL))stop("if eQTL weighting is needed, a transeQTL data.frame must be provided")
    
    #extracting list info if relevant
    if(class(transeQTL)=="list"){
      if(simulate){
        snp<-sample(names(transeQTL)[sapply(transeQTL,class)=="data.frame"],1)	
      }else{
        snp<-input[i,"newsnp"]
      }
      if(!snp%in%names(transeQTL))stop(paste("transeQTL was a list, which is allowed - but the snp",snp,"was not found in it"))
      transeQTL<-transeQTL[[snp]]
    }
    
    #handling NA cases by simply removing the checks
    originalChecks<-checks
    if((class(transeQTL)=="logical" && is.na(transeQTL)) | (class(transeQTL)=="character" && transeQTL[1]=="No pvalues")){
      print(paste("transeQTL was NA. Removing all eQTL weighted checks from setup and moving on"))
      checks<-checks[!checks%in%c("eQTL-P cutoff","eQTL-P raw")]
      if(length(checks)==0)stop("No checks left after removal")
    }else{
      
      if(class(transeQTL)!="data.frame")stop(paste("transeQTL must be a data.frame, not",class(transeQTL)))
      if(!"best eQTL"%in%colnames(transeQTL))stop("transeQTL must have colname 'best eQTL'")
      
      networkGenes<-unique(c(master_network[,node1name],master_network[,node2name]))
      transeQTLGenes<-rownames(transeQTL)
      print(paste(length(networkGenes),"network-genes and",length(transeQTLGenes),"trans-eQTL genes had",length(intersect(networkGenes,transeQTLGenes)),"overlaps"))
      if(length(intersect(networkGenes,transeQTLGenes))/length(union(networkGenes,transeQTLGenes))<0.2)stop("Too little overlap between trans-eQTL genes and network genes")
    }
  }
  
  
  if(any(c("Edge-confidence cutoff","Edge-confidence")%in%checks)){
    if(!"Confidence"%in%colnames(master_network))stop("If any confidence cheks are needed, the master_network needs a 'Confidence' column")
    if(class(master_network[,"Confidence"])!="numeric")stop("Confidence column must be numeric")
    if(max(master_network[,"Confidence"],na.rm=T) == min(master_network[,"Confidence"],na.rm=T))stop("Confidence column range must be greater than 0")
  }
  
  
  
  if(!is.null(hgnc_to_ensembl)){
    if(class(hgnc_to_ensembl)!="data.frame")stop("hgnc_to_ensembl must be class data.frame")
    if(!all(c("ensembl_gene_id","hgnc_symbol")%in%colnames(hgnc_to_ensembl)))stop("hgnc_to_ensembl must have columns ensembl_gene_id and hgnc_symbol")
    if(!all(rownames(hgnc_to_ensembl) == hgnc_to_ensembl[,"ensembl_gene_id"]))stop("hgnc_to_ensembl must have a column ensembl_gene_id that is also set as rownames")
    
    annotation_genes<-rownames(hgnc_to_ensembl) 
    networkGenes<-unique(c(master_network[,node1name],master_network[,node2name]))
    print(paste(length(networkGenes),"network-genes and",length(annotation_genes),"annotation-genes had",length(intersect(networkGenes,annotation_genes)),"overlaps"))
    if(length(intersect(networkGenes,annotation_genes))/length(union(networkGenes,annotation_genes))<0.2)stop("Too little overlap between trans-eQTL genes and annotation_genes genes")
    
  }
  
  #Getting start and end points as well as output container
  if(simulate){
    z<-sample(c(master_network[,1],master_network[,2]),2)
    cis_id<-z[1]
    trans_id<-z[2]
    output<-matrix(nrow=1,ncol=length(checks),dimnames=list(i,checks))
  }else{
    output<-input[i,,drop=FALSE]
    cis_id<-input[i,cis_gene_col]
    trans_id<-input[i,trans_gene_col]
  }
  
  
  #check they are present in the master network
  networkGenes<-unique(c(master_network[,node1name],master_network[,node2name]))
  if(!all(c(cis_id,trans_id)%in%networkGenes))stop("The start and end point was not present in the master network")
  
  
  #set NA eQTL values to 1
  if(any(c("eQTL-P cutoff","eQTL-P raw")%in%checks)){
    transeQTL[is.na(transeQTL[,"best eQTL"]),"best eQTL"]<-1
  }
  
  
  networks<-list()
  
  ######################################
  #This is basic length
  # W = 1
  if("basic length"%in%checks){
    print(paste("network: basic length"))
    network<-graph.edgelist(as.matrix(master_network[,c(node1name,node2name)]), directed=F)
    if(simulate)network<-rewire(network, with = keeping_degseq(niter = niter)) 
    networks[["basic length"]]<-network
  }
  
  ######################################
  #This is the eQTL-P raw
  #W = 0.5 + trans-eQTL P-value	/ 2
  if("eQTL-P raw"%in%checks){
    stop("Not implemented yet!")
    print(paste("network: eQTL-P raw"))
    p1<-transeQTL[master_network[,node1name],"best eQTL"]
    p1[is.na(p1)]<-1
    p2<-transeQTL[master_network[,node2name],"best eQTL"]
    p2[is.na(p2)]<-1
    master_network[,"edge P"]<-apply(cbind(p1, p2),1,min,na.rm=F)
    master_network[,"score"] <- master_network[,"edge P"]/2 + 0.5
    network<-graph.edgelist(as.matrix(master_network[,c(node1name,node2name)]), directed=F)
    E(network)$weight<-master_network[,"score"]
    if(simulate)network<-rewire(network, with = keeping_degseq(niter = niter)) 
    networks[["eQTL-P raw"]]<-network
    
  }
  
  
  ######################################
  #This is the eQTL-P cutoff
  #W = 1, except 
  # 	trans-eQTL P<0.05 -> W = 0.8
  # 	trans-eQTL P<0.005 -> W = 0.6
  if("eQTL-P cutoff"%in%checks){
    print(paste("network: eQTL-P cutoff"))
    p1<-transeQTL[master_network[,node1name],"best eQTL"]
    p1[is.na(p1)]<-1
    p2<-transeQTL[master_network[,node2name],"best eQTL"]
    p2[is.na(p2)]<-1
    master_network[,"edge P"]<-apply(cbind(p1, p2),1,min,na.rm=F)
    master_network[,"score"] <- 1
    master_network[master_network[,"edge P"] < 0.05, "score"]<-0.8
    master_network[master_network[,"edge P"] < 0.005, "score"]<-0.6
    network<-graph.edgelist(as.matrix(master_network[,c(node1name,node2name)]), directed=F)
    E(network)$weight<-master_network[,"score"]
    if(simulate)network<-rewire(network, with = keeping_degseq(niter = niter)) 
    networks[["eQTL-P cutoff"]]<-network
  }
  
  
  
  ######################################
  #This is the "eQTL-P cutoff unweighted-ends"
  #W = 1, except 
  # 	trans-eQTL P<0.05 -> W = 0.8
  # 	trans-eQTL P<0.005 -> W = 0.6
  #     --- and CIS and TRANS is set to P=1
  if("eQTL-P cutoff unweighted-ends"%in%checks){
    print(paste("network: eQTL-P cutoff unweighted-ends"))
    t<-transeQTL
    if(cis_id%in%rownames(t))t[cis_id,"best eQTL"]<-1
    if(trans_id%in%rownames(t))t[trans_id,"best eQTL"]<-1
    p1<-t[master_network[,node1name],"best eQTL"]
    p1[is.na(p1)]<-1
    p2<-t[master_network[,node2name],"best eQTL"]
    p2[is.na(p2)]<-1
    master_network[,"edge P"]<-apply(cbind(p1, p2),1,min,na.rm=F)
    master_network[,"score"] <- 1
    master_network[master_network[,"edge P"] < 0.05, "score"]<-0.8
    master_network[master_network[,"edge P"] < 0.005, "score"]<-0.6
    network<-graph.edgelist(as.matrix(master_network[,c(node1name,node2name)]), directed=F)
    E(network)$weight<-master_network[,"score"]
    if(simulate)network<-rewire(network, with = keeping_degseq(niter = niter)) 
    networks[["eQTL-P cutoff unweighted-ends"]]<-network
  }
  
  
  ######################################
  #This is the "eQTL-P cutoff mean-P"
  #W = 1, except 
  # 	trans-eQTL P<0.05 -> W = 0.8
  # 	trans-eQTL P<0.005 -> W = 0.6
  #     --- edge-P values are calculated as means of the P-values of the adjacent nodes
  if("eQTL-P cutoff mean-P"%in%checks){
    print(paste("network: eQTL-P cutoff mean-P"))
    master_network[,"edge P"]<-(transeQTL[master_network[,node1name],"best eQTL"] + transeQTL[master_network[,node2name],"best eQTL"]) / 2
    master_network[,"score"] <- 1
    master_network[!is.na(master_network[,"edge P"]) & master_network[,"edge P"] < 0.05, "score"]<-0.5
    master_network[!is.na(master_network[,"edge P"]) & master_network[,"edge P"] < 0.001, "score"]<-0.2
    network<-graph.edgelist(as.matrix(master_network[,c(node1name,node2name)]), directed=F)
    E(network)$weight<-master_network[,"score"]
    if(simulate)network<-rewire(network, with = keeping_degseq(niter = niter)) 
    networks[["eQTL-P cutoff mean-P"]]<-network
  }
  
  
  
  ######################################
  #This is the Edge-confidence
  #W = 1 - edge-confidence x 0.5
  
  if("Edge-confidence"%in%checks){
    print(paste("network: Edge-confidence"))
    master_network[,"score"] <- master_network[,"Confidence"] * 0.5 + 0.5
    network<-graph.edgelist(as.matrix(master_network[,c(node1name,node2name)]), directed=F)
    E(network)$weight<-master_network[,"score"]
    if(simulate)network<-rewire(network, with = keeping_degseq(niter = niter)) 
    networks[["Edge-confidence"]]<-network
  }
  
  
  ######################################
  #This is the Edge-confidence cutoff
  # 	W = 1, except 
  # 	edge-confidence > 0.9 -> W = 0.5
  # 	
  if("Edge-confidence cutoff"%in%checks){
    print(paste("network: Edge-confidence cutoff"))
    p1<-transeQTL[master_network[,node1name],"best eQTL"]
    p1[is.na(p1)]<-1
    p2<-transeQTL[master_network[,node2name],"best eQTL"]
    p2[is.na(p2)]<-1
    master_network[,"edge P"]<-apply(cbind(p1, p2),1,min,na.rm=F)
    master_network[,"score"] <- 1
    master_network[master_network[,"edge P"] < 0.05, "score"]<-0.8
    master_network[master_network[,"edge P"] < 0.005, "score"]<-0.6
    network<-graph.edgelist(as.matrix(master_network[,c(node1name,node2name)]), directed=F)
    E(network)$weight<-master_network[,"score"]
    if(simulate)network<-rewire(network, with = keeping_degseq(niter = niter)) 
    networks[["Edge-confidence cutoff"]]<-network
  }
  
  
  
  ######################################
  # This is the Degree centrality	
  # W = degree centrality
  if("Degree centrality"%in%checks){
    print(paste("network: Degree centrality"))
    network<-graph.edgelist(as.matrix(master_network[,c(node1name,node2name)]), directed=F)
    centrality<-data.frame(row.names=V(network)$name, Centrality=centralization.degree(network)[["res"]])
    p1<-centrality[master_network[,node1name],"Centrality"]
    p1[is.na(p1)]<-1
    p2<-centrality[master_network[,node2name],"Centrality"]
    p2[is.na(p2)]<-1
    master_network[,"Centrality"]<-apply(cbind(p1, p2),1,max,na.rm=F)
    E(network)$weight<-master_network[,"Centrality"]
    if(simulate)network<-rewire(network, with = keeping_degseq(niter = niter)) 
    networks[["Degree centrality"]]<-network
  }
  
  
  
  if(simulate){
    for(check in checks){
      network<-networks[[check]]
      lengthOfPath<-shortest.paths(network,cis_id,trans_id)
      output[as.character(i),check]<-lengthOfPath[1,1]
    }
  }else{
    for(check in originalChecks){
      #handling the NA transeQTL case
      if(!check %in% checks){
        output[1,paste(check,"Path length",sep=" - ")]<-NA
        output[1,paste(check,"Path-ENSG",sep=" - ")]<-NA
        if(!is.null(hgnc_to_ensembl))output[1,paste(check,"Path-HGNC",sep=" - ")]<-NA
      }else{
        print(paste("Analysing network",check))
        network<-networks[[check]]
        
        pathLength<-as.numeric(shortest.paths(network,cis_id,trans_id))
        path<-get.shortest.paths(network,cis_id,trans_id)[["vpath"]][[1]]
        path_ensembl<-V(network)$name[path]
        output[1,paste(check,"Path length",sep=" - ")]<-pathLength
        output[1,paste(check,"Path-ENSG",sep=" - ")]<-paste(path_ensembl,collapse=" // ")
        if(!is.null(hgnc_to_ensembl)){
          path_hgnc<-hgnc_to_ensembl[path_ensembl,"hgnc_symbol"]
          path_hgnc[which(is.na(path_hgnc))]<-path_ensembl[which(is.na(path_hgnc))]
          output[1,paste(check,"Path-HGNC",sep=" - ")]<-paste(path_hgnc,collapse=" // ")
        }
      }
    }
  }
  return(output)
}

