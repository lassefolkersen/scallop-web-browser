


rm(list=ls())
library(openxlsx)
#read the systems biology data (this is basically the supplementary data file to the manuscript, so should be openly available elsewhere)
d<-read.xlsx("C:/Users/FOLK/Documents/Work/Bioinformatics/2017-03-31_scallop/2018-11-23_systems_biology/network_analysis_combined.sorted.xlsx",startRow=2)


#also read S1, for translating to rs-ids
S1<-read.xlsx("C:/Users/FOLK/Documents/Work/Bioinformatics/2017-03-31_scallop/article_cvdI/supplementary/S01_table_hit_list.xlsx")

proteins<-c('ADM','AGRP','Beta-NGF','CA-125','CASP-8','CCL20','CCL3','CCL4','CD40','CD40-L','CHI3L1','CSF-1','CSTB','CTSD','CTSL1','CX3CL1','CXCL1','CXCL16','CXCL6','Dkk-1','ECP','EGF','EN-RAGE','ESM-1','FABP4','FAS','FGF-23','FS','GAL','Gal-3','GDF-15','GH','HB-EGF','HGF','hK11','HSP_27','IL-18','IL-1ra','IL-27','IL-6','IL-6RA','IL-8','IL16','ITGB1BP2','KIM-1','KLK6','LEP','LOX-1','mAmP','MB','MCP-1','MMP-1','MMP-10','MMP-12','MMP-3','MMP-7','MPO','NEMO','NT-pro_BNP','OPG','PAPPA','PAR-1','PDGF_subunit_B','PECAM-1','PlGF','PSGL-1','PTX3','RAGE','REN','RETN','SCF','SELE','SIRT2','SPON1','ST2','t-PA','TF','TIE2','TM','TNF-R1','TNF-R2','TNFSF14','TRAIL','TRAIL-R2','TRANCE','U-PAR','VEGF-A','VEGF-D')



for(target_protein in proteins){
  
  
  if(!target_protein %in% d[,"trait_protein_olink_name.(network_analysis_01)"])stop("Target protein not found")
  d1<-d[d[,"trait_protein_olink_name.(network_analysis_01)"] %in% target_protein,]
  
  
  
  # creating edge-pair list for cytoscape. First step - show all paths, for the one strongest hit for each markername in a single trait protein
  output <-data.frame(source=vector(),target=vector(),type=vector(),strength=vector())
  for(snp in unique(d1[,"marker_name.(network_analysis_01)"])){
    
    #get the block of connections present for this trait protein AND this SNP
    d2<-d1[d1[,"marker_name.(network_analysis_01)"]%in%snp,]
    
    
    #get the RS-id (this was the only purpose of why we loaded S1)
    markerName <- d2[1,1]
    rsid <- unique(S1[S1[,"MarkerName"] %in% markerName,"rs-id"])
    if(length(rsid)!=1)stop("!!!")
    if(d2[1,"path_pvalue.(network_analysis_01)"  ] > d2[2,"path_pvalue.(network_analysis_01)"  ])stop("odd, these data should be sorted")
    
    
    #take either path with the best p-value - or all significant (whichever is most paths)
    significant_count <- sum(d2[,"path_pvalue.(network_analysis_01)"  ]<0.05,na.rm=T)
    if(significant_count>1){
      d3<-d2[1:significant_count,]
      # stop("This - that there is more than one significant mediator-gene - could be implement but it is currently not necessary because there is not")
    }else{
      d3<-d2[1,]
    }
    
    
    for(cis_gene_i in 1:nrow(d3)){ 
      
      
      paths <- strsplit(d3[cis_gene_i,"all_shortest_paths.(network_analysis_01)"],";")[[1]]
      best_pvalue <- d3[cis_gene_i,"path_pvalue.(network_analysis_01)" ]
      
      
      #iterate over all these paths 
      for(path in paths){
        path_split <- strsplit(path," - ")[[1]]
        #add in each gene-to-gene edge
        for(j in 1:(length(path_split)-1)){
          o1 <- data.frame(
            source=gsub(" ","",path_split[j+1]),
            target=gsub(" ","",path_split[j]),
            type="PPI",
            strength= -log10(best_pvalue))
          output <- rbind(output,o1)
        }
        
        #add in the SNP to source-gene edge
        o2 <- data.frame(
          source=rsid,
          target=gsub(" ","",path_split[length(path_split)]),
          type="SNP_to_PPI",
          strength= -log10(best_pvalue))
        output <- rbind(output,o2)
        
        
        
      }
    }
  }
  
  
  #Add in the CIS-connections
  S1a<-S1[S1[,"Protein"] %in% target_protein,]
  cis_snps <- S1a[S1a[,"Cis/Trans.(1.MB)"]%in%"cis","rs-id"]
  for(cis_snp in cis_snps){
    o2 <- data.frame(
      source=cis_snp,
      target=sub("_HUMAN","",gsub(" ","",d1[1,"trait_protein_upid.(network_analysis_01)"])),
      type="cis_SNP",
      strength= 4) #completely arbitrary - yes
    output <- rbind(output,o2)
  }
  
  
  #remove _HUMAN suffix
  output[,"target"]<-sub("_HUMAN","",output[,"target"])
  output[,"source"]<-sub("_HUMAN","",output[,"source"])
  
  #adding some example size edges
  # o3 <- data.frame(source="A-1-0",target="B-1-0",type="legend",strength= -log10(1.0))
  # o4 <- data.frame(source="A-0-05",target="B-0-05",type="legend",strength= -log10(0.05))
  # o5 <- data.frame(source="A-0-005",target="B-0-005",type="legend",strength= -log10(0.005))
  # o6 <- data.frame(source="A-0-0005",target="B-0-0005",type="legend",strength= -log10(0.0005))
  # output<-rbind(output,o3,o4,o5,o6)
  
  
  #remove double edges
  output<-output[order(output[,"strength"],decreasing=T),]
  output <- output[!duplicated(paste(output[,"source"],output[,"target"],sep="-")),]
  
  
  
  write.table(output,file=paste0("scallop_pathways/2019-01-31_",target_protein,"_pathway.txt"),sep="\t",col.names=T,row.names=F,quote=F)
  

  
}