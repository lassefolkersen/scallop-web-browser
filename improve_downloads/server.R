library("shiny")



#gene positions
load("~/srv/scallop-web-browser/2014-07-16 gene locations.rdata")
protein_pos_file<-"~/srv/scallop-web-browser/2017-04-07_protein_pos_data.rdata"
load(protein_pos_file)
p<-data.frame(
  row.names=data[,"hgnc_symbol"],
  pheno_id=data[,"No_in_GWAS_files"],
  CHR = data[,"trait_chr"],
  BP = data[,"trait_pos"],
  stringsAsFactors=F)
p<-p[order(rownames(p)),]
phenotypes_vector<-p[,"pheno_id"]
names(phenotypes_vector) <- rownames(p)






shinyServer(function(input, output) {
  

  
  output$explanatoryText <- renderText({
    if(input$goButton > 0){
      # email <- isolate(input$email)
      link<-paste("http://www.scallop-consortium.com/www/dl/Olink_panel_IMPROVE_May_28th.txt",sep="")
      o<-paste("Users of <i>wget</i> may wish to batch-download data and refer to <u><a href='",link,"'>this file for protein-number to protein-name conversion</a></u>.<br><br>")
    }else{
      o<-""
    }
    return(o)
  })
  
  
  
  output$mainTable <- renderDataTable({ 
    if(input$goButton > 0){
      # email <- isolate(input$email)
      phenotype <- isolate(input$phenotype)
      
      # get_data()
      
      allFiles<-list.files("~/srv/scallop-web-browser/www/dl/")
      chr<-sub("_pheno.+$","",sub("^.+chr","",allFiles)	)
      phenoNumber<-sub("\\.txt.+$","",sub("^.+pheno","",allFiles)	)
      phenoName<-rownames(p)[p[,"pheno_id"]%in%phenotype]
      data<-data.frame(
        chr=chr,
        phenotypeName=phenoName,
        phenotypeNumber=phenoNumber,
        link=paste("http://www.scallop-consortium.com/www/dl/",allFiles,sep=""),
        stringsAsFactors = 
      )
      data<-data[data[,"phenotypeNumber"]%in%phenotype,]
      data<-data[order(as.numeric(data[,"chr"])),]
      data[,"phenotypeNumber"]<-NULL
      colnames(data)[2]<-"phenotype"
      return(data)
    }
  })
})







