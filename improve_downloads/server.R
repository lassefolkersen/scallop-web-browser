library("shiny")


#permissions
# accepted_users<-tolower(read.table("/home/ubuntu/misc/accepted_emails.txt",sep="\t",header=F,stringsAsFactors=F)[,1])


#gene positions
load("~/srv/olink-scallop/2014-07-16 gene locations.rdata")
protein_pos_file<-"~/srv/olink-scallop/2017-04-07_protein_pos_data.rdata"
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
  
  
  
  #main data gathering function	
  # get_data <- reactive({
    # if(input$goButton == 0){
      # return(NULL)
    # }else{
      
      ##################################
      #input-variables, log and register	
      ##################################
      # email <- isolate(input$email)
      # gene <- isolate(input$gene)
      # distance <- isolate(input$distance)
      # top_label_count<-isolate(input$top_label_count)
      # phenotype <- isolate(input$phenotype)
      
      # if(!tolower(email) %in% accepted_users ){
      #   m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email)
      #   m<-paste(m,collapse="\t")
      #   write(m,file="/home/ubuntu/logs/illegal_access_log.txt",append=TRUE)
      #   Sys.sleep(2)
      #   stop("In the test-phase non-privileged users are not allowed")
      # }else{
      #   m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),email,"download",phenotype)
      #   m<-paste(m,collapse="\t")
      #   write(m,file="/home/ubuntu/logs/log.txt",append=TRUE)
      # }
    # }
  # })
  
  
  
  output$explanatoryText <- renderText({
    if(input$goButton > 0){
      # email <- isolate(input$email)
      link<-paste("http://www.scallop-consortium.com/www/dl/Olink_panel_IMPROVE_May_28th.txt",sep="")
      # if(tolower(email) %in% accepted_users ){
      o<-paste("Users of <i>wget</i> may wish to batch-download data and refer to <u><a href='",link,"'>this file for protein-number to protein-name conversion</a></u>.<br><br>")
      # }else{
      # o<-""
      # }
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
      
      allFiles<-list.files("~/srv/olink-scallop/www/dl/")
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



