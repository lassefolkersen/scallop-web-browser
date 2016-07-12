library("shiny")


#saving global variables
#chr lengths #from here https://support.bioconductor.org/p/14766/
# cl<-c(247197891, 242713278, 199439629, 191246650, 180727832,
#       170735623, 158630410, 146252219, 140191642, 135347681,
#       134361903, 132289533, 114110907, 106354309, 100334282,
#       88771793,  78646005,  76106388,  63802660,  62429769,
#       46935585, 49396972, 154908521,  57767721)
# chrLengths<-data.frame(row.names=c(1:22,"X","Y"),sum=cl,number=1:length(cl),endsum=cumsum(cl))
# chrLengths[,"startsum"]<- chrLengths[,"endsum"] - chrLengths[,"sum"]
# 
# #Misc small variables
# p_range<-1
# maxPos<-sum(cl)#max(data[,"snp_abs_pos"],na.rm=T)
# base<-2
# ip<-read.table("/home/ubuntu/misc/current_address.txt",stringsAsFactors=F)[1,1]
# 
# #permissions
# accepted_users<-tolower(read.table("/home/ubuntu/misc/accepted_emails.txt",sep="\t",header=F,stringsAsFactors=F)[,1])
# 
# #colouring scheme
# #Setting colours
# # sum(cl)/1e7 this is 307 (meaning we could divide the genome into 307 1e7 bp chunks. Let's.
# cols<-rainbow(307, s = 1, v = 1, start = 0, end = 1, alpha = 1)
# names(cols) <- as.character(1:length(cols))


#gene positions
load("/srv/shiny-server/olink-improve/2014-07-16 gene locations.rdata")
protein_pos_file<-"/srv/shiny-server/olink-improve/2016-02-22_protein_pos_data.rdata"
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
  get_data <- reactive({
    if(input$goButton == 0){
      return(NULL)
    }else{
      
      ##################################
      #input-variables, log and register	
      ##################################
      email <- isolate(input$email)
      # gene <- isolate(input$gene)
      # distance <- isolate(input$distance)
      # top_label_count<-isolate(input$top_label_count)
      phenotype <- isolate(input$phenotype)
      
      if(!tolower(email) %in% accepted_users ){
        m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email)
        m<-paste(m,collapse="\t")
        write(m,file="/home/ubuntu/logs/illegal_access_log.txt",append=TRUE)
        Sys.sleep(2)
        stop("In the test-phase non-privileged users are not allowed")
      }else{
        m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),email,"download",phenotype)
        m<-paste(m,collapse="\t")
        write(m,file="/home/ubuntu/logs/log.txt",append=TRUE)
      }
    }
  })
  
  
  
  output$explanatoryText <- renderText({
    if(input$goButton > 0){
      email <- isolate(input$email)
      link<-paste(ip,"www/","Olink_panel_IMPROVE_May_28th.txt",sep="")
      if(tolower(email) %in% accepted_users ){
        o<-paste("Users of <i>wget</i> may wish to batch-download data and refer to <u><a href='",link,"'>this file for protein-number to protein-name conversion</a></u>.<br><br>")
      }else{
        o<-""
      }
    }else{
      o<-""
    }
    return(o)
  })
  
  
  
  output$mainTable <- renderDataTable({ 
    email <- isolate(input$email)
    phenotype <- isolate(input$phenotype)
    
    get_data()
    
    allFiles<-list.files("/srv/shiny-server/www")
    chr<-sub("_pheno.+$","",sub("^.+chr","",allFiles)	)
    pheno<-sub("\\.txt.+$","",sub("^.+pheno","",allFiles)	)
    data<-data.frame(
      chr=chr,
      phenotype=pheno,
      link=paste("http://www.olink-improve.com/www/",allFiles,sep="")
    )
    data<-data[data[,"phenotype"]%in%phenotype,]
    data<-data[order(data[,"chr"]),]
    return(data)
    
  })
})



