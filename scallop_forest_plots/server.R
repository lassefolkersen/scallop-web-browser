library("shiny")
library("forestplot")
library("openxlsx")

#data path
hit_list_path <- "/home/ubuntu/data/2019-04-23_forest_plot_data//S01_table_hit_list.xlsx"
per_study_path <- "/home/ubuntu/data/2019-04-23_forest_plot_data/2019-02-24_per_study_data.txt.gz"


#read S1
data <- read.xlsx(hit_list_path)
data[,"data_line"] <- 1:nrow(data)
rownames(data) <- apply(data[,c("Protein","MarkerName")],1,paste,collapse=" / ")


#read new per study data and adding colnames (per SCALLOP standards)
per_study <- read.table(per_study_path,sep="\t",stringsAsFactors = F)
colnames(per_study) <- c("MARKERID","SNPID","CHR","POS","STRAND","N","EFFECT_ALLELE","REFERENCE_ALLELE","CODE_ALL_FQ","BETA","SE","PVAL","RSQ","INFO","IMP","Protein","Study")



#for matching with the direction column in S1
study_order <- c("EpiHealth","Estonian_Biobank","IMPROVE","INTERVAL","LifeLinesDeep","MPP_RES","NSPHS","ORCADES","PIVUS","STABILITY","STANLEY_lah1","STANLEY_swe6","ULSAM","VIS")



#permissions
accepted_users<-tolower(read.table("/home/ubuntu/misc/accepted_emails.txt",sep="\t",header=F,stringsAsFactors=F)[,1])


# Define server logic for a template
shinyServer(function(input, output) {
  
  
  #Get the pre-calculated genetic data for this user
  output$mainPlot <- renderPlot({ 
    if(input$goButton == 0){
      return(NULL)
    }
    
    
    s1_entry<-isolate(input$s1_entry)
    protein<-data[s1_entry,"Protein"]
    chr<-sub(":.+$","",data[s1_entry,"chr:pos.(MB)"])
    MarkerName<-data[s1_entry,"MarkerName"]
    
    email <- NA #isolate(input$email)

    #the mail checker and logger - first part is not needed when we don't check mail, but left in anyway
    if(!tolower(email) %in% accepted_users &  FALSE){
      m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email)
      m<-paste(m,collapse="\t")
      write(m,file="/home/ubuntu/logs/illegal_access_log.txt",append=TRUE)
      Sys.sleep(2)
      stop(safeError("In the test-phase non-privileged users are not allowed"))
    }else{
      m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),email,"scallop_forrest_plots",s1_entry)
      m<-paste(m,collapse="\t")
      write(m,file="/home/ubuntu/logs/log.txt",append=TRUE)
    }
    
    
    
    d <- per_study[per_study[,"MARKERID"] %in% MarkerName & per_study[,"Protein"] %in% protein,]
    
    #re-implement the IMPROVE genotype/impute filter
    if(sum(d[,"Study"]%in%"IMPROVE")>1){
      d<-d[!(d[,"Study"] %in% "IMPROVE" & d[,"IMP"] %in% 1),]  
    }
    
    
    
    #In the step_4 deletions, there is a few bug fixes that are *not* included in the step_3 files. And the step_4 half-way 
    # calculations were deleted because of space-limitations on the SCALLOP server
    direction <- strsplit(data[s1_entry,"Direction"],"")[[1]]
    names(direction) <- study_order
    direction<-direction[!direction%in%"?"]
    
    #then check it's consistent with re-loaded step-3 data
    betas <- d[,"BETA"]
    names(betas)<- d[,"Study"]
    
    
    #allele flipped check-up
    d[,"allele_flipped"] <- suppressWarnings(!as.numeric(paste0(direction[names(betas)],1)) != sign(betas))
    d<-d[!is.na(d[,"allele_flipped"]),]
    
    #Flipping alleles as needed - the benefit of taking from S1 directly is that
    #it ensures that all step-4 and step-5 filters are in.
    d[!d[,"allele_flipped"],"BETA"] <- -d[!d[,"allele_flipped"],"BETA"]
    
    
    
    #optional: order according to alternative study-order
    # study_order_2 <- c("STANLEY_lah1","STANLEY_swe6","ULSAM","PIVUS","NSPHS","IMPROVE","Estonian_Biobank","VIS","ORCADES","EpiHealth","INTERVAL","LifeLinesDeep","STABILITY","MPP_RES")
    # d<-d[order(factor(d[,"Study"],levels=study_order_2)),]
    # re_naming <- study_order_2
    # names(re_naming) <- study_order_2
    # re_naming[1:6] <- paste(re_naming[1:6],"(CVD1)")
    # re_naming[7:12] <- paste(re_naming[7:12],"(non-CVD1)")
    # d[,"Study"]<-re_naming[d[,"Study"]]
    
    
    
    #plotting
    forestplot(
      labeltext=d[,"Study"],
      mean=d[,"BETA"], 
      lower= d[,"BETA"]-d[,"SE"]*1.96, 
      upper= d[,"BETA"]+d[,"SE"]*1.96,
      title=paste0(MarkerName," and ",protein)
    )
    
    
    
    
    
  })
  

})








