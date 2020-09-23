library("shiny")



#saving global variables
#chr lengths #from here https://support.bioconductor.org/p/14766/
cl<-c(247197891, 242713278, 199439629, 191246650, 180727832,
      170735623, 158630410, 146252219, 140191642, 135347681,
      134361903, 132289533, 114110907, 106354309, 100334282,
      88771793,  78646005,  76106388,  63802660,  62429769,
      46935585, 49396972, 154908521,  57767721)
chrLengths<-data.frame(row.names=c(1:22,"X","Y"),sum=cl,number=1:length(cl),endsum=cumsum(cl))
chrLengths[,"startsum"]<- chrLengths[,"endsum"] - chrLengths[,"sum"]

#Misc small variables
p_range<-1
maxPos<-sum(cl)#max(data[,"snp_abs_pos"],na.rm=T)
base<-2
ip<-read.table("/home/ubuntu/misc/current_address.txt",stringsAsFactors=F)[1,1]

#permissions
# accepted_users<-tolower(read.table("/home/ubuntu/misc/accepted_emails.txt",sep="\t",header=F,stringsAsFactors=F)[,1])

#colouring scheme
#Setting colours
# sum(cl)/1e7 this is 307 (meaning we could divide the genome into 307 1e7 bp chunks. Let's.
cols<-rainbow(307, s = 1, v = 1, start = 0, end = 1, alpha = 1)
names(cols) <- as.character(1:length(cols))


#gene positions
load("~/srv/scallop-web-browser/location_files/2014-07-16 gene locations.rdata")
protein_pos_file<-"~/srv/scallop-web-browser/location_files/2017-04-07_protein_pos_data.rdata"
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
      # email <- isolate(input$email)
      gene <- isolate(input$gene)
      distance <- isolate(input$distance)
      top_label_count<-isolate(input$top_label_count)
      phenotype <- isolate(input$phenotype)
      
      # if(!tolower(email) %in% accepted_users ){
      #   m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email)
      #   m<-paste(m,collapse="\t")
      #   write(m,file="/home/ubuntu/logs/illegal_access_log.txt",append=TRUE)
      #   Sys.sleep(2)
      #   stop("In the test-phase non-privileged users are not allowed")
      # }else{
      #   
      m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),NA,"improve_dna_region",phenotype, gene, distance, top_label_count)
      m<-paste(m,collapse="\t")
      write(m,file="/home/ubuntu/logs/log.txt",append=TRUE)
      # }
      
      
      
      ##################################
      #DNA-aspect --- starting at a 
      #genetic location
      ##################################
      in_dir<-"~/data/2016-02-19_splits/"
      window<-1000000
      if(!gene%in%rownames(geneLocations)){stop(safeError(paste(gene,"not found. Please only use human genesymbols (all upper-case letters).")))}
      chr<-sub("^chr","",geneLocations[gene,"chr_name"])
      start<-geneLocations[gene,"start"] - distance
      end<-geneLocations[gene,"end"] + distance
      p1<-floor(start/window)
      p2<-floor(end/window)
      filename1<-paste(in_dir,"split_pheno_",phenotype,"_chr",chr,"_",p1,"_",p1+1,".txt.gz",sep="")
      if(!file.exists(filename1))stop(safeError(paste("Could not find file",filename1)))
      dh<-read.table(filename1,sep="\t",header=T,stringsAsFactors=F)
      
      #in case we are at a window breakpoint
      if(p1!=p2){
        dh1<-dh
        filename2<-paste(in_dir,"split_pheno_",phenotype,"_chr",chr,"_",p2,"_",p2+1,".txt.gz",sep="")
        dh2<-read.table(filename2,sep="\t",header=T,stringsAsFactors=F)
        dh<-rbind(dh1,dh2)
      }
      
      data<-dh[dh[,"BP"]>start & dh[,"BP"]<end,]
      if(nrow(data)==0){stop(safeError(paste("No SNPs found around gene",gene)))}
      
      data[,"-log10(P)"] <- -log10(data[,"P" ])
      rownames(data)<-data[,"SNP"]
      
      return(data)				
      
    }
  })
  
  output$mainPlot <- renderPlot({ 
    # email <- isolate(input$email)
    gene <- isolate(input$gene)
    distance <- isolate(input$distance)
    # p_value_cutoff <- isolate(input$p_value_cutoff)
    top_label_count<-isolate(input$top_label_count)
    phenotype <- isolate(input$phenotype)
    
    data<-get_data()
    if(is.null(data))return(NULL)
    
    
    chr<-sub("^chr","",geneLocations[gene,"chr_name"])
    start<-geneLocations[gene,"start"] - distance
    end<-geneLocations[gene,"end"] + distance
    
    
    ylim<-c(0,max(c(5,data[,"-log10(P)"])))
    xlim<-c(start,end)
    plot(NULL,
         xlim=xlim,
         ylim=ylim,
         xlab=paste("chr",chr),
         ylab="-log10(P)",
         main=paste0(rownames(p)[p[,"pheno_id"]%in%phenotype],"-affecting SNPs around ",gene)
    )
    
    points(
      x=data[,"BP"],
      y=data[,"-log10(P)"],
      pch=19,
      col="#821556"
    )
    
    #highlight top-3 SNPs
    tooClose<-vector()
    tooCloseDist<-4000
    count<-min(c(nrow(data),top_label_count))
    for(i in 1:count){
      snps<-rownames(data)[order(data[,"-log10(P)"],decreasing=T)]
      if(sum(!snps%in%tooClose)==0)break
      snp<-snps[!snps%in%tooClose][1]
      text(x=data[snp,"BP"],y=data[snp,"-log10(P)"],label=snp,adj=0,cex=0.9)
      tooCloseHere<-data[snp,"BP"] - tooCloseDist < data[,"BP"] & data[snp,"BP"] + tooCloseDist > data[,"BP"]
      tooClose<-c(tooClose,rownames(data)[tooCloseHere])
      
    }
    lines(x=c(geneLocations[gene,"start"],geneLocations[gene,"end"]),y=c(0,0),lwd=4,col="black")
    
    
    
  })
  
  output$mainTable <- renderDataTable({ 
    # email <- isolate(input$email)
    
    gene <- isolate(input$gene)
    distance <- isolate(input$distance)
    # p_value_cutoff <- isolate(input$p_value_cutoff)
    top_label_count<-isolate(input$top_label_count)
    phenotype <- isolate(input$phenotype)
    
    data<-get_data()
    if( is.null(data))return(NULL)
    
    
      data<-data[,c("SNP","CHR","BP","P")]
      return(data)
    
    
  })
})



