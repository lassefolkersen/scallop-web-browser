library("shiny")


#permissions
accepted_users<-tolower(read.table("/home/ubuntu/misc/accepted_emails.txt",sep="\t",header=F,stringsAsFactors=F)[,1])



#colours (could speed by saving as fixed)
# library(RColorBrewer)
# colours <- c(brewer.pal(9,"Set1"), brewer.pal(12,"Set3") ,brewer.pal(8,"Pastel2"),brewer.pal(8,"Pastel1"),brewer.pal(8,"Pastel2"),brewer.pal(8,"Pastel1"),brewer.pal(8,"Pastel2"),brewer.pal(8,"Pastel1"),brewer.pal(8,"Pastel2"),brewer.pal(8,"Pastel1") ,brewer.pal(8,"Pastel2"),brewer.pal(8,"Pastel1") )
colours<-c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999','#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F','#B3E2CD','#FDCDAC','#CBD5E8','#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC','#FBB4AE','#B3CDE3','#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC','#B3E2CD','#FDCDAC','#CBD5E8','#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC','#FBB4AE','#B3CDE3','#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC','#B3E2CD','#FDCDAC','#CBD5E8','#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC','#FBB4AE','#B3CDE3','#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC','#B3E2CD','#FDCDAC','#CBD5E8','#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC','#FBB4AE','#B3CDE3','#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC','#B3E2CD','#FDCDAC','#CBD5E8','#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC','#FBB4AE','#B3CDE3','#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC')


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
phenotypes_vector<-c("none",rownames(p))




proteins<-c('ADM','AGRP','Beta-NGF','CA-125','CASP-8','CCL20','CCL3','CCL4','CD40','CD40-L','CHI3L1','CSF-1','CSTB','CTSD','CTSL1','CX3CL1','CXCL1','CXCL16','CXCL6','Dkk-1','ECP','EGF','EN-RAGE','ESM-1','FABP4','FAS','FGF-23','FS','GAL','Gal-3','GDF-15','GH','HB-EGF','HGF','hK11','HSP_27','IL-18','IL-1ra','IL-27','IL-6','IL-6RA','IL-8','IL16','ITGB1BP2','KIM-1','KLK6','LEP','LOX-1','mAmP','MB','MCP-1','MMP-1','MMP-10','MMP-12','MMP-3','MMP-7','MPO','NEMO','NT-pro_BNP','OPG','PAPPA','PAR-1','PDGF_subunit_B','PECAM-1','PlGF','PSGL-1','PTX3','RAGE','REN','RETN','SCF','SELE','SIRT2','SPON1','ST2','t-PA','TF','TIE2','TM','TNF-R1','TNF-R2','TNFSF14','TRAIL','TRAIL-R2','TRANCE','U-PAR','VEGF-A','VEGF-D')



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
      gene <- isolate(input$gene)
      distance <- isolate(input$distance)
      top_label_count<-isolate(input$top_label_count)
      protein_to_highlight <- isolate(input$protein_to_highlight)

      
            
      if(!tolower(email) %in% accepted_users ){
        m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email)
        m<-paste(m,collapse="\t")
        write(m,file="/home/ubuntu/logs/illegal_access_log.txt",append=TRUE)
        Sys.sleep(2)
        stop(safeError("In the test-phase non-privileged users are not allowed"))
      }else{
        
        m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),NA,"scallop_regional",protein_to_highlight, gene, distance, top_label_count)
        m<-paste(m,collapse="\t")
        write(m,file="/home/ubuntu/logs/log.txt",append=TRUE)
      }
      
      
      
      ##################################
      #DNA-aspect --- starting at a 
      #genetic location
      ##################################
      data_dir<-"~/data/2019-01-23_regional/"
      window<-1000000
      distance <- distance * 1.5
      
      #position or gene
      if(length(grep(":",gene))>0){
        gene<-sub("^chr","",gene)
        gene<-gsub(" ","",gene)
        s <- strsplit(gene,":")[[1]]
        if(length(s)!=2)stop(safeError("If giving position as e.g. chr4:43254 - there must exactly one :"))
        chr <- as.numeric(s[1])
        end <- start <- as.numeric(s[2])
        if(is.na(chr)|is.na(end)| is.na(start))stop(safeError(paste("couldn't recognize", gene, "as a chr5:423432 position. Make sure there's no non-numeric characters")))

      }else{
        if(!gene%in%rownames(geneLocations)){stop(safeError(paste(gene,"not found. Please only use human genesymbols (all upper-case letters).")))}
        chr<-sub("^chr","",geneLocations[gene,"chr_name"])
        start<-geneLocations[gene,"start"] - distance
        end<-geneLocations[gene,"end"] + distance
      }

      p1<-floor(start/window)
      p2<-floor(end/window)
      

      
      
      dh <- list()
      for(protein in proteins){
        filename1<-paste(data_dir,"2019-01-22_",protein,"_",chr,"_",p1,"_region.txt",sep="")
        if(!file.exists(filename1))stop(safeError(paste("Could not find file",filename1)))
        if(!file.exists(filename1))next
        dh[[protein]]<-read.table(filename1,stringsAsFactors = F, header=T)
      }
      d<-do.call(rbind,dh)      
      
      
      
      
      
      
      
      # #in case we are at a window breakpoint
      # if(p1!=p2){
      #   print(paste("breakpoint",p1,p2))
      #   d1 <- d
      #   dh2 <- list()
      #   for(protein in proteins){
      #     filename1<-paste(data_dir,"2019-01-22_",protein,"_",chr,"_",p2,"_region.txt",sep="")
      #     if(!file.exists(filename1))stop(safeError(paste("Could not find file",filename1)))
      #     if(!file.exists(filename1))next
      #     dh2[[protein]]<-read.table(filename1,stringsAsFactors = F, header=T)
      #     
      #   }
      #   d2<-do.call(rbind,dh2)
      #   d<-rbind(d1,d2)
      # }
      #some duplicates because of double load in
      # d <- d[!duplicated(apply(d[,c("MarkerName","protein")],1,paste,collapse="_")),]
      
      
      
      d<-d[d[,"pos"]>start & d[,"pos"]<end,]
      if(nrow(d)==0){stop(safeError(paste("No SNPs found around gene",gene)))}
      
      
      
      
      return(d)				
      
    }
  })
  
  output$mainPlot <- renderPlot({ 
    gene <- isolate(input$gene)
    distance <- isolate(input$distance)
    top_label_count<-isolate(input$top_label_count)
    protein_to_highlight <- isolate(input$protein_to_highlight)
    
    d<-get_data()
    if(is.null(d) || nrow(d)==0){
      print("no data ready")
      return(NULL)
    }
    
    #calculate in mb (because X-axis becomes nicer then)
    d[,"pos_mb"] <- d[,"pos"] / 1000000
    d<-d[order(d[,"logP"],decreasing=T),]
    
    #then get colours - strongest colours for strongest P-values
    colours<-colours[c(1:3,5:length(colours))] #remove colour #4, it's almost the same as #2
    if(length(unique(d[,"protein"])) > length(colours)){stop(safeError("add more colours"))}
    colours<-colours[1:(length(unique(d[,"protein"])))]
    names(colours) <- unique(d[,"protein"])
    
    #set tuning-parameters
    cutoff <- 5e-8
    hit_per_kb <- 50
    artifical_max <- 50
    artifical_min <- 20
    placement <- "topleft"
    
    #auto-extract n_legend (=all stronger then cutoff)
    n_legend <- length(unique(d[d[,"logP"] > -log10(cutoff),"protein"]))

    #set ylim
    ylim <- c(0, max(d[,"logP"],na.rm=T))
    if(ylim[2]>artifical_max){ylim[2]<-artifical_max  }
    if(ylim[2]<artifical_min){ylim[2]<-artifical_min  }
    
    
    #set xlim and chr
    start<-geneLocations[gene,"start"] - distance
    end<-geneLocations[gene,"end"] + distance
    xlim <- range(start/ 1000000,end/ 1000000)
    xlab <- paste0("Chr ",sub(":.+$","",d[1,"MarkerName"])," (MB)")
    
    #start plot
    plot(NULL,ylim=ylim,xlim=xlim,xlab=xlab,ylab="-log10(P)")
    
    #prepare for saving max P-values (for legend)
    max_p_values <- vector()
    
    for(protein in rev(names(colours))){
      #extract relevant protein data
      d1<-d[d[,"protein"] %in% protein,]
      
      
      #extract max P-value
      max_p_value <- d1[1,"P.value.character"]
      names(max_p_value) <- protein
      max_p_values <- c(max_p_values, max_p_value)
      
      #skip if too low
      if(max(d1[,"logP"])< 4){
        if(protein_to_highlight == protein){
          stop(safeError(paste("Highlighted protein",protein_to_highlight,"not shown, because association strength too low")))
        }
        next
      }
      
      #remove any but strongest hit per n kp
      # keep <- which(!duplicated(round(d1[,"pos_mb"] / (hit_per_kb/1000)   )))
      # d2<-d1[keep,]
      d2 <- d1
      d2<-d2[order(d2[,"pos_mb"]),]
      
      
      #plot top-line
      #set lwd 2 if within the top hits in legend and 10/black if highlighted
      if(protein_to_highlight == protein){
        lines(d2[,"pos_mb"],d2[,"logP"],col="black",lwd=10)
      }else{
        if(which(names(colours)%in%protein) <= n_legend){
          lwd <- 2
        }else{
          lwd <- 1
        }
        lines(d2[,"pos_mb"],d2[,"logP"],col=colours[protein],lwd=lwd)
      }
      
      
      
      
      #plot filling: first do nice transparent colour, then polygonize it
      solid_col <- col2rgb(colours[protein])
      transparent_col <-  rgb(solid_col[1],solid_col[2],solid_col[3],255*0.3,maxColorValue=256)
      x_lines <- c(min(d2[,"pos_mb"]),min(d2[,"pos_mb"]),d2[,"pos_mb"],max(d2[,"pos_mb"]),max(d2[,"pos_mb"]),min(d2[,"pos_mb"]))
      y_lines <- c(0,d2[1,"logP"],d2[,"logP"],d2[nrow(d2),"logP"],0,0)
      polygon(x=x_lines, y = y_lines, density = NULL, angle = 45,border = NA, col = transparent_col)
      
    }
    
    
    
    
    legend(
      placement,
      legend=paste0(
        names(colours)[1:n_legend],
        " (P<",
        sub("\\.[0-9]+e-","e-",max_p_values[names(colours)[1:n_legend]]),
        ")"
      ),
      col=colours[1:n_legend],
      lwd=2,
      cex=0.7)
    
    
    abline(h=-log10(5e-8),lty=2,lwd=2)
    
    
    
  })
  
  output$mainTable <- renderDataTable({

    gene <- isolate(input$gene)
    distance <- isolate(input$distance)
    # p_value_cutoff <- isolate(input$p_value_cutoff)
    top_label_count<-isolate(input$top_label_count)
    protein_to_highlight <- isolate(input$protein_to_highlight)

    d<-get_data()
    if(is.null(d) || nrow(d)==0){
      print("no data ready")
      return(NULL)
    }
    d[,"Direction"] <- d[,"chr"]<- d[,"pos_mb"] <- d[,"P.value.character"] <- d[,"logP"] <- NULL
    
    #order by P-value
    d <- d[order(d[,"P.value"]),]
    
    #only get top-list
    if(top_label_count>nrow(d))top_label_count<-nrow(d)
    d<-d[1:top_label_count,]
    
    
    return(d) 
     
  })
})



