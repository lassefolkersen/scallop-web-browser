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
accepted_users<-tolower(read.table("/home/ubuntu/misc/accepted_emails.txt",sep="\t",header=F,stringsAsFactors=F)[,1])

#colouring scheme
#Setting colours
# sum(cl)/1e7 this is 307 (meaning we could divide the genome into 307 1e7 bp chunks. Let's.
cols<-rainbow(308, s = 1, v = 1, start = 0, end = 1, alpha = 1)
names(cols) <- as.character(1:(length(cols)-1))


#gene positions
# load("~/srv/olink-scallop/2014-07-16 gene locations.rdata")
#load these when needed instead


#protein locations
load("~/srv/olink-scallop/scallop_protein_centric/2019-01-24_cvd1_protein_positions.rdata")

#significant hits - loaded at run-time instead
# all_data<-read.table("/home/ubuntu/data/2019-01-24_significant_bits_scallop/2019-01-24_top_hits_pruned.txt.gz",sep="\t",header=T,stringsAsFactors = F)
# all_data<-all_data[nchar(gsub("\\?","",all_data[,"Direction"]))>1,]
# 
# for(protein in unique(all_data[,"protein"])){
#   d1 <-all_data[all_data[,"protein"] %in% protein,]
#   gz1 <- gzfile(paste0("/home/ubuntu/data/2019-01-24_significant_bits_scallop/2019-01-24_",protein,"_top_hits_pruned.txt.gz"), "w",compression = 9)
#   write.table(d1, gz1,sep="\t",col.names=T,row.names=F,quote=F)
#   close(gz1)
# }






shinyServer(function(input, output) {
  
  
  
  #main data gathering function	
  get_data <- reactive({
    if(input$goButton == 0){
      return(NULL)
    }else{
      
      ##################################
      #input-variables, log and register	
      ##################################
      email <- NA #isolate(input$email)
      p_value_cutoff <- isolate(input$p_value_cutoff)
      top_label_count<-isolate(input$top_label_count)
      protein <- isolate(input$protein)
      include_closest_genes <- isolate(input$include_closest_genes)

      #the mail checker and logger - first part is not needed when we don't check mail, but left in anyway
      if(!tolower(email) %in% accepted_users &  FALSE){
        m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email)
        m<-paste(m,collapse="\t")
        write(m,file="/home/ubuntu/logs/illegal_access_log.txt",append=TRUE)
        Sys.sleep(2)
        stop(safeError("In the test-phase non-privileged users are not allowed"))
      }else{
        m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),email,"scallop_protein_centric",protein, p_value_cutoff, top_label_count,include_closest_genes)
        m<-paste(m,collapse="\t")
        write(m,file="/home/ubuntu/logs/log.txt",append=TRUE)
      }
      
      
      ##################################
      #Protein-aspect --- starting with 
      #a specific plasma protein
      ##################################
      
      in_dir <- "~/data/2019-01-24_significant_bits_scallop/"
      file_name<-paste(in_dir,"2019-01-24_",protein,"_top_hits_pruned.txt.gz",sep="")
      
      data<-read.table(file_name, header=T,sep="\t",stringsAsFactors=FALSE)
      data<-data[data[,"logP"] > p_value_cutoff,]
      if(nrow(data)==0)stop(safeError("No SNPs were significant at this p_value_cutoff"))
      data[,"neglogp"] <- data[,"logP"]
      
      
      data[,"trait_chr"]<-protein_pos[protein,"Chr"]
      data[,"trait_pos"]<-protein_pos[protein,"CDS_start"]
      
      
      #calculate the "universal-BP" --- i.e. the BP-sum on chr + all-other-chr before it (=were on the circle to plot)
      data[,"snp_abs_pos"]<-chrLengths[data[,"chr"],"startsum"] + data[,"pos"]
      data[,"trait_abs_pos"]<-chrLengths[data[,"trait_chr"],"startsum"] + data[,"trait_pos"]
      
      
      #This block transforms all the (chr-pos, -logP) positions from cartesian coordinates to polar coordinates (=this is the circularizing step)
      for(i in 1:nrow(data)){
        data[i,"r"] <- max(c(0,log10(data[i,"neglogp"]) / p_range))
        data[i,"theta"] <- (data[i,"snp_abs_pos"] / maxPos) * 2 * pi
        data[i,"x"] = (base+data[i,"r"]) * cos( data[i,"theta"] )
        data[i,"y"] = (base+data[i,"r"]) * sin( data[i,"theta"] )
        data[i,"x0"] = (base) * cos( data[i,"theta"] )
        data[i,"y0"] = (base) * sin( data[i,"theta"] )
      }
      data<-data[!is.na(data[,"x"]),]
      data<-data[order(data[,"snp_abs_pos"]),]
      
      data[,"colours"]<-cols[as.character(round(data[,"snp_abs_pos"]/1e7))]
      return(data)
    }
  })
  
  output$mainPlot <- renderPlot({ 
    p_value_cutoff <- isolate(input$p_value_cutoff)
    top_label_count<-isolate(input$top_label_count)
    protein <- isolate(input$protein)
    include_closest_genes <- isolate(input$include_closest_genes)
    
    #load the gene locations only if necessary
    if(include_closest_genes){
      if(!exists("geneLocations")){
        load("~/srv/olink-scallop/2014-07-16 gene locations.rdata")  
      }
    }
    
    
    data<-get_data()
    if(is.null(data))return(NULL)
    
    #Open up the actual plotting mechanism with an empty canvas
    plot(NULL,xlim=c(-4,4),ylim=c(-4,4),ylab="",xlab="",xaxt="n",yaxt="n")
    
    #Drawing the lines for the manhattan plot (using polar coordinates from above)
    for(i in 1:nrow(data)){
      lines(
        x=c(data[i,"x"],data[i,"x0"]),
        y=c(data[i,"y"],data[i,"y0"]),
        col=data[i,"colours"])
    }
    
    
    #drawing the frame-work circles (we want a P 0, 10^-8 and a 'current-cutoff-circle)
    ch<-data.frame(offset=c(0,log10(p_value_cutoff),log10(7.30103)),lty=c(1,2,1),col=c("black","grey70","grey50"),stringsAsFactors=F)
    for(k in 1:nrow(ch)){
      res<-50
      b<-seq(0,2*pi,length.out=res)
      for(i in 1:(res-1)){
        x1 = (ch[k,"offset"] + base) * cos( b[i] )
        x2 = (ch[k,"offset"] + base) * cos( b[i+1] )
        y1 = (ch[k,"offset"] + base) * sin( b[i] )
        y2 = (ch[k,"offset"] + base) * sin( b[i+1] )
        lines(x=c(x1,x2),y=c(y1,y2),lty=ch[k,"lty"],col=ch[k,"col"])
      }
    }
    
    #drawing the internal lines
    for(j in 1:nrow(data)){
      orig<-data[j,"snp_abs_pos"]
      dest<-data[j,"trait_abs_pos"]
      # colour<-rainbow(nrow(data), s = 1, v = 1, start = 0, end = 1, alpha = 1)[sum(data[,"snp_abs_pos"]<orig)+1]
      colour <- data[j,"colours"]
      x1 = (base) * cos( 2*pi*(orig/maxPos) )
      y1 = (base) * sin( 2*pi*(orig/maxPos) )
      x2 = (base) * cos( 2*pi*(dest/maxPos) )
      y2 = (base) * sin( 2*pi*(dest/maxPos) )
      lines(x=c(x1,x2),y=c(y1,y2),col=colour,lwd=1)
      # print(paste(x1,y1,x2,y2))
    }
    
    
    
    #drawing the chromosome sizes
    for(i in 1:nrow(chrLengths)){
      start<-chrLengths[i,"startsum"]
      # end<-chrLengths[i,"cumsum"]
      x1 = (base) * cos( 2*pi*(start/maxPos) )
      x2 = (base+0.1) * cos( 2*pi*(start/maxPos) )
      y1 = (base) * sin( 2*pi*(start/maxPos)  )
      y2 = (base+0.1) * sin( 2*pi*(start/maxPos)  )
      
      lines(x=c(x1,x2),y=c(y1,y2),lty=1,col="black")
      text(x1,y1,label=rownames(chrLengths)[i],cex=0.5,adj=-0.2,col="grey50")
    }
    
    
    
    #labelling top-X hits - with gene-labels if requested
    if(top_label_count > 0){
      tooClose<-vector()
      tooCloseDist<-20000000
      count<-min(c(nrow(data),top_label_count))
      for(i in 1:count){
        MarkerNames<-data[order(data[,"neglogp"],decreasing=T),"MarkerName"]
        if(sum(!MarkerNames%in%tooClose)==0)break
        MarkerName<-MarkerNames[!MarkerNames%in%tooClose][1]
        x<-data[data[,"MarkerName"]%in%MarkerName,"x"]
        y<-data[data[,"MarkerName"]%in%MarkerName,"y"]
        text(x=x,y=y,label=MarkerName,adj=0,cex=0.8)
        pos<-data[data[,"MarkerName"]%in%MarkerName,"snp_abs_pos"]
        tooCloseHere<-data[abs(data[,"snp_abs_pos"] - pos) <tooCloseDist,"MarkerName"]
        tooClose<-c(tooClose,tooCloseHere)
        
        
        #optionally also include closest genes
        if(include_closest_genes){
          chr <- paste0("chr",sub(":.+$","",MarkerName ))
          pos <- as.numeric(sub(":.+$","",sub("[0-9]+:","",MarkerName)))
          g1 <- geneLocations[geneLocations[,"chr_name"]%in%chr,]
          g1[,"start_to_pos"] <- abs(g1[,"start"] - pos)
          g1[,"end_to_pos"] <- abs(g1[,"end"] - pos)
          g1[,"distance"] <- apply(g1[,c("start_to_pos","end_to_pos")],1,min)
          g1[g1[,"start"] < pos & g1[,"end"] > pos,"distance"] <- 0 #set distance to 0 if within gene
          g1<-g1[order(g1[,"distance"]),]
          gene_label_count <- 4
          gene_label <- paste(rownames(g1)[1:gene_label_count],collapse=", ")
          text(x=x,y=y-0.15,label=gene_label,adj=0,cex=0.8)
        }
      }
    }
    
    
    #labelling the trait
    # trait<-protein
    dest<-data[i,"trait_abs_pos"]
    x3 = (base) * cos( 2*pi*(dest/maxPos) )
    y3 = (base) * sin( 2*pi*(dest/maxPos) )
    text(x3,y3,protein,adj=0,font=4)
    
    
    
    #drawing a legend
    cols<-cols[seq(1,length(cols),by=3)]
    min<-3
    max<-4
    scale = (length(cols)-1)/(max-min)
    # plot(NULL,xlim=c(-4,4),ylim=c(-4,4),ylab="",xlab="",xaxt="n",yaxt="n")
    for (i in 1:(length(cols)-1)) {
      y = (i-1)/scale + min
      rect(-4,y,-3.5,y+1/scale, col=cols[i], border=NA)
    }
    text(x=-3.4,y=4.05,"Chr",cex=0.7)
    
    for(chr in c(1,2,3,5,7,9,11,13,15,18,22)){
      f <- chrLengths[chr,"startsum"] / chrLengths[nrow(chrLengths),"endsum"]
      text(x=-3.4, y=f+min * (max-min),chr,cex=0.5)
    }
    legend(
      x=-3,y=4,
      legend=c("P=5e-8",paste("P=1e-",signif(p_value_cutoff,2),sep="")),
      lty=c(1,2),
      lwd=c(1,1),
      col=c("grey50","grey70"),
      cex=0.9,bty="n"
    )
  })
  
  output$mainTable <- renderDataTable({ 
    # email <- isolate(input$email)
    # gene <- isolate(input$gene)
    p_value_cutoff <- isolate(input$p_value_cutoff)
    top_label_count<-isolate(input$top_label_count)
    protein <- isolate(input$protein)
    
    data<-get_data()
    if( is.null(data))return(NULL)
    d<-data[data[,"logP"]>p_value_cutoff,]
    d <- d[order(d[,"logP"],decreasing=T),]
    d<-d[,c("MarkerName","Freq1","Effect","P.value.character","logP","Direction","TotalSampleSize")]
    colnames(d) <- c("MarkerName","Frequency","Effect","P-value","logP","Direction","SampleSize")
    
    
    #round sample size
    d[,"SampleSize"] <- round(d[,"SampleSize"])
    
    #shorten other numbers
    d[,"Frequency"] <- signif(d[,"Frequency"],2)
    d[,"Effect"] <- signif(d[,"Effect"],3)
    d[,"P-value"] <- signif(d[,"P-value"],2)
    d[,"logP"] <- signif(d[,"logP"],3)
    
    return(d)
    
    
  })
  
  
  
  
  
  
  #The long methods section written after the SNPs-table
  output$text_3 <- renderText({ 
    
    if(input$goButton == 0){
      methodsToReturn<-""
      
      
    }else{
      o<-get_data()
      methodsToReturn<-o[["methodsToReturn"]]
    }
    return(methodsToReturn)
  })
  
  
  
  #The methods text box
  output$text_1 <- renderText({ 
    
    
    if(input$goButton == 0){
      methodsToReturn
    }else{
      methodsToReturn<-paste0("<small><br><b>Methods</b><br>The circular manhattan plot shows the association strength of all observed pQTL effects, both cis- and trans-effects. The highlighted protein-name shows the location of the protein-encoding gene and any proximal peaks are cis-effects. Trans-effects are indicated with connecting lines towards the protein-encoding gene. The fully drawn grey circle indicates a significance level of P=5e-8. The dashed line indicates the customizable treshold, defaulting to P=1e-6. The markerNames of the 3 strongest associations are labelled on plot along with the names of their 4 most proximal genes. These numbers can be altered in advanced options. 

                              <br><br>The table shows extended information for each of the pQTL shown in the circular manhattan plot. The markername in the format of chr:pos:A1_A2, where A1/A2 is by alphabetical sorting.  The frequency-column shows the observed A1 alelle frequency. The effect-column shows the effect of the A1 allele in standardized units, meaning 1 equals one SD of protein level change. P-value and logP indicates significance, the logP-column is included because values less than 1e-300 often results in failure of analysis software (including R) and logP-units are recommended. Direction is given for the participating SCALLOP studies, alphabetically: EpiHealth, Estonian_Biobank, IMPROVE, INTERVAL, LifeLinesDeep, MPP_RES, NSPHS, ORCADES, PIVUS, STABILITY, STANLEY_lah1, STANLEY_swe6, ULSAM, VIS. A question mark indicates that either the protein or the SNP was not available in this study. Sample size is the effective sample size for the indicated pQTL. Only the strongest effect in each 50kb window is shown in the table.<br></small>")
    }
  })
  
  
  
})



