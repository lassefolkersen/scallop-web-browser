library("shiny")
library("openxlsx")
library("jsonlite")
library("igraph")
library("visNetwork")


#data path
input_path <- "/home/ubuntu/data/2019-01-31_pathways/"

#permissions
accepted_users<-tolower(read.table("/home/ubuntu/misc/accepted_emails.txt",sep="\t",header=F,stringsAsFactors=F)[,1])


# Define server logic for a template
shinyServer(function(input, output) {
  
  
  #Get the pre-calculated genetic data for this user
  get_data <- reactive({
    if(input$goButton == 0){
      return(NULL)
    }
    protein<-isolate(input$protein)
    email <- NA #isolate(input$email)
    
    #the mail checker and logger - first part is not needed when we don't check mail, but left in anyway
    if(!tolower(email) %in% accepted_users &  FALSE){
      m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"plot",email)
      m<-paste(m,collapse="\t")
      write(m,file="/home/ubuntu/logs/illegal_access_log.txt",append=TRUE)
      Sys.sleep(2)
      stop(safeError("In the test-phase non-privileged users are not allowed"))
    }else{
      m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),email,"scallop_pathways",protein)
      m<-paste(m,collapse="\t")
      write(m,file="/home/ubuntu/logs/log.txt",append=TRUE)
    }
    
    
    
    path_to_read<-paste0(input_path,"2019-01-31_",protein,"_pathway.txt")
    if(!file.exists(path_to_read))stop(safeError(paste("Didn't find the network for protein",protein,"likely because no paths were possible.")))
    edges <- read.table(path_to_read,stringsAsFactors = F,header=T)    
    
    e<-graph_from_edgelist(as.matrix(edges[,c("source","target")]), directed = TRUE)
    
    #set edge types and strength
    e<-set_edge_attr(e,"value",value=edges[,"strength"])
    e<-set_edge_attr(e,"type",value=edges[,"type"])
    
    #set node type
    V(e)$type <- "gene"
    V(e)$type[grep("^[0-9r]",V(e)$name)]<-"snp"
    
    
    
    
    return(e)
  })
  
  
  
  
  
  #create the network
  output$plot1 <- renderVisNetwork({
    
    
    e1 <- get_data()
    
    if(is.null(e1))return(NULL)
    
    
    #set shape
    V(e1)$shape <- "circle"
    
    #big and grey for genes/proteins
    V(e1)$label.cex <- 1
    V(e1)$color<-"#CCCCCC"
    
    #big and light for snps
    V(e1)$label.cex[V(e1)$type=="snp"] <- 0.7
    V(e1)$color[V(e1)$type=="snp"]<-"#F5F5F5"
    
    # #dash the non-significant
    E(e1)$dashes <- FALSE
    E(e1)$dashes[E(e1)$value < -log10(0.05)] <- TRUE
    
    
    # E(e1)$dash <- E(e1)$dashes
    # E(e1)$lty <- 2
    
    # For legend
    ledges <- data.frame(color = c("gray","gray","gray"),
                         label = c("P>0.05", "P<0.05","P<0.005"), 
                         width=c(0.5,1,3),
                         dashes=c(T,F,F),
                         arrows=c("none","none","none"))
    
    
    #then create the visNetwork from this igraph object    
    a<-visIgraph(e1)%>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = FALSE) %>%
      visIgraphLayout(layout = "layout_as_tree",flip.y=TRUE, mode="in" ) %>%
      visLegend(addEdges = ledges)
    
    
    
    return(a)
  })
  
  
  
  
  
  #The methods text box
  output$text_1 <- renderText({ 
    if(input$goButton == 0){
      methodsToReturn<-""
    }else{
      methodsToReturn<-paste0("<small><br><b>Methods</b><br>
This plot shows an overview of all starting points for observed pQTL effects for a given protein-encoding gene (bottom dark-grey node). Any cis-pQTL effects are shown immediately above this (white nodes). If any trans-pQTL effects were observed, the graph provides a calculation of the most likely molecular path (dark nodes) from the trans-pQTL variants (white nodes) to the protein-encoding gene (bottom node). The paths are calculated using the <u><a href='https://www.ncbi.nlm.nih.gov/pubmed/?term=27892958'>protein-protein interaction (PPI)</a></u> data from inBio Map™ using <u><a href='https://www.intomics.com/inbio/map.html'>Intomics pipelines</a></u>. The overall network had 13033 nodes and 147882 edges that were used as input for the pathway calculation. We define a <i>cis-gene intermediary</i> as the immediate gene(s) through which a trans-pQTL variant exerts its effect, and each analysis consequently tries to identify the shortest path between any of the cis-gene intermediaries to the plasma protein. The width of each edge indicates the likelihood of a path arising by chance, as calculated by permutation sampling, using randomly sampled networks with a conserved node degree distribution.<br> ")
    }
  })
  
  
  
  
})








