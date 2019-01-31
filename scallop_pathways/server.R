library("shiny")
library("openxlsx")
library("jsonlite")
library("igraph")
library("visNetwork")


#data path
input_path <- "/home/ubuntu/data/2019-01-31_pathways/"

# Define server logic for a template
shinyServer(function(input, output) {
  
  
  #Get the pre-calculated genetic data for this user
  get_data <- reactive({
    if(input$goButton == 0){
      return(NULL)
    }
    protein<-isolate(input$protein)

    
    path_to_read<-paste0(input_path,"2019-01-31_",protein,"_pathway.txt")
    if(!file.exists(path_to_read))stop(safeError(paste("Didn't find the network for protein",protein,"likely because no paths were possible.")))
    edges <- read.table(path_to_read,stringsAsFactors = F,header=T)    
    
    e<-graph_from_edgelist(as.matrix(edges[,c("source","target")]), directed = TRUE)
    e<-set_edge_attr(e,"strength",value=edges[,"strength"])
    e<-set_edge_attr(e,"type",value=edges[,"type"])
    
    
    
    
    return(e)
  })
  
  
  
  
  
  #create the network
  network_proxy_select <- reactive({

    
    e1 <- get_data()
    
    #get the tooltip
    # V(e1)$title <- V(e1)$niceName #to get the tooltip
    
    #get the sizes - large for center and close proximity - smaller for farther out. This has to be set within the igraph object
    # dr<-range(V(e1)$distance)
    # V(e1)$norm_dist <- (dr[2]-V(e1)$distance)/(dr[2]-dr[1])
    # V(e1)$size <- (V(e1)$norm_dist + 0.5)*200
    
    
    #get the shape etc.     
    # show_small <- which(V(e1)$distance == max(V(e1)$distance))
    
    V(e1)$shape <- "circle"
    V(e1)$label.cex <- 1.5
    # V(e1)$label.cex[show_small] <- 0.4
    
    V(e1)$color<-"#F5F5F5"

    
    
    #then create the visNetwork from this igraph object    
    a<-visIgraph(e1)%>%
      visInteraction(tooltipStyle = 'position: fixed;visibility:hidden;padding: 5px;white-space: wrap;font-family: arial;font-size:18px;font-color:black;') %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = FALSE) %>%
      visEvents(select = "function(nodes) {
            Shiny.onInputChange('focus_node', nodes.nodes);
            ;}") %>%
      visIgraphLayout(layout = "layout_as_tree",flip.y=TRUE, mode="in" )
      

    
    return(a)
  })
  
  
  
  # 
  # #function to get ID-hover working (I think)
  # observe({
  #   nodes_selection <- input$selnodes
  #   visNetworkProxy("network_proxy_select") %>%
  #     visSelectNodes(id = nodes_selection) 
  # })
  # 
  
  
  #getting a table of hits in a bubble
  output$table1 <- renderTable({ 
    uniqueID <- gsub(" ","",input$uniqueID)
    o<-get_colour_code()
    
    if(is.null(o) | input$goButton == 0 | is.null(input$focus_node)){
      return(NULL)
    }
    focus_node <- input$focus_node
    
    #write the score to the log file
    log_function<-function(uniqueID,focus_node){
      user_log_file<-paste("/home/ubuntu/data/",uniqueID,"/user_log_file.txt",sep="")
      m<-c(format(Sys.time(),"%Y-%m-%d-%H-%M-%S"),"scallop_pathways",uniqueID,focus_node)
      m<-paste(m,collapse="\t")
      if(file.exists(user_log_file)){
        write(m,file=user_log_file,append=TRUE)
      }else{
        write(m,file=user_log_file,append=FALSE)
      }
    }
    try(log_function(uniqueID,focus_node))
    
    
    
    o1<-o[["scores"]]
    
    o1<- o1[o1[,"ICD_code"] %in% focus_node,]
    
    if(nrow(o1)==0)return(NULL)
    
    
    #rename genetics study (only in AllDiseases module)
    w1 <- which(o1[,"module"]%in%"AllDiseases")
    n <- o1[w1,"study_code"]
    n <- sub("([0-9]+)$","(PMID \\1)",gsub("_"," ",n)) #remove underscore and add PMID
    for(i in 1:length(n)){
      n[i] <- paste(toupper(substring(n[i], 1,1)), substring(n[i], 2),sep="", collapse=" ") #capitalize first letter  
    }
    o1[w1,"study_code"] <- n #re-insert
    
    

    #insert disease name
    o1[,"disease"]<- paste0(V(e)[o1[,"ICD_code"]]$niceName," (",o1[,"ICD_code"],")")
    
    #remove disease code if found in "Feeling fine" (doesn't make sense there)
    o1[o1[,"disease"]%in%"Feeling fine (Feeling fine)","disease"] <- ""
    
    
    #remove Z-score if found in BRCA or rareDiseases (doesn't make sense there)
    w2 <- which(o1[,"module"]%in%c("rareDiseases","BRCA"))
    o1[w2,"score"] <- "+"
    
    # Translate module names
    niceNames <- c("GWAS calculator","Drug response","UK-biobank","Rare Diseases","BRCA")
    names(niceNames) <- c("AllDiseases","drugResponse","ukbiobank","rareDiseases","BRCA")
    o1[,"module"] <- niceNames[o1[,"module"]]
    
    
    #rename and re-order columns
    select <- c("disease","study_code","module","score")
    names(select) <- c("Disease (code)","Genetic study","Further details in this Module","Z-score")
    if(!all(select%in%colnames(o1)))stop("Not all columns found")
    
    
    

    
    
    
    
    
    o1 <- o1[,select]
    colnames(o1)<-names(select)

    
    return(o1)
  })
  
  #show the network with focus based on input
  output$plot1 <- renderVisNetwork({
    network <- network_proxy_select()
    return(network)
  })
  
  
  output$text_1 <- renderText({
    if(input$goButton == 0){
      m <- "Before"

    }else{
      m<-"After"
    }
    return(m)
    
  })
  
  

})








