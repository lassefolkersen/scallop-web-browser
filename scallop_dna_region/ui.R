library("shiny")

source("../uifunctions.R")
initialize('con',TRUE)



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
phenotypes_vector<-c("none",rownames(p))




shinyUI(bootstrapPage(
  head(),
  navigation(),
  titlePanel("SNP-centric view - by a specific gene"),
  beginPage(),
  beginPanel('1/3'),
  HTML("This page can be used to browse the pQTL effects of all SNPs proximal to any gene in the genome:<br><br>"),
  textInput(inputId="email", label = "E-mail", value = ""),
  textInput(inputId="gene", label = "Gene (HGNC) or position (chr:bp)", value = ""),
  checkboxInput("advanced", "Advanced options", value = FALSE),
  conditionalPanel(
    condition = "input.advanced",
    checkboxInput("show_gene_map", "Show gene map", value = FALSE)
  ),
  conditionalPanel(
    condition = "input.advanced",
    selectInput("protein_to_highlight", "Highlight protein", choices = phenotypes_vector,selected="none")
  ),
  conditionalPanel(
    condition = "input.advanced",
    sliderInput("distance","Distance from gene (bp)",min=100000,max=500000,value=200000)
  ),
  conditionalPanel(
    condition = "input.advanced",
    sliderInput("top_label_count","#SNPs in table",min=10,max=300,value=50,step=1)
  ),
  conditionalPanel(
    condition = "input.advanced",
    HTML("Documentation for analysis is available at <u><a href='http://github.com/lassefolkersen/olink-improve'>github</a></u>. All positions are in GRCh37/hg19 coordinates.")
  ),
  actionButton("goButton","Run analysis"),
  endPanel(),
  beginPanel('2/3'),
  
  plotOutput("mainPlot",height = "700px"),
  # conditionalPanel(
    # condition = "input.show_gene_map == TRUE",
    # 
  # ),
  # 
  # conditionalPanel(
    # condition = "input.show_gene_map == FALSE",
    # plotOutput("mainPlot",height = "400px")
  # ),
  
  
  
  
  HTML("<br><br>"),
  dataTableOutput("mainTable"),
  endPanel(),
  
  endPage(),
  footer()
  
)	
)







