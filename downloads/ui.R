source("../uifunctions.R")
initialize('con',TRUE)



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


# 
# ui <- fluidPage(
#   titlePanel("Improve-OLINK browser"),
#   sidebarLayout(
#     sidebarPanel(
#       
#       radioButtons("type", "Focus", 
#                    c("Cis-effects"="dna",
#                      "Trans-effects" = "protein",
#                      "Data download"="download"), 
#                    selected = "dna", inline = TRUE),
#       
#       textInput(inputId="email", label = "E-mail", value = ""),
#       conditionalPanel(
#         condition = "input.type == 'dna'",
#         textInput(inputId="gene", label = "Gene", value = "")
#       ),
#       selectInput("phenotype", "Protein", choices = phenotypes_vector),
#       checkboxInput("advanced", "Advanced options", value = FALSE),
#       conditionalPanel(
#         condition = "input.type == 'dna'  & input.advanced",
#         sliderInput("distance","Distance from gene (bp)",min=100000,max=500000,value=200000)
#       ),
#       conditionalPanel(
#         condition = "input.type == 'protein'  & input.advanced",
#         sliderInput("p_value_cutoff","P-value cutoff (log)",min=4,max=12,value=6,step=0.1)
#       ),
#       conditionalPanel(
#         condition = "input.advanced & (input.type == 'protein' | input.type == 'dna')",
#         sliderInput("top_label_count","#SNPs to label",min=3,max=30,value=3,step=1)
#       ),
#       conditionalPanel(
#         condition = "input.advanced",
#         HTML("In the review phase access to this browser is available only through the reviewer-provided pass-email<br><br>Documentation for analysis is available at <u><a href='http://github.com/lassefolkersen/olink-improve'>github</a></u>..")
#       ),
#       
#       
#       
#       actionButton("goButton","Start"),
#       width=4
#     ),
#     mainPanel(
#       conditionalPanel(
#         condition = "input.type == 'download'",
#         htmlOutput("explanatoryText")
#       ),
#       conditionalPanel(
#         condition = "input.type == 'protein' | input.type == 'dna'",
#         plotOutput("mainPlot",width = "800px", height = "800px")
#       ),
#       dataTableOutput("mainTable")
#       
#     )
#   )
# )
# 
# shinyApp(ui = ui, server = server)
# 


shinyUI(bootstrapPage(
  head(),
  navigation(),
  titlePanel("Per-SNP summary downloads"),
  beginPage(),
  beginPanel('1/3'),
  HTML("This page can be used to download per-SNP summary data:<br><br>"),
  textInput(inputId="email", label = "E-mail", value = ""),
  textInput(inputId="gene", label = "Gene", value = ""),
  selectInput("phenotype", "Protein", choices = phenotypes_vector),
  checkboxInput("advanced", "Advanced options", value = FALSE),
  conditionalPanel(
    condition = "input.advanced",
    HTML("In the review phase access to this browser is available only through the reviewer-provided pass-email<br><br>Documentation for analysis is available at <u><a href='http://github.com/lassefolkersen/olink-improve'>github</a></u>..")
  ),
  actionButton("goButton","Run analysis"),
  endPanel(),
  beginPanel('2/3'),
  plotOutput("mainPlot"),
  
  dataTableOutput("mainTable"),
  endPanel(),
  
  endPage(),
  footer()
  
)	
)







