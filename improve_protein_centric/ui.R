source("../uifunctions.R")
initialize('con',TRUE)



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




shinyUI(bootstrapPage(
  head(),
  navigation(),
  titlePanel("IMPROVE: Protein-centric browser"),
  beginPage(),
  beginPanel('1/3'),
  HTML("<i>Note: This module queries the data presented in our <u><a href='http://dx.plos.org/10.1371/journal.pgen.1006706'>2017-paper on pQTL</a></u>. They are conserved here for consistency-reasons. You may however wish to start analysis with the larger SCALLOP analysis in other modules.</i><br><br>This page can be used to browse the significant pQTL effects of SNPs anywhere in the genome:<br><br>"),
  # textInput(inputId="email", label = "E-mail", value = ""),
  selectInput("phenotype", "Protein", choices = phenotypes_vector),
  checkboxInput("advanced", "Advanced options", value = FALSE),
  conditionalPanel(
    condition = "input.advanced",
    sliderInput("p_value_cutoff","P-value cutoff (log)",min=4,max=12,value=6,step=0.1)
  ),
  conditionalPanel(
    condition = "input.advanced",
    sliderInput("top_label_count","#SNPs to label",min=3,max=30,value=3,step=1)
  ),
  conditionalPanel(
    condition = "input.advanced",
    HTML("Documentation for analysis is available at <u><a href='http://github.com/lassefolkersen/scallop-web-browser'>github</a></u>.")
  ),
  actionButton("goButton","Run analysis"),
  endPanel(),
  beginPanel('2/3'),
  plotOutput("mainPlot",height="800px",width="800px"),
  
  dataTableOutput("mainTable"),
  endPanel(),
  
  endPage(),
  footer()
  
)	
)







