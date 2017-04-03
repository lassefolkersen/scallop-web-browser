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



shinyUI(bootstrapPage(
  head(),
  navigation(),
  titlePanel("Per-SNP summary downloads"),
  beginPage(),
  beginPanel('1/3'),
  HTML("This page can be used to download per-SNP summary data:<br><br>"),
  # textInput(inputId="email", label = "E-mail", value = ""),
  # textInput(inputId="gene", label = "Gene", value = ""),
  selectInput("phenotype", "Protein", choices = phenotypes_vector),
  checkboxInput("advanced", "Advanced options", value = FALSE),
  conditionalPanel(
    condition = "input.advanced",
    HTML("Documentation for analysis is available at <u><a href='http://github.com/lassefolkersen/olink-improve'>github</a></u>.")
  ),
  actionButton("goButton","Run analysis"),
  endPanel(),
  beginPanel('2/3'),
  
  htmlOutput("explanatoryText"),
  dataTableOutput("mainTable"),
  endPanel(),
  
  endPage(),
  footer()
  
)	
)







