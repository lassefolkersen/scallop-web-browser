library("shiny")

source("../uifunctions.R")
initialize('con',TRUE)



#gene positions
load("~/srv/olink-scallop/2014-07-16 gene locations.rdata")
proteins<-c('ADM','AGRP','Beta-NGF','CA-125','CASP-8','CCL20','CCL3','CCL4','CD40','CD40-L','CHI3L1','CSF-1','CSTB','CTSD','CTSL1','CX3CL1','CXCL1','CXCL16','CXCL6','Dkk-1','ECP','EGF','EN-RAGE','ESM-1','FABP4','FAS','FGF-23','FS','GAL','Gal-3','GDF-15','GH','HB-EGF','HGF','hK11','HSP_27','IL-18','IL-1ra','IL-27','IL-6','IL-6RA','IL-8','IL16','ITGB1BP2','KIM-1','KLK6','LEP','LOX-1','mAmP','MB','MCP-1','MMP-1','MMP-10','MMP-12','MMP-3','MMP-7','MPO','NEMO','NT-pro_BNP','OPG','PAPPA','PAR-1','PDGF_subunit_B','PECAM-1','PlGF','PSGL-1','PTX3','RAGE','REN','RETN','SCF','SELE','SIRT2','SPON1','ST2','t-PA','TF','TIE2','TM','TNF-R1','TNF-R2','TNFSF14','TRAIL','TRAIL-R2','TRANCE','U-PAR','VEGF-A','VEGF-D')





shinyUI(bootstrapPage(
  head(),
  navigation(),
  titlePanel("SCALLOP: Genome-region browser"),
  beginPage(),
  beginPanel('1/3'),
  HTML("This page can be used to browse the pQTL effects of all SNPs proximal to any gene in the genome:<br><br>"),
  # textInput(inputId="email", label = "E-mail", value = ""),
  textInput(inputId="gene", label = "Gene (HGNC) or position (chr:bp)", value = ""),
  checkboxInput("advanced", "Advanced options", value = FALSE),
  conditionalPanel(
    condition = "input.advanced",
    checkboxInput("show_gene_map", "Show gene map", value = FALSE)
  ),
  conditionalPanel(
    condition = "input.advanced",
    selectInput("protein_to_highlight", "Highlight protein", choices = c("none",proteins),selected="none")
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
  HTML("<br><br>"),
  dataTableOutput("mainTable"),
  endPanel(),
  
  endPage(),
  footer()
  
)	
)







