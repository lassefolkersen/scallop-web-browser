source("../uifunctions.R")
initialize('con',TRUE)




phenotypes_vector<-c('ADM','AGRP','Beta-NGF','CA-125','CASP-8','CCL20','CCL3','CCL4','CD40','CD40-L','CHI3L1','CSF-1','CSTB','CTSD','CTSL1','CX3CL1','CXCL1','CXCL16','CXCL6','Dkk-1','ECP','EGF','EN-RAGE','ESM-1','FABP4','FAS','FGF-23','FS','GAL','Gal-3','GDF-15','GH','HB-EGF','HGF','hK11','HSP_27','IL-18','IL-1ra','IL-27','IL-6','IL-6RA','IL-8','IL16','ITGB1BP2','KIM-1','KLK6','LEP','LOX-1','mAmP','MB','MCP-1','MMP-1','MMP-10','MMP-12','MMP-3','MMP-7','MPO','NEMO','NT-pro_BNP','OPG','PAPPA','PAR-1','PDGF_subunit_B','PECAM-1','PlGF','PSGL-1','PTX3','RAGE','REN','RETN','SCF','SELE','SIRT2','SPON1','ST2','t-PA','TF','TIE2','TM','TNF-R1','TNF-R2','TNFSF14','TRAIL','TRAIL-R2','TRANCE','U-PAR','VEGF-A','VEGF-D')



shinyUI(bootstrapPage(
  head(),
  navigation(),
  titlePanel("SCALLOP: Protein-centric browser"),
  beginPage(),
  beginPanel('1/3'),
  HTML("This page can be used to browse the significant pQTL effects of SNPs anywhere in the genome:<br><br>"),
  selectInput("protein", "Protein", choices = phenotypes_vector),
  checkboxInput("advanced", "Advanced options", value = FALSE),
  conditionalPanel(
    condition = "input.advanced",
    sliderInput("p_value_cutoff","P-value cutoff (log)",min=5,max=12,value=6,step=0.1)
  ),
  conditionalPanel(
    condition = "input.advanced",
    sliderInput("top_label_count","Top-SNPs to label",min=0,max=30,value=3,step=1)
  ),
  conditionalPanel(
    condition = "input.advanced",
    checkboxInput("include_closest_genes","Also label closest 4 genes",value = TRUE)
  ),
  
  actionButton("goButton","Run analysis"),
  endPanel(),
  beginPanel('2/3'),
  plotOutput("mainPlot",height="800px",width="800px"),
  
  dataTableOutput("mainTable"),
  htmlOutput("text_1"),
  endPanel(),
  
  endPage(),
  footer()
  
)	
)







