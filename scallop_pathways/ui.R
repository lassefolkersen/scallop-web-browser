
source("../uifunctions.R")
initialize('sti',TRUE)

library("igraph")
library("visNetwork")


proteins<-c('ADM','AGRP','Beta-NGF','CA-125','CASP-8','CCL20','CCL3','CCL4','CD40','CD40-L','CHI3L1','CSF-1','CSTB','CTSD','CTSL1','CX3CL1','CXCL1','CXCL16','CXCL6','Dkk-1','ECP','EGF','EN-RAGE','ESM-1','FABP4','FAS','FGF-23','FS','GAL','Gal-3','GDF-15','GH','HB-EGF','HGF','hK11','HSP_27','IL-18','IL-1ra','IL-27','IL-6','IL-6RA','IL-8','IL16','ITGB1BP2','KIM-1','KLK6','LEP','LOX-1','mAmP','MB','MCP-1','MMP-1','MMP-10','MMP-12','MMP-3','MMP-7','MPO','NEMO','NT-pro_BNP','OPG','PAPPA','PAR-1','PDGF_subunit_B','PECAM-1','PlGF','PSGL-1','PTX3','RAGE','REN','RETN','SCF','SELE','SIRT2','SPON1','ST2','t-PA','TF','TIE2','TM','TNF-R1','TNF-R2','TNFSF14','TRAIL','TRAIL-R2','TRANCE','U-PAR','VEGF-A','VEGF-D')


proteins_to_omit <- c('Beta-NGF','CSF-1','ECP','EGF','GH','mAmp','NEMO','NT-pro_BNP')
#mAmp gave a weird error - the others are just not possible to show

proteins<-proteins[!proteins%in%proteins_to_omit]


shinyUI(bootstrapPage(
	head(),
	navigation(),

	titlePanel("Trans-pQTL pathways"),
	beginPage(),
	beginPanel('1/3'),
	HTML("Visualization of trans-pathways. For each of the target proteins, this module shows likely paths from trans-pQTL SNPs through possible gene intermediaries. The width of each edge indicates the significance of the finding, as calculated using permutation analysis in randomly generated networks."),

	selectInput("protein", "Protein", choices = proteins),
	
	actionButton("goButton","Go"),

	endPanel(),
	beginPanel('2/3'),
	visNetworkOutput("plot1"),
	# tableOutput("table1"),	
	htmlOutput("text_1"),
	endPanel(),
	
	
	endPage(),
	footer()
))






