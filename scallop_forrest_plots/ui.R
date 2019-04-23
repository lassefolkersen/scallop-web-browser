library("openxlsx")
source("../uifunctions.R")
initialize('sti',TRUE)

#read S1
hit_list_path <- "/home/ubuntu/data/2019-04-23_forest_plot_data//S01_table_hit_list.xlsx"
data <- read.xlsx(hit_list_path)
rownames(data) <- apply(data[,c("Protein","MarkerName")],1,paste,collapse=" / ")


shinyUI(bootstrapPage(
	head(),
	navigation(),

	titlePanel("Forest plots"),
	beginPage(),
	beginPanel('1/3'),
	selectInput("s1_entry", "pQTLs", choices = rownames(data)),
	
	actionButton("goButton","Go"),

	endPanel(),
	beginPanel('2/3'),
	plotOutput("mainPlot"),
	endPanel(),
	
	
	endPage(),
	footer()
))






