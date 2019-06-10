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
	HTML("This page contains forest plots for all pQTL associations that were significant in the CVD1 meta-analysis. As further described in the main paper, all values are given in standard deviations of protein expression. Error bars indicate 95% confidence intervals.<br><br>"),
	
	selectInput("s1_entry", "pQTLs", choices = rownames(data)),
	
	actionButton("goButton","Go"),

	endPanel(),
	beginPanel('2/3'),
	plotOutput("mainPlot"),
	endPanel(),
	
	
	endPage(),
	footer()
))






