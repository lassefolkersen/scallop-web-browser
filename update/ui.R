source("../uifunctions.R")
initialize('con',TRUE)

shinyUI(bootstrapPage(
	head(),
	navigation(),
	titlePanel("Update"),
	beginPage(),
	actionButton("goButton","update"),
	htmlOutput("text1"),
	footer()
))












