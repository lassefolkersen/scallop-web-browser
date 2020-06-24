source("../uifunctions.R")
initialize('con',TRUE)

shinyUI(bootstrapPage(
	head(),
	navigation(),
	titlePanel("Update"),
	beginPage(),
	HTML("Control button for pulling a new github commit into the server"),
	actionButton("goButton","update"),
	htmlOutput("text1"),
	footer()
))












