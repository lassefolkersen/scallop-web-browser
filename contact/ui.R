source("../uifunctions.R")
initialize('con',TRUE)

shinyUI(bootstrapPage(
	head(),
	navigation(),
	titlePanel("Contact"),
	beginPage(),
	HTML(	"
<b>Web-app design:</b><br><a href='http://orcid.org/0000-0003-0708-9530'>Lasse Folkersen</a>.<br><br>
				
<b>Further details:</b><br><a href='https://www.olink.com/scallop/'>SCALLOP consortium</a>.<br><br>



	      				"),
	footer()
))












