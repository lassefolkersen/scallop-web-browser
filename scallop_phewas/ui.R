source("../uifunctions.R")
initialize('con',TRUE)



shinyUI(bootstrapPage(
  head(),
  navigation(),
  titlePanel("PheWas"),
  beginPage(),
  beginPanel('1/3'),
  HTML("This is intended to hold the phewas code:<br><br>"),
  actionButton("goButton","Run analysis"),
  endPanel(),
  beginPanel('2/3'),
  
  htmlOutput("text1"),
  endPanel(),
  
  endPage(),
  footer()
  
)	
)







