library("shiny")


shinyServer(function(input, output) {
	
	

	output$text1 <- renderText({
	  if(input$goButton == 0){
	    return(NULL)
	  }
		
		if(input$goButton > 0){
	    setwd("/home/ubuntu/srv/scallop-web-browser")
	    m1<-system("git pull",intern=T)
	    setwd("/home/ubuntu/srv/scallop-web-browser/update")
	    Sys.sleep(1)
	    m2<-paste("<br>Executed update command with this output:<br>",m1)
	    return(m2)  
	  }
	  
		
		
		
	})


	
	 
	  
	
})


