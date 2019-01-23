library("shiny")


shinyServer(function(input, output) {
  
  
  
  output$text1 <- renderText({
    stop(safeError("Not implemented yet"))
  })
  
  
})



