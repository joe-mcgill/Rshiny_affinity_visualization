#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)

ui <- function(){
  navbarPage(title = 'Hello', 
             tabPanel("title2"),
             tabPanel("title3"),
             tags$head(
               tags$style(type = 'text/css', 
                          HTML('.navbar { background-color: red;}
                          .navbar-default .navbar-brand{color: white;}
                          .tab-panel{ background-color: red; color: white}
                          .navbar-default .navbar-nav > .active > a, 
                           .navbar-default .navbar-nav > .active > a:focus, 
                           .navbar-default .navbar-nav > .active > a:hover {
                                color: #555;
                                background-color: white;
                            }')
               )
             )
  )
  
}

server <- function(input, output, session){
}


shinyApp(ui = ui, server = server)