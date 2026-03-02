library(shiny)
library(ggplot2)

# The SIR Model
SIR <- function(S,I,R,beta,gamma,T){
  i = 1
  N = S+I+R
  
  while(i<T){
    
    dS = -1*beta*S[i]*I[i]/N
    dI = beta*S[i]*I[i]/N - gamma*I[i]
    dR = gamma*I[i]
    
    S[i+1] = S[i] + dS
    I[i+1] = I[i] + dI
    R[i+1] = R[i] + dR
    
    i=i+1
  }
  
  return(data.frame(Time=1:T,S,I,R))
}


ui <- fluidPage(
  titlePanel("SIR Model"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("binput", 
                  label = "beta",
                  min=0, max=1, step=0.01, value=0.9),
                  
      
      sliderInput("ginput", 
                  label = "gamma",
                  min=0, max=1, step=0.01, value=0.3),
      
      numericInput("N", label="Population Size", value=10000),
      
      numericInput("I", label="Number of initial infections", value=10),
      
      numericInput("Time", label="Number of days to run", value=100)
      
      ),
    
    mainPanel(
      "The dynamics of infectious disease through a population.",
      plotOutput("plot")
    )
  )
  )




server <- function(input,output){
  
  
  
  output$plot <- renderPlot({ 
    
    data <- SIR(input$N-input$I,input$I,0,input$binput, input$ginput,input$Time)
    
    ggplot(data, aes(x=Time)) + geom_path(aes(y=I, color="B"), size=1) + geom_path(aes(y=S, color="A")) + geom_path(aes(y=R, color="C")) + xlab("Day") +ylab("People") + scale_colour_manual(name = '', values =c('A'='Orange','B'='Red', "C"="Green4"), labels = c("S","I","R")) + theme_bw()
  })
}
  
shinyApp(ui, server)
