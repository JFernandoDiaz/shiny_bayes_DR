library(shiny)
library(shinydashboard)
library("clinDR")

exdat <- read.csv(file="dosedat2.csv")

mcmc<-mcmc.control(chains=3)


ui <- dashboardPage(
  dashboardHeader(title="Bayesian Dose Response",
                  tags$li(class="dropdown",
                          tags$style(".main-header {max-height: 50px}")
                          )
                  ),
  dashboardSidebar(
  
    numericInput("epmu","Mean for E0 in a t-prior distribution",value=0),
    numericInput("epsca","Scale parameter in a t-prior distribution",value=4),
    numericInput("difTargetmu","Mean Target difference",value=0),
    numericInput("difTargetsca","Scale Target difference",value=4),
    numericInput("dTarget","Target dose for prior dose",value=50),
    numericInput("p50","Projected ED50",value=10),
    numericInput("sigmalow","Lower bound, prior distribution, SD",value=0.01),
    numericInput("sigmaup","Upper bound, prior distribution, SD",value=3)
  ),
  dashboardBody(
    fluidRow (
      box(   title = "Dose Response", width = 20,
             plotOutput("plots1",height = 500) 
      )
    )
  )	
)


server <- function(input, output, session) {

    output$plots1 <- renderPlot({

      prior<-emaxPrior.control(epmu=input$epmu,epsca=input$epsca,difTargetmu=input$difTargetmu,difTargetsca=input$difTargetsca,dTarget=input$dTarget,
                               #                         p50=(2+5)/2,
                               p50=input$p50, sigmalow=input$sigmalow,sigmaup=input$sigmaup)
      
      fitout<-fitEmaxB(exdat$RESP,exdat$DOSE,prior,modType=4,prot=rep(1, length(exdat$RESP) ),
                       count=rep(1, length(exdat$RESP) ),msSat=NULL,mcmc=mcmc)
      
      parms <- coef(fitout)[,1:4] #use first intercept
      
      plotB(exdat$RESP, exdat$DOSE, parms, sigma2 = (sigma(fitout) )^2, ylab = "Effect" , xlab = "Dose (mg)", xat = c(0,5,10,20,30,40,50), 
            dgrid = sort(unique(c(seq(0,max(exdat$DOSE), length = 50), exdat$DOSE) ) ) ) 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
