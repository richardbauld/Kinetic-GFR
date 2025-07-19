#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyTime)

# Helpers for eGFR calculation
kappa <- function(sex)
  if (sex == "Female"){
  return(0.7)
} else {
  return(0.9)
}

alpha <- function(sex)
  if (sex == "Female"){
    return(-0.241)
  } else {
    return(-0.302)
  }

sex_coeff <- function(sex)
  if (sex == "Female"){
    return(1.012)
  } else {
    return(1)
  }


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Kinetic GFR Visualiser"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30),
            numericInput("age", "Age:", min = 18, max = 100, value = 65),
            numericInput("weight", "Weight (kg):", min = 30, max = 220, value = 70),
            selectInput("sex", "Sex:", choices = c("Male", "Female")),
            numericInput("base_creat", "Baseline creatinine (umol/L):", min = 30, max = 1000, value = 80),
            dateInput("date", "Date of baseline creatinine:", value = Sys.Date(), format = "dd/mm/yyyy"),
            timeInput("time_input", "Enter time of baseline creatinine:", value = Sys.time(), seconds = FALSE),
            numericInput("creat_2", "New creatinine (umol/L):", min = 30, max = 3000, value = 80),
            dateInput("date_2", "Date of new creatinine:", value = Sys.Date(), format = "dd/mm/yyyy"),
            timeInput("time_input_2", "Enter time of 2nd creatinine:", value = Sys.time(), seconds = FALSE),
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("kineticPlot"),
           textOutput("dynamicText"),
           textOutput("debugText")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Time difference in hours
  time_diff <- reactive({
    datetime1 <- as.POSIXct(paste(input$date, format(input$time_input, "%H:%M:%S")))
    datetime2 <- as.POSIXct(paste(input$date_2, format(input$time_input_2, "%H:%M:%S")))
  as.numeric(difftime(datetime2, datetime1, units = "hours"))
    
    })

  #calc baseline GFR via CKD-EPI 2021 without race
  base_egfr <- reactive({
    cr_over_kappa <- (input$base_creat / 88.4) / kappa(input$sex)
    142 * 
      (pmin(cr_over_kappa, 1) ^ alpha(input$sex)) *
      (pmax(cr_over_kappa, 1) ^ (-1.200)) *
      (0.9938 ^ input$age) *
      sex_coeff(input$sex)
  })
  
  kinetic_gfr <- reactive({
    req(input$base_creat, input$creat_2, time_diff())
    avg_creat <- (input$base_creat + input$creat_2) / 2
    delta_creat <- input$creat_2 - input$base_creat
    Vd <- if (input$sex == "Female") 
      0.5 * input$weight
    else
      0.6 * input$weight
    adjustment <- max(0, 1 - ((24 * delta_creat) / (time_diff() * Vd)))
  
    input$base_creat * base_egfr() / avg_creat * adjustment
    
  })
  
    output$kineticPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
    
    output$dynamicText <- renderText({
      paste("Baseline eGFR by CKD-EPI:", round(base_egfr(), 0), "ml/min/1.73m²\n",
      "Kinetic GFR:", round(kinetic_gfr(),0))
      
    }
    )
      
      output$debugText <- renderText({
        cr1 <- input$base_creat
        cr2 <- input$creat_2
        delta <- cr2 - cr1
        avg <- (cr1 + cr2) / 2
        egfr <- base_egfr()
        td <- time_diff()
        vd <- if (input$sex == "Female") 0.5 * input$weight else 0.6 * input$weight
        adj <- max(0, 1 - ((24 * delta) / (td * vd)))
        kgfr <- cr1 * egfr / avg * adj
        
        paste(
          "Creat1:", cr1,
          "Creat2:", cr2,
          "Delta:", delta,
          "Time diff (hrs):", round(td, 2),
          "Vd:", round(vd, 1),
          "Adj factor:", round(adj, 3),
          "Base eGFR:", round(egfr, 1),
          "Kinetic GFR:", round(kgfr, 1)
        )
      })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
