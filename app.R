library(shiny)
library(DT)

# eGFR helper functions
kappa <- function(sex) if (sex == "Female") 0.7 else 0.9
alpha <- function(sex) if (sex == "Female") -0.241 else -0.302
sex_coeff <- function(sex) if (sex == "Female") 1.012 else 1

compute_egfr <- function(creat, age, sex) {
  cr_over_kappa <- (creat / 88.4) / kappa(sex)
  142 *
    (pmin(cr_over_kappa, 1) ^ alpha(sex)) *
    (pmax(cr_over_kappa, 1) ^ (-1.200)) *
    (0.9938 ^ age) *
    sex_coeff(sex)
}

ui <- fluidPage(
  titlePanel("Kinetic GFR Visualiser (with DT + Modal Input)"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("age", "Age:", min = 18, max = 100, value = 65),
      numericInput("weight", "Weight (kg):", min = 30, max = 220, value = 70),
      selectInput("sex", "Sex:", choices = c("Male", "Female")),
      actionButton("add_entry", "Add Creatinine Entry")
    ),
    mainPanel(
      DTOutput("creat_table"),
      plotOutput("gfr_plot")
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    data = data.frame(
      Date = as.Date(character()),
      Time = character(),
      Creatinine = numeric(),
      stringsAsFactors = FALSE
    )
  )
  
  # Modal dialog for adding entries
  observeEvent(input$add_entry, {
    showModal(modalDialog(
      title = "Add Creatinine Entry",
      dateInput("modal_date", "Date:", value = Sys.Date(), format = "dd/mm/yyyy"),
      textInput("modal_time", "Time (HH:MM):", value = format(Sys.time(), "%H:%M")),
      numericInput("modal_creat", "Creatinine (µmol/L):", value = 80, min = 0),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("submit_entry", "Submit")
      )
    ))
  })
  
  # Handle submit from modal
  observeEvent(input$submit_entry, {
    new_row <- data.frame(
      Date = input$modal_date,
      Time = input$modal_time,
      Creatinine = input$modal_creat,
      stringsAsFactors = FALSE
    )
    rv$data <- rbind(rv$data, new_row)
    removeModal()
  })
  
  # Render DT table
  output$creat_table <- renderDT({
    datatable(rv$data, editable = TRUE, options = list(pageLength = 10), rownames = FALSE) %>%
      formatDate(columns = "Date", method = "toLocaleDateString", params = list("en-GB"))
  })
  
  observeEvent(input$creat_table_cell_edit, {
    info <- input$creat_table_cell_edit
    str(info)  # optional: view in R console
    i <- info$row
    j <- info$col + 1  # DT uses 0-based indexing
    v <- info$value
    
    # Coerce value based on column
    if (names(rv$data)[j] == "Date") {
      rv$data[i, j] <- as.Date(v)
    } else if (names(rv$data)[j] == "Creatinine") {
      rv$data[i, j] <- as.numeric(v)
    } else {
      rv$data[i, j] <- v
    }
  })
  
  
  # GFR Plot
  output$gfr_plot <- renderPlot({
    df <- rv$data
    req(nrow(df) >= 2)
    
    # Convert to datetime
    df$datetime <- as.POSIXct(paste(df$Date, df$Time), format = "%Y-%m-%d %H:%M", tz = "UTC")
    df <- df[order(df$datetime), ]
    validate(
      need(!any(is.na(df$datetime)), "Date/time formatting issue — check entries.")
    )
    
    df$eGFR <- NA_real_
    df$eGFR[1] <- compute_egfr(df$Creatinine[1], input$age, input$sex)
    
    for (i in 2:nrow(df)) {
      cr1 <- df$Creatinine[i - 1]
      cr2 <- df$Creatinine[i]
      egfr1 <- df$eGFR[i - 1]
      
      delta_creat <- cr2 - cr1
      time_diff <- as.numeric(difftime(df$datetime[i], df$datetime[i - 1], units = "hours"))
      avg_creat <- mean(c(cr1, cr2))
      Vd <- if (input$sex == "Female") 0.5 * input$weight else 0.6 * input$weight
      adj <- min(1, max(0, 1 - ((24 * delta_creat) / (time_diff * Vd))))
      
      df$eGFR[i] <- cr1 * egfr1 / avg_creat * adj
    }
    
    plot(df$datetime, df$eGFR, type = "b", lwd = 2,
         xlab = "Time", ylab = "Kinetic eGFR (ml/min/1.73m²)",
         main = "Kinetic eGFR Over Time")
    grid()
  })
}

shinyApp(ui, server)
