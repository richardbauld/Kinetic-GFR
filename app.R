library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library(ggpubr)
library(rmarkdown)
library(knitr)

# CKD-EPI 2021 helper functions
kappa <- function(sex) if (sex == "Female") 0.7 else 0.9
alpha <- function(sex) if (sex == "Female") -0.241 else -0.302
sex_coeff <- function(sex) if (sex == "Female") 1.012 else 1

compute_egfr <- function(creat_umol, age, sex) {
  cr_over_kappa <- (creat_umol / 88.4) / kappa(sex)
  142 *
    (pmin(cr_over_kappa, 1) ^ alpha(sex)) *
    (pmax(cr_over_kappa, 1) ^ (-1.200)) *
    (0.9938 ^ age) *
    sex_coeff(sex)
}

# Example patient (with empty Events)
example_patients <- list(
  "Patient A" = data.frame(
    Date = as.Date(c("2025-04-01", "2025-05-05", "2025-05-05", "2025-05-05",
                     "2025-05-06", "2025-05-06", "2025-05-06", "2025-05-07",
                     "2025-05-07", "2025-05-08", "2025-05-09", "2025-05-10",
                     "2025-05-11", "2025-05-12", "2025-05-13", "2025-05-13",
                     "2025-05-14", "2025-05-15", "2025-05-16", "2025-05-17",
                     "2025-05-18", "2025-05-19")),
    Time = c("08:00", "02:00", "04:45", "12:00", "04:00", "12:00", "22:00",
             "06:00", "14:00", "04:00", "04:00", "04:00", "04:00", "04:00", "04:00",
             "18:00", "04:00", "02:00", "04:00", "04:00", "04:00", "04:00"),
    Creatinine = c(3.0*88.4, 2.93*88.4, 3.1*88.4, 3.6*88.4, 4.6*88.4, 4.93*88.4, 3.7*88.4, 2.3*88.4, 1.9*88.4, 1.2*88.4, 1.1*88.4,
                   1.1*88.4, 1.15*88.4, 1.1*88.4, 1.0*88.4, 1.2*88.4, 1.8*88.4, 3.2*88.4, 3.6*88.4, 3.15*88.4, 3*88.4, 4.7*88.4),
    Event = rep("", 22),
    stringsAsFactors = FALSE
  ),
  "Patient B" = data.frame(
    Date = as.Date(c("2022-10-20", "2022-11-04", "2022-11-05", "2022-11-06",
                     "2022-11-06", "2022-11-06", "2022-11-07", "2022-11-08",
                     "2022-11-08", "2022-11-09", "2022-11-10", "2022-11-11")),
    Time = c("06:00", "07:27", "15:12", "05:19", "13:52", "22:41", "05:14",
             "05:44", "15:48", "04:41", "04:52", "04:59"),
    Creatinine = c(102, 171, 324, 251, 236, 233, 236, 199, 179, 156, 104, 89),
    Event = c("Baseline", "RIE Admission", "CRRT started", "Theatre", "", "", "", "", "",
              "CRRT Stopped", "", ""),
    stringsAsFactors = FALSE
  ),
  "Patient C" = data.frame(
    Date = as.Date(c(
      "2025-05-05", "2025-05-06", "2025-05-07", "2025-05-08", "2025-05-09",
      "2025-05-10", "2025-05-11", "2025-05-12", "2025-05-13", "2025-05-14",
      "2025-05-15", "2025-05-16", "2025-05-17", "2025-05-18", "2025-05-19"
    )),
    Time = rep("08:00", 15),
    Creatinine = c(
      3.9, 4.8, 6.7, 5.0, 3.2,
      2.1, 1.9, 1.8, 1.8, 2.0,
      3.0, 3.8, 3.5, 3.0, 5.2
    ) * 88.4,
    Event = rep("", 15),
    stringsAsFactors = FALSE
  )
)

ui <- fluidPage(
  titlePanel("Kinetic GFR Visualiser (Chen 2017, µmol/L Input)"),
  sidebarLayout(
    sidebarPanel(
      numericInput("age", "Age:", min = 18, max = 100, value = 65),
      numericInput("weight", "Weight (kg):", min = 30, max = 220, value = 70),
      numericInput("base_creat", "Baseline Creatinine (umol/L):", min = 30, max = 1000, value = 80),
      selectInput("sex", "Sex:", choices = c("Male", "Female")),
      selectInput("example_patient", "Load Example Patient:", choices = c("None", names(example_patients))),
      actionButton("add_entry", "Add Creatinine Entry"),
      actionButton("reset", "Reset All Data"),
      downloadButton("export_csv", "Export CSV"),
      downloadButton("export_pdf", "Export PDF"),
      downloadButton("save_session", "Save Session"),
      fileInput("load_session", "Load Session"),
      checkboxInput("show_creat_plot", "Show Creatinine Plot", value = FALSE)
    ),
    mainPanel(
      tableOutput(("debug_table")),
      DTOutput("creat_table"),
      h4("Kinetic GFR Over Time:"),
      plotlyOutput("gfr_plot"),
      textOutput("fixed_vales"),
      conditionalPanel(
        condition = "input.show_creat_plot == true",
        h4("Serum Creatinine Over Time:"),
        plotlyOutput("creat_plot")
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    data = data.frame(
      Date = as.Date(character()),
      Time = character(),
      Creatinine = numeric(),
      Event = character(),
      stringsAsFactors = FALSE
    )
  )
  
  observeEvent(input$example_patient, {
    if (input$example_patient != "None") {
      rv$data <- example_patients[[input$example_patient]]
    }
  })
  
  observeEvent(input$add_entry, {
    showModal(modalDialog(
      title = "Add Creatinine Entry",
      dateInput("modal_date", "Date:", value = Sys.Date(), format = "dd/mm/yyyy"),
      textInput("modal_time", "Time (HH:MM):", value = format(Sys.time(), "%H:%M")),
      numericInput("modal_creat", "Creatinine (µmol/L):", value = 80, min = 0),
      textInput("modal_event", "Event (optional):", value = ""),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("submit_entry", "Submit")
      )
    ))
  })
  
  observeEvent(input$submit_entry, {
    new_row <- data.frame(
      Date = input$modal_date,
      Time = input$modal_time,
      Creatinine = input$modal_creat,
      Event = input$modal_event,
      stringsAsFactors = FALSE
    )
    rv$data <- rbind(rv$data, new_row)
    removeModal()
  })
  
  observeEvent(input$reset, {
    rv$data <- rv$data[0, ]
  })
  
  output$creat_table <- renderDT({
    datatable(rv$data, editable = TRUE, rownames = FALSE)
  })
  
  observeEvent(input$creat_table_cell_edit, {
    info <- input$creat_table_cell_edit
    i <- info$row
    j <- info$col + 1
    v <- info$value
    rv$data[i, j] <- v
  })
  
  processed_data <- reactive({
    df <- rv$data
    req(nrow(df) >= 2)
    df$Time <- format(strptime(df$Time, "%H:%M"), "%H:%M:%S")
    df$datetime <- as.POSIXct(paste(df$Date, df$Time), tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
    df <- df[order(df$datetime), ]
    
    cr_umol <- df$Creatinine
    cr_mgdl <- cr_umol / 88.4
    weight <- input$weight
    sex <- input$sex
    age <- input$age
    
 #   TBW <- if (input$sex == "Female") 0.6 * weight else 0.7 * weight
    TBW <- 0.6 * weight
    
 #   cr_ss <- cr_mgdl[1]  # steady-state creatinine
    cr_ss <- input$base_creat / 88.4
 #   egfr_ss <- compute_egfr(df$Creatinine[1], input$age, input$sex)
    egfr_ss <- compute_egfr(input$base_creat, input$age, input$sex)
    
    kGFR <- rep(NA, length(cr_mgdl))
    kGFR[1] <- NA
    
    eGFR <- rep(NA, length(cr_mgdl))
    eGFR[1] <- egfr_ss

    for (i in 2:length(cr_mgdl)) {
      cr1 <- cr_mgdl[i - 1]
      cr2 <- cr_mgdl[i]
      cr_mean <- (cr1 + cr2) / 2
      delta_cr <- cr2 - cr1
      delta_t <- as.numeric(difftime(df$datetime[i], df$datetime[i - 1], units = "hours"))
      max_delta_cr <- (cr_ss * egfr_ss) / TBW
      
      numerator <- cr_ss * egfr_ss
      denominator <- cr_mean
      correction <- (24 * delta_cr) / (delta_t * max_delta_cr)
      adj <- 1 - correction
      adj <- max(0, min(1, adj))
      kGFR[i] <- (numerator / denominator) * adj
    }
    
    
    for(i in 2:length(cr_mgdl)) {
      cr_egfr <- cr_mgdl[i]
      cr_over_kappa <- cr_egfr / kappa(sex)
      
      eGFR[i] <- 142 *
        (pmin(cr_over_kappa, 1) ^ alpha(sex)) *
        (pmax(cr_over_kappa, 1) ^ (-1.200)) *
        (0.9938 ^ age) *
        sex_coeff(sex)
    }
    
    df$eGFR <- eGFR
    df$kGFR <- kGFR
    df$Label <- paste0(
      "Date: ", format(df$datetime, "%d/%m/%Y %H:%M"), "<br>",
      "Creatinine: ", round(df$Creatinine, 1), " µmol/L<br>",
      "eGFR: ", round(df$eGFR, 1), " µmol/L<br>",
      "kGFR: ", ifelse(is.na(df$kGFR), "", paste0(round(df$kGFR, 1), " ml/min/1.73m²")),
      ifelse(df$Event != "", paste0("<br><b>Event:</b> ", df$Event), "")
    )
    list(data = df, egfr_ss = egfr_ss)
  })
  
  
  
  clean_data <- reactive({
    df <- processed_data()$data
    df[!is.na(df$kGFR), ]
  })

 
    
  output$debug_table <- renderTable({
    head(clean_data())
  })
  
  output$gfr_plot <- renderPlotly({
    df <- clean_data()
    p <- ggplot(df, aes(x = datetime)) +
      geom_line(aes(y = kGFR), color = "blue") +
      geom_point(aes(y = kGFR, text = Label), size = 2) +
      geom_line(aes(y = eGFR), color = "darkgreen") +
      geom_point(aes(y = eGFR, text = Label), size = 2) +
      labs(x = "Time", y = "Kinetic eGFR (ml/min/1.73m²)", title = "Kinetic GFR Over Time") +
      theme_minimal()
    event_df <- df[df$Event != "", ]
    if (nrow(event_df) > 0) {
      p <- p + geom_vline(data = event_df, aes(xintercept = as.numeric(datetime)), linetype = "dashed", color = "red")
    }
    ggplotly(p, tooltip = "text")
  })
  
  
  output$creat_plot <- renderPlotly({
    df <- clean_data()
    p <- ggplot(df, aes(x = datetime, y = Creatinine)) +
      geom_line(color = "darkgreen") +
      geom_point(aes(text = Label), size = 2) +
      labs(x = "Time", y = "Creatinine (µmol/L)", title = "Serum Creatinine Over Time") +
      theme_minimal()
    event_df <- df[df$Event != "", ]
    if (nrow(event_df) > 0) {
      p <- p + geom_vline(data = event_df, aes(xintercept = as.numeric(datetime)), linetype = "dashed", color = "red")
    }
    ggplotly(p, tooltip = "text")
  })
  
  output$fixed_vales <- renderText({
    req(processed_data()$egfr_ss)
    paste("Baseline eGFR (CKD-EPI 2021):", round(processed_data()$egfr_ss, 1), "ml/min/1.73m²")
  })
}

shinyApp(ui, server)
