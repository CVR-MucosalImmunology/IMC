library(shiny)
library(SpatialExperiment)
library(tidyverse)
library(shinythemes)
library(data.table)

spe_subset = readRDS("../analysis/RDSFiles/shinyapp_spe.rds")

# Define UI for application with a navbar
ui <- navbarPage(
  theme = shinytheme("united"),
  # Application title
  "Adjusting Celltype Populations",
  
  # First tab for marker expression
  tabPanel("Marker Expression",
           # Sidebar layout for marker expression
           sidebarLayout(
             sidebarPanel(
               # Select input for target (Cluster)
               selectInput("target", 
                           "Select Target Celltype:", 
                           choices = sort(unique(as.character(spe_subset$Cluster))), 
                           selected = sort(unique(as.character(spe_subset$Cluster)))[1]),
               
               # Select input for marker
               selectInput("marker", 
                           "Select Marker:", 
                           choices = sort(rownames(spe_subset)), 
                           selected = sort(rownames(spe_subset))[1]),
               
               # Numeric input for value threshold
               numericInput("value", 
                            "Set Value Threshold:", 
                            value = 0.5, 
                            min = NA, 
                            max = NA, 
                            step = 0.05),
               
               # Conditional donor selection checkbox group
               uiOutput("donor_selection")
             ),
             
             # Main panel with tabs for overall marker expression and donor-specific expression
             mainPanel(
               tabsetPanel(
                 id = "tabs",  # Add an id to track selected tab
                 # Tab for overall marker expression
                 tabPanel(
                   "Overall Marker Expression",
                   br(),
                   plotOutput("Plot")
                 ),
                 # Tab for viewing by donor
                 tabPanel(
                   "Marker Expression by Donor ID",
                   br(),
                   plotOutput("DonorPlot", height = "auto")  # Set dynamic height
                 )
               )
             )
           )
  ),
  
  # New tab at the top of the page for "Move Cells"
  # New tab at the top of the page for "Move Cells"
  tabPanel("Move Cells",
           sidebarLayout(
             sidebarPanel(
               fluidRow(
                 column(6, selectInput("from", "Move From:", 
                                       choices = sort(unique(as.character(spe_subset$Cluster))),
                                       selected = sort(unique(as.character(spe_subset$Cluster)))[1])),
                 column(6, selectInput("to", "To:", 
                                       choices = c(sort(unique(as.character(spe_subset$Cluster))), "Other"),
                                       selected = sort(unique(as.character(spe_subset$Cluster)))[2]))
               ),
               conditionalPanel(
                 condition = "input.to == 'Other'",
                 textInput("custom_cluster_to", "Please specify the population:")
               ),
               fluidRow(
                 column(4, selectInput("marker_move", "If Marker:", choices = sort(rownames(spe_subset)))),
                 column(4, selectInput("direction", "Direction:", choices = c("<", ">"))),
                 column(4, numericInput("value_move", "Cutoff:", value = 0.5, min = NA, max = NA, step = 0.05))
               ),
               fluidRow(
                 column(6, actionButton("move_cells", "Reassign Cells")),
                 column(6, actionButton("undo_last", "Undo Last Move"))
               ),
               # Text input for the user to input the filename
               br(),
               textInput("filename_input", "Enter Filename to Save RDS:", value = "new_spe_subset"),
               # Save button with a save icon
               actionButton("save_spe", "Save", icon = icon("save"))
             ),
             mainPanel(
               fluidRow(
                 # Left column: Most Recent Cell Movement
                 column(
                   width = 6,  # Use 6 out of 12 units for the left column
                   h4("Most Recent Cell Movement"),
                   tableOutput("cell_moved_summary")  # Display only the latest cell moved info
                 ),
                 
                 # Right column: Edit History
                 column(
                   width = 6,  # Use 6 out of 12 units for the right column
                   h4("Edit History"),
                   tableOutput("edit_history")  # Display only the accumulated edit history
                 )
               )
             )
           )
  )
)

# Define server logic required to draw the plots and handle inputs
server <- function(input, output, session) {
  
  # Render the donor selection UI conditionally
  output$donor_selection <- renderUI({
    if (input$tabs == "Marker Expression by Donor ID") {
      checkboxGroupInput("donor", 
                         "Select Donor ID(s):", 
                         choices = sort(unique(spe_subset$DonorID)),
                         selected = sort(unique(spe_subset$DonorID)))
    }
  })
  
  # Plot for overall marker expression
  output$Plot <- renderPlot({
    spe_subset = spe_subset_rv()
    
    exprs_vector = assay(spe_subset, "norm_exprs")[input$marker, colData(spe_subset)[["Cluster"]] == input$target]
    data.frame(x=exprs_vector) %>%
      ggplot(aes(x=x)) +
      geom_density(fill="#69b3a2", colour="black", alpha=0.8) +
      theme_bw() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      geom_vline(xintercept=input$value, colour="black", lwd=1, linetype="dashed") +
      xlim(0, 1) +
      labs(
        x = sprintf("%s Expression ('norm_exprs')", input$marker),
        y = "Density",
        title = sprintf("Expression of %s in '%s'", input$marker, input$target)
      ) +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  # Dynamically set the height of the donor plot based on the number of selected donors
  output$DonorPlot <- renderPlot({
    req(input$donor)  # Make sure donors are selected before rendering the plot
    spe_subset = spe_subset_rv()
    
    # Filter expression data for selected donors
    exprs_vector = assay(spe_subset, "norm_exprs")[input$marker, 
                                                   colData(spe_subset)[["Cluster"]] == input$target & 
                                                     colData(spe_subset)[["DonorID"]] %in% input$donor]
    
    # Create data frame for plotting, including donor information
    donor_data <- data.frame(
      x = exprs_vector, 
      donor = factor(colData(spe_subset)[["DonorID"]][colData(spe_subset)[["Cluster"]] == input$target & 
                                                        colData(spe_subset)[["DonorID"]] %in% input$donor])
    )
    
    # Plot multiple donor expression data in faceted rows
    donor_data %>%
      ggplot(aes(x=x)) +
      geom_density(aes(fill=donor), colour="black", alpha=0.8) +
      facet_wrap(~donor, ncol = 1, scales = "free_y") +  # Facet by donor, one per row
      theme_bw() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      geom_vline(xintercept=input$value, colour="black", lwd=1, linetype="dashed") +
      xlim(0, 1) +
      labs(
        x = sprintf("%s Expression ('norm_exprs')", input$marker),
        y = "Density",
        title = sprintf("Expression of %s in '%s' by Donor ID", input$marker, input$target),
        fill = "Donor ID"
      ) +
      theme(plot.title = element_text(hjust = 0.5))
  }, height = reactive({
    max(400, length(input$donor) * 150)  # Ensure minimum height and scale with the number of donors
  }))
  
  # ------------------------------ Move Cells --------------------------------
  # Initialize reactive values to store the dataset and movement history
  spe_subset_rv <- reactiveVal()  # Make spe_subset reactive
  
  # Load spe_subset when the app starts
  observe({
    spe_subset <- readRDS("../analysis/RDSFiles/shinyapp_spe.rds")
    spe_subset_rv(spe_subset)  # Store it in reactive value
  })
  
  # Initialize a reactive data frame to store movement history (for edit history)
  movement_history <- reactiveVal(data.frame(
    From = character(0),
    To = character(0),
    Marker = character(0),
    Direction = character(0),
    Value = numeric(0),
    stringsAsFactors = FALSE
  ))
  
  # Initialize a reactive value to store the most recent movement summary
  recent_move_summary <- reactiveVal(data.frame(
    From = character(0),
    To = character(0),
    Donor = character(0),
    CellsMoved = integer(0),
    stringsAsFactors = FALSE
  ))
  
  # Function to move cells based on the provided criteria
  moveCells <- function(from, to, marker, direction, value) {
    spe_subset <- spe_subset_rv()  # Get the current spe_subset
    
    # Extract expression data for the specified marker and from-cluster
    cells_in_from <- which(colData(spe_subset)[["Cluster"]] == from)  # Cells still in 'from' cluster
    if (length(cells_in_from) == 0) {
      return(table())  # Return an empty table if no cells are in the 'from' cluster
    }
    
    exprs_vector <- assay(spe_subset, "norm_exprs")[marker, cells_in_from]  # Expression data for those cells
    
    # Get donor IDs for cells in the from-cluster
    donor_ids <- colData(spe_subset)[["DonorID"]][cells_in_from]
    
    # Identify cells that meet the condition (e.g., based on marker expression and direction)
    if (direction == "<") {
      cells_to_move <- which(exprs_vector < value)
    } else {
      cells_to_move <- which(exprs_vector > value)
    }
    
    # Filter out cells that have already been moved (already in 'to' cluster)
    cells_to_move <- cells_to_move[colData(spe_subset)$Cluster[cells_in_from[cells_to_move]] != to]
    
    # Move only the cells that haven't already been moved
    if (length(cells_to_move) > 0) {
      # Update the cluster assignment for the selected cells in the spe_subset object
      colData(spe_subset)$Cluster[cells_in_from[cells_to_move]] <- to
      spe_subset_rv(spe_subset)  # Update the reactive value with the modified spe_subset
    }
    
    # Count how many cells were moved for each donor
    moved_donor_counts <- table(donor_ids[cells_to_move])
    
    if (length(moved_donor_counts) == 0) {
      moved_donor_counts <- table(NA)  # Return a table with NA if no cells were moved
    }
    
    return(moved_donor_counts)
  }
  
  # Function to apply all changes from the movement history
  apply_changes <- function(movements) {
    spe_subset <- readRDS("../analysis/RDSFiles/shinyapp_spe.rds")  # Reset to the original state
    spe_subset_rv(spe_subset)  # Store it in the reactive value
    
    for (i in 1:nrow(movements)) {
      moveCells(
        from = movements$From[i],
        to = movements$To[i],
        marker = movements$Marker[i],
        direction = movements$Direction[i],
        value = movements$Value[i]
      )
    }
  }
  
  # Observe when the "Reassign Cells" button is clicked
  observeEvent(input$move_cells, {
    # Get the current spe_subset
    spe_subset <- spe_subset_rv()
    
    # Get all unique donor IDs and sort them in ascending order
    all_donor_ids <- sort(unique(colData(spe_subset)[["DonorID"]]))
    
    # Perform the cell movement
    moved_counts <- moveCells(input$from, ifelse(input$to == "Other", input$custom_cluster_to, input$to), 
                              input$marker_move, input$direction, input$value_move)
    
    # Add this move to the edit history
    new_entry <- data.frame(
      From = input$from,
      To = ifelse(input$to == "Other", input$custom_cluster_to, input$to),
      Marker = input$marker_move,
      Direction = input$direction,
      Value = input$value_move,
      stringsAsFactors = FALSE
    )
    history <- rbind(movement_history(), new_entry)
    movement_history(history)
    
    # Save the updated history to CSV
    fwrite(history, "ClusterModTracker.csv")
    
    # Handle case where no cells were moved (i.e., empty moved_counts)
    if (length(moved_counts) == 0) {
      # No cells moved, create a data.frame with Donor IDs and CellsMoved = 0
      recent_summary <- data.frame(
        From = input$from,
        To = ifelse(input$to == "Other", input$custom_cluster_to, input$to),
        Donor = all_donor_ids,  # All Donor IDs included, sorted in ascending order
        CellsMoved = 0,  # CellsMoved set to 0
        stringsAsFactors = FALSE
      )
    } else {
      # Cells were moved, create the data frame with the actual values
      recent_summary <- data.frame(
        From = input$from,
        To = ifelse(input$to == "Other", input$custom_cluster_to, input$to),
        Donor = all_donor_ids,
        CellsMoved = sapply(all_donor_ids, function(donor) {
          if (donor %in% names(moved_counts)) {
            return(moved_counts[donor])
          } else {
            return(0)
          }
        }),
        stringsAsFactors = FALSE
      )
    }
    
    # Update the most recent move summary for display
    recent_move_summary(recent_summary)
  })
  
  # Output the edit history (accumulative)
  output$edit_history <- renderTable({
    history <- movement_history()
    
    # Handle case when history is empty
    if (nrow(history) > 0) {
      history[, c("From", "To", "Marker", "Direction", "Value")]
    } else {
      data.frame(From = character(0), To = character(0), Marker = character(0), Direction = character(0), Value = numeric(0))
    }
  })
  
  # Output the most recent cell movement summary
  output$cell_moved_summary <- renderTable({
    recent_summary <- recent_move_summary()
    
    # Handle case when there are no recent moves
    if (nrow(recent_summary) > 0) {
      recent_summary[, c("From", "To", "Donor", "CellsMoved")]
    } else {
      data.frame(From = character(0), To = character(0), Donor = character(0), CellsMoved = integer(0))
    }
  })
  
  # Function to reset the data to its original state (before any cell movements)
  reset_to_original_state <- function() {
    spe_subset <- readRDS("../analysis/RDSFiles/shinyapp_spe.rds")  # Reload the original data
    spe_subset_rv(spe_subset)  # Store it in the reactive value
  }
  
  # Add an undo functionality
  observeEvent(input$undo_last, {
    history <- movement_history()
    
    # Only proceed if there is at least one entry in the history
    if (nrow(history) > 0) {
      # Remove the last entry
      history <- head(history, -1)  # Remove the last row
      movement_history(history)
      
      # Check if the history is empty after removing the last row
      if (nrow(history) == 0) {
        # History is empty, reset the spe_subset to the original state
        reset_to_original_state()
      } else {
        # Reapply all remaining changes from history after resetting to the original state
        apply_changes(history)
      }
      
      # Save the updated history to CSV
      fwrite(history, "ClusterModTracker.csv")
      
      # Clear the recent move summary (since it was undone)
      recent_move_summary(data.frame(From=character(0), To=character(0), Donor=character(0), CellsMoved=integer(0)))
    }
  })
  
  # Function to reset the data to its original state (before any cell movements)
  reset_to_original_state <- function() {
    spe_subset <- readRDS("../analysis/RDSFiles/shinyapp_spe.rds")  # Reload the original data
    spe_subset_rv(spe_subset)  # Store it in the reactive value
  }
  
  # Function to apply all changes from the movement history
  apply_changes <- function(movements) {
    if (nrow(movements) == 0) {
      reset_to_original_state()  # Reset to the original state if no moves exist
    } else {
      reset_to_original_state()  # Reset to the original state
      for (i in 1:nrow(movements)) {
        moveCells(
          from = movements$From[i],
          to = movements$To[i],
          marker = movements$Marker[i],
          direction = movements$Direction[i],
          value = movements$Value[i]
        )
      }
    }
  }
  
  # Save the spe_subset as an RDS file when the save button is clicked
  observeEvent(input$save_spe, {
    # Get the current spe_subset
    spe_subset <- spe_subset_rv()
    
    # Get the filename from the input field
    filename <- input$filename_input
    
    # Ensure that the filename is valid
    if (nchar(filename) > 0) {
      # Construct the full path with ".rds" extension
      filepath <- paste0(filename, ".rds")
      
      # Save the spe_subset as an RDS file
      saveRDS(spe_subset, filepath)
      
      # Show a notification or message that the file has been saved
      showNotification(paste("File saved as:", filepath), type = "message")
    } else {
      # Show a warning notification if the filename is invalid
      showNotification("Please enter a valid filename", type = "warning")
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)