# Load necessary libraries
library(ShortRead)   # For reading FASTQ files
library(Biostrings)  # For handling sequences
library(dplyr)       # For data manipulation
library(readr)       # For reading/writing files
library(shiny)       # For the Shiny app
library(ggplot2)     # For plotting
library(UpSetR)      # For UpSet plots
library(tidyr)       # For data reshaping
library(DT)          # For interactive tables
library(gridExtra)   # For PDF table generation

# Step 1: Automatically detect all FASTQ files in the directory
fastq_files <- list.files(pattern = "\\.fastq$", full.names = TRUE)

if (length(fastq_files) == 0) {
  stop("No FASTQ files found in the directory!")
}

# Step 2: Define BLAST database parameters
db_file <- "unique_list.fasta"  # Your reference FASTA file
db_name <- "Unique-CC-DB"  # BLAST database name (without .fasta)

# Function to convert FASTQ to FASTA
convert_fastq_to_fasta <- function(fastq_file) {
  cat("Converting", fastq_file, "to FASTA...\n")
  
  # Define output FASTA file name
  fasta_file <- sub("\\.fastq$", ".fasta", fastq_file)
  
  # Read FASTQ file
  fastq_data <- readFastq(fastq_file)
  
  # Extract sequences and headers
  fasta_data <- sread(fastq_data)
  names(fasta_data) <- sub("^@", ">", as.character(ShortRead::id(fastq_data)))  # Convert headers
  
  # Write to FASTA file
  writeXStringSet(DNAStringSet(fasta_data), fasta_file, format = "fasta")
  
  cat("Conversion completed: Saved as", fasta_file, "\n")
  return(fasta_file)  # Return new FASTA file path
}

# Function to check if BLAST database exists, and create it if necessary
create_blast_db <- function(db_file, db_name) {
  if (!file.exists(paste0(db_name, ".nin"))) {  # Check if BLAST DB index exists
    cat("BLAST database not found. Creating database...\n")
    
    # Run makeblastdb command
    system(paste("makeblastdb -in", db_file, "-dbtype nucl -out", db_name))
    
    cat("BLAST database created.\n")
  } else {
    cat("BLAST database already exists. Skipping database creation.\n")
  }
}

# Function to run BLASTN
run_blastn <- function(fasta_file, db_name) {
  cat("Running BLASTN for", fasta_file, "...\n")
  
  # Define output file name
  output_file <- sub("\\.fasta$", "_blast_results.txt", fasta_file)
  
  # Construct BLAST command
  blast_cmd <- paste("blastn -query", fasta_file, "-db", db_name, 
                     '-outfmt "6 qseqid sseqid evalue"', "-out", output_file)
  
  # Run BLAST
  system(blast_cmd)
  
  cat("BLAST completed. Results saved in", output_file, "\n")
  return(output_file)  # Return BLAST results file path
}

# Step 3: Ensure BLAST database is ready
create_blast_db(db_file, db_name)

# Step 4: Loop through all FASTQ files
blast_results_files <- c()  # To store paths of BLAST results files
for (fastq_file in fastq_files) {
  # Extract the GeneID from the file name (e.g., "Tb927.1.1500" from "Tb927.1.1500.fastq")
  gene_id <- sub("\\.fastq$", "", basename(fastq_file))
  
  # Convert FASTQ to FASTA
  fasta_file <- convert_fastq_to_fasta(fastq_file)
  
  # Run BLASTN on the converted file
  blast_results_file <- run_blastn(fasta_file, db_name)
  
  # Store the BLAST results file
  blast_results_files <- c(blast_results_files, blast_results_file)
}

cat("Processing completed for all FASTQ files!\n")

# Step 5: Process BLAST results
# Step 5.1: Parse BLAST tabular results (-outfmt 6)
parse_blast_tabular <- function(blast_tabular_file) {
  # Read the tabular file and assign column names
  blast_results <- read_tsv(
    blast_tabular_file,
    col_names = c("Query", "Prey", "E_value"),  # Explicitly assign column names
    col_types = cols(
      Query = col_character(),
      Prey = col_character(),
      E_value = col_double()
    )
  )
  
  # Convert E_value to numeric (in case it's read as character)
  blast_results <- blast_results %>%
    mutate(E_value = as.numeric(E_value))
  
  return(blast_results)
}

# Step 5.2: Retain only the prey with the lowest E-value for each bait
filter_lowest_evalue <- function(blast_df) {
  blast_df %>%
    group_by(Query) %>%
    filter(E_value == min(E_value)) %>%  # Filter for the row with the lowest E-value per Query
    ungroup()
}

# Step 5.3: Load contaminants list (if present)
contaminants_file <- list.files(path = ".", pattern = "contaminants.txt", full.names = TRUE)
if (length(contaminants_file) > 0) {
  contaminants <- read_lines(contaminants_file)
} else {
  contaminants <- character(0)  # Empty vector if no contaminants file is found
}

# Step 5.4: Filter out contaminants
filter_contaminants <- function(blast_df, contaminants) {
  blast_df %>%
    filter(!Prey %in% contaminants)
}

# Step 5.5: Process BLAST results
final_interactors_list <- list()  # To store final interactors for each bait
for (blast_results_file in blast_results_files) {
  print(paste("Processing file:", blast_results_file))
  
  # Extract the GeneID from the file name (e.g., "Tb927.1.1500" from "Tb927.1.1500_blast_results.txt")
  gene_id <- sub("_blast_results\\.txt$", "", basename(blast_results_file))
  
  # Parse the BLAST results for the current file
  blast_results <- parse_blast_tabular(blast_results_file)
  
  # Add the GeneID to the results
  blast_results <- blast_results %>%
    mutate(Bait = gene_id)
  
  # Filter out contaminants (if contaminants file is present)
  if (length(contaminants) > 0) {
    blast_results <- filter_contaminants(blast_results, contaminants)
  }
  
  # Retain only the prey with the lowest E-value for each Query-Prey pair
  final_interactors <- filter_lowest_evalue(blast_results)
  
  # Store the final interactors for this bait
  final_interactors_list[[blast_results_file]] <- final_interactors
  
  print(paste("Processed bait:", gene_id))
}

# Step 5.6: Combine all final interactors into a single data frame and ensure only lowest E-value preys are kept
all_final_interactors <- bind_rows(final_interactors_list) %>%
  group_by(Bait, Prey) %>%  # Group by both Bait and Prey
  filter(E_value == min(E_value)) %>%  # Keep only the lowest E-value for each Bait-Prey pair
  ungroup() %>%
  distinct(Bait, Prey, .keep_all = TRUE)  # Ensure no duplicates remain

# Step 5.7: Print the combined data frame for debugging
print("Combined final interactors:")
print(head(all_final_interactors))

# Step 5.8: Save the combined results to a CSV file
write_csv(all_final_interactors, "all_final_interactors.csv")

print("All final interactors saved to all_final_interactors.csv")

# Step 6: Shiny App for Visualization
ui <- fluidPage(
  titlePanel("InteractR"),
  sidebarLayout(
    sidebarPanel(
      uiOutput("baitSelector"),  # Dynamic UI for bait selection
      hr(),
      helpText("Select one or more baits to visualize their interactions."),
      uiOutput("uniquePreyDownload"),  # Download button for unique preys
      actionButton("showScatterPlot", "Show Scatter Plot")  # Button to show scatter plot
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Upset Plot", 
                 plotOutput("upsetPlot"),
                 downloadButton("downloadUpsetPlot", "Download UpSet Plot")),
        tabPanel("Protein Level", 
                 plotOutput("scatterPlot"),
                 downloadButton("downloadScatterPlotPDF", "Download as PDF")),
        tabPanel("Subcellular Localization", 
                 plotOutput("preyLocalizationsBarPlot"),
                 downloadButton("downloadLocalizationPlotPDF", "Download as PDF")),
        tabPanel("InterProScan",
                 plotOutput("interproBarPlot"),
                 downloadButton("downloadInterproPlotPDF", "Download as PDF"),
                 DTOutput("interproTable"),
                 downloadButton("downloadInterproTableCSV", "Download Table as CSV"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive function to read and process the data directly from files
  data <- reactive({
    # Read the BLAST results (generated earlier in the script)
    blast_results <- read_csv("all_final_interactors.csv", show_col_types = FALSE) %>%
      group_by(Bait, Prey) %>%
      filter(E_value == min(E_value)) %>%
      ungroup() %>%
      distinct(Bait, Prey, .keep_all = TRUE)
    
    # Read the localization data from localization2.csv
    localizations <- read_csv("localization2.csv", col_names = TRUE, show_col_types = FALSE)
    
    # Combine genes with NA or Unknown location
    localizations <- localizations %>%
      mutate(Localization = ifelse(is.na(Localization) | Localization == "Unknown", "Unknown", Localization))
    
    # Define the total localizations from the provided list
    total_localizations <- c(
      "Unknown" = 6087,
      "secretory/endocytic 2" = 34,
      "nucleus" = 378,
      "mitochondrion - IM" = 252,
      "endoplasmic reticulum" = 304,
      "secretory/endocytic 1" = 219,
      "secretory/endocytic 3" = 124,
      "mitochondrion - matrix 1" = 153,
      "mitochondrion - matrix 2" = 294,
      "microtubule structures 1" = 235,
      "nucleus - chromatin" = 80,
      "flagellum 1" = 489,
      "glycosomes" = 62,
      "microtubule structures 2" = 132,
      "cytosol" = 424,
      "acidocalcisomes" = 17,
      "ribonucleoproteins 1" = 201,
      "mitochondrion - OM" = 41,
      "proteasome regulatory particle" = 18,
      "intraflagellar transport" = 11,
      "proteasome" = 14
    )
    
    # Merge localization data with BLAST results
    merged_data <- blast_results %>%
      left_join(localizations, by = c("Prey" = "Gene ID"), relationship = "many-to-many") %>%
      mutate(Localization = ifelse(is.na(Localization), "Unknown", Localization))
    
    # Create a binary matrix for UpSet plot
    binary_matrix <- blast_results %>%
      select(Bait, Prey) %>%
      distinct() %>%
      mutate(Value = 1) %>%
      pivot_wider(names_from = Bait, values_from = Value, values_fill = list(Value = 0)) %>%
      as.data.frame()
    
    # Set row names for the binary matrix
    row.names(binary_matrix) <- binary_matrix$Prey
    binary_matrix <- binary_matrix %>% select(-Prey)
    
    # Return the binary matrix, merged_data, list of baits, and total localizations
    list(binary_matrix = binary_matrix, merged_data = merged_data, baits = colnames(binary_matrix), total_localizations = total_localizations)
  })
  
  # Read the InterproScan data with empty value filtering
  interpro_data <- reactive({
    req(file.exists("filtered_results.csv"))
    read_csv("filtered_results.csv", show_col_types = FALSE) %>%
      select(`Gene ID`, Definition) %>%
      filter(!is.na(Definition),
             Definition != "",
             !Definition %in% c("NA", "-", "empty")) %>%
      group_by(`Gene ID`) %>%
      summarise(Definitions = paste(unique(Definition), collapse = "; ")) %>%
      ungroup()
  })
  
  # Dynamic UI for bait selection
  output$baitSelector <- renderUI({
    req(data())
    baits <- data()$baits
    selectInput("selectedBaits", "Select Baits", choices = baits, multiple = TRUE)
  })
  
  # Generate the UpSet plot with single-bait handling
  output$upsetPlot <- renderPlot({
    req(data(), input$selectedBaits)
    binary_matrix <- data()$binary_matrix
    if (!is.null(input$selectedBaits)) {
      binary_matrix <- binary_matrix %>% select(all_of(input$selectedBaits))
    }
    if (ncol(binary_matrix) < 2) {
      # Plot a placeholder if only one bait is selected
      ggplot() +
        annotate("text", x = 1, y = 1, label = "Select at least two baits for UpSet plot", size = 6) +
        theme_void()
    } else {
      upset(binary_matrix, nsets = ncol(binary_matrix), nintersects = 20, mb.ratio = c(0.6, 0.4),
            order.by = "freq", decreasing = TRUE, sets = colnames(binary_matrix),
            keep.order = TRUE, text.scale = 1.5)
    }
  })
  
  # Updated download handler for UpSet plot
  output$downloadUpsetPlot <- downloadHandler(
    filename = function() {
      "upset_plot.pdf"
    },
    content = function(file) {
      req(data(), input$selectedBaits)
      binary_matrix <- data()$binary_matrix
      
      if (!is.null(input$selectedBaits)) {
        binary_matrix <- binary_matrix %>% select(all_of(input$selectedBaits))
      }
      
      if (ncol(binary_matrix) < 2) {
        # Save placeholder ggplot for single bait
        placeholder_plot <- ggplot() +
          annotate("text", x = 1, y = 1, label = "Select at least two baits for UpSet plot", size = 6) +
          theme_void()
        ggsave(file, plot = placeholder_plot, device = "pdf", width = 8, height = 6)
      } else {
        # Create a temporary file for the plot
        temp_plot <- tempfile(fileext = ".pdf")
        
        # Save the UpSet plot properly
        pdf(temp_plot, width = 8, height = 6)
        print(upset(binary_matrix, 
                    nsets = ncol(binary_matrix), 
                    nintersects = 20, 
                    mb.ratio = c(0.6, 0.4),
                    order.by = "freq", 
                    decreasing = TRUE, 
                    sets = colnames(binary_matrix),
                    keep.order = TRUE, 
                    text.scale = 1.5))
        dev.off()
        
        # Copy the temporary file to the download location
        file.copy(temp_plot, file)
        
        # Remove the temporary file
        unlink(temp_plot)
      }
    }
  )
  
  # Generate the bar plot for prey localizations with y-axis in percentage
  output$preyLocalizationsBarPlot <- renderPlot({
    req(data(), input$selectedBaits)
    merged_data <- data()$merged_data
    
    # Filter data for selected baits and remove duplicate preys
    filtered_data <- merged_data %>%
      filter(Bait %in% input$selectedBaits) %>%
      distinct(Bait, Prey, Localization)  # Ensure each prey is counted only once per bait
    
    # Count the frequency of each localization
    localization_counts <- filtered_data %>%
      group_by(Localization) %>%
      summarise(Count = n(), .groups = 'drop')
    
    # Calculate the percentage for each localization category
    localization_counts <- localization_counts %>%
      mutate(Percentage = (Count / data()$total_localizations[Localization]) * 100)
    
    # Filter out NA and Unknown for the plot (optional)
    plot_data <- localization_counts %>%
      filter(!is.na(Localization) & Localization != "Unknown")
    
    # Generate the bar plot with y-axis as percentage
    ggplot(plot_data, aes(x = reorder(Localization, -Percentage), y = Percentage, fill = Localization)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = paste0(round(Percentage, 2), "%")), vjust = -0.5, size = 4) +
      labs(title = "Prey Subcellular Localizations", 
           x = "Location", 
           y = "Percentage (%)",
           caption = "Data from Moloney et al. Nat Commun 14, 4401 (2023)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.caption = element_text(size = 8)) +
      scale_fill_discrete()
  })
  
  # Create a bar plot of top enriched InterProScan features
  output$interproBarPlot <- renderPlot({
    req(data(), input$selectedBaits, interpro_data())
    
    # Get preys for selected baits
    merged_data <- data()$merged_data
    selected_preys <- merged_data %>%
      filter(Bait %in% input$selectedBaits) %>%
      pull(Prey) %>%
      unique()
    
    # Filter InterproScan data and remove empty/NA values
    domain_counts <- interpro_data() %>%
      filter(`Gene ID` %in% selected_preys) %>%
      separate_rows(Definitions, sep = "; ") %>%
      filter(Definitions != "") %>%
      group_by(Definitions) %>%
      summarise(Unique_Preys = n_distinct(`Gene ID`)) %>%
      arrange(desc(Unique_Preys)) %>%
      head(20)
    
    # Plot top domains by unique prey count
    if(nrow(domain_counts) > 0) {
      ggplot(domain_counts, aes(x = reorder(Definitions, Unique_Preys), y = Unique_Preys)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = "Top InterProScan Domains by Unique Prey Count",
             x = "Domain",
             y = "Number of Unique Preys") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 10))
    } else {
      ggplot() + 
        annotate("text", x = 1, y = 1, label = "No significant domains found", size = 6) + 
        theme_void()
    }
  })
  
  # Create the InterProScan definitions table
  output$interproTable <- renderDT({
    req(data(), input$selectedBaits, interpro_data())
    
    # Get preys for selected baits
    merged_data <- data()$merged_data
    selected_preys <- merged_data %>%
      filter(Bait %in% input$selectedBaits) %>%
      pull(Prey) %>%
      unique()
    
    # Process InterproScan data for table
    interpro_table_data <- interpro_data() %>%
      filter(`Gene ID` %in% selected_preys) %>%
      rename(Prey = `Gene ID`, Domains = Definitions) %>%
      left_join(merged_data %>% select(Prey, Bait, E_value), by = "Prey") %>%
      select(Bait, Prey, E_value, Domains) %>%
      arrange(Bait, Prey)
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # Add download handler for Subcellular Localization plot
  output$downloadLocalizationPlotPDF <- downloadHandler(
    filename = function() {
      "localization_plot.pdf"
    },
    content = function(file) {
      merged_data <- data()$merged_data
      
      # Filter data for selected baits and remove duplicate preys
      filtered_data <- merged_data %>%
        filter(Bait %in% input$selectedBaits) %>%
        distinct(Bait, Prey, Localization)
      
      # Count the frequency of each localization
      localization_counts <- filtered_data %>%
        group_by(Localization) %>%
        summarise(Count = n(), .groups = 'drop')
      
      # Calculate the percentage for each localization category
      localization_counts <- localization_counts %>%
        mutate(Percentage = (Count / data()$total_localizations[Localization]) * 100)
      
      # Filter out NA and Unknown for the plot
      plot_data <- localization_counts %>%
        filter(!is.na(Localization) & Localization != "Unknown")
      
      # Generate the bar plot
      localization_plot <- ggplot(plot_data, aes(x = reorder(Localization, -Percentage), y = Percentage, fill = Localization)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(Percentage, 2), "%")), vjust = -0.5, size = 4) +
        labs(title = "Prey Subcellular Localizations", 
             x = "Location", 
             y = "Percentage (%)",
             caption = "Data from Moloney et al. Nat Commun 14, 4401 (2023)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none",
              plot.caption = element_text(size = 8)) +
        scale_fill_discrete()
      
      ggsave(file, plot = localization_plot, device = "pdf", width = 8, height = 6)
    }
  )
  
  # Add download handler for InterProScan plot
  output$downloadInterproPlotPDF <- downloadHandler(
    filename = function() {
      "interpro_plot.pdf"
    },
    content = function(file) {
      merged_data <- data()$merged_data
      selected_preys <- merged_data %>%
        filter(Bait %in% input$selectedBaits) %>%
        pull(Prey) %>%
        unique()
      
      # Filter InterproScan data and remove empty/NA values
      domain_counts <- interpro_data() %>%
        filter(`Gene ID` %in% selected_preys) %>%
        separate_rows(Definitions, sep = "; ") %>%
        filter(Definitions != "") %>%
        group_by(Definitions) %>%
        summarise(Unique_Preys = n_distinct(`Gene ID`)) %>%
        arrange(desc(Unique_Preys)) %>%
        head(20)
      
      # Generate the bar plot
      interpro_plot <- if(nrow(domain_counts) > 0) {
        ggplot(domain_counts, aes(x = reorder(Definitions, Unique_Preys), y = Unique_Preys)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          labs(title = "Top InterProScan Domains by Unique Prey Count",
               x = "Domain",
               y = "Number of Unique Preys") +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 10))
      } else {
        ggplot() + 
          annotate("text", x = 1, y = 1, label = "No significant domains found", size = 6) + 
          theme_void()
      }
      
      ggsave(file, plot = interpro_plot, device = "pdf", width = 8, height = 6)
    }
  )
  
  # Add download handler for InterProScan table
  output$downloadInterproTableCSV <- downloadHandler(
    filename = function() {
      "interpro_table.csv"
    },
    content = function(file) {
      merged_data <- data()$merged_data
      selected_preys <- merged_data %>%
        filter(Bait %in% input$selectedBaits) %>%
        pull(Prey) %>%
        unique()
      
      # Process InterproScan data for table
      table_data <- interpro_data() %>%
        filter(`Gene ID` %in% selected_preys) %>%
        rename(Prey = `Gene ID`, Domains = Definitions) %>%
        left_join(merged_data %>% select(Prey, Bait, E_value), by = "Prey") %>%
        select(Bait, Prey, E_value, Domains) %>%
        arrange(Bait, Prey)
      
      write_csv(table_data, file)
    }
  )
  
  # Observe the scatter plot button click
  observeEvent(input$showScatterPlot, {
    req(input$selectedBaits)  # Ensure a bait is selected
    
    tryCatch({
      # Read the protein level data
      if (!file.exists("protein_level.csv")) {
        stop("The file 'protein_level.csv' does not exist in the working directory.")
      }
      
      protein_data <- read.csv("protein_level.csv", header = TRUE)
      
      # Check if the required columns exist
      if (!all(c("GeneID", "bsf", "pcf") %in% colnames(protein_data))) {
        stop("The file 'protein_level.csv' must contain the columns 'GeneID', 'bsf', and 'pcf'.")
      }
      
      # Get the selected bait
      selected_bait <- input$selectedBaits[1]  # Get the selected bait
      
      # Get the specific preys for the selected bait
      binary_matrix <- data()$binary_matrix
      specific_preys <- if (ncol(binary_matrix) >= 2) {
        rownames(binary_matrix)[binary_matrix[[selected_bait]] == 1 & rowSums(binary_matrix) == 1]
      } else {
        rownames(binary_matrix)[binary_matrix[[selected_bait]] == 1]
      }
      
      # Add a new column to indicate the color group
      protein_data <- protein_data %>%
        mutate(ColorGroup = case_when(
          GeneID == selected_bait ~ "Bait",
          GeneID %in% specific_preys ~ "Specific Preys",
          TRUE ~ "Rest of Population"
        ))
      
      # Create the scatter plot
      scatter_plot <- ggplot() +
        geom_point(
          data = filter(protein_data, ColorGroup == "Rest of Population"),
          aes(x = bsf, y = pcf),
          color = "grey",
          size = 3,
          alpha = 0.7
        ) +
        geom_point(
          data = filter(protein_data, ColorGroup == "Specific Preys"),
          aes(x = bsf, y = pcf),
          color = "blue",
          size = 3,
          alpha = 0.7
        ) +
        geom_point(
          data = filter(protein_data, ColorGroup == "Bait"),
          aes(x = bsf, y = pcf),
          color = "red",
          size = 5,
          alpha = 1
        ) +
        scale_x_continuous(limits = c(2.5, 10)) +
        scale_y_continuous(limits = c(2.5, 10)) +
        labs(
          title = paste("Scatter Plot: PCF vs. BSF for", selected_bait),
          x = "BSF",
          y = "PCF",
          caption = "Data from Tinti et al. Wellcome Open Res (2023)"
        ) +
        theme_minimal() +
        theme(plot.caption = element_text(size = 8))
      
      # Render the scatter plot
      output$scatterPlot <- renderPlot(scatter_plot)
      
      # Download handler for the scatter plot as PDF
      output$downloadScatterPlotPDF <- downloadHandler(
        filename = function() {
          paste0("scatter_plot_", selected_bait, ".pdf")
        },
        content = function(file) {
          ggsave(file, plot = scatter_plot, device = "pdf", width = 8, height = 6)
        }
      )
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred while generating the scatter plot:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # Download handler for unique preys with domain information (as CSV)
  output$uniquePreyDownload <- renderUI({
    req(data(), input$selectedBaits)
    
    lapply(input$selectedBaits, function(bait) {
      downloadButton(paste0("download_unique_preys_", bait), label = paste("Download Unique Preys for", bait))
    })
  })
  
  observe({
    req(data(), input$selectedBaits, interpro_data())
    
    # Extract the binary matrix
    binary_matrix <- data()$binary_matrix
    
    # Filter the binary matrix to include only selected baits
    if (!is.null(input$selectedBaits)) {
      binary_matrix <- binary_matrix %>% select(all_of(input$selectedBaits))
    }
    
    # Create download handlers for unique preys per bait
    lapply(input$selectedBaits, function(bait) {
      output[[paste0("download_unique_preys_", bait)]] <- downloadHandler(
        filename = function() {
          paste0(bait, "_unique_preys.csv")
        },
        content = function(file) {
          # Identify unique preys for the current bait
          unique_preys <- if (ncol(binary_matrix) >= 2) {
            rownames(binary_matrix)[binary_matrix[[bait]] == 1 & rowSums(binary_matrix) == 1]
          } else {
            rownames(binary_matrix)[binary_matrix[[bait]] == 1]
          }
          
          # Combine all data for the unique preys of this bait
          table_data <- data()$merged_data %>%
            filter(Bait == bait, Prey %in% unique_preys) %>%
            left_join(
              read_csv("localization2.csv", col_names = TRUE, show_col_types = FALSE) %>% 
                rename(Prey = `Gene ID`),
              by = "Prey",
              relationship = "many-to-many",
              suffix = c("", ".loc")  # Avoid column overwrite
            ) %>%
            left_join(
              interpro_data() %>% rename(Prey = `Gene ID`),
              by = "Prey",
              suffix = c("", ".ipro")  # Avoid column overwrite
            ) %>%
            select(Bait, Prey, E_value, Localization, `Product Description`, Category, CtermTryptag, NtermTryptag, Definitions) %>%
            rename(InterProScan_Domains = Definitions) %>%
            distinct() %>%
            arrange(Prey)
          
          # Write to CSV
          write_csv(table_data, file)
        }
      )
    })
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)