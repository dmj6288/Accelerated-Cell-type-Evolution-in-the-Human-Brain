rm(list = ls())

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(readxl)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressMessages(suppressWarnings(library(rstatix)))
suppressMessages(suppressWarnings(library(ggprism)))
suppressMessages(suppressWarnings(library(ggVennDiagram)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(VennDiagram)))

# Your vector
my_vector <- c("A", "B", "C", "D") # Example vector

# Function to get all combinations of size greater than 2
get_combinations <- function(vec) {
  results <- list()
  n <- length(vec)
  
  for (i in 2:n) {
    combos <- combn(vec, i, simplify = FALSE)
    results <- c(results, combos)
  }
  
  return(results)
}

assign_names_with_serial <- function(lst) {
  # Using seq_along to iterate with indexes and sapply for operation
  names(lst) <- sapply(seq_along(lst), function(i) {
    vec <- lst[[i]]
    paste("Serial", i, "Length", length(vec), sep = "_")
  })
  return(lst)
}

# Function to merge dataframes by name across sublists
aggregate_dataframes <- function(sublists) {
  # Identifying unique dataframe names across all sublists
  df_names <- unique(unlist(lapply(sublists, names)))
  
  # Aggregating dataframes by name
  aggregated_dfs <- lapply(df_names, function(df_name) {
    # Extracting the dataframe with the same name across all sublists
    dfs_to_merge <- lapply(sublists, function(sublist) sublist[[df_name]])
    # Combine all dataframes row-wise assuming they have the same columns
    combined_df <- do.call(rbind, dfs_to_merge)
    return(combined_df)
  })
  names(aggregated_dfs) <- df_names # Assigning the correct names to the aggregated dataframes
  return(aggregated_dfs)
}

# Function to calculate averages across numeric columns in a list of dataframes
calculate_averages_numeric_only <- function(df_list) {
  lapply(df_list, function(df) {
    # Identify numeric columns
    numeric_cols <- sapply(df, is.numeric)
    # Calculate averages only for numeric columns
    if (any(numeric_cols)) {
      numeric_df <- df[, numeric_cols, drop = FALSE]
      averages <- colMeans(df[, numeric_cols], na.rm = TRUE)
      return(averages)
    } else {                                     
      return(list())
    }
  })
}

calculate_average_df <- function(df_name, complex_list) {
  
  
  # Extracting the dataframes by name from each sublist
  dfs <- lapply(complex_list, function(x) x[[df_name]])
  dfs <- dfs[!sapply(dfs, is.null)]  # Remove NULL entries if any sublist doesn't have the df
  
  # Assuming all dataframes have the same structure
  # Initialize an empty dataframe with the same dimensions and column names
  average_df <- dfs[[1]]
  average_df[,] <- 0  # Set all values to 0 as a starting point
  
  # Calculate averages
  for (i in seq_along(average_df)) {
    for (j in seq_len(ncol(average_df))) {
      # Extracting the i,jth element across all dfs, calculating mean
      values_to_average <- sapply(dfs, function(df) as.numeric(df[i, j]))
      average_df[i, j] <- mean(values_to_average, na.rm = TRUE)
    }
  }
  
  return(average_df)
}

# Helper function to calculate mean for common columns in dataframes with the same name
calculate_mean_for_common_columns <- function(dfs) {
  # Identify common columns
  common_cols <- Reduce(intersect, lapply(dfs, colnames))
  
  # Extract these columns and convert to numeric if possible
  dfs_common <- lapply(dfs, function(df) {
    sapply(df[common_cols], function(column) {
      suppressWarnings(as.numeric(as.character(column)))
    })
  })
  
  # Calculate mean for each common column
  mean_df <- data.frame(lapply(seq_along(dfs_common[[1]]), function(i) {
    col_means <- sapply(dfs_common, function(df) df[, i])
    mean(col_means, na.rm = TRUE)
  }))
  
  # Set the column names to the common columns
  names(mean_df) <- common_cols
  
  return(mean_df)
}

# Step 1: Extract dataframes with the same name from each sublist
extract_dfs <- function(main_list, df_name) {
  do.call(rbind, lapply(main_list, function(sublist) {
    if (!is.null(sublist[[df_name]])) {
      # Assuming the dataframes may need to be transformed for consistent structure
      as.data.frame(lapply(sublist[[df_name]], function(x) as.numeric(as.character(x))))
    }
  }))
}


calculate_averages_by_length <- function(df) {
  # Extracting the unique 'Length_#num' parts from the column names
  length_groups <- unique(gsub(".*_(Length_\\d+)_.*", "\\1", names(df)))
  
  # Initialize an empty list to store the average results
  averages_by_length <- list()
  
  # Loop over each length group to calculate averages
  for (length_group in length_groups) {
    # Find columns that belong to the current length group
    cols_in_group <- grep(length_group, names(df), value = TRUE)
    
    # Subset the dataframe to only include columns from the current length group
    df_subset <- df[, cols_in_group, drop = FALSE]
    
    # Calculate the average for each row in the subset
    averages_by_length[[length_group]] <- rowMeans(df_subset, na.rm = TRUE)
  }
  
  # Convert the list to a dataframe for easier viewing
  averages_df <- data.frame(averages_by_length)
  
  # Include the row names as a column in the averages dataframe
  averages_df$Description <- rownames(df)
  
  # Return the result
  return(averages_df)
}

# Example usage:
# Assuming 'df' is your original dataframe
# averages_df <- calculate_averages_by_length(df)
# View the result
# print(averages_df)


studies            <- c("Caglayan", "Ma")
species            <- c("HUMAN", "CHIMP")
regions            <- c("PCC", "dlPFC")
names(regions)     <- studies
regulation         <- c("UR", "DR")
FDR                <- 0.05

cell_names_of_interest_Excite  <- c("L2-3 IT",      "L3-5 IT-1",    "L3-5 IT-2",    "L3-5 IT-3",    "L5 ET",        
                                    "L5-6 NP",      "L6 CT",        "L6 IT-1",      "L6 IT-2",      "L6B")
excite_cell_names_combinations <- get_combinations(cell_names_of_interest_Excite)
named_excite_cell_names_combinations <- assign_names_with_serial(excite_cell_names_combinations)

cell_names_of_interest_Inhibit <- c("LAMP5 LHX6",   "LAMP5 RELN", "PVALB ChC",    "PVALB",  "SST HGF",      "SST NPY",      "SST")
inhibit_cell_names_combinations <- get_combinations(cell_names_of_interest_Inhibit)
named_inhibit_cell_names_combinations <- assign_names_with_serial(inhibit_cell_names_combinations)

cell_names_of_interest         <- list(cell_names_of_interest_Excite,
                                       cell_names_of_interest_Inhibit)

combination_list               <- list(named_excite_cell_names_combinations, 
                                       named_inhibit_cell_names_combinations)

names(cell_names_of_interest)  <- c("Excite", "Inhibit")
names(combination_list)        <- c("Excite", "Inhibit")

subtype_folder_name   <- "../../Cloud Data/subtype/"
majortype_folder_name <- "../../Cloud Data/major/"



