#List all the packages available
library()

#Install the Diagnostics package
install.packages("/Users/marcellovendruscolo/Documents/rstudio-workspace/Diagnostics/Diagnostics_1.2.0.tar.gz", repos = NULL, type = "source")

#Load Diagnostics package for this session
library(Diagnostics)

#List all the packages currently loaded
search()

#Function to create and learn the Bayesian network. The only parameter is the historical medical data. It returns the learned network as a list of matrices
learn <- function(historical_data) {
  total_observations <- nrow(historical_data)
  network <- list()

  pneumonia_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(pneumonia_matrix) <- c('Pn == 0', 'Pn == 1')
  
  pneumonia_df <- historical_data['Pn']
  pneumonia_counts <- nrow(subset(pneumonia_df, Pn == 1)) + 1
  p_pneumonia <- (pneumonia_counts / total_observations)
  
  pneumonia_matrix[1,1] <- (1 - p_pneumonia)
  pneumonia_matrix[1,2] <- (p_pneumonia)
  
  network$pneumonia <- as.data.frame(pneumonia_matrix)
  
  #print('network$pneumonia')
  #print(network$pneumonia)
  
  temperature_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(temperature_matrix) <- c('Pn', 'Mean', 'sd')
  temperature_matrix[1,1] <- 0
  temperature_matrix[2,1] <- 1
  
  temperature_df <- historical_data['Te']
  temperature_pneumonia_df <- data.frame(pneumonia_df, temperature_df)
  temperature_vector_pneumonia0 <- c()
  temperature_vector_pneumonia1 <- c()
  
  for (counter in 1:total_observations) {
    if (temperature_pneumonia_df$Pn[[counter]] == 0) {
      temperature_vector_pneumonia0[length(temperature_vector_pneumonia0) + 1] <- temperature_pneumonia_df$Te[[counter]]
    } else {
      temperature_vector_pneumonia1[length(temperature_vector_pneumonia1) + 1] <- temperature_pneumonia_df$Te[[counter]]
    }
  }
  
  mean_temperature_pneumonia0 <- mean(temperature_vector_pneumonia0)
  sd_temperature_pneumonia0 <- sd(temperature_vector_pneumonia0)
  mean_temperature_pneumonia1 <- mean(temperature_vector_pneumonia1)
  sd_temperature_pneumonia1 <- sd(temperature_vector_pneumonia1)
  
  temperature_matrix[1,2] <- mean_temperature_pneumonia0
  temperature_matrix[1,3] <- sd_temperature_pneumonia0
  temperature_matrix[2,2] <- mean_temperature_pneumonia1
  temperature_matrix[2,3] <- sd_temperature_pneumonia1
  
  network$temperature <- as.data.frame(temperature_matrix)
  
  #print('network$temperature')
  #print(network$temperature)
  
  visitedTBSpot_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(visitedTBSpot_matrix) <- c('VTB == 0', 'VTB == 1')
  
  visitedTBSpot_df <- historical_data['VTB']
  visitedTBSpot_counts <- nrow(subset(visitedTBSpot_df, VTB == 1)) + 1
  p_visitedTBSpot <- (visitedTBSpot_counts / total_observations)
  
  visitedTBSpot_matrix[1,1] <- (1 - p_visitedTBSpot)
  visitedTBSpot_matrix[1,2] <- (p_visitedTBSpot)
  
  network$visitedTBSpot <- as.data.frame(visitedTBSpot_matrix)
  
  #print('network$visitedTBSpot')
  #print(network$visitedTBSpot)
  
  tuberculosis_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(tuberculosis_matrix) <- c('VTB', 'TB == 0', 'TB == 1')
  tuberculosis_matrix[1,1] <- 0
  tuberculosis_matrix[2,1] <- 1
  
  tuberculosis_df <- historical_data['TB']
  tuberculosis_visitedTBSpot_df <- data.frame(visitedTBSpot_df, tuberculosis_df)
  tuberculosis_VTBSpot0_counts <- nrow(subset(tuberculosis_visitedTBSpot_df, VTB == 0 & TB == 1)) + 1
  tuberculosis_VTBSpot1_counts <- nrow(subset(tuberculosis_visitedTBSpot_df, VTB == 1 & TB == 1)) + 1
  p_tuberculosis_VTBSpot0 <- (tuberculosis_VTBSpot0_counts / ((total_observations - visitedTBSpot_counts) + 2))
  p_tuberculosis_VTBSpot1 <- (tuberculosis_VTBSpot1_counts / (visitedTBSpot_counts + 2))
  
  tuberculosis_matrix[1,2] <- (1 - p_tuberculosis_VTBSpot0)
  tuberculosis_matrix[1,3] <- p_tuberculosis_VTBSpot0
  tuberculosis_matrix[2,2] <- (1 - p_tuberculosis_VTBSpot1)
  tuberculosis_matrix[2,3] <- p_tuberculosis_VTBSpot1
  
  network$tuberculosis <- as.data.frame(tuberculosis_matrix)
  
  #print('network$tuberculosis')
  #print(network$tuberculosis)
  
  smokes_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(smokes_matrix) <- c('Sm == 0', 'Sm == 1')
  
  smokes_df <- historical_data['Sm']
  smokes_counts <- nrow(subset(smokes_df, Sm == 1)) + 1
  p_smokes <- (smokes_counts / total_observations)
  
  smokes_matrix[1,1] <- (1 - p_smokes)
  smokes_matrix[1,2] <- (p_smokes)
  
  network$smokes <- as.data.frame(smokes_matrix)
  
  #print('network$smokes')
  #print(network$smokes)
  
  lungCancer_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(lungCancer_matrix) <- c('Sm', 'LC == 0', 'LC == 1')
  lungCancer_matrix[1,1] <- 0
  lungCancer_matrix[2,1] <- 1
  
  lungCancer_df <- historical_data['LC']
  lungCancer_smokes_df <- data.frame(smokes_df, lungCancer_df)
  lungCancer_Sm0_counts <- nrow(subset(lungCancer_smokes_df, Sm == 0 & LC == 1)) + 1
  lungCancer_Sm1_counts <- nrow(subset(lungCancer_smokes_df, Sm == 1 & LC == 1)) + 1
  p_lungCancer_Sm0 <- (lungCancer_Sm0_counts / ((total_observations - smokes_counts) + 2))
  p_lungCancer_Sm1 <- (lungCancer_Sm1_counts / (smokes_counts + 2))
  
  lungCancer_matrix[1,2] <- (1 - p_lungCancer_Sm0)
  lungCancer_matrix[1,3] <- p_lungCancer_Sm0
  lungCancer_matrix[2,2] <- (1 - p_lungCancer_Sm1)
  lungCancer_matrix[2,3] <- p_lungCancer_Sm1
  
  network$lungCancer <- as.data.frame(lungCancer_matrix)
  
  #print('network$lungCancer')
  #print(network$lungCancer)
  
  bronchitis_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(bronchitis_matrix) <- c('Sm', 'Br == 0', 'Br == 1')
  bronchitis_matrix[1,1] <- 0
  bronchitis_matrix[2,1] <- 1
  
  bronchitis_df <- historical_data['Br']
  bronchitis_smokes_df <- data.frame(smokes_df, bronchitis_df)
  bronchitis_Sm0_counts <- nrow(subset(bronchitis_smokes_df, Sm == 0 & Br == 1)) + 1
  bronchitis_Sm1_counts <- nrow(subset(bronchitis_smokes_df, Sm == 1 & Br == 1)) + 1
  p_bronchitis_Sm0 <- (bronchitis_Sm0_counts / ((total_observations - smokes_counts) + 2))
  p_bronchitis_Sm1 <- (bronchitis_Sm1_counts / (smokes_counts + 2))
  
  bronchitis_matrix[1,2] <- (1 - p_bronchitis_Sm0)
  bronchitis_matrix[1,3] <- p_bronchitis_Sm0
  bronchitis_matrix[2,2] <- (1 - p_bronchitis_Sm1)
  bronchitis_matrix[2,3] <- p_bronchitis_Sm1
  
  network$bronchitis <- as.data.frame(bronchitis_matrix)
  
  #print('network$bronchitis')
  #print(network$bronchitis)
  
  dyspnea_matrix <- matrix(NA, nrow = 4, ncol = 4, byrow = TRUE)
  colnames(dyspnea_matrix) <- c('LC', 'Br', 'Dy == 0', 'Dy == 1')
  dyspnea_matrix[1,1] <- 0
  dyspnea_matrix[2,1] <- 0
  dyspnea_matrix[1,2] <- 0
  dyspnea_matrix[3,2] <- 0
  dyspnea_matrix[3,1] <- 1
  dyspnea_matrix[4,1] <- 1
  dyspnea_matrix[2,2] <- 1
  dyspnea_matrix[4,2] <- 1
  
  dyspnea_df <- historical_data['Dy']
  dyspnea_lungCancer_Bronchitis_df <- data.frame(lungCancer_df, bronchitis_df, dyspnea_df)
  dyspnea_LC0_Br0_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 0 & Br == 0 & Dy == 1)) + 1
  dyspnea_LC0_Br1_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 0 & Br == 1 & Dy == 1)) + 1
  dyspnea_LC1_Br0_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 1 & Br == 0 & Dy == 1)) + 1
  dyspnea_LC1_Br1_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 1 & Br == 1 & Dy == 1)) + 1
  lungCancer0_bronchitis0_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 0 & Br == 0))
  lungCancer0_bronchitis1_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 0 & Br == 1))
  lungCancer1_bronchitis0_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 1 & Br == 0))
  lungCancer1_bronchitis1_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 1 & Br == 1))
  p_dyspnea_LC0Br0 <- (dyspnea_LC0_Br0_counts / (lungCancer0_bronchitis0_counts + 2))
  p_dyspnea_LC0Br1 <- (dyspnea_LC0_Br1_counts / (lungCancer0_bronchitis1_counts + 2))
  p_dyspnea_LC1Br0 <- (dyspnea_LC1_Br0_counts / (lungCancer1_bronchitis0_counts + 2))
  p_dyspnea_LC1Br1 <- (dyspnea_LC1_Br1_counts / (lungCancer1_bronchitis1_counts + 2))
  
  dyspnea_matrix[1,3] <- (1 - p_dyspnea_LC0Br0)
  dyspnea_matrix[1,4] <- p_dyspnea_LC0Br0
  dyspnea_matrix[2,3] <- (1 - p_dyspnea_LC0Br1)
  dyspnea_matrix[2,4] <- p_dyspnea_LC0Br1
  dyspnea_matrix[3,3] <- (1 - p_dyspnea_LC1Br0)
  dyspnea_matrix[3,4] <- p_dyspnea_LC1Br0
  dyspnea_matrix[4,3] <- (1 - p_dyspnea_LC1Br1)
  dyspnea_matrix[4,4] <- p_dyspnea_LC1Br1
  
  network$dyspnea <- as.data.frame(dyspnea_matrix)
  
  #print('network$dyspnea')
  #print(network$dyspnea)
  
  xRay_matrix <- matrix(NA, nrow = 8, ncol = 5, byrow = TRUE)
  colnames(xRay_matrix) <- c('Pn', 'TB', 'LC', 'XR == 0', 'XR == 1')
  
  reset = FALSE
  for (counter in 1:8) {
    if (counter <= 4) {
      xRay_matrix[counter,1] <- 0
    } else {
      xRay_matrix[counter,1] <- 1
    }
    if (reset == FALSE) {
      xRay_matrix[counter,2] <- 0
    } else {
      xRay_matrix[counter,2] <- 1
    }
    if (counter %% 2 == 0) {
      reset <- !reset
      xRay_matrix[counter,3] <- 1
    } else {
      xRay_matrix[counter,3] <- 0
    }
  }
  
  xRay_df <- historical_data['XR']
  xRay_pneumonia_tuberculosis_lungCancer_df <- data.frame(pneumonia_df, tuberculosis_df, lungCancer_df, xRay_df)
  xRay_Pn0_TB0_LC0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 0 & LC == 0 & XR == 1)) + 1
  xRay_Pn0_TB0_LC1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 0 & LC == 1 & XR == 1)) + 1
  xRay_Pn0_TB1_LC0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 1 & LC == 0 & XR == 1)) + 1
  xRay_Pn0_TB1_LC1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 1 & LC == 1 & XR == 1)) + 1
  xRay_Pn1_TB0_LC0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 0 & LC == 0 & XR == 1)) + 1
  xRay_Pn1_TB0_LC1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 0 & LC == 1 & XR == 1)) + 1
  xRay_Pn1_TB1_LC0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 1 & LC == 0 & XR == 1)) + 1
  xRay_Pn1_TB1_LC1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 1 & LC == 1 & XR == 1)) + 1
  pneumonia0_tuberculosis0_lungCancer0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 0 & LC == 0))
  pneumonia0_tuberculosis0_lungCancer1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 0 & LC == 1))
  pneumonia0_tuberculosis1_lungCancer0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 1 & LC == 0))
  pneumonia0_tuberculosis1_lungCancer1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 1 & LC == 1))
  pneumonia1_tuberculosis0_lungCancer0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 0 & LC == 0))
  pneumonia1_tuberculosis0_lungCancer1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 0 & LC == 1))
  pneumonia1_tuberculosis1_lungCancer0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 1 & LC == 0))
  pneumonia1_tuberculosis1_lungCancer1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 1 & LC == 1))
  p_xRay_Pn0TB0LC0 <- (xRay_Pn0_TB0_LC0_counts / (pneumonia0_tuberculosis0_lungCancer0_counts + 2))
  p_xRay_Pn0TB0LC1 <- (xRay_Pn0_TB0_LC1_counts / (pneumonia0_tuberculosis0_lungCancer1_counts + 2))
  p_xRay_Pn0TB1LC0 <- (xRay_Pn0_TB1_LC0_counts / (pneumonia0_tuberculosis1_lungCancer0_counts + 2))
  p_xRay_Pn0TB1LC1 <- (xRay_Pn0_TB1_LC1_counts / (pneumonia0_tuberculosis1_lungCancer1_counts + 2))
  p_xRay_Pn1TB0LC0 <- (xRay_Pn1_TB0_LC0_counts / (pneumonia1_tuberculosis0_lungCancer0_counts + 2))
  p_xRay_Pn1TB0LC1 <- (xRay_Pn1_TB0_LC1_counts / (pneumonia1_tuberculosis0_lungCancer1_counts + 2))
  p_xRay_Pn1TB1LC0 <- (xRay_Pn1_TB1_LC0_counts / (pneumonia1_tuberculosis1_lungCancer0_counts + 2))
  p_xRay_Pn1TB1LC1 <- (xRay_Pn1_TB1_LC1_counts / (pneumonia1_tuberculosis1_lungCancer1_counts + 2))
  
  xRay_matrix[1,4] <- (1 - p_xRay_Pn0TB0LC0)
  xRay_matrix[1,5] <- p_xRay_Pn0TB0LC0
  xRay_matrix[2,4] <- (1 - p_xRay_Pn0TB0LC1)
  xRay_matrix[2,5] <- p_xRay_Pn0TB0LC1
  xRay_matrix[3,4] <- (1 - p_xRay_Pn0TB1LC0)
  xRay_matrix[3,5] <- p_xRay_Pn0TB1LC0
  xRay_matrix[4,4] <- (1 - p_xRay_Pn0TB1LC1)
  xRay_matrix[4,5] <- p_xRay_Pn0TB1LC1
  xRay_matrix[5,4] <- (1 - p_xRay_Pn1TB0LC0)
  xRay_matrix[5,5] <- p_xRay_Pn1TB0LC0
  xRay_matrix[6,4] <- (1 - p_xRay_Pn1TB0LC1)
  xRay_matrix[6,5] <- p_xRay_Pn1TB0LC1
  xRay_matrix[7,4] <- (1 - p_xRay_Pn1TB1LC0)
  xRay_matrix[7,5] <- p_xRay_Pn1TB1LC0
  xRay_matrix[8,4] <- (1 - p_xRay_Pn1TB1LC1)
  xRay_matrix[8,5] <- p_xRay_Pn1TB1LC1
  
  network$xRay <- as.data.frame(xRay_matrix)
  
  #print('network$xRay')
  #print(network$xRay)
  
  return(network)
}

#Function to estimate probabilities of unknown variables in a set of medical cases. 
diagnose <- function(network, cases) {
  number_cases <- nrow(cases)
  unknown_variables <- c('Pn', 'TB', 'LC', 'Br')
  
  diagnose_matrix <- matrix(NA, nrow = number_cases, ncol = 4, byrow = TRUE)
  colnames(diagnose_matrix) <- unknown_variables
  
  storage_randomBinaries <- rbinom(n = 40000, size = 1, prob = 0.5 )
  storage_randomProbabilites <- runif(n = 40000, min = 0, max = 1)
  counter_binaryStorage_index <- 1
  counter_randomProbabilites_index <- 1
  
  for (case_id in 1:number_cases) {
    
    sampling_results <- matrix(NA, nrow = 1, ncol = 9, byrow = TRUE)
    colnames(sampling_results) <- c('Pn', 'Te', 'VTB', 'TB', 'Sm', 'LC', 'Br', 'XR', 'Dy')
    
    sample_matrix <- matrix(NA, nrow = 1, ncol = 9, byrow = TRUE)
    colnames(sample_matrix) <- c('Pn', 'Te', 'VTB', 'TB', 'Sm', 'LC', 'Br', 'XR', 'Dy')
    
    for (column_index in 1:9) {
      sample_matrix[1, column_index] <- cases[case_id, column_index]
    }
    
    for (sampling in 1:1000) {
      
      for (unknown_var in 1:4) {
        sample_matrix[1, unknown_variables[unknown_var]] <- storage_randomBinaries[counter_binaryStorage_index]
        counter_binaryStorage_index <- counter_binaryStorage_index + 1
      }
      
      for (unknown_var in 1:4) {
        current_unknownVariable <- unknown_variables[unknown_var]
        
        current_probability <- network$pneumonia[[1, paste("Pn == ", sample_matrix[1, "Pn"], sep = "")]] *
          dnorm(sample_matrix[1, "Te"], subset(network$temperature, Pn == sample_matrix[1, "Pn"])[[1,'Mean']], subset(network$temperature, Pn == sample_matrix[1, "Pn"])[[1,'sd']])[[1]] *
          network$visitedTBSpot[[1, paste("VTB == ", sample_matrix[1, "VTB"], sep = "")]] *
          subset(network$tuberculosis, VTB == sample_matrix[1, "VTB"])[[1, paste("TB == ", sample_matrix[1, "TB"], sep = "")]] *
          network$smokes[[1, paste("Sm == ", sample_matrix[1, "Sm"], sep = "")]] *
          subset(network$lungCancer, Sm == sample_matrix[1, "Sm"])[[1, paste("LC == ", sample_matrix[1, "LC"], sep = "")]] *
          subset(network$bronchitis, Sm == sample_matrix[1, "Sm"])[[1, paste("Br == ", sample_matrix[1, "Br"], sep = "")]] *
          subset(network$dyspnea, LC == sample_matrix[1, "LC"] & Br == sample_matrix[1, "Br"])[[1, paste("Dy == ", sample_matrix[1, "Dy"], sep = "")]] *
          subset(network$xRay, Pn == sample_matrix[1, "Pn"] & TB == sample_matrix[1, "TB"] & LC == sample_matrix[1, "LC"])[[1, paste("XR == ", sample_matrix[1, "XR"], sep = "")]]
        
        current_unknownVariable_value <- sample_matrix[1, current_unknownVariable]
        
        if (current_unknownVariable_value == 0) {
          proposed_value <- 1
        } else {
          proposed_value <- 0
        }
        
        proposed_matrix <- sample_matrix
        proposed_matrix[1, current_unknownVariable] <- proposed_value
        
        proposed_probability <- network$pneumonia[[1, paste("Pn == ", proposed_matrix[1, "Pn"], sep = "")]] *
          dnorm(proposed_matrix[1, "Te"], subset(network$temperature, Pn == proposed_matrix[1, "Pn"])[[1,'Mean']], subset(network$temperature, Pn == proposed_matrix[1, "Pn"])[[1,'sd']])[[1]] *
          network$visitedTBSpot[[1, paste("VTB == ", proposed_matrix[1, "VTB"], sep = "")]] *
          subset(network$tuberculosis, VTB == proposed_matrix[1, "VTB"])[[1, paste("TB == ", proposed_matrix[1, "TB"], sep = "")]] *
          network$smokes[[1, paste("Sm == ", proposed_matrix[1, "Sm"], sep = "")]] *
          subset(network$lungCance, Sm == proposed_matrix[1, "Sm"])[[1, paste("LC == ", proposed_matrix[1, "LC"], sep = "")]] *
          subset(network$bronchitis, Sm == proposed_matrix[1, "Sm"])[[1, paste("Br == ", proposed_matrix[1, "Br"], sep = "")]] *
          subset(network$dyspnea, LC == proposed_matrix[1, "LC"] & Br == proposed_matrix[1, "Br"])[[1, paste("Dy == ", proposed_matrix[1, "Dy"], sep = "")]] *
          subset(network$xRay, Pn == proposed_matrix[1, "Pn"] & TB == proposed_matrix[1, "TB"] & LC == proposed_matrix[1, "LC"])[[1, paste("XR == ", proposed_matrix[1, "XR"], sep = "")]]
        
        if (proposed_probability > current_probability) {
          sample_matrix <- proposed_matrix
        } else {
          ratio <- (proposed_probability / current_probability)
          random_number <- storage_randomProbabilites[counter_randomProbabilites_index]
          counter_randomProbabilites_index <- counter_randomProbabilites_index + 1
          if (random_number < ratio) {
            sample_matrix <- proposed_matrix
          }
        }
      }
      
      if (sampling > 100) {
        sampling_results <- rbind(sampling_results, sample_matrix)
      }
    }
    #Analyse sampling_results table to get new probabilities
    #print(sampling_results)
  }
  
  #print(diagnose_matrix)
  return(diagnose_matrix)
}

runDiagnostics(learn, diagnose, verbose = 2)
