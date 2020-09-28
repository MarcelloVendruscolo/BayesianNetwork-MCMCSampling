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

  network$pneumonia_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(network$pneumonia_matrix) <- c('Pn == 0', 'Pn == 1')
  
  pneumonia_df <- historical_data['Pn']
  pneumonia_counts <- nrow(subset(pneumonia_df, Pn == 1)) + 1
  p_pneumonia <- (pneumonia_counts / total_observations)
  
  network$pneumonia_matrix[1,1] <- (1 - p_pneumonia)
  network$pneumonia_matrix[1,2] <- (p_pneumonia)
  
  print('network$pneumonia_matrix')
  print(network$pneumonia_matrix)
  
  network$temperature_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(network$temperature_matrix) <- c('Pn', 'Mean', 'sd')
  network$temperature_matrix[1,1] <- 0
  network$temperature_matrix[2,1] <- 1
  
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
  
  network$temperature_matrix[1,2] <- mean_temperature_pneumonia0
  network$temperature_matrix[1,3] <- sd_temperature_pneumonia0
  network$temperature_matrix[2,2] <- mean_temperature_pneumonia1
  network$temperature_matrix[2,3] <- sd_temperature_pneumonia1
  
  print('network$temperature_matrix')
  print(network$temperature_matrix)
  
  network$visitedTBSpot_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(network$visitedTBSpot_matrix) <- c('VTB == 0', 'VTB == 1')
  
  visitedTBSpot_df <- historical_data['VTB']
  visitedTBSpot_counts <- nrow(subset(visitedTBSpot_df, VTB == 1)) + 1
  p_visitedTBSpot <- (visitedTBSpot_counts / total_observations)
  
  network$visitedTBSpot_matrix[1,1] <- (1 - p_visitedTBSpot)
  network$visitedTBSpot_matrix[1,2] <- (p_visitedTBSpot)
  
  print('network$visitedTBSpot_matrix')
  print(network$visitedTBSpot_matrix)
  
  network$tuberculosis_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(network$tuberculosis_matrix) <- c('VTB', 'TB == 0', 'TB == 1')
  network$tuberculosis_matrix[1,1] <- 0
  network$tuberculosis_matrix[2,1] <- 1
  
  tuberculosis_df <- historical_data['TB']
  tuberculosis_visitedTBSpot_df <- data.frame(visitedTBSpot_df, tuberculosis_df)
  tuberculosis_VTBSpot0_counts <- nrow(subset(tuberculosis_visitedTBSpot_df, VTB == 0 & TB == 1)) + 1
  tuberculosis_VTBSpot1_counts <- nrow(subset(tuberculosis_visitedTBSpot_df, VTB == 1 & TB == 1)) + 1
  p_tuberculosis_VTBSpot0 <- (tuberculosis_VTBSpot0_counts / ((total_observations - visitedTBSpot_counts) + 2))
  p_tuberculosis_VTBSpot1 <- (tuberculosis_VTBSpot1_counts / (visitedTBSpot_counts + 2))
  
  network$tuberculosis_matrix[1,2] <- (1 - p_tuberculosis_VTBSpot0)
  network$tuberculosis_matrix[1,3] <- p_tuberculosis_VTBSpot0
  network$tuberculosis_matrix[2,2] <- (1 - p_tuberculosis_VTBSpot1)
  network$tuberculosis_matrix[2,3] <- p_tuberculosis_VTBSpot1
  
  print('network$tuberculosis_matrix')
  print(network$tuberculosis_matrix)
  
  network$smokes_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(network$smokes_matrix) <- c('Sm == 0', 'Sm == 1')
  
  smokes_df <- historical_data['Sm']
  smokes_counts <- nrow(subset(smokes_df, Sm == 1)) + 1
  p_smokes <- (smokes_counts / total_observations)
  
  network$smokes_matrix[1,1] <- (1 - p_smokes)
  network$smokes_matrix[1,2] <- (p_smokes)
  
  print('network$smokes_matrix')
  print(network$smokes_matrix)
  
  network$lungCancer_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(network$lungCancer_matrix) <- c('Sm', 'LC == 0', 'LC == 1')
  network$lungCancer_matrix[1,1] <- 0
  network$lungCancer_matrix[2,1] <- 1
  
  lungCancer_df <- historical_data['LC']
  lungCancer_smokes_df <- data.frame(smokes_df, lungCancer_df)
  lungCancer_Sm0_counts <- nrow(subset(lungCancer_smokes_df, Sm == 0 & LC == 1)) + 1
  lungCancer_Sm1_counts <- nrow(subset(lungCancer_smokes_df, Sm == 1 & LC == 1)) + 1
  p_lungCancer_Sm0 <- (lungCancer_Sm0_counts / ((total_observations - smokes_counts) + 2))
  p_lungCancer_Sm1 <- (lungCancer_Sm1_counts / (smokes_counts + 2))
  
  network$lungCancer_matrix[1,2] <- (1 - p_lungCancer_Sm0)
  network$lungCancer_matrix[1,3] <- p_lungCancer_Sm0
  network$lungCancer_matrix[2,2] <- (1 - p_lungCancer_Sm1)
  network$lungCancer_matrix[2,3] <- p_lungCancer_Sm1
  
  print('network$lungCancer_matrix')
  print(network$lungCancer)
  
  network$bronchitis_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(network$bronchitis_matrix) <- c('Sm', 'Br == 0', 'Br == 1')
  network$bronchitis_matrix[1,1] <- 0
  network$bronchitis_matrix[2,1] <- 1
  
  bronchitis_df <- historical_data['Br']
  bronchitis_smokes_df <- data.frame(smokes_df, bronchitis_df)
  bronchitis_Sm0_counts <- nrow(subset(bronchitis_smokes_df, Sm == 0 & Br == 1)) + 1
  bronchitis_Sm1_counts <- nrow(subset(bronchitis_smokes_df, Sm == 1 & Br == 1)) + 1
  p_bronchitis_Sm0 <- (bronchitis_Sm0_counts / ((total_observations - smokes_counts) + 2))
  p_bronchitis_Sm1 <- (bronchitis_Sm1_counts / (smokes_counts + 2))
  
  network$bronchitis_matrix[1,2] <- (1 - p_bronchitis_Sm0)
  network$bronchitis_matrix[1,3] <- p_bronchitis_Sm0
  network$bronchitis_matrix[2,2] <- (1 - p_bronchitis_Sm1)
  network$bronchitis_matrix[2,3] <- p_bronchitis_Sm1
  
  print('network$bronchitis_matrix')
  print(network$bronchitis_matrix)
  
  network$dyspnea_matrix <- matrix(NA, nrow = 4, ncol = 4, byrow = TRUE)
  colnames(network$dyspnea_matrix) <- c('LC', 'Br', 'Dy == 0', 'Dy == 1')
  network$dyspnea_matrix[1,1] <- 0
  network$dyspnea_matrix[2,1] <- 0
  network$dyspnea_matrix[1,2] <- 0
  network$dyspnea_matrix[3,2] <- 0
  network$dyspnea_matrix[3,1] <- 1
  network$dyspnea_matrix[4,1] <- 1
  network$dyspnea_matrix[2,2] <- 1
  network$dyspnea_matrix[4,2] <- 1
  
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
  
  network$dyspnea_matrix[1,3] <- (1 - p_dyspnea_LC0Br0)
  network$dyspnea_matrix[1,4] <- p_dyspnea_LC0Br0
  network$dyspnea_matrix[2,3] <- (1 - p_dyspnea_LC0Br1)
  network$dyspnea_matrix[2,4] <- p_dyspnea_LC0Br1
  network$dyspnea_matrix[3,3] <- (1 - p_dyspnea_LC1Br0)
  network$dyspnea_matrix[3,4] <- p_dyspnea_LC1Br0
  network$dyspnea_matrix[4,3] <- (1 - p_dyspnea_LC1Br1)
  network$dyspnea_matrix[4,4] <- p_dyspnea_LC1Br1
  
  print('network$dyspnea_matrix')
  print(network$dyspnea_matrix)
  
  network$xRay_matrix <- matrix(NA, nrow = 8, ncol = 5, byrow = TRUE)
  colnames(network$xRay_matrix) <- c('Pn', 'TB', 'LC', 'XR == 0', 'XR == 1')
  
  reset = FALSE
  for (counter in 1:8) {
    if (counter <= 4) {
      network$xRay_matrix[counter,1] <- 0
    } else {
      network$xRay_matrix[counter,1] <- 1
    }
    if (reset == FALSE) {
      network$xRay_matrix[counter,2] <- 0
    } else {
      network$xRay_matrix[counter,2] <- 1
    }
    if (counter %% 2 == 0) {
      reset <- !reset
      network$xRay_matrix[counter,3] <- 1
    } else {
      network$xRay_matrix[counter,3] <- 0
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
  
  network$xRay_matrix[1,4] <- (1 - p_xRay_Pn0TB0LC0)
  network$xRay_matrix[1,5] <- p_xRay_Pn0TB0LC0
  network$xRay_matrix[2,4] <- (1 - p_xRay_Pn0TB0LC1)
  network$xRay_matrix[2,5] <- p_xRay_Pn0TB0LC1
  network$xRay_matrix[3,4] <- (1 - p_xRay_Pn0TB1LC0)
  network$xRay_matrix[3,5] <- p_xRay_Pn0TB1LC0
  network$xRay_matrix[4,4] <- (1 - p_xRay_Pn0TB1LC1)
  network$xRay_matrix[4,5] <- p_xRay_Pn0TB1LC1
  network$xRay_matrix[5,4] <- (1 - p_xRay_Pn1TB0LC0)
  network$xRay_matrix[5,5] <- p_xRay_Pn1TB0LC0
  network$xRay_matrix[6,4] <- (1 - p_xRay_Pn1TB0LC1)
  network$xRay_matrix[6,5] <- p_xRay_Pn1TB0LC1
  network$xRay_matrix[7,4] <- (1 - p_xRay_Pn1TB1LC0)
  network$xRay_matrix[7,5] <- p_xRay_Pn1TB1LC0
  network$xRay_matrix[8,4] <- (1 - p_xRay_Pn1TB1LC1)
  network$xRay_matrix[8,5] <- p_xRay_Pn1TB1LC1
  
  print('network$xRay_matrix')
  print(network$xRay_matrix)
  
  return(network)
}

#Function to estimate probabilities of unknown variables in a set of medical cases. 
diagnose <- function(network, cases) {
  
}

runDiagnostics(learn, diagnose, verbose = 0)
