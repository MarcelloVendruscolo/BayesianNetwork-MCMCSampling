#List all the packages available
library()

#Install the Diagnostics package
install.packages("/Users/marcellovendruscolo/Documents/rstudio-workspace/Diagnostics/Diagnostics_1.2.0.tar.gz", repos = NULL, type = "source")

#Load Diagnostics package for this session
library(Diagnostics)

#List all the packages currently loaded
search()

#Function to create and learn the Bayesian network. The only parameter is the historical medical data. It returns the learned network as a list of dataframe
learn <- function(historical_data) {
  total_observations <- nrow(historical_data) #Variable holding the number of observations in the historical data
  network <- list() #Variable holding the learned network to be passed onto diagnose function
  
  #This section fits the pneumonia categorical distribution to data
  #Conservatism is added by adding 1 sample of each case ('Pn == 0' and 'Pn == 1') to the observed data
  pneumonia_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(pneumonia_matrix) <- c('Pn == 0', 'Pn == 1')
  
  pneumonia_df <- historical_data['Pn']
  pneumonia_counts <- nrow(subset(pneumonia_df, Pn == 1))
  p_pneumonia <- (pneumonia_counts + 1) / (total_observations + 2)
  
  pneumonia_matrix[1,'Pn == 0'] <- 1 - p_pneumonia
  pneumonia_matrix[1,'Pn == 1'] <- p_pneumonia
  
  network$pneumonia <- as.data.frame(pneumonia_matrix)
  
  #This section fits the temperature conditional normal distribution to data
  temperature_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(temperature_matrix) <- c('Pn', 'Mean', 'sd')
  temperature_matrix[1,'Pn'] <- 0
  temperature_matrix[2,'Pn'] <- 1
  
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
  
  temperature_matrix[1,'Mean'] <- mean_temperature_pneumonia0
  temperature_matrix[1,'sd'] <- sd_temperature_pneumonia0
  temperature_matrix[2,'Mean'] <- mean_temperature_pneumonia1
  temperature_matrix[2,'sd'] <- sd_temperature_pneumonia1
  
  network$temperature <- as.data.frame(temperature_matrix)
  
  #This section fits the Visited TB Spot categorical distribution to data
  #Conservatism is added by adding 1 sample of each case ('VTB == 0' and 'VTB == 1') to the observed data
  visitedTBSpot_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(visitedTBSpot_matrix) <- c('VTB == 0', 'VTB == 1')
  
  visitedTBSpot_df <- historical_data['VTB']
  visitedTBSpot_counts <- nrow(subset(visitedTBSpot_df, VTB == 1))
  p_visitedTBSpot <- (visitedTBSpot_counts + 1) / (total_observations + 2)
  
  visitedTBSpot_matrix[1,1] <- 1 - p_visitedTBSpot
  visitedTBSpot_matrix[1,2] <- p_visitedTBSpot
  
  network$visitedTBSpot <- as.data.frame(visitedTBSpot_matrix)
  
  #This section fits the Tuberculosis conditional categorical distribution to data
  #Conservatism is added by adding 1 sample of each case ('TB == 0' and 'TB == 1') for each conditional to the observed data
  tuberculosis_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(tuberculosis_matrix) <- c('VTB', 'TB == 0', 'TB == 1')
  tuberculosis_matrix[1,'VTB'] <- 0
  tuberculosis_matrix[2,'VTB'] <- 1
  
  tuberculosis_df <- historical_data['TB']
  tuberculosis_visitedTBSpot_df <- data.frame(visitedTBSpot_df, tuberculosis_df)
  tuberculosis_VTBSpot0_counts <- nrow(subset(tuberculosis_visitedTBSpot_df, VTB == 0 & TB == 1))
  tuberculosis_VTBSpot1_counts <- nrow(subset(tuberculosis_visitedTBSpot_df, VTB == 1 & TB == 1))
  p_tuberculosis_VTBSpot0 <- (tuberculosis_VTBSpot0_counts + 1) / (total_observations - visitedTBSpot_counts + 2)
  p_tuberculosis_VTBSpot1 <- (tuberculosis_VTBSpot1_counts + 1) / (visitedTBSpot_counts + 2)
  
  tuberculosis_matrix[1,'TB == 0'] <- (1 - p_tuberculosis_VTBSpot0)
  tuberculosis_matrix[1,'TB == 1'] <- p_tuberculosis_VTBSpot0
  tuberculosis_matrix[2,'TB == 0'] <- (1 - p_tuberculosis_VTBSpot1)
  tuberculosis_matrix[2,'TB == 1'] <- p_tuberculosis_VTBSpot1
  
  network$tuberculosis <- as.data.frame(tuberculosis_matrix)
  
  #This section fits the Smokes categorical distribution to data
  #Conservatism is added by adding 1 sample of each case ('Sm == 0' and 'Sm == 1') to the observed data
  smokes_matrix <- matrix(NA, nrow = 1, ncol = 2, byrow = TRUE)
  colnames(smokes_matrix) <- c('Sm == 0', 'Sm == 1')
  
  smokes_df <- historical_data['Sm']
  smokes_counts <- nrow(subset(smokes_df, Sm == 1))
  p_smokes <- (smokes_counts + 1) / (total_observations + 2)
  
  smokes_matrix[1,'Sm == 0'] <- 1 - p_smokes
  smokes_matrix[1,'Sm == 1'] <- p_smokes
  
  network$smokes <- as.data.frame(smokes_matrix)
  
  #This section fits the Lung Cancer conditional categorical distribution to data
  #Conservatism is added by adding 1 sample of each case ('LC == 0' and 'LC == 1') for each conditional to the observed data
  lungCancer_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(lungCancer_matrix) <- c('Sm', 'LC == 0', 'LC == 1')
  lungCancer_matrix[1,'Sm'] <- 0
  lungCancer_matrix[2,'Sm'] <- 1
  
  lungCancer_df <- historical_data['LC']
  lungCancer_smokes_df <- data.frame(smokes_df, lungCancer_df)
  lungCancer_Sm0_counts <- nrow(subset(lungCancer_smokes_df, Sm == 0 & LC == 1))
  lungCancer_Sm1_counts <- nrow(subset(lungCancer_smokes_df, Sm == 1 & LC == 1))
  p_lungCancer_Sm0 <- (lungCancer_Sm0_counts + 1) / (total_observations - smokes_counts + 2)
  p_lungCancer_Sm1 <- (lungCancer_Sm1_counts + 1) / (smokes_counts + 2)
  
  lungCancer_matrix[1,'LC == 0'] <- 1 - p_lungCancer_Sm0
  lungCancer_matrix[1,'LC == 1'] <- p_lungCancer_Sm0
  lungCancer_matrix[2,'LC == 0'] <- 1 - p_lungCancer_Sm1
  lungCancer_matrix[2,'LC == 1'] <- p_lungCancer_Sm1
  
  network$lungCancer <- as.data.frame(lungCancer_matrix)
  
  #This section fits the Bronchitis conditional categorical distribution to data
  #Conservatism is added by adding 1 sample of each case ('Br == 0' and 'Br == 1') for each conditional to the observed data
  bronchitis_matrix <- matrix(NA, nrow = 2, ncol = 3, byrow = TRUE)
  colnames(bronchitis_matrix) <- c('Sm', 'Br == 0', 'Br == 1')
  bronchitis_matrix[1,'Sm'] <- 0
  bronchitis_matrix[2,'Sm'] <- 1
  
  bronchitis_df <- historical_data['Br']
  bronchitis_smokes_df <- data.frame(smokes_df, bronchitis_df)
  bronchitis_Sm0_counts <- nrow(subset(bronchitis_smokes_df, Sm == 0 & Br == 1))
  bronchitis_Sm1_counts <- nrow(subset(bronchitis_smokes_df, Sm == 1 & Br == 1))
  p_bronchitis_Sm0 <- (bronchitis_Sm0_counts + 1) / (total_observations - smokes_counts + 2)
  p_bronchitis_Sm1 <- (bronchitis_Sm1_counts + 1) / (smokes_counts + 2)
  
  bronchitis_matrix[1,'Br == 0'] <- (1 - p_bronchitis_Sm0)
  bronchitis_matrix[1,'Br == 1'] <- p_bronchitis_Sm0
  bronchitis_matrix[2,'Br == 0'] <- (1 - p_bronchitis_Sm1)
  bronchitis_matrix[2,'Br == 1'] <- p_bronchitis_Sm1
  
  network$bronchitis <- as.data.frame(bronchitis_matrix)
  
  #This section fits the Dyspnea conditional categorical distribution to data
  #Conservatism is added by adding 1 sample of each case ('Dy == 0' and 'Dy == 1') for each conditional to the observed data
  dyspnea_matrix <- matrix(NA, nrow = 4, ncol = 4, byrow = TRUE)
  colnames(dyspnea_matrix) <- c('LC', 'Br', 'Dy == 0', 'Dy == 1')
  dyspnea_matrix[1,'LC'] <- 0
  dyspnea_matrix[2,'LC'] <- 0
  dyspnea_matrix[3,'LC'] <- 1
  dyspnea_matrix[4,'LC'] <- 1
  dyspnea_matrix[1,'Br'] <- 0
  dyspnea_matrix[2,'Br'] <- 1
  dyspnea_matrix[3,'Br'] <- 0
  dyspnea_matrix[4,'Br'] <- 1
  
  dyspnea_df <- historical_data['Dy']
  dyspnea_lungCancer_Bronchitis_df <- data.frame(lungCancer_df, bronchitis_df, dyspnea_df)
  dyspnea_LC0_Br0_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 0 & Br == 0 & Dy == 1))
  dyspnea_LC0_Br1_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 0 & Br == 1 & Dy == 1))
  dyspnea_LC1_Br0_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 1 & Br == 0 & Dy == 1))
  dyspnea_LC1_Br1_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 1 & Br == 1 & Dy == 1))
  lungCancer0_bronchitis0_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 0 & Br == 0))
  lungCancer0_bronchitis1_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 0 & Br == 1))
  lungCancer1_bronchitis0_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 1 & Br == 0))
  lungCancer1_bronchitis1_counts <- nrow(subset(dyspnea_lungCancer_Bronchitis_df, LC == 1 & Br == 1))
  p_dyspnea_LC0Br0 <- (dyspnea_LC0_Br0_counts + 1) / (lungCancer0_bronchitis0_counts + 2)
  p_dyspnea_LC0Br1 <- (dyspnea_LC0_Br1_counts + 1) / (lungCancer0_bronchitis1_counts + 2)
  p_dyspnea_LC1Br0 <- (dyspnea_LC1_Br0_counts + 1) / (lungCancer1_bronchitis0_counts + 2)
  p_dyspnea_LC1Br1 <- (dyspnea_LC1_Br1_counts + 1) / (lungCancer1_bronchitis1_counts + 2)
  
  dyspnea_matrix[1,'Dy == 0'] <- 1 - p_dyspnea_LC0Br0
  dyspnea_matrix[1,'Dy == 1'] <- p_dyspnea_LC0Br0
  dyspnea_matrix[2,'Dy == 0'] <- 1 - p_dyspnea_LC0Br1
  dyspnea_matrix[2,'Dy == 1'] <- p_dyspnea_LC0Br1
  dyspnea_matrix[3,'Dy == 0'] <- 1 - p_dyspnea_LC1Br0
  dyspnea_matrix[3,'Dy == 1'] <- p_dyspnea_LC1Br0
  dyspnea_matrix[4,'Dy == 0'] <- 1 - p_dyspnea_LC1Br1
  dyspnea_matrix[4,'Dy == 1'] <- p_dyspnea_LC1Br1
  
  network$dyspnea <- as.data.frame(dyspnea_matrix)
  
  #This section fits the X-Ray Result conditional categorical distribution to data
  #Conservatism is added by adding 1 sample of each case ('XR == 0' and 'XR == 1') for each conditional to the observed data
  xRay_matrix <- matrix(NA, nrow = 8, ncol = 5, byrow = TRUE)
  colnames(xRay_matrix) <- c('Pn', 'TB', 'LC', 'XR == 0', 'XR == 1')
  
  #This loop fills up the first three columns ('Pn', 'TB', 'LC') of the matrix in a manner to get all possible combinations (2^3 = 8 rows)
  reset = FALSE
  for (counter in 1:8) {
    if (counter <= 4) {
      xRay_matrix[counter,'Pn'] <- 0
    } else {
      xRay_matrix[counter,'Pn'] <- 1
    }
    if (reset == FALSE) {
      xRay_matrix[counter,'TB'] <- 0
    } else {
      xRay_matrix[counter,'TB'] <- 1
    }
    if (counter %% 2 == 0) {
      reset <- !reset
      xRay_matrix[counter,'LC'] <- 1
    } else {
      xRay_matrix[counter,'LC'] <- 0
    }
  }
  
  xRay_df <- historical_data['XR']
  xRay_pneumonia_tuberculosis_lungCancer_df <- data.frame(pneumonia_df, tuberculosis_df, lungCancer_df, xRay_df)
  xRay_Pn0_TB0_LC0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 0 & LC == 0 & XR == 1))
  xRay_Pn0_TB0_LC1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 0 & LC == 1 & XR == 1))
  xRay_Pn0_TB1_LC0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 1 & LC == 0 & XR == 1))
  xRay_Pn0_TB1_LC1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 1 & LC == 1 & XR == 1))
  xRay_Pn1_TB0_LC0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 0 & LC == 0 & XR == 1))
  xRay_Pn1_TB0_LC1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 0 & LC == 1 & XR == 1))
  xRay_Pn1_TB1_LC0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 1 & LC == 0 & XR == 1))
  xRay_Pn1_TB1_LC1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 1 & LC == 1 & XR == 1))
  pneumonia0_tuberculosis0_lungCancer0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 0 & LC == 0))
  pneumonia0_tuberculosis0_lungCancer1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 0 & LC == 1))
  pneumonia0_tuberculosis1_lungCancer0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 1 & LC == 0))
  pneumonia0_tuberculosis1_lungCancer1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 0 & TB == 1 & LC == 1))
  pneumonia1_tuberculosis0_lungCancer0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 0 & LC == 0))
  pneumonia1_tuberculosis0_lungCancer1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 0 & LC == 1))
  pneumonia1_tuberculosis1_lungCancer0_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 1 & LC == 0))
  pneumonia1_tuberculosis1_lungCancer1_counts <- nrow(subset(xRay_pneumonia_tuberculosis_lungCancer_df, Pn == 1 & TB == 1 & LC == 1))
  p_xRay_Pn0TB0LC0 <- (xRay_Pn0_TB0_LC0_counts + 1) / (pneumonia0_tuberculosis0_lungCancer0_counts + 2)
  p_xRay_Pn0TB0LC1 <- (xRay_Pn0_TB0_LC1_counts + 1) / (pneumonia0_tuberculosis0_lungCancer1_counts + 2)
  p_xRay_Pn0TB1LC0 <- (xRay_Pn0_TB1_LC0_counts + 1) / (pneumonia0_tuberculosis1_lungCancer0_counts + 2)
  p_xRay_Pn0TB1LC1 <- (xRay_Pn0_TB1_LC1_counts + 1) / (pneumonia0_tuberculosis1_lungCancer1_counts + 2)
  p_xRay_Pn1TB0LC0 <- (xRay_Pn1_TB0_LC0_counts + 1) / (pneumonia1_tuberculosis0_lungCancer0_counts + 2)
  p_xRay_Pn1TB0LC1 <- (xRay_Pn1_TB0_LC1_counts + 1) / (pneumonia1_tuberculosis0_lungCancer1_counts + 2)
  p_xRay_Pn1TB1LC0 <- (xRay_Pn1_TB1_LC0_counts + 1) / (pneumonia1_tuberculosis1_lungCancer0_counts + 2)
  p_xRay_Pn1TB1LC1 <- (xRay_Pn1_TB1_LC1_counts + 1) / (pneumonia1_tuberculosis1_lungCancer1_counts + 2)
  
  xRay_matrix[1,'XR == 0'] <- 1 - p_xRay_Pn0TB0LC0
  xRay_matrix[1,'XR == 1'] <- p_xRay_Pn0TB0LC0
  xRay_matrix[2,'XR == 0'] <- 1 - p_xRay_Pn0TB0LC1
  xRay_matrix[2,'XR == 1'] <- p_xRay_Pn0TB0LC1
  xRay_matrix[3,'XR == 0'] <- 1 - p_xRay_Pn0TB1LC0
  xRay_matrix[3,'XR == 1'] <- p_xRay_Pn0TB1LC0
  xRay_matrix[4,'XR == 0'] <- 1 - p_xRay_Pn0TB1LC1
  xRay_matrix[4,'XR == 1'] <- p_xRay_Pn0TB1LC1
  xRay_matrix[5,'XR == 0'] <- 1 - p_xRay_Pn1TB0LC0
  xRay_matrix[5,'XR == 1'] <- p_xRay_Pn1TB0LC0
  xRay_matrix[6,'XR == 0'] <- 1 - p_xRay_Pn1TB0LC1
  xRay_matrix[6,'XR == 1'] <- p_xRay_Pn1TB0LC1
  xRay_matrix[7,'XR == 0'] <- 1 - p_xRay_Pn1TB1LC0
  xRay_matrix[7,'XR == 1'] <- p_xRay_Pn1TB1LC0
  xRay_matrix[8,'XR == 0'] <- 1 - p_xRay_Pn1TB1LC1
  xRay_matrix[8,'XR == 1'] <- p_xRay_Pn1TB1LC1
  
  network$xRay <- as.data.frame(xRay_matrix)
  
  return(network)
}

#Function to estimate probabilities of unknown variables in a set of medical cases
diagnose <- function(network, cases) {
  
  sampling_size <- 1000 #Hard-coded sampling size
  cat('Metropolis in Gibbs sampling size:', sampling_size, '\n')
  unknown_variables <- c('Pn', 'TB', 'LC', 'Br') #Hard-coded unknown variables
  
  total_cases <- nrow(cases) #Variable holding the number of observations in the medical cases
  total_unknown_var <- length(unknown_variables)
  
  #Matrix to be returned by the function
  diagnose_matrix <- matrix(NA, nrow = total_cases, ncol = 4, byrow = TRUE)
  colnames(diagnose_matrix) <- unknown_variables
  
  #Preparing in advance a set of random binary numbers and probabilities for sampling.
  storage_randomBinaries <- rbinom(n = (total_cases * sampling_size * length(unknown_variables)), size = 1, prob = 0.5)
  storage_randomProbabilites <- runif(n = (total_cases * sampling_size * length(unknown_variables)), min = 0, max = 1)
  counter_binaryStorage_index <- 1
  counter_randomProbabilites_index <- 1
  
  for (case_id in 1:total_cases) {
    
    #Matrix variable that holds the outcomes of each sampling for analyse when each medical test case is completed
    sampling_results <- matrix(NA, nrow = 1, ncol = 9, byrow = TRUE)
    colnames(sampling_results) <- c('Pn', 'Te', 'VTB', 'TB', 'Sm', 'LC', 'Br', 'XR', 'Dy')
    
    #Temporary matrix variable that holds the data from medical test cases to avoid modifying original data
    sample_matrix <- matrix(NA, nrow = 1, ncol = 9, byrow = TRUE)
    colnames(sample_matrix) <- c('Pn', 'Te', 'VTB', 'TB', 'Sm', 'LC', 'Br', 'XR', 'Dy')
    
    #Prepare the temporary matrix copying the corresponding medical test data row
    for (column_index in 1:9) {
      sample_matrix[1, column_index] <- cases[case_id, column_index]
    }
    
    #Allocate randomly generated binary values (0 or 1) to the unknown variables in the temporary matrix
    for (unknown_var in 1:total_unknown_var) {
      sample_matrix[1, unknown_variables[unknown_var]] <- storage_randomBinaries[counter_binaryStorage_index]
      counter_binaryStorage_index <- counter_binaryStorage_index + 1
    }
    
    old_prob <- 0 #Variable that holds the probability achieved in the previous sample
    first_sample <- TRUE #Variable that indicates whether it is the first sample of a medical test case
    
    #The following for block represents the remaining sampling iterations. Metropolis in Gibbs algorithm applied
    for (sampling in 1:sampling_size) {
      
      for (unknown_var in 1:total_unknown_var) {
        current_unknownVariable <- unknown_variables[unknown_var]
        
        if (first_sample) {
          current_probability <- network$pneumonia[[1, sample_matrix[1, "Pn"] + 1]] *
            dnorm(sample_matrix[1, "Te"], subset(network$temperature, Pn == sample_matrix[1, "Pn"])[[1,'Mean']], subset(network$temperature, Pn == sample_matrix[1, "Pn"])[[1,'sd']])[[1]] *
            network$visitedTBSpot[[1, sample_matrix[1, "VTB"] + 1]] *
            subset(network$tuberculosis, VTB == sample_matrix[1, "VTB"])[[1, sample_matrix[1, "TB"] + 2]] *
            network$smokes[[1, sample_matrix[1, "Sm"] + 1]] *
            subset(network$lungCancer, Sm == sample_matrix[1, "Sm"])[[1, sample_matrix[1, "LC"] + 2]] *
            subset(network$bronchitis, Sm == sample_matrix[1, "Sm"])[[1, sample_matrix[1, "Br"] + 2]] *
            subset(network$dyspnea, LC == sample_matrix[1, "LC"] & Br == sample_matrix[1, "Br"])[[1, sample_matrix[1, "Dy"] + 3]] *
            subset(network$xRay, Pn == sample_matrix[1, "Pn"] & TB == sample_matrix[1, "TB"] & LC == sample_matrix[1, "LC"])[[1, sample_matrix[1, "XR"] + 4]]
          old_prob <- current_probability
          first_sample <- FALSE
        }
        
        current_unknownVariable_value <- sample_matrix[1, current_unknownVariable]
        proposed_value <- 1 - current_unknownVariable_value
        proposed_matrix <- sample_matrix
        proposed_matrix[1, current_unknownVariable] <- proposed_value
        
        proposed_probability <- network$pneumonia[[1, proposed_matrix[1, "Pn"] + 1]] *
          dnorm(proposed_matrix[1, "Te"], subset(network$temperature, Pn == proposed_matrix[1, "Pn"])[[1,'Mean']], subset(network$temperature, Pn == proposed_matrix[1, "Pn"])[[1,'sd']])[[1]] *
          network$visitedTBSpot[[1, proposed_matrix[1, "VTB"] + 1]] *
          subset(network$tuberculosis, VTB == proposed_matrix[1, "VTB"])[[1, proposed_matrix[1, "TB"] + 2]] *
          network$smokes[[1, proposed_matrix[1, "Sm"] + 1]] *
          subset(network$lungCancer, Sm == proposed_matrix[1, "Sm"])[[1, proposed_matrix[1, "LC"] + 2]] *
          subset(network$bronchitis, Sm == proposed_matrix[1, "Sm"])[[1, proposed_matrix[1, "Br"] + 2]] *
          subset(network$dyspnea, LC == proposed_matrix[1, "LC"] & Br == proposed_matrix[1, "Br"])[[1, proposed_matrix[1, "Dy"] + 3]] *
          subset(network$xRay, Pn == proposed_matrix[1, "Pn"] & TB == proposed_matrix[1, "TB"] & LC == proposed_matrix[1, "LC"])[[1, proposed_matrix[1, "XR"] + 4]]
        
        if (proposed_probability > old_prob) {
          sample_matrix <- proposed_matrix
          old_prob <- proposed_probability
        } else {
          ratio <- (proposed_probability / old_prob)
          random_number <- storage_randomProbabilites[counter_randomProbabilites_index]
          counter_randomProbabilites_index <- counter_randomProbabilites_index + 1
          if (random_number < ratio) {
            sample_matrix <- proposed_matrix
            old_prob <- proposed_probability
          }
        }
      }
      
      if (sampling > sampling_size %/% 10) {
        sampling_results <- rbind(sampling_results, sample_matrix)
      }
    }
    
    #Analyse sampling_results table to get new probabilities
    sampling_results <- as.data.frame(sampling_results[-1, ])
    diagnose_matrix[case_id, 'Pn'] <- (nrow(subset(sampling_results, Pn == 1)) / 900)
    diagnose_matrix[case_id, 'TB'] <- (nrow(subset(sampling_results, TB == 1)) / 900)
    diagnose_matrix[case_id, 'LC'] <- (nrow(subset(sampling_results, LC == 1)) / 900)
    diagnose_matrix[case_id, 'Br'] <- (nrow(subset(sampling_results, Br == 1)) / 900)
    cat('Test case', case_id, '/', total_cases, 'completed.\n')
  }
  return(diagnose_matrix)
}

runDiagnostics(learn, diagnose, verbose = 2)

