#List all the packages available
library()

#Install the Diagnostics package
install.packages("/Users/marcellovendruscolo/Documents/rstudio-workspace/Diagnostics/Diagnostics_1.2.0.tar.gz", repos = NULL, type = "source")

#Load Diagnostics package for this session
library(Diagnostics)

#List all the packages currently loaded
search()


#Function to create and learn the Bayesian network. The only parameter is the historical medical data.
learn <- function(historical_data) {
  print(historical_data)
}

diagnose <- function(network, cases) {
  print("diagnose function")
}

runDiagnostics(learn, diagnose, verbose = 0)
