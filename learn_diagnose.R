#List all the packages available
library()

#Install the Diagnostics package
install.packages("/Users/marcellovendruscolo/Documents/rstudio-workspace/Diagnostics/Diagnostics_1.2.0.tar.gz", repos = NULL, type = "source")

#Load Diagnostics package for this session
library(DeliveryMan)

#List all the packages currently loaded
search()

learn <- function(historical_data) {
  
}