# This script reads the information from every "minimal_info.csv" file 
# in the directoreies of each subject and creates a file with data tables
# that will be crucial for the rest of the analysis.

# Save path to directory in a variable to enhance portability.
project.directory <- "/d/projects/u/cd002"

# Packages used.
library(data.table)
library(beepr)

setwd(project.directory)
directories <- list.dirs(path = ".", full.names = F, recursive = F)

# Queries are stored in a .csv file
queries <- data.table(read.csv("queries.csv",
                               header = T, stringsAsFactors = F))

# There are 8 directories, one per individual.

for (n in 1:8){
  
  start.time <- Sys.time()
  
  current.directory <- directories[n]
  setwd(paste(project.directory, current.directory, sep = "/"))
  
  # There are 18 .txt files in each directory and 1 .csv file. We will use this
  # last one, named "minimal_info.csv".
  
  # Read data ----
  
  # Data tables (dt) are enhanced data frames objects. Create three data
  # tables:
  # - All information in minimal_info.csv
  # - dt for IgM
  # - dt for IgG
  dt <- paste("dt", current.directory, sep = ".")
  assign(dt, data.table(read.csv("minimal_info.csv",
                                 header = T, stringsAsFactors = F)))
  
  dt.IgM <- paste("dt.IgM", current.directory, sep = ".")
  assign(dt.IgM, subset(get(dt), isotype == "IgM",
                        select = colnames(get(dt))
                        [-which(colnames(get(dt)) == "isotype")]))
  
  dt.IgG <- paste("dt.IgG", current.directory, sep = ".")
  assign(dt.IgG, subset(get(dt),
                        isotype %in% c("IgG1", "IgG2", "IgG3", "IgG4"),
                        select = colnames(get(dt))
                        [-which(colnames(get(dt)) == "isotype")]))
  end.time <- Sys.time()
  print(end.time - start.time)
}

# Remove useless objects.
rm(dt, dt.IgM, dt.IgG, start.time, end.time, directories, current.directory, n)

# Save all objects.
setwd(project.directory)
save(list = ls(), file = "tables", compress = FALSE)