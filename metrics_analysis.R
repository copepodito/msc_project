# Packages used.
library(stringdist)
library(data.table)
library(beepr)  # Optional.

# Save path to directory in a variable to enhance portability.
project.directory <- "/d/projects/u/cd002"
setwd(project.directory)
directories <- list.dirs(path = ".", full.names = F, recursive = F)

# In this version of the code all tables are attached from an external file
# named "tables". Attach tables to the R session before running anyything else.
attach(file.choose())

# Queries are stored in a .csv file.
queries <- data.table(read.csv("queries.csv",
                               header = T, stringsAsFactors = F))

# Use SingleQuery function to analyse Hamming distances of a single query.
# Function value (query.pdb) must be a character object (e.g. "3fku").
# Results are stored in a file named "data_[query.pdb]".
SingleQuery <- function(query.pdb){
  
  query <- queries$cdr3_aa[match(query.pdb, queries$pdb)]
  query.length <- queries$length[match(query.pdb, queries$pdb)]
  
  # There are 8 directories, one per individual.
  
  for (n in 1:8){
    
    start.time <- Sys.time()
    
    current.directory <- directories[n]
    setwd(paste(project.directory, current.directory, sep = "/"))
    
    dt <- paste("dt", current.directory, sep = ".")
    dt.IgM <- paste("dt.IgM", current.directory, sep = ".")
    dt.IgG <- paste("dt.IgG", current.directory, sep = ".")
    
    # Query analysis ----
    
    # Condition to be addressed: query length == CDR3 length.
    indexes.IgM <- paste("indexes.IgM", current.directory, sep = ".")
    assign(indexes.IgM, which(query.length == get(dt.IgM)$cdr3_length))
    
    indexes.IgG <- paste("indexes.IgG", current.directory, sep = ".")
    assign(indexes.IgG, which(query.length == get(dt.IgG)$cdr3_length))
    
    # Previous vectorisation speeds up the value allocation (instead of
    # creating a NULL vector).
    hamming.distances.IgM <- paste("hamming.distances.IgM", current.directory,
                                   sep = ".")
    assign(hamming.distances.IgM, numeric(length = length(get(indexes.IgM))))
    
    hamming.distances.IgG <- paste("hamming.distances.IgG", current.directory,
                                   sep = ".")
    assign(hamming.distances.IgG, numeric(length = length(get(indexes.IgG))))
    
    # stringdist() is case sensitive, but normalising it within the loop would 
    # hinder performance. Make sure query is in capital letters!
    # A temporary copy of the vectors containing Hamming distances is used to
    # allocate values.
    temp <- get(hamming.distances.IgM)
    for (i in 1:length(get(indexes.IgM))){
      temp[i] <- stringdist(query,
                            get(dt.IgM)$cdr3_aa[get(indexes.IgM)[i]],
                            "hamming")
    }
    assign(hamming.distances.IgM, temp)
    
    temp <- get(hamming.distances.IgG)
    for (i in 1:length(get(indexes.IgG))){
      temp[i] <- stringdist(query,
                            get(dt.IgG)$cdr3_aa[get(indexes.IgG)[i]],
                            "hamming")
    }
    assign(hamming.distances.IgG, temp)
    
    rm(temp, i)
    
    # Subset for a given similarity value or percentage (deprecated)
    # cut.off <- [a given percentage]
    # threshold <- query.length * (1 - cut.off)
    # threshold.indexes.IgM <- paste("threshold.indexes.IgM", current.directory,
    #                                sep = ".")
    # assign(threshold.indexes.IgM,
    #        get(indexes.IgM)[which(get(hamming.distances.IgM) < threshold)])
    # threshold.indexes.IgG <- paste("threshold.indexes.IgG", current.directory,
    #                                sep = ".")
    # assign(threshold.indexes.IgG,
    #        get(indexes.IgG)[which(get(hamming.distances.IgG) < threshold)])
    
    end.time <- Sys.time()
    print(end.time - start.time)
  }
  
  # Remove variables with irrelevant information.
  rm(current.directory, dt, dt.IgG, dt.IgM,
     hamming.distances.IgG, hamming.distances.IgM, indexes.IgG, indexes.IgM,
     n, start.time, end.time)
  
  # Save information in a file.
  setwd(project.directory)
  save(list = ls(pattern = "query|hamming|indexes"),
       file = paste("data", query.pdb, sep = "_"), compress = F)
}; beep()

# Use AllQueries function to analyse Hamming distances of all queries within
# the the queries.csv file. WARNING: It takes 1-2 hours per query.
AllQueries <- function(){
  for (i in queries$pdb){
    
    query.pdb <- i
    query <- queries$cdr3_aa[match(query.pdb, queries$pdb)]
    query.length <- queries$length[match(query.pdb, queries$pdb)]
    
    # There are 8 directories, one per individual.
    
    for (n in 1:8){
      
      start.time <- Sys.time()
      
      current.directory <- directories[n]
      setwd(paste(project.directory, current.directory, sep = "/"))
      
      dt <- paste("dt", current.directory, sep = ".")
      dt.IgM <- paste("dt.IgM", current.directory, sep = ".")
      dt.IgG <- paste("dt.IgG", current.directory, sep = ".")
      
      # Query analysis ----
      
      # Condition to be addressed: query length == CDR3 length.
      indexes.IgM <- paste("indexes.IgM", current.directory, sep = ".")
      assign(indexes.IgM, which(query.length == get(dt.IgM)$cdr3_length))
      
      indexes.IgG <- paste("indexes.IgG", current.directory, sep = ".")
      assign(indexes.IgG, which(query.length == get(dt.IgG)$cdr3_length))
      
      # Previous vectorisation speeds up the value allocation (instead of
      # creating a NULL vector).
      hamming.distances.IgM <- paste("hamming.distances.IgM", current.directory,
                                     sep = ".")
      assign(hamming.distances.IgM, numeric(length = length(get(indexes.IgM))))
      
      hamming.distances.IgG <- paste("hamming.distances.IgG", current.directory,
                                     sep = ".")
      assign(hamming.distances.IgG, numeric(length = length(get(indexes.IgG))))
      
      # stringdist() is case sensitive, but normalising it within the loop would 
      # hinder performance. Make sure query is in capital letters!
      # A temporary copy of the vectors containing Hamming distances is used to
      # allocate values.
      temp <- get(hamming.distances.IgM)
      for (i in 1:length(get(indexes.IgM))){
        temp[i] <- stringdist(query,
                              get(dt.IgM)$cdr3_aa[get(indexes.IgM)[i]],
                              "hamming")
      }
      assign(hamming.distances.IgM, temp)
      
      temp <- get(hamming.distances.IgG)
      for (i in 1:length(get(indexes.IgG))){
        temp[i] <- stringdist(query,
                              get(dt.IgG)$cdr3_aa[get(indexes.IgG)[i]],
                              "hamming")
      }
      assign(hamming.distances.IgG, temp)
      
      rm(temp, i)
      
      # Subset for a given similarity value or percentage (deprecated)
      # cut.off <- [a given percentage]
      # threshold <- query.length * (1 - cut.off)
      # threshold.indexes.IgM <- paste("threshold.indexes.IgM", current.directory,
      #                                sep = ".")
      # assign(threshold.indexes.IgM,
      #        get(indexes.IgM)[which(get(hamming.distances.IgM) < threshold)])
      # threshold.indexes.IgG <- paste("threshold.indexes.IgG", current.directory,
      #                                sep = ".")
      # assign(threshold.indexes.IgG,
      #        get(indexes.IgG)[which(get(hamming.distances.IgG) < threshold)])
      
      end.time <- Sys.time()
      print(end.time - start.time)
    }
    
    # Remove variables with irrelevant information.
    rm(current.directory, dt, dt.IgG, dt.IgM,
       hamming.distances.IgG, hamming.distances.IgM, indexes.IgG, indexes.IgM,
       n, start.time, end.time)
    
    # Save information in a file.
    setwd(project.directory)
    save(list = ls(pattern = "query|hamming|indexes"),
         file = paste("data", query.pdb, sep = "_"), compress = F)
  }
}