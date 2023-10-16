library(igraph)
library(ggplot2)

setwd("D:/Daniel_Jimenez/Data_Science_master_Sapienza_University/SDS/SDS_Homework3")
load("hw3_data.RData")

dim(asd_sel[[1]])
length(asd_sel)
td_sel[[1]]
dim(td_sel[[1]])

corrs_asd_sel[[1]]
dim(corrs_asd_sel[[1]])
length(td_sel)


get_estimates <-function(data,alpha = 0.05, appl_bonf = T, categories = 116){
  # data: Initial list of raw data to process
  # alpha: Significance level to calculate confidence intervals
  # appl_bonf: Boolean to select if the Bonferroni correction should be applied or not
  # categories: Number of categories considered for the Bonferroni correction
  
  # Calculate correlation matrix inside of each element of the list of data
  correlations <- lapply(data, cor)
  
  # Create matrix to store in each position a list with the correlations of each individual
  aux <- matrix(data = list(),nrow = dim(correlations[[1]])[1],ncol = dim(correlations[[1]])[2])
  
  # Iterate over each element of the list
  for (i in 1:length(data)){
    # Iterate over each dimension of the matrix of correlation
    for (j in 1:dim(correlations[[i]])[1]){
      for (k in 1:dim(correlations[[i]])[2]){
        
        # Store in the aux matrix 12 values (as a list) corresponding to 12 individuals
        aux[j, k] <- list(append(unlist(aux[j, k]),correlations[[i]][j, k],after = length(aux[j, k])))
        
      }
    }

  }
  
  # Initialize matrix of 80th percentiles
  percentiles <- matrix(data = NA,nrow = dim(correlations[[1]])[1],ncol = dim(correlations[[1]])[2])
  
  # Initialize matrix of confidence intervals
  conf_int <- matrix(data = list(),nrow = dim(correlations[[1]])[1],ncol = dim(correlations[[1]])[2])
  
  # Initialize matrix of adjacency for the final graph
  adj_mat <- matrix(data = NA,nrow = dim(correlations[[1]])[1],ncol = dim(correlations[[1]])[2])
  
  # Iterate over each row and column (each positon) of the matrix created before
  for (i in 1:dim(aux)[1]){
    for (j in 1:dim(aux)[2]){
      
      # Calculate the 80th percentile of the 12 individuals
      percentiles[i, j] <-  quantile(unlist(aux[i, j]), p = 0.8)
      
      # Perform Fisher Z-transformation for normality
      aux_fisher <-  0.5*(log(1 + unlist(aux[i, j])) - log(1 - unlist(aux[i, j])))
      
      # Correct infinite values caused by transformation
      aux_fisher[is.infinite(aux_fisher)] <- 1
      
      # Check if Bonferroni correction is declared
      if(appl_bonf == T){
        
        # Calculate the alpha corrected by Bonferroni
        alpha_bonf <- alpha / categories
        
        # Calculate the quantile of the normal distribution with the previous corrected significance level
        z_bonf <- qnorm(1 - alpha_bonf / 2)
        
        # Calculate the sample size
        n <- length(data)
        
        # Compute the confidence intervals intervals and revert the Fisher Z-transformation 
        conf_int[i, j] <- list(tanh(mean(aux_fisher) + c(-1,1)*z_bonf*sqrt((1/(n - 3))/n)))
        
      } else{
        # Calculate the quantile of the normal distribution with the original significance level
        z <- qnorm(1 - alpha / 2)
        
        # Calculate the sample size
        n <- length(data)
        
        # Compute the confidence intervals intervals and revert the Fisher Z-transformation 
        conf_int[i, j] <- list(tanh(mean(aux_fisher) + c(-1,1)*z*sqrt((1/(n - 3))/n)))
        
      }
      
      # Extract the current lower bound of the confidence interval 
      low <- unlist(conf_int[i, j])[1]
      
      # Extract the current upper bound of the confidence interval
      upp <- unlist(conf_int[i, j])[2]
      
      # Extract the current threshold
      t <- percentiles[i, j]
      
      # Evaluate whether the threshold is contained on the confidence intervals or not
      if((t >= low & t <= upp) & (-t >= low & -t <= upp)){
        # Do not place an edge
        adj_mat[i, j] <- 0
      } else {
        # Place an edge
        adj_mat[i, j] <- 1
      }

    }
      
  }
  
  # Update the main diagonal of the adjacency matrix with 0 to avoid nodes connected to itself 
  diag(adj_mat) <- 0  
  
  return(adj_mat)
}

# Calculate adjacency matrix applying Bonferroni correction
adjacency_matrix <- get_estimates(asd_sel, categories = 116)

# Build the graph
G <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
# Calculate number of edges
E <- length(E(G))
# Calculate number of vertexes
V <- length(V(G))
# Plot graph
plot(G, vertex.color="dodgerblue4", main=paste("Undirected Graph with ", V," vertexes and ", E, " edges",sep=""), vertex.label.color="white")

# Store the degrees of the graph when appyting Bonferroni
degrees_with_bonf <- data.frame(bonferroni = "With Bonferroni", degrees = degree(G))

# Calculate density of the graph
D <- 2*E / (V*(V-1))
D
# Calculate adjacency matrix without applying Bonferroni correction
adjacency_matrix <- get_estimates(asd_sel, categories = 116, appl_bonf = F)

# Build the graph
G <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
# Calculate number of edges
E <- length(E(G))
# Calculate number of vertexes
V <- length(V(G))
# Plot graph
plot(G, vertex.color="dodgerblue4", main=paste("Undirected Graph with ", V," vertexes and ", E, " edges",sep=""), vertex.label.color="white")

# Store the degrees of the graph when avoiding Bonferroni
degrees_without_bonf <- data.frame(bonferroni = "Without Bonferroni", degrees = degree(G))

# Calculate density of the graph
D <- 2*E / (V*(V-1))
D
# Plot densities of the degrees when applying/avoiding Bonferroni correction
df_to_plot <- rbind(degrees_with_bonf,degrees_without_bonf)
ggplot(df_to_plot,aes(x=degrees, fill=bonferroni)) + geom_density(alpha=0.7) + scale_fill_manual(values=c("#2e598f","#940f15"))




# Calculate adjacency matrix applying Bonferroni correction
adjacency_matrix <- get_estimates(td_sel, categories = 116)

# Build the graph
G <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
# Calculate number of edges
E <- length(E(G))
# Calculate number of vertexes
V <- length(V(G))
# Plot graph
plot(G, vertex.color="dodgerblue4", main=paste("Undirected Graph with ", V," vertexes and ", E, " edges",sep=""), vertex.label.color="white")

# Store the degrees of the graph when appyting Bonferroni
degrees_with_bonf <- data.frame(bonferroni = "With Bonferroni", degrees = degree(G))

# Calculate density of the graph
D <- 2*E / (V*(V-1))
D
# Calculate adjacency matrix without applying Bonferroni correction
adjacency_matrix <- get_estimates(td_sel, categories = 116, appl_bonf = F)

# Build the graph
G <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
# Calculate number of edges
E <- length(E(G))
# Calculate number of vertexes
V <- length(V(G))
# Plot graph
plot(G, vertex.color="dodgerblue4", main=paste("Undirected Graph with ", V," vertexes and ", E, " edges",sep=""), vertex.label.color="white")

# Store the degrees of the graph when avoiding Bonferroni
degrees_without_bonf <- data.frame(bonferroni = "Without Bonferroni", degrees = degree(G))

# Calculate density of the graph
D <- 2*E / (V*(V-1))
D
# Plot densities of the degrees when applying/avoiding Bonferroni correction
df_to_plot <- rbind(degrees_with_bonf,degrees_without_bonf)
ggplot(df_to_plot,aes(x=degrees, fill=bonferroni)) + geom_density(alpha=0.7) + scale_fill_manual(values=c("#2e598f","#940f15"))






reshape_data <- function(corr_mat){
  # corr_mat: List of correlation matrix to reshape
  
  # Create matrix to store in each position a list with the correlations of each individual
  aux <- matrix(data = list(),nrow = dim(corr_mat[[1]])[1],ncol = dim(corr_mat[[1]])[2])
  
  # Iterate over each element of the list
  for (i in 1:length(corr_mat)){
    # Iterate over each dimension of the matrix of correlation
    for (j in 1:dim(corr_mat[[i]])[1]){
      for (k in 1:dim(corr_mat[[i]])[2]){
        
        # Store in the aux matrix n values (as a list) corresponding to n individuals
        aux[j, k] <- list(append(unlist(aux[j, k]),corr_mat[[i]][j, k],after = length(aux[j, k])))
        
      }
    }
    
  }
  
  return(aux)
}
# Bootstrap
# Set number of iterations
M <- 5

# Set sample size in each iteration
n <- 10

# Define the seed for replicability
seed <- 1234

# Set threshold
t <- 0.5

# Set number of categories for Bonferroni correction
categories <- 116

# Significance level to calculate confidence intervals
alpha <- 0.05

# Pre-allocate matrix to store mean of the correlations
mean_corr_boot <- vector(mode = "list", length = M)

# Loop over desired Bootstrap size
for(l in 1:M){
  # Get sample from ASD
  set.seed(seed)
  sample_asd <-sample(asd_sel, n, replace=TRUE)
  
  # Get sample from TD
  set.seed(seed)
  sample_td <-sample(td_sel, n, replace=TRUE)
  
  # Calculate correlations inside of list for each sample
  correlations_asd <- lapply(sample_asd, cor)
  correlations_td  <- lapply(sample_td, cor)
  
  # Calculate difference of correlation in absolute value
  delta_correlation <- lapply(mapply('-', correlations_asd, correlations_td, SIMPLIFY = FALSE),abs)
  
  # Reshape data to get a matrix with list inside of each position
  aux <- reshape_data(delta_correlation)
  
  # Calculate the mean of the correlation values for each position
  mean_corr <- matrix(data = NA,nrow = dim(delta_correlation[[1]])[1],ncol = dim(delta_correlation[[1]])[2])
  
  # Loop over the matrix of means of correlation
  for (i in 1:dim(aux)[1]){
    for (j in 1:dim(aux)[2]){
      
      # Calculate the mean over correlations and store in matrix of Boostrap means
      mean_corr[i, j] <-  mean(unlist(aux[i, j]), na.rm = T)
    }
  }
  
  # Save the matrix of means
  mean_corr_boot[[l]] <- mean_corr
  
  # Modify seed
  seed <- seed + 100
}

# Reshape data to get a matrix with list of means inside of each position
aux_means_boot <- reshape_data(mean_corr_boot)

# Pre-allocate matrix to store confidence intervals
conf_int_boot <- matrix(data = list(),nrow = dim(mean_corr_boot[[1]])[1],ncol = dim(mean_corr_boot[[1]])[2])
# Pre-allocate matrix of adjacency 
adj_mat_boot <- matrix(data = NA,nrow = dim(mean_corr_boot[[1]])[1],ncol = dim(mean_corr_boot[[1]])[2])

# Lopp over means of Bootstrap
for (i in 1:dim(aux_means_boot)[1]){
  for (j in 1:dim(aux_means_boot)[2]){
    # Calculate the alpha corrected by Bonferroni
    alpha_bonf <- alpha / categories
  
    # Calculate the quantile of the normal distribution with the previous corrected significance level    
    z_bonf <- qnorm(1 - alpha_bonf / 2)
    
    # Compute the confidence intervals intervals
    conf_int_boot[i, j] <- list(mean(unlist(aux_means_boot[i, j])) + c(-1,1)*z_bonf*sd(unlist(aux_means_boot[i, j])))
    
    # Extract the current lower bound of the confidence interval   
    low <- unlist(conf_int_boot[i, j])[1]
    
    # Extract the current lower bound of the confidence interval 
    upp <- unlist(conf_int_boot[i, j])[2]
    
    # Evaluate whether the threshold is contained on the confidence intervals or not
    if((t >= low & t <= upp)){
      # Do not place an edge
      adj_mat_boot[i, j] <- 0
    } else {
      # Place an edge
      adj_mat_boot[i, j] <- 1
    }
    
  }
  
}

# Update the main diagonal of the adjacency matrix with 0 to avoid nodes connected to itself 
diag(adj_mat_boot) <- 0  

# Plot graph
G <- graph_from_adjacency_matrix(adj_mat_boot, mode = "undirected")
# Calculate number of edges
E <- length(E(G))
# Calculate number of vertexes
V <- length(V(G))
# Plot graph
plot(G, vertex.color="dodgerblue4", main=paste("Undirected Graph with ", V," vertexes and ", E, " edges",sep=""), vertex.label.color="white")

# Calculate density of the graph
D <- 2*E / (V*(V-1))

# summary(adj_mat_boot)
# length(sample_asd)
# dim(sample_asd[[1]])
# 
# 
# aux_means_boot[[1]]
# 
# aux[[1,2]]
# # View(delta_correlation[[1]])
# length(mean_corr_boot)
# mean_corr_boot[[1]]