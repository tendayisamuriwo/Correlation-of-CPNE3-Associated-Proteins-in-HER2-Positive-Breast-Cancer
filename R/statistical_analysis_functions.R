# Function to plot correlation matrix with clustering
plot_correlation_with_clusters <- function(cor_matrix, plot_title) {
  corrplot(cor_matrix,
           method="circle",
           type="full",
           order="hclust",
           hclust.method="complete",
           addrect=3,
           rect.col="black",
           tl.col="black",
           tl.cex=0.6,
           main=plot_title,
           mar=c(0,0,3,0))  # Increased top margin for title
}

# Function to get Pearson correlation analysis for CPNE3
get_cpne3_correlations <- function(data, cohort_name) {
  data_numeric <- data[,-1]
  n <- nrow(data_numeric)
  
  # Get CPNE3 correlations
  cors <- cor(data_numeric)["CPNE3",]
  
  # Calculate p-values
  t_stat <- cors * sqrt((n-2)/(1-cors^2))
  p_values <- 2 * pt(abs(t_stat), n-2, lower.tail=FALSE)
  
  # Filter for significant correlations (excluding CPNE3 itself)
  significant_idx <- p_values < 0.05 & names(cors) != "CPNE3"
  
  cat("\nSignificant CPNE3 correlations in", cohort_name, "cohort:\n")
  cat("----------------------------------------\n")
  cat("Protein\tCorrelation\tP-Value\n")
  
  # Sort by absolute correlation value and print
  results_order <- order(abs(cors[significant_idx]), decreasing=TRUE)
  sig_proteins <- names(cors)[significant_idx][results_order]
  sig_cors <- cors[significant_idx][results_order]
  sig_pvals <- p_values[significant_idx][results_order]
  
  for(i in seq_along(sig_proteins)) {
    cat(sprintf("%s\t%.3f\t%.3e\n", sig_proteins[i], sig_cors[i], sig_pvals[i]))
  }
  
  return(list(
    proteins = sig_proteins,
    correlations = sig_cors,
    p_values = sig_pvals
  ))
}

# Function to get Spearman correlation analysis for CPNE3
get_spearman_cpne3_correlations <- function(data, cohort_name) {
  data_numeric <- data[,-1]
  n <- nrow(data_numeric)
  
  # Get CPNE3 correlations
  cors <- cor(data_numeric, method="spearman")["CPNE3",]
  
  # Calculate p-values using correlation test
  sig_results <- sapply(colnames(data_numeric), function(col) {
    if(col != "CPNE3") {
      test <- cor.test(data_numeric[,"CPNE3"], data_numeric[,col], 
                       method="spearman", exact=FALSE)
      return(c(cors[col], test$p.value))
    } else {
      return(c(1, 1))
    }
  })
  
  p_values <- sig_results[2,]
  
  # Filter for significant correlations (excluding CPNE3 itself)
  significant_idx <- p_values < 0.05 & names(cors) != "CPNE3"
  
  # Create output message
  output <- sprintf("\nSignificant CPNE3 Spearman correlations in %s cohort:\n", cohort_name)
  output <- paste0(output, "----------------------------------------\n")
  output <- paste0(output, "Protein\tCorrelation\tP-Value\n")
  
  # Sort by absolute correlation value and add to output
  results_order <- order(abs(cors[significant_idx]), decreasing=TRUE)
  sig_proteins <- names(cors)[significant_idx][results_order]
  sig_cors <- cors[significant_idx][results_order]
  sig_pvals <- p_values[significant_idx][results_order]
  
  for(i in seq_along(sig_proteins)) {
    output <- paste0(output, 
                     sprintf("%s\t%.3f\t%.3e\n", 
                             sig_proteins[i], sig_cors[i], sig_pvals[i]))
  }
  
  return(output)
}

# Function to summarise statistical assumptions
summarise_assumptions <- function(data, cohort_name) {
  cat("\nAnalysing", cohort_name, "\n")
  cat("----------------------------------------\n")
  
  # Remove first column (ID)
  data_numeric <- data[,-1]
  
  # Calculate normality pass rate
  shapiro_results <- apply(data_numeric, 2, function(x) {
    shapiro.test(x)$p.value > 0.05
  })
  normality_rate <- mean(shapiro_results) * 100
  cat(sprintf("Normality: %.1f%% passed\n", normality_rate))
  
  # Count strong correlations
  cor_matrix <- cor(data_numeric)
  strong_cors <- sum(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1) / 2
  cat(sprintf("Linearity: %d strong correlations\n", strong_cors))
  
  # Calculate homoscedasticity pass rate
  bp_tests <- 0
  bp_passed <- 0
  for(i in 1:(ncol(data_numeric)-1)) {
    for(j in (i+1):ncol(data_numeric)) {
      bp_tests <- bp_tests + 1
      bp <- ncvTest(lm(data_numeric[,i] ~ data_numeric[,j]))
      if(bp$p > 0.05) bp_passed <- bp_passed + 1
    }
  }
  homo_rate <- (bp_passed/bp_tests) * 100
  cat(sprintf("Homoscedasticity: %.1f%% passed\n", homo_rate))
  
  return(list(
    normality = round(normality_rate, 1),
    strong_correlations = strong_cors,
    homoscedasticity = round(homo_rate, 1)
  ))
}

# Function to perform detailed assumption testing
check_assumptions <- function(data, cohort_name) {
  # Remove first column (ID)
  data_numeric <- data[,-1]
  
  # 1. Normality Testing
  cat("\nNormality Tests for", cohort_name, "cohort:\n")
  cat("----------------------------------------\n")
  shapiro_results <- apply(data_numeric, 2, function(x) {
    test <- shapiro.test(x)
    return(c(W = test$statistic, p.value = test$p.value))
  })
  print(round(t(shapiro_results), 4))
  
  # 2. Linearity Check
  cat("\nLinearity Assessment:\n")
  cat("----------------------------------------\n")
  # Calculate correlation matrix
  cor_matrix <- cor(data_numeric)
  # Print strongest correlations
  strong_cors <- which(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1, arr.ind = TRUE)
  if(length(strong_cors) > 0) {
    cat("Strong correlations (|r| > 0.7):\n")
    for(i in 1:nrow(strong_cors)) {
      if(strong_cors[i,1] < strong_cors[i,2]) {
        cat(sprintf("%s - %s: r = %.3f\n", 
                    colnames(data_numeric)[strong_cors[i,1]], 
                    colnames(data_numeric)[strong_cors[i,2]], 
                    cor_matrix[strong_cors[i,1], strong_cors[i,2]]))
      }
    }
  }
  
  # 3. Homoscedasticity Testing
  cat("\nHomoscedasticity Tests:\n")
  cat("----------------------------------------\n")
  bp_results <- matrix(NA, ncol(data_numeric), ncol(data_numeric))
  colnames(bp_results) <- colnames(data_numeric)
  rownames(bp_results) <- colnames(data_numeric)
  
  for(i in 1:(ncol(data_numeric)-1)) {
    for(j in (i+1):ncol(data_numeric)) {
      bp <- ncvTest(lm(data_numeric[,i] ~ data_numeric[,j]))
      bp_results[i,j] <- bp$p
    }
  }
  print("Breusch-Pagan test p-values for selected pairs:")
  print(round(bp_results[1:5, 1:5], 4))  # Print first 5x5 subset for brevity
  
  return(list(shapiro = shapiro_results,
              homoscedasticity = bp_results))
}