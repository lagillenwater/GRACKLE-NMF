calculate_eigengap <- function(data) {
  # Calculate the covariance matrix
  covariance_matrix <- cov(data)
  
  # Calculate the eigenvalues
  eigenvalues <- eigen(covariance_matrix)$values
  
  # Sort the eigenvalues in descending order
  sorted_eigenvalues <- sort(eigenvalues, decreasing = TRUE)
  
  # Calculate the eigengap (difference between the largest and second largest eigenvalues)
  eigengap <- sorted_eigenvalues[1] - sorted_eigenvalues[2]
  
  return(eigengap)
}