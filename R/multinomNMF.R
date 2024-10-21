
multinomNNMF <- function(
    A,
    exposure,
    condition,
    metadata,
    lv_num = 2,
    learning_rate = .01,
    iterations =100, 
    alpha = .5) {
  
    # Create a sample non-negative matrix
    A <- as.matrix(A)
    
    # Initialize W and H with non-negative values
    set.seed(123)
    W <- matrix(runif(nrow(A) * lv_num, min = 0, max = 1), nrow = nrow(A), ncol = lv_num)
    H <- matrix(runif(lv_num * ncol(A), min = 0, max = 1), nrow = lv_num, ncol = ncol(A))
    
    # create vectors for loss values
    loss_values <- numeric(iterations)
    loss_rmse <- numeric(iterations)
    loss_causal <- numeric(iterations)
    
    # calculate the direct correlation
    direct_cor <- suppressWarnings({cor.test(meta_W[,condition],meta_W[,exposure], method = "spearman")$estimate})
    
    # progress bar
    pb = txtProgressBar(min=0, max=iterations, style = 3, initial = 0)
    
    # Perform NMF using multiplicative update rules
    for (i in 1:iterations) {
      
      # reconstructed
      WH <- W %*% H
      
      # gradient descent
      gradient_W_list = gradient_W(A, W, H, WH,exposure,condition,metadata,lv_num,alpha)
                  
      
      W <- W - learning_rate * gradient_W_list$combined_gradient
      H <- H - learning_rate * gradient_H(A, W, H,WH)
      
      W[W < 0] <- 0  # Ensure non-negativity
      H[H < 0] <- 0  # Ensure non-negativity
    

      
      loss_rmse[i] <- rmse(A,WH)
      loss_causal[i] <- causal_loss( exposure,condition,W, metadata, lv_num, direct_cor, gradient_W_list$greatest_mediators)
      
    #  print(gradient_W_list$greatest_mediators)
      
      setTxtProgressBar(pb,i)
    }
    
  close(pb)
  
  loss_rmse <- loss_rmse/loss_rmse[1]
  loss_causal <- loss_causal/loss_causal[1]
  loss_values <- alpha *loss_rmse+ (1-alpha) * loss_causal
  
  plot(loss_values, main = "total loss")
  plot(loss_rmse, main = "rmse")
  plot(loss_causal, main = "Partial Cor Loss")
  
  return(list(W = W, H = H, greatest_mediators = gradient_W_list$greatest_mediators))
}


rmse <- function(A, WH) {
  return(sqrt(mean((A - WH)^2)))
}



gradient_W <- function(A, W, H, WH, exposure,condition, metadata,lv_num, alpha) {
#  tmp <- A/WH
 # tmp[is.infinite(tmp)|is.na(tmp)] <- 0
  #gradient_1 <-  W * (tmp %*% t(H)) / rowSums(H)  # some infinite values based on dividing by zero
  # Update W and H
  gradient_1 <- -2 * (A - WH) %*% t(H)
  gradient_1 <- apply(gradient_1,2,min_max_scale)
  causal_gradient <- causalGradient(exposure,condition,W, metadata,lv_num)
  gradient_2 <- apply(causal_gradient$gradient,2,min_max_scale)
  combined_gradient <- alpha * gradient_1 + (1-alpha) * gradient_2
  return(list(combined_gradient = combined_gradient, greatest_mediators = causal_gradient$greatest_mediators))
}

gradient_H <- function(A, W, H, WH) {
  # tmp <- A/WH
  # tmp[is.infinite(tmp)|is.na(tmp)] <- 0
  # gradient_1<- H * (t(W) %*% tmp) / colSums(W)
  gradient_1 =  -2 * t(W) %*% (A - WH)
  gradient_1 <- t(apply(gradient_1,1,min_max_scale))
  
  
  return(gradient_1)
}

causal_loss <- function(exposure,condition,W, metadata, lv_num, direct_cor, mediators) {
  
  
  meta_W <- cbind(W[,mediators],metadata[, c(exposure,condition)])
  colnames(meta_W)[1:length(mediators)] <- paste0("LV",mediators)

  
#  meta_W <- cbind(W,metadata[, c(exposure,condition)])
 # colnames(meta_W)[1:lv_num] <- paste0("LV",lv_num)
  
  meta_W <- meta_W %>%
    na.omit()
  
  partial_cor <- spcor(meta_W, method = "spearman")$estimate[condition, exposure]
  loss = 1-(direct_cor-partial_cor)/direct_cor
  
  return(loss)
}


causalGradient <- function(exposure,condition,W, metadata,lv_num) {
  
  meta_W <- cbind(W,metadata[, c(exposure,condition)])
  colnames(meta_W)[1:lv_num] <- paste0("LV",1:lv_num)
  
  residual_differences <- list()
  percent_change <- numeric()
  
  for(i in 1:lv_num) {

    mod0 <- sprintf("%s ~%s", condition,paste0("LV",i))
    m0 <-  glm(mod0, data = meta_W, family = binomial(link = 'logit'))
      
    mod1 <- sprintf("%s ~%s+%s", condition,paste0("LV",i),exposure)
    m1 <-  glm(mod1, data = meta_W, family = binomial(link = 'logit'))
        
    residual_differences[[i]] <- m1$residuals-m0$residuals
    m0_estimate = summary(m0)$coefficients[2,4]
    m1_estimate = summary(m1)$coefficients[2,4]
    percent_change[i] = (m0_estimate - m1_estimate)/m0_estimate
    
    }
    
  # greatest_mediator 
  greatest_mediators <- order(-percent_change)[1]
  
  
    # combine residual differences to be added to the LV's to encourage divergence  
    dependent_residual_changes <- do.call(cbind, residual_differences[greatest_mediators])

    not_NA <-which(!(is.na(meta_W[,exposure])) & !(is.na(meta_W[, condition])))
    
    
    ## update LV's with significant mediators
    meta_W[not_NA,greatest_mediators] <- meta_W[not_NA,greatest_mediators]+dependent_residual_changes
    
    
      
  return(list(gradient = meta_W[,1:lv_num], greatest_mediators = greatest_mediators))
  
}


min_max_scale <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}




######### OLD CODE ##########
# And loss function is:

mse <- function(original, reconstructed){
  residuals = original - reconstructed
  return(mean(residuals ^ 2))
}



# KL Divergence Loss Function
kl_divergence_loss <- function(A, W, H,WH) {
  tmp <-  log(A + 1e-9 / WH+1e-9)
  tmp[is.infinite(tmp)] <- log(1e-9)
  loss <- sum(A * tmp - A + WH)
  return(loss)
}

