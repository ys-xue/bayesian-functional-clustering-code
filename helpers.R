### helper functions

### data generation
genData <- function(seed, p1, signal = c("weak", "strong")) {
  set.seed(seed)
  income <- matrix(0, 10000, 51)
  if (signal == "weak") {
    income[, stateinfo$abbr %in% c1] <-
      rgamma(length(c1) * 10000, shape = 1.15, scale = 50000) + 
      rbinom(length(c1) * 10000, 1, p1) *
      rgamma(length(c1) * 10000, shape = 0.3, scale = 50000)
    income[, stateinfo$abbr %in% c2] <-
      rgamma(length(c2) * 10000, shape = 1.2, scale = 50000) + 
      rbinom(length(c2) * 10000, 1, p1) *
      rgamma(length(c2) * 10000, shape = 0.3, scale = 50000)
    income[, stateinfo$abbr %in% c3] <-
      rgamma(length(c3) * 10000, shape = 1.25, scale = 50000) + 
      rbinom(length(c3) * 10000, 1, p1) *
      rgamma(length(c3) * 10000, shape = 0.3, scale = 50000)
  } else if (signal == "strong") {
    income[, stateinfo$abbr %in% c1] <-
      rgamma(length(c1) * 10000, shape = 1.1, scale = 50000) + 
      rbinom(length(c1) * 10000, 1, p1) *
      rgamma(length(c1) * 10000, shape = 0.5, scale = 50000)
    income[, stateinfo$abbr %in% c2] <-
      rgamma(length(c2) * 10000, shape = 1.2, scale = 50000) + 
      rbinom(length(c2) * 10000, 1, p1) *
      rgamma(length(c2) * 10000, shape = 0.5, scale = 50000)
    income[, stateinfo$abbr %in% c3] <-
      rgamma(length(c3) * 10000, shape = 1.3, scale = 50000) + 
      rbinom(length(c3) * 10000, 1, p1) *
      rgamma(length(c3) * 10000, shape = 0.5, scale = 50000)
  } else {
    stop("Incorrect signal type.")
  }
  lorenzcurve <- matrix(0, 1000, 51)
  gini <- rep(0, 51)
  for (i in 1:51) {
    lorenzcurve[, i] <-
      quantile(Lc(income[, i])$L, seq(0.001, 1, by = 0.001))
    gini[i] <- Gini(income[, i])
  }
  
  srsf <- f_to_srvf(lorenzcurve, seq(0.001, 1, by = 0.001))
  
  inner_prod_matrix <- matrix(0, 51, 51)
  for (i in 1:51) {
    for (j in 1:51) {
      inner_prod_matrix[i, j] = 
        trapz(seq(0.001, 1, by = 0.001), srsf[, i] * srsf[, j])
    }
  }
  
  return(list(
    lorenzcurve = lorenzcurve,
    inner_prod_matrix = inner_prod_matrix,
    gini = gini,
    srvf = srsf
  ))
}






### four output lists
### 

summary_output1 <- function(outputList) {
  result1 <- purrr::map(outputList, ~getDahl(.x$fit1, 200))
  result2 <- purrr::map(outputList, ~getDahl(.x$fit2, 200))
  result3 <- purrr::map(outputList, ~getDahl(.x$fit3, 200))
  result4 <- purrr::map(outputList, ~getDahl(.x$fit4, 200))
  result5 <- purrr::map(outputList, ~getDahl(.x$fit5, 200))
  
  RI1 <- purrr::map_dbl(result1, ~fossil::rand.index(.x$zout, IND1))
  RI2 <- purrr::map_dbl(result2, ~fossil::rand.index(.x$zout, IND1))
  RI3 <- purrr::map_dbl(result2, ~fossil::rand.index(.x$zout, IND1))
  RI4 <- purrr::map_dbl(result3, ~fossil::rand.index(.x$zout, IND1))
  RI5 <- purrr::map_dbl(result5, ~fossil::rand.index(.x$zout, IND1))
  
  K1 <- data.frame(count = purrr::map_int(result1, ~length(unique(.x$zout))))
  K2 <- data.frame(count = purrr::map_int(result2, ~length(unique(.x$zout))))
  K3 <- data.frame(count = purrr::map_int(result3, ~length(unique(.x$zout))))
  K4 <- data.frame(count = purrr::map_int(result4, ~length(unique(.x$zout))))
  K5 <- data.frame(count = purrr::map_int(result5, ~length(unique(.x$zout))))
  
  finalK <- rbind(K1, K2, K3, K4, K5)
  finalK$smooth <- rep(c("0.5", "1", "1.5", "2", "MFM"), each = length(outputList))
  Klevels <- c(expression(paste(lambda, " = 0.5")),
               expression(paste(lambda, " = 1")),
               expression(paste(lambda, " = 1.5")),
               expression(paste(lambda, " = 2")),
               "MFM")
  finalK$smooth <- as.factor(finalK$smooth)
  levels(finalK$smooth) <- Klevels
  
  RImeans <- round(c(mean(RI1), mean(RI2), mean(RI3), mean(RI4), mean(RI5)), 3)
  correctK <- c(sum(K1 == 3), sum(K2 == 3), sum(K3 == 3), sum(K4 == 3), sum(K5 == 3))
  return(list(RImeans = RImeans, correctK = correctK, finalK = finalK))
}


summary_output2 <- function(outputList) {
  result1 <- purrr::map(outputList, ~getDahl(.x$fit1, 200))
  result2 <- purrr::map(outputList, ~getDahl(.x$fit2, 200))
  result3 <- purrr::map(outputList, ~getDahl(.x$fit3, 200))
  result4 <- purrr::map(outputList, ~getDahl(.x$fit4, 200))
  result5 <- purrr::map(outputList, ~getDahl(.x$fit5, 200))
  
  RI1 <- purrr::map_dbl(result1, ~fossil::rand.index(.x$zout, IND2))
  RI2 <- purrr::map_dbl(result2, ~fossil::rand.index(.x$zout, IND2))
  RI3 <- purrr::map_dbl(result2, ~fossil::rand.index(.x$zout, IND2))
  RI4 <- purrr::map_dbl(result3, ~fossil::rand.index(.x$zout, IND2))
  RI5 <- purrr::map_dbl(result5, ~fossil::rand.index(.x$zout, IND2))
  
  K1 <- data.frame(count = purrr::map_int(result1, ~length(unique(.x$zout))))
  K2 <- data.frame(count = purrr::map_int(result2, ~length(unique(.x$zout))))
  K3 <- data.frame(count = purrr::map_int(result3, ~length(unique(.x$zout))))
  K4 <- data.frame(count = purrr::map_int(result4, ~length(unique(.x$zout))))
  K5 <- data.frame(count = purrr::map_int(result5, ~length(unique(.x$zout))))
  
  finalK <- rbind(K1, K2, K3, K4, K5)
  finalK$smooth <- rep(c("1.5", "2", "2.5", "3", "MFM"), each = length(outputList))
  Klevels <- c(expression(paste(lambda, " = 1.5")),
               expression(paste(lambda, " = 2")),
               expression(paste(lambda, " = 2.5")),
               expression(paste(lambda, " = 3")),
               "MFM")
  finalK$smooth <- as.factor(finalK$smooth)
  levels(finalK$smooth) <- Klevels
  
  RImeans <- round(c(mean(RI1), mean(RI2), mean(RI3), mean(RI4), mean(RI5)), 3)
  correctK <- c(sum(K1 == 3), sum(K2 == 3), sum(K3 == 3), sum(K4 == 3), sum(K5 == 3))
  return(list(RImeans = RImeans, correctK = correctK, finalK = finalK))
}
