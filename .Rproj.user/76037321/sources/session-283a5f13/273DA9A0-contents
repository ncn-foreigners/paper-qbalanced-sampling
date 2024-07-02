library("data.table")
library("jointCalib") ## joint_calib_create_matrix
library("sampling") ## data and inclusionprobabilities
library("BalancedSampling") ## BalancedSampling::cube (faster than sampling)
library("laeken") ## weightedMean and weightedQuantile
library("xtable") ## xtable
library("doSNOW") ## parallel
library("progress") ## progress for parallel computing
library("foreach") ## foreach for parallel computing

seed_number <- 2024-7-22
set.seed(seed_number)

cores <- 8
R <- 10000
N <- 100000
x1 <- rexp(N,1)
x2 <- rexp(N,1)
alp <- rnorm(N)
epsilon <- rnorm(N)

y11 <- 1 + x1 + x2 + alp + epsilon
y12 <- 0.5*(x1-1.5)^2 + x2 + alp + epsilon
y21 <- rbinom(N,1,plogis(-1 + x1 + x2 + alp))
y22 <- rbinom(N,1,plogis(0.5*(x1-1.5)^2 + x2 + alp))
pop_data <- data.table(x1,x2,y11,y12,y21,y22)
pop_data[, pik_500:=500/N]
pop_data[, pik_1000:=1000/N]
#p_quantiles1 <- seq(0.25, 0.75, 0.25)
p_quantiles1 <- seq(0.10, 0.90, 0.10)
probs_est <- 0.9
vars_est <- c("y11", "y12", "y21", "y22")

X <- as.matrix(pop_data[, .(x1, x2)])

pop_quantiles <- list(x1=quantile(pop_data$x1, p_quantiles1), 
                      x2=quantile(pop_data$x2, p_quantiles1))

Xs <- joint_calib_create_matrix(X_q = X, 
                                N = N, 
                                pop_quantiles = pop_quantiles)

# population means and medians
trues <- pop_data[, lapply(.SD, function(x) c(mean = mean(x),
                                           median = median(x),
                                           quantile = quantile(x, probs = probs_est),
                                           sum = sum(x))), 
               .SDcols = vars_est][, type:=c("mean", "median", "quantile", "total")] |> 
  melt(id.vars = "type", value.name = "true")


# sample size 500 ---------------------------------------------------------

cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = R)
opts <- list(progress = \(n) pb$tick())

results500 <- foreach(
  k=1:R, .combine = rbind, 
  .packages = c("BalancedSampling", "jointCalib", "data.table", "laeken", "sampling"), 
  .options.snow = opts) %dopar% {

  ## standard balancing sampling
  sample1 <- BalancedSampling::cube(prob = pop_data$pik_500,  x = X)
  ## without calibration
  pop_data_sample1 <- copy(pop_data[sample1])
  pop_data_sample1[, g:=1]
  sample_cal <- pop_data_sample1[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g/pik_500),
                                                              median = weightedMedian(x = x, weights = g/pik_500),
                                                              quantile = weightedQuantile(x, weights = g/pik_500, probs = probs_est),
                                                              sum = sum(x*g/pik_500))), 
                                  .SDcols = vars_est][
                                    , type:=c("mean", "median", "quantile", "total")][ 
                                      , sample:="cal-no"] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")
  
  ## with calibration
  gcal_x <- calib(Xs = model.matrix(~ x1 + x2, pop_data_sample1),
             d = 1/pop_data_sample1$pik_500,
             total = c(N, colSums(X)),
             method = "raking")
  
  pop_data_sample1[, g1:=gcal_x]
  
  sample_cal2 <- pop_data_sample1[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g1/pik_500),
                                                             median = weightedMedian(x = x, weights = g1/pik_500),
                                                             quantile = weightedQuantile(x, weights = g1/pik_500, probs = probs_est),
                                                             sum = sum(x*g1/pik_500))), 
                      .SDcols = vars_est][
                        , type:=c("mean", "median", "quantile", "total")][ 
                          , sample:="cal-yes"] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")
  
  ## quantile balancing
  sample2 <- BalancedSampling::cube(prob = pop_data$pik_500, x = cbind(X,Xs))
  ## without calibration
  pop_data_sample2 <- copy(pop_data[sample2])
  pop_data_sample2[, g:=1]
  sample_qcal <- pop_data_sample2[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g/pik_500),
                                                              median = weightedMedian(x = x, weights = g/pik_500),
                                                              quantile = weightedQuantile(x, weights = g/pik_500, probs = probs_est),
                                                              sum = sum(x*g/pik_500))), 
                                  .SDcols = vars_est][
                                    , type:=c("mean", "median", "quantile", "total")][ 
                                      , sample:="qcal-no"] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")
  
  ## with calibration
  gcal <- joint_calib(formula_totals =  ~ x1 + x2,
                   formula_quantiles = ~ x1 + x2,
                   data = pop_data_sample2,
                   dweights = 1/pop_data_sample2$pik_500,
                   N = N,
                   pop_totals = colSums(X),
                   pop_quantiles = pop_quantiles,
                   method = "raking")
  
  pop_data_sample2[, g1:=gcal$g]
  
  sample_qcal2 <- pop_data_sample2[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g1/pik_500),
                                                              median = weightedMedian(x = x, weights = g1/pik_500),
                                                              quantile = weightedQuantile(x, weights = g1/pik_500, probs = probs_est),
                                                              sum = sum(x*g1/pik_500))), 
                       .SDcols = vars_est][
                         , type:=c("mean", "median", "quantile", "total")][ 
                           , sample:="qcal-yes"] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")
  
  df <- rbind(sample_cal, sample_cal2, sample_qcal, sample_qcal2)
  df$k <- k
  df
}

stopCluster(cl)


# sample size 1000 --------------------------------------------------------

cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = R)
opts <- list(progress = \(n) pb$tick())

results1000 <- foreach(
  k=1:R, .combine = rbind, 
  .packages = c("BalancedSampling", "jointCalib", "data.table", "laeken", "sampling"), 
  .options.snow = opts) %dopar% {
    
    ## standard balancing sampling
    sample1 <- BalancedSampling::cube(prob = pop_data$pik_1000,  x = X)
    ## without calibration
    pop_data_sample1 <- copy(pop_data[sample1])
    pop_data_sample1[, g:=1]
    sample_cal <- pop_data_sample1[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g/pik_1000),
                                                               median = weightedMedian(x = x, weights = g/pik_1000),
                                                               quantile = weightedQuantile(x, weights = g/pik_1000, probs = probs_est),
                                                               sum = sum(x*g/pik_1000))), 
                                   .SDcols = vars_est][
                                     , type:=c("mean", "median", "quantile", "total")][ 
                                       , sample:="cal-no"] |> 
      melt(id.vars = c("type", "sample"), value.name = "value")
    
    ## with calibration
    gcal_x <- calib(Xs = model.matrix(~ x1 + x2, pop_data_sample1),
                    d = 1/pop_data_sample1$pik_1000,
                    total = c(N, colSums(X)),
                    method = "raking")
    
    pop_data_sample1[, g1:=gcal_x]
    
    sample_cal2 <- pop_data_sample1[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g1/pik_1000),
                                                                median = weightedMedian(x = x, weights = g1/pik_1000),
                                                                quantile = weightedQuantile(x, weights = g1/pik_1000, probs = probs_est),
                                                                sum = sum(x*g1/pik_1000))), 
                                    .SDcols = vars_est][
                                      , type:=c("mean", "median", "quantile", "total")][ 
                                        , sample:="cal-yes"] |> 
      melt(id.vars = c("type", "sample"), value.name = "value")
    
    ## quantile balancing
    sample2 <- BalancedSampling::cube(prob = pop_data$pik_1000, x = cbind(X,Xs))
    ## without calibration
    pop_data_sample2 <- copy(pop_data[sample2])
    pop_data_sample2[, g:=1]
    sample_qcal <- pop_data_sample2[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g/pik_1000),
                                                                median = weightedMedian(x = x, weights = g/pik_1000),
                                                                quantile = weightedQuantile(x, weights = g/pik_1000, probs = probs_est),
                                                                sum = sum(x*g/pik_1000))), 
                                    .SDcols = vars_est][
                                      , type:=c("mean", "median", "quantile", "total")][ 
                                        , sample:="qcal-no"] |> 
      melt(id.vars = c("type", "sample"), value.name = "value")
    
    ## with calibration
    gcal <- joint_calib(formula_totals =  ~ x1 + x2,
                        formula_quantiles = ~ x1 + x2,
                        data = pop_data_sample2,
                        dweights = 1/pop_data_sample2$pik_1000,
                        N = N,
                        pop_totals = colSums(X),
                        pop_quantiles = pop_quantiles,
                        method = "raking")
    
    pop_data_sample2[, g1:=gcal$g]
    
    sample_qcal2 <- pop_data_sample2[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g1/pik_1000),
                                                                 median = weightedMedian(x = x, weights = g1/pik_1000),
                                                                 quantile = weightedQuantile(x, weights = g1/pik_1000, probs = probs_est),
                                                                 sum = sum(x*g1/pik_1000))), 
                                     .SDcols = vars_est][
                                       , type:=c("mean", "median", "quantile", "total")][ 
                                         , sample:="qcal-yes"] |> 
      melt(id.vars = c("type", "sample"), value.name = "value")
    
    df <- rbind(sample_cal, sample_cal2, sample_qcal, sample_qcal2)
    df$k <- k
    df
  }

stopCluster(cl)


# reporting ---------------------------------------------------------------

results500[trues, on = c("variable", "type"), trues := true]
results500 <- results500[!(type %in% c("quantile", "median") & variable %in% c("y21", "y22"))]
results500[, size := 500]
results1000[trues, on = c("variable", "type"), trues := true]
results1000 <- results1000[!(type %in% c("quantile", "median") & variable %in% c("y21", "y22"))]
results1000[, size := 1000]
results <- rbind(results500, results1000)

results[type != "total", .(m = (mean(value) - mean(trues))*100, 
                           var = (var(value))*100,
                           rmse = (sqrt( (mean(value) - mean(trues))^2 + var(value)))*100), 
        keyby=.(size, type, variable, sample)] |>
  melt(id.vars = c("size", "type", "sample", "variable")) |>
  transform(variable=paste(variable, variable.1, sep = "_"),
            variable.1 = NULL) |>
  transform(variable = factor(variable, 
                              levels = c("y11_m", "y11_var", "y11_rmse", 
                                         "y12_m", "y12_var", "y12_rmse",
                                         "y21_m", "y21_var", "y21_rmse",
                                         "y22_m", "y22_var", "y22_rmse"))) |>
  dcast(size + type + sample ~ variable, value.var = "value") |>
  {\(x) x[ , size := NULL]}() |>
  xtable() |>
  print.xtable(include.rownames = F)


results[type == "total", 
        .(m = mean(value) - mean(trues), 
          var = var(value)/1000,
          rmse = sqrt( (mean(value) - mean(trues))^2 + var(value))), 
        keyby=.(size, variable, sample)][, size := NULL] |>
  xtable() |>
  print.xtable(include.rownames = F, format.args = list(big.mark = ","))

saveRDS(results, "results/appendix-simulation.RDS")
