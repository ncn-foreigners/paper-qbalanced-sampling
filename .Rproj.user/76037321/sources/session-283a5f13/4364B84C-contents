library("data.table")
library("jointCalib") ## joint_calib_create_matrix
library("sampling") ## data and inclusionprobabilities
library("BalancedSampling") ## BalancedSampling::cube (faster than sampling)
library("laeken") ## weightedMean and weightedQuantile
library("xtable") ## xtable
library("ggplot2")
library("scales") 

data(MU284, package = "sampling")
setDT(MU284)
## outliers
#MU284 <- MU284[!LABEL %in% c(137,16,114)]
MU284[, pik:= inclusionprobabilities(P75, 50)]

X <- as.matrix(MU284[, .(P75,CS82,SS82,S82,ME84,REV84)])
vars  <- c("P75", "CS82", "SS82", "S82", "ME84", "REV84", "RMT85")
p_quantiles1 <- 0.5
p_quantiles2 <- c(0.25,0.5,0.75)
p_quantiles3 <- seq(0.10, 0.90, 0.10)
probs_est <- 0.9

pop_quantiles1 <- list(quantile(X[,1], p_quantiles1),
                       quantile(X[,2], p_quantiles1),
                       quantile(X[,3], p_quantiles1),
                       quantile(X[,4], p_quantiles1),
                       quantile(X[,5], p_quantiles1),
                       quantile(X[,6], p_quantiles1))

pop_quantiles2 <- list(quantile(X[,1], p_quantiles2),
                       quantile(X[,2], p_quantiles2),
                       quantile(X[,3], p_quantiles2),
                       quantile(X[,4], p_quantiles2),
                       quantile(X[,5], p_quantiles2),
                       quantile(X[,6], p_quantiles2))

names(pop_quantiles2) <- names(pop_quantiles1) <- vars[1:6]

Xs_median <- joint_calib_create_matrix(X_q = X, N = nrow(X), pop_quantiles = pop_quantiles1)
Xs_quartiles <- joint_calib_create_matrix(X_q = X, N = nrow(X), pop_quantiles = pop_quantiles2)

trues <- MU284[, lapply(.SD, function(x) c(mean = mean(x),
                                           median = median(x),
                                           quantile = quantile(x, probs = probs_est),
                                           sum = sum(x))), 
                  .SDcols = vars][, type:=c("mean", "median", "quantile", "total")] |> 
  melt(id.vars = "type", value.name = "true")

trues_quants <- MU284[, lapply(.SD, quantile, probs = p_quantiles3), .SDcols = vars][
  , ":="(type=p_quantiles3)] |> 
  melt(id.vars = c("type"), value.name = "true")

R <- 10000
results <- list()
results_quants <- list()
for (r in 1:R) {
  set.seed(2024+r)
  if (r %% 100 == 0) print(r)
  
  # standard balanced sampling ----------------------------------------------
  sample1 <- BalancedSampling::cube(prob = MU284$pik,  x = X)
  MU284_sample <- copy(MU284[sample1])

  # g <- sampling::calib(Xs = model.matrix(~ P75 + CS82 + SS82 + S82 + ME84 + REV84, MU284_sample), 
  #                      d = 1/MU284_sample$pik,
  #                      total = c(nrow(MU284), colSums(X)), method = "raking")
  # MU284_sample[, g:=g]
  MU284_sample[, g:=1]
  
  sample_cal <- MU284_sample[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g/pik),
                                                         median = weightedMedian(x = x, weights = g/pik),
                                                         quantile = weightedQuantile(x, weights = g/pik, probs = probs_est),
                                                         sum = sum(x*g/pik),
                                                         var_tot = sampling::varest(Y=x, pik = pik))), 
                    .SDcols = vars][
                      , type:=c("mean", "median", "quantile", "total", "var_tot")][ 
                        , sample:="cal"] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")
  
  sample_cal_quants <- MU284_sample[, lapply(.SD, weightedQuantile, weights = g/pik, probs = p_quantiles3), .SDcols = vars][
    , ":="(type=p_quantiles3,sample="cal")] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")

  # quantile balanced sampling with medians ----------------------------------------------
  sample2 <- BalancedSampling::cube(prob = MU284$pik, x = cbind(X,Xs_median))
  MU284_sample2 <- copy(MU284[sample2])
  
  # g <- sampling::calib(Xs = model.matrix(~ P75 + CS82 + SS82 + S82 + ME84 + REV84, MU284_sample2), 
  #                      d = 1/MU284_sample2$pik,
  #                      total = c(nrow(MU284), colSums(X)), method = "raking")
  # MU284_sample2[, g:=g]
  MU284_sample2[, g:=1]
  sample_qcal <- MU284_sample2[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g/pik),
                                                          median = weightedMedian(x = x, weights = g/pik),
                                                          quantile = weightedQuantile(x, weights = g/pik, probs = probs_est),
                                                          sum = sum(x*g/pik),
                                                          var_tot = sampling::varest(Y=x, pik = pik))), 
                          .SDcols = vars][
                            , type:=c("mean", "median", "quantile", "total", "var_tot")][ 
                              , sample:="qcal"] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")
  
  sample_qcal_quants <- MU284_sample2[, lapply(.SD, weightedQuantile, weights = g/pik, probs = p_quantiles3), .SDcols = vars][
    , ":="(type=p_quantiles3,sample="qcal")] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")

  # quantile balanced sampling with deciles ----------------------------------------------
  sample3 <- BalancedSampling::cube(prob = MU284$pik, x = cbind(X,Xs_quartiles))
  MU284_sample3 <- copy(MU284[sample3])

  # g <- sampling::calib(Xs = model.matrix(~ P75 + CS82 + SS82 + S82 + ME84 + REV84, MU284_sample3), 
  #                      d = 1/MU284_sample3$pik,
  #                      total = c(nrow(MU284), colSums(X)), method = "raking")
  # MU284_sample3[, g:=g]
  MU284_sample3[, g:=1]
  
  sample_qcal2 <- MU284_sample3[, lapply(.SD, function(x) c(mean = weighted.mean(x=x, w = g/pik),
                                                           median = weightedMedian(x = x, weights = g/pik),
                                                           quantile = weightedQuantile(x, weights = g/pik, probs = probs_est),
                                                           sum = sum(x*g/pik),
                                                           var_tot = sampling::varest(Y=x, pik = pik))), 
                               .SDcols = vars][
                                 , type:=c("mean", "median", "quantile", "total", "var_tot")][ 
                                   , sample:="qcal2"] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")
  
  sample_qcal_quants2 <- MU284_sample3[, lapply(.SD, weightedQuantile, weights = g/pik, probs = p_quantiles3), .SDcols = vars][
    , ":="(type=p_quantiles3,sample="qcal2")] |> 
    melt(id.vars = c("type", "sample"), value.name = "value")
  
  results[[r]] <- rbind(sample_cal, sample_qcal, sample_qcal2)
  results_quants[[r]] <- rbind(sample_cal_quants, sample_qcal_quants,sample_qcal_quants2)
  
}


# main results ------------------------------------------------------------


results_all <- rbindlist(results, idcol = 'rep')
results_all[trues, on = c("variable", "type"), trues := true]

results_all[variable == "RMT85" & type != "var_tot", .(m = mean(value) - mean(trues), 
                               var = var(value),
                               rmse = sqrt( (mean(value) - mean(trues))^2 + var(value))), 
            keyby=.(type, sample)] |>
  xtable(digits = 2) |>
  print.xtable(include.rownames = F,
               format.args = list(big.mark = ","))

# quantiles ----------------------------------------------------------------

results_quants_all <- rbindlist(results_quants, idcol = 'rep')

saveRDS(results_quants_all, file = "results/main-paper-simulation.RDS")
# plot --------------------------------------------------------------------


results_quants_all[variable !="RMT85", 
                   .(cv=sd(value)/mean(value)), .(variable, type, sample)] |>
  ggplot(data = _, aes(x = type, y = cv, 
                       group = sample, 
                       linetype = sample, shape = sample)) +
  geom_line() +
  geom_point() + 
  facet_wrap(~variable) +
  scale_x_continuous(breaks = p_quantiles3) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Quantile", y = "Monte Carlo CV", linetype = 'Method', shape = "Method") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), 
                        labels = c("Balanced", "QB medians", "QB quartiles")) +
  scale_shape_manual(values = c(1, 2, 3), 
                     labels = c("Balanced", "QB medians", "QB quartiles")) -> p1

ggsave(plot = p1, filename = "figures/figure-sim.pdf", width = 10, height = 5)
