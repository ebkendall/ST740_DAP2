library(matrixStats)
library(mvtnorm)
library(MASS)
library(latex2exp)
library(LaplacesDemon, quietly=T)
library(evd)

args <- commandArgs(TRUE)
mod_num = as.numeric(args[1])

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at which the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 10000
# Matrix row indices for the posterior sample to use for GFF computation
index_post = (steps - burnin - n_post + 1):(steps - burnin)

m=92; n = 112

labels = NULL
if(mod_num == 2) {
    labels = c(TeX(r'($\mu$)'), TeX(r'($\sigma^2$)'), 
               paste0("theta(", 1:m, ", 1)")) 
} else {
    labels = c(TeX(r'($\mu$)'), TeX(r'($\sigma^2$)'), paste0("M(", 1:n, ")"),
               paste0("theta(", 1:m, ", 1)")) 
}

index_seeds = c(1:3)

true_par = NULL
par_index = list('mu' = 1, 'sigma2' = 2, 'M_i' = (2+1):(2+n),
                 'theta' = (2+n+1):(2+n+m*n))

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))

ind = 0

for(seed in index_seeds){
  file_name = paste0('Model_out/mcmc_out_mod', mod_num, '_', toString(seed),'.rda')
  if (file.exists(file_name)) {
    load(file_name)
    ind = ind + 1
    print(paste0(ind, ": ", file_name))
    
    # # Geweke diagnostic information ------------------------------------
    # diagnostic_res = Geweke.Diagnostic(mcmc_out$chain[18000:19000,])
    # print(paste0(length(which(is.nan(diagnostic_res))), " of the ", ncol(mcmc_out$chain), 
    #              " parameters lead to NaN"))
    # diagnostic_res = diagnostic_res[!is.nan(diagnostic_res)]
    # print(mean(abs(diagnostic_res) < 2))
    
    
    chain_list[[ind]] = mcmc_out$chain[index_post,]
    rm(mcmc_out)
  }
}

stacked_chains = do.call( rbind, chain_list)
# png(paste0('Plots/trace_plot_dap2_Mi0.png'), width = 1200, height = 1200)
# par(mfcol=c(4, 3))
pdf(paste0('Plots/trace_plot_mod', mod_num, '.pdf'))
par(mfrow=c(3, 2))
lab_ind = 0
fo_loop = NULL
if(mod_num == 2) {
    fo_loop = c(1:2, par_index$theta[1:m])
} else {
    fo_loop = c(1:2, par_index$M_i, par_index$theta[1:m])
}
for(r in fo_loop) {
  lab_ind = lab_ind + 1
  # stacked_chains[,r]
  plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,length(index_post)),
        ylim=range(stacked_chains[,r]),
        xlab="") #,  cex.lab=3, cex.axis=3, cex.main=4)
  
  for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)
  
  parMean = round( mean(stacked_chains[,r]), 4)
  parMedian = round( median(stacked_chains[,r]), 4)
  upper = quantile( stacked_chains[,r], prob=.975)
  lower = quantile( stacked_chains[,r], prob=.025)
  
  abline( h=parMean, col='purple', lwd=2, lty=2)
}
dev.off()

# # Summary plots for the posterior ---------------------------------------------
# library(ggplot2)
# library(viridis)
# library(latex2exp)
# file <- "https://www4.stat.ncsu.edu/~bjreich/ST740/HCDN_annual_max.RData"
# load(url(file))
# 
# keeper_rows = rowSums(is.na(Y))
# Y_sub = Y[which(keeper_rows == 0), ]
# 
# n = nrow(Y_sub)
# 
# load(paste0('Model_out/post_mcmc_out_dev_FINAL', 2, '_2.rda'))
# theta_2 = mcmc_out$chain[15000:20000, 257:492]
# 
# cred_set = matrix(nrow = ncol(theta_2), ncol = 2)
# for(i in 1:ncol(theta_2)) {
#   upper = quantile( theta_2[,i], prob=.95)
#   lower = quantile( theta_2[,i], prob=.05)
#   cred_set[i, ] = c(lower, upper)
# }
# 
# less_or_great = rep(0,nrow(cred_set))
# for(i in 1:nrow(cred_set)) {
#   if(cred_set[i,2] < 0) {
#     less_or_great[i] = -1
#   }
#   if(cred_set[i,1] > 0) {
#     less_or_great[i] = 1
#   }
# }
# 
# theta_2_mean = colMeans(theta_2)
# theta_2_prob = colMeans(theta_2 > 0)
# 
# s_1 = s[which(keeper_rows == 0), ]
# 
# # posterior mean
# legend_title <- TeX(r'($\theta_{2,s}$ mean)')
# dat <- data.frame(long=s_1[,1],lat=s_1[,2],b=theta_2_mean)
# ggplot(dat, aes(long, lat)) +
#   borders("state") +
#   geom_point(aes(colour = b)) +
#   scale_colour_gradientn(legend_title, colours = viridis(10)) +
#   coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
#   xlab("")+ylab("")+labs(title=TeX(r'(Posterior Means of $\theta_{2,s}$)'))
# 
# # posterior probability
# legend_title <- TeX(r'($P(\theta_{2,s} > 0)$)')
# dat <- data.frame(long=s_1[,1],lat=s_1[,2],b=theta_2_prob)
# ggplot(dat, aes(long, lat)) +
#   borders("state") +
#   geom_point(aes(colour = b)) +
#   scale_colour_gradientn(legend_title, colours = viridis(10)) +
#   coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
#   xlab("")+ylab("")+labs(title=TeX(r'(Posterior Probability of $\theta_{2,s} > 0$)'))
# 
# # significant slope
# legend_title <- TeX(r'(Sign of $\theta_{2,s}$)')
# dat <- data.frame(long=s_1[less_or_great!=0,1],lat=s_1[less_or_great!=0,2],b=less_or_great[less_or_great!=0])
# ggplot(dat, aes(long, lat)) +
#   borders("state") +
#   geom_point(aes(colour = as.factor(b))) +
#   scale_color_manual(name = legend_title,
#                      values = c("-1" = "red",
#                                 "1" = "blue"),
#                      labels = c("< 0", "> 0")) +
#   coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
#   xlab("")+ylab("")+labs(title=TeX(r'(Is 0 outside the 90\% Credible Set?)'))
