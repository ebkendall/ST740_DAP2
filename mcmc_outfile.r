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
    
    # Geweke diagnostic information ------------------------------------
    diagnostic_res = Geweke.Diagnostic(mcmc_out$chain[4000:5000,])
    print(paste0(length(which(is.nan(diagnostic_res))), " of the ", ncol(mcmc_out$chain),
                 " parameters lead to NaN"))
    diagnostic_res = diagnostic_res[!is.nan(diagnostic_res)]
    print(mean(abs(diagnostic_res) < 2))
    
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

if (mod_num == 1) {
    png(paste0('Plots/trace_plot_mod2_final.png'), width = 1200, height = 600)
    par(mfcol=c(2, 4))
    lab_ind = 0
    fo_loop = c(1:2, par_index$M_i[1], par_index$theta[1], par_index$M_i[2], par_index$theta[2],
                par_index$M_i[3], par_index$theta[3])
    new_labs = c(TeX(r'($\mu$)'), TeX(r'($\sigma^2$)'), TeX(r'($M_1$)'), TeX(r'($\theta_{1,1}$)'),
                 TeX(r'($M_2$)'), TeX(r'($\theta_{2,1}$)'), TeX(r'($M_3$)'), TeX(r'($\theta_{3,1}$)')) 
    for(r in fo_loop) {
        lab_ind = lab_ind + 1
        # stacked_chains[,r]
        plot( NULL, ylab=NA, main=new_labs[lab_ind], xlim=c(1,length(index_post)),
              ylim=range(stacked_chains[,r]),
              xlab="",  cex.lab=3, cex.axis=3, cex.main=4)
        
        for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)
        
        parMean = round( mean(stacked_chains[,r]), 4)
        parMedian = round( median(stacked_chains[,r]), 4)
        upper = quantile( stacked_chains[,r], prob=.975)
        lower = quantile( stacked_chains[,r], prob=.025)
        
        abline( h=parMean, col='purple', lwd=2, lty=2)
    }
    dev.off()
}

