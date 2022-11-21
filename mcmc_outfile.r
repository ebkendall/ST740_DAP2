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
dic = NULL
waic = NULL

for(seed in index_seeds){
  file_name = paste0('Model_out/mcmc_out_mod', mod_num, '_', toString(seed),'.rda')
  if (file.exists(file_name)) {
    load(file_name)
    ind = ind + 1
    print(paste0(ind, ": ", file_name))
    
    # DIC and WAIC
    dic = c(dic, mcmc_out$DIC$DIC)
    waic = c(waic, mcmc_out$WAIC$WAIC)
    
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

print("DIC"); print(mean(dic))
print("WAIC"); print(mean(waic))

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


# Answering the model comparison question -------------------------------------
if (mod_num == 1) {
    library(devtools)
    # install_github("microbiome/microbiome")
    library(microbiome)
    data(dietswap)
    Y_all <- otu_table(dietswap)    # all data
    X_all <- sample_data(dietswap)
    Y <- Y_all[,X_all[,7]==1]       # Extract baseline data
    X <- X_all[X_all[,7]==1,3]
    
    delta_j = matrix(nrow = nrow(stacked_chains), ncol = m)
    x_1 = which(X == "AAM")
    x_2 = which(X == "AFR")
    print(nrow(stacked_chains))
    for(i in 1:nrow(stacked_chains)) {
        theta_temp = matrix(stacked_chains[i, par_index$theta], ncol = n, nrow = m)
        
        theta_1 = theta_temp[,x_1]
        theta_2 = theta_temp[,x_2]
        
        dj_vec_1 = rowSums(theta_1) / length(x_1)
        dj_vec_2 = rowSums(theta_2) / length(x_2)
        
        delta_j[i, ] = c(dj_vec_1 - dj_vec_2)
    }

    load('Model_out/taxa_names.rda')
    cred_set = data.frame("taxa" = taxa_names, "lower" = rep(NA, m), "upper" = rep(NA, m))
    pdf(paste0('Plots/delta_', mod_num, '.pdf'))
    par(mfrow=c(3, 2))
    for(r in 1:m) {
        
        parMean = round( mean(delta_j[,r]), 4)
        parMedian = round( median(delta_j[,r]), 4)
        upper = quantile( delta_j[,r], prob=.975)
        lower = quantile( delta_j[,r], prob=.025)
        
        hist( delta_j[,r], breaks=sqrt(nrow(delta_j)), ylab=NA, main=NA,
              freq=F, xlab=paste0('Mean = ',toString(parMean),
                                  ' Median = ',toString(parMedian)))
        abline( v=upper, col='red', lwd=2, lty=2)
        abline( v=lower, col='purple', lwd=2, lty=2)
        
        cred_set$lower[r] = lower
        cred_set$upper[r] = upper
        
        print(paste0(r, ": Contains 0? ", (0 > lower & 0 < upper)))
    }
    dev.off()
    
    png("Plots/taxa_cred_set.png", width = 1600, height = 700)
    col_seq = rep(1, m)
    col_seq[which(cred_set$lower < 0 & cred_set$upper > 0)] = 2
    ggplot(cred_set, aes(x = taxa)) + 
        geom_errorbar(aes(ymin = lower, ymax = upper), color = factor(col_seq), position = position_dodge(1), size = 0.5, width=0.8) +
        scale_x_discrete(labels = taxa_names) +
        scale_colour_manual(values=col_seq) + 
        labs(title="95% Credible Sets", x ="", y="") +
        ylim(-0.03, 0.03) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, colour = col_seq))
    dev.off()
    
    
}



