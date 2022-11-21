source('mcmc_routine.r')

args <- commandArgs(TRUE)

mod_num = as.numeric(args[1])

ind = args[2]
set.seed(ind)

print(ind)


# Loading Data   ---------------------------------------------------------------
library(devtools)
# install_github("microbiome/microbiome")
library(microbiome)
data(dietswap)
Y_all <- otu_table(dietswap)    # all data
X_all <- sample_data(dietswap)
Y <- Y_all[,X_all[,7]==1]       # Extract baseline data
X <- X_all[X_all[,7]==1,3]

steps  = 10000
burnin = 5000

r_name = rownames(Y)
c_name = colnames(Y)
Y_c = c(Y)

Y_mat = matrix(Y_c, nrow=130, ncol=112)
rownames(Y_mat) = r_name
colnames(Y_mat) = c_name

# Determining importance of variables in the model
relative_importance = rowSums(Y_mat); names(relative_importance) = NULL
relative_importance = relative_importance / sum(relative_importance)

important_taxa = which(relative_importance > 0.0001)
Y_mat = Y_mat[important_taxa, ]
min_list = NULL
for(i in 1:ncol(Y_mat)) {
     min_list = c(min_list, which.min(Y_mat[,i]))
}
min_list = unique(min_list)
Y_mat = Y_mat[-min_list, ]

n = ncol(Y_mat)
m = nrow(Y_mat)

taxa_names = rownames(Y_mat)
save(taxa_names, file = 'Model_out/taxa_names.rda')

N_i = colSums(Y_mat)

# alpha calc
sum_Y_i = rowSums(Y_mat)
which_zer0 = which(sum_Y_i == 0)
sum_Y_i[which_zer0] = 1
temp_N_i = N_i
temp_N_i[which_zer0] = temp_N_i[which_zer0] + 1

alpha = sum_Y_i / sum(temp_N_i)
names(alpha) = NULL


par_index = list('mu' = 1, 'sigma2' = 2, 'M_i' = (2+1):(2+n),
                 'theta' = (2+n+1):(2+n+m*n))
par = NULL
# initial values ---------------------------------------------------------------
if(mod_num != 2) {
    if(mod_num == 1) {
        load('Model_out/mcmc_out_mod2_3.rda')
        sub_chain = mcmc_out$chain[4000:5000, ]
        theta_est = colMeans(sub_chain[,par_index$theta])
        theta_est_mat = matrix(theta_est, nrow = m, ncol = n)
        
        M_i = rep(0, n)
        for(i in 1:n) {
            E_theta = mean(theta_est_mat[,i])
            V_theta = sd(theta_est_mat[,i])^2
            hold1 = E_theta * (1-E_theta) / V_theta
            M_i[i] = log((hold1 - 1) / sum(alpha))
        }
        
        mu = mean(M_i); sigma2 = sd(M_i)^2
        par = c(mu, sigma2, M_i, theta_est)
    } else { # Model (3)
        mu = 0; sigma2 = 0.1;
        theta = rep(alpha, n);
        M_i = rep(Inf, 112)
        par = c(mu, sigma2, M_i, theta)
        
        print("AIC")
        dic_calc = 0
        for(www in 1:n) {
            x_temp = Y_mat[,www]
            if(sum(x_temp == 0) > 0) x_temp = x_temp + 1
            dic_calc = dic_calc + log(dmultinom(x = x_temp, prob = alpha))
        }
        aic = -2 * dic_calc + 2*length(par_index$theta)
        print(aic)
        return(0)
    }
    
    
} else {
    E_theta = alpha; names(E_theta) = NULL
    M_i = rep(0, 112)
    mu = 0; sigma2 = 0.1;
    theta = rep(E_theta, n)
    par = c(mu, sigma2, M_i, theta)
}

names(par) = NULL
# -----------------------------------------------------------------------------

s_time = Sys.time()
mcmc_out = mcmc_routine(Y_mat, X, par, par_index, steps, burnin, n, m, alpha, mod_num)
e_time = Sys.time() - s_time; print(e_time)

save( mcmc_out, file=paste0('Model_out/mcmc_out_mod', mod_num, '_',ind,'.rda'))


# Checking Model Fit ----------------------------------------------------------
summary_stats = function(Y) {
    stats        <- c(rowMeans(Y), apply(Y, 1, sd), rowMeans(Y / colSums(Y)))
    return(stats)
}

if (mod_num == 1) {
    D0 = summary_stats(Y_mat)
    
    D <- NULL
    for(iter in 1:5000){
        print(iter)
        Y = matrix(nrow = nrow(Y_mat), ncol = ncol(Y_mat))
        theta = matrix(mcmc_out$chain[iter, par_index$theta], nrow = nrow(Y_mat), ncol = ncol(Y_mat))
        for(i in 1:n) {
            Y[,i] = rmultinom(n = 1, size = sum(Y_mat[,i]), prob = theta[,i])
        }
        D <- rbind(D,summary_stats(Y))
    }
    significance = rep(NA, ncol(D))
    for(j in 1:ncol(D)){
        pval <- round(mean(D0[j]>D[,j]),3)
        xlim <- quantile(D[,j],c(0.02,0.98))
        ex   <- (D[,j] < xlim[1]) | (D[,j]>xlim[2])
        hist(D[!ex,j],breaks=50,
             xlab=names(D0)[j],
             main=paste("p-value = ",pval))
        abline(v=D0[j],col=2,lwd=2)
        significance[j] = pval
    }
}
