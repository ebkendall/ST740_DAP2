source('mcmc_routine.r')

args <- commandArgs(TRUE)
ind = args[1]
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

steps  = 20000
burnin = 1000

r_name = rownames(Y)
c_name = colnames(Y)
Y_c = c(Y)

n = ncol(Y)
m = nrow(Y)

Y_mat = matrix(Y_c, nrow=m, ncol=n)
rownames(Y_mat) = r_name
colnames(Y_mat) = c_name

N_i = colSums(Y_mat)

# alpha calc
sum_Y_i = rowSums(Y_mat)
which_zer0 = which(sum_Y_i == 0)
sum_Y_i[which_zer0] = 1
temp_N_i = N_i
temp_N_i[which_zer0] = temp_N_i[which_zer0] + 1

alpha = sum_Y_i / sum(temp_N_i)

# initial values ---------------------------------------------------------------
E_theta = alpha / sum(alpha); names(E_theta) = NULL

# prop_Y = Y_mat / colSums(Y_mat)
# M_i = apply(prop_Y, 2, sd); names(M_i) = NULL
# M_i = apply(Y_mat, 2, sd); names(M_i) = NULL
# M_i = log(M_i)
# mu = mean(M_i); sigma2 = sd(M_i)^2

M_i = rep(0, 112)
mu = 0; sigma2 = 0.1;

theta = rep(E_theta, n)

par = c(mu, sigma2, M_i, theta)
names(par) = NULL
# -----------------------------------------------------------------------------

par_index = list('mu' = 1, 'sigma2' = 2, 'M_i' = (2+1):(2+n),
                 'theta' = (2+n+1):(2+n+m*n))

s_time = Sys.time()
mcmc_out = mcmc_routine(Y_mat, X, par, par_index, steps, burnin, n, m, alpha)
e_time = Sys.time() - s_time; print(e_time)

save( mcmc_out, file=paste0('Model_out/mcmc_out_mod2_',ind,'.rda'))
