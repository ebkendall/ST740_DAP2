library(gtools)
library(invgamma)
library(mvtnorm)

# -----------------------------------------------------------------------------
# The mcmc routine for sampling the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function(Y_mat, X, par, par_index, steps, burnin, n, m, alpha, mod_num){
    
    chain = matrix(0, nrow = steps, ncol = length(par))
    
    group = as.list(par_index$M_i)
    n_group = length(group)
    
    # Model selection values
    dev = rep(0, steps)
    logpdf1 = matrix(ncol = n, nrow = steps - burnin)
    curll = NULL
    
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
    pscale = rep( 0.0001, n_group)
    accept = rep( 0, n_group)
    
    # Hyperparameter initialization -------------------------------------------
    a = 2; b = par[par_index$sigma2]
    mean_mu = par[par_index$mu]
    s2 = 100
    
    # Begin the MCMC algorithm -------------------------------------------------
    chain[1,] = par
    for(ttt in 2:steps){
        
        # Gibbs steps ----------------------------------------------------------
        # Update theta
        theta_mat = matrix(par[par_index$theta], nrow = m, ncol = n)
        M_i = par[par_index$M_i]
        for(i in 1:n) {
            theta_mat[,i] = rdirichlet(1, Y_mat[,i] + exp(M_i[i])*alpha)
        }
        par[par_index$theta] = c(theta_mat)
        
        # Update mu
        par[par_index$mu] = update_mu(par, par_index, mean_mu, s2, n)
        
        # Update sigma2
        par[par_index$sigma2] = update_sigma2(par, par_index, a, b, n)
        
        chain[ttt, ] = par
        # Metropolis step ------------------------------------------------------
        if(mod_num != 2) {
            for(j in 1:n_group){
                
                # Propose an update
                ind_j = group[[j]]
                proposal = par
                if(length(ind_j > 1)) {
                    proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],sigma=pcov[[j]]*pscale[j])
                } else {
                    proposal[ind_j] = rnorm( n=1, mean=par[ind_j],sigma=pcov[[j]]*pscale[j])
                }
                
                # Evaluate the log posterior of the initial parameters
                ind_j_adj = which(par_index$M_i %in% ind_j)
                log_post_prev = fn_log_post(par, par_index, m, n, alpha, ind_j_adj)
                
                # Compute the log density for the proposal
                log_post = fn_log_post(proposal, par_index, m, n, alpha, ind_j_adj)
                
                # Only propose valid parameters during the burnin period
                if(ttt < burnin){
                    while(!is.finite(log_post)){
                        print('bad proposal')
                        proposal = par
                        if(length(ind_j > 1)) {
                            proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],sigma=pcov[[j]]*pscale[j])
                        } else {
                            proposal[ind_j] = rnorm( n=1, mean=par[ind_j],sigma=pcov[[j]]*pscale[j])
                        }
                        log_post = fn_log_post(proposal, par_index, m, n, alpha, ind_j_adj)
                    }
                }
                
                # Evaluate the Metropolis-Hastings ratio
                if( log_post - log_post_prev > log(runif(1,0,1)) ){
                    log_post_prev = log_post
                    par[ind_j] = proposal[ind_j]
                    accept[j] = accept[j] +1
                }
                chain[ttt,ind_j] = par[ind_j]
                
                # Proposal tuning scheme ------------------------------------------------
                if(ttt < burnin){
                    # During the burnin period, update the proposal covariance in each step
                    # to capture the relationships within the parameters vectors for each
                    # transition.  This helps with mixing.
                    if(ttt == 100)  pscale[j] = 1
                    
                    if (length(ind_j) > 1) {
                        if(100 <= ttt & ttt <= 2000){
                            temp_chain = chain[1:ttt,ind_j]
                            pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                            
                        } else if(2000 < ttt){
                            temp_chain = chain[(ttt-2000):ttt,ind_j]
                            pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                        }
                    }
                    
                    if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )
                    
                    # Tune the proposal covariance for each transition to achieve
                    # reasonable acceptance ratios.
                    if(ttt %% 30 == 0){
                        if(ttt %% 480 == 0){
                            accept[j] = 0
                            
                        } else if( accept[j] / (ttt %% 480) < .4 ){
                            pscale[j] = (.75^2)*pscale[j]
                            
                        } else if( accept[j] / (ttt %% 480) > .5 ){
                            pscale[j] = (1.25^2)*pscale[j]
                        }
                    }
                }
                # -----------------------------------------------------------------------
            }
        }
        
        
        # Model selection and evaluation
        curll = full_log_post(chain[ttt, ], par_index, m, n, alpha, Y_mat, 1:n, mod_num)
        dev[ttt]      <- -2*curll
        if(ttt>burnin){
            for(www in 1:n_group) {
                logpdf1[ttt-burnin, www] = full_log_post(chain[ttt, ], par_index, m, n, alpha, Y_mat, www, mod_num)
            }
        }
        
        # Restart the acceptance ratio at burnin.
        if(ttt == burnin)  accept = rep( 0, n_group)
        
        if(ttt%%1==0)  cat('--->',ttt,'\n')
    }
    # ---------------------------------------------------------------------------
    
    # DIC computations
    dbar <- mean(dev[burnin:steps])
    dhat <- -2*full_log_post(colMeans(chain[burnin:steps,]), par_index, m, n, alpha, Y_mat, 1:n, mod_num)
    pD   <- dbar-dhat
    DIC  <- list(dbar=dbar,dhat=dhat,pD=pD,DIC=dbar+pD) 
    
    # WAIC computations
    mn_logpdf  <- colMeans(logpdf1)
    var_logpdf <- apply(logpdf1, 2, sd)^2
    pW         <- sum(var_logpdf)
    WAIC       <- list(WAIC=-2*sum(mn_logpdf)+2*pW,pW=pW)
    
    print(accept/(steps-burnin))
    return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
                 pscale=pscale, pcov=pcov, DIC = DIC, WAIC = WAIC))
}
# -----------------------------------------------------------------------------



update_mu = function(par, par_index, mean_mu, s2, n) {
    sigma2 = par[par_index$sigma2]
    M_i = par[par_index$M_i]
    
    V = (sigma2 * s2) / (n*s2 + sigma2)
    M = V * ( (sum(M_i) / sigma2) + (mean_mu / s2) )
    
    return(rnorm(n = 1, mean = M, sd = sqrt(V)))
}

update_sigma2 = function(par, par_index, a, b, n) {
    M_i = par[par_index$M_i]
    mu = par[par_index$mu]
    temp = (M_i - mu)^2
    a1 = (n/2) + a
    b1 = b + 0.5 * sum(temp)
    
    return(rinvgamma(n = 1, shape = a1, rate = b1))
}

fn_log_post = function(par, par_index, m, n, alpha, ind_j_adj) {
    log_post_temp = 0
    theta_mat = matrix(par[par_index$theta], nrow = m, ncol = n)
    M_i = par[par_index$M_i]
    mu = par[par_index$mu]
    sigma2 = par[par_index$sigma2]
    for(i in ind_j_adj) {
        dir_temp = ddirichlet(x = theta_mat[,i], alpha = exp(M_i[i])*alpha)
        norm_temp = dnorm(x = M_i[i], mean = mu, sd = sqrt(sigma2))
        log_post_temp = log_post_temp + log(dir_temp) + log(norm_temp)
    }
    return(log_post_temp)
}

full_log_post = function(par, par_index, m, n, alpha, Y_mat, loop_it, mod_num) {
    log_post_temp = 0
    theta_mat = matrix(par[par_index$theta], nrow = m, ncol = n)
    M_i = par[par_index$M_i]
    mu = par[par_index$mu]
    sigma2 = par[par_index$sigma2]
    for(i in loop_it) {
        y_temp = dmultinom(x = Y_mat[,i], prob = theta_mat[,i])
        dir_temp = ddirichlet(x = theta_mat[,i], alpha = exp(M_i[i])*alpha)
        
        if (mod_num == 2) {
            log_post_temp = log_post_temp + log(y_temp) + log(dir_temp)
        } else {
            norm_temp = dnorm(x = M_i[i], mean = mu, sd = sqrt(sigma2))
            log_post_temp = log_post_temp + log(y_temp) + log(dir_temp) + log(norm_temp)
        }
        
        
    }
    return(log_post_temp)
}
