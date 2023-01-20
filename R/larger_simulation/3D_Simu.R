library(tidyverse)
library(tidybayes)
library(rstan)
library(cmdstanr)
library(car)
library(ggpubr)
source('R/to_docs/mcmc_utils.R')
options(mc.cores = parallel::detectCores())
rstan_options(javascript=FALSE)

model.odr = cmdstan_model('stan_models/calibration_ODR.stan')
model.odr$code()
#library(hash)

init_fun <- function() {
  list(
    s = c(1,1,1),
    b = c(0,0,0),
    sigma=0.01
  )
}

fit_mc = function(model, data, warmup = 10000){
  iter_sampling = 2000
  fit2 = model$sample(data = data, init=init_fun, iter_warmup = warmup, iter_sampling = iter_sampling, refresh = 2000)
  res_mcmc = extract_info(fit = fit2, data = data, num_draws = iter_sampling)
  
  mle = model$optimize(data = data) 
  mle
  
  fit_vb = model$variational(data = data, iter = 40000, algorithm='fullrank') 
  fit_vb$summary() 
  fit_vb$time()
 
  return(list(fit_vb=fit_vb, mle=mle, res_mcmc=res_mcmc))
}

b = c(0.1,-0.2, 0.3)
S = diag(c(0.9,1.0, 1.1))
if (FALSE){
  set.seed(42)
  df_hmc = NULL
  res = list()
  for (n in seq(1,20,1)) {
    acc = get_3D_Data(b=b, N = n, S = S)
    data = list(N = nrow(acc), X = acc, D=3, s_mu=0, s_sigma=0.5)
    res[[n]] = fit_mc(model.odr, data)
    save(list = ls(), file='R/larger_simulation/run_ODR_13Jan_1_20.rda')
  }
  print("----- FINISHED SIMULATIONS ---------")
} else{
  load('R/larger_simulation/run_ODR_13Jan_1_20.rda')
}

#Comparison ODR vs Full 3D
if (FALSE){
  load('R/larger_simulation/run_3D_16Jan.rda')
  res3d = res
  #Some 3D examples
  d = res3d[[4]]
  df_3d = data.frame(
    var = d$res_mcmc$summary$variable,
    q5 = d$res_mcmc$summary$q5,
    q50 = d$res_mcmc$summary$median,
    q95 = d$res_mcmc$summary$q95,
    type = '3D'
    )
  
  load('R/larger_simulation/run_ODR_13Jan_1_20.rda')
  dd = res[[5]]
  df_ODR = data.frame(
    var = dd$res_mcmc$summary$variable,
    q5 = dd$res_mcmc$summary$q5,
    q50 = dd$res_mcmc$summary$median,
    q95 = dd$res_mcmc$summary$q95,
    type = 'ODR'
  )
  
  df_3d[1:10,]
  df_ODR
  
}
    



df = NULL
for (n in 1:30) {
  if (is.null(res[[n]])){
    next
  }
  d = res[[n]]
  dd = as.data.frame(t(as.data.frame(d$mle$mle())))
  dd$n = n
  dd$method = 'MLE'
  df = bind_rows(df,dd)
  
  ####### Processing the VI Run ########
  #Checking if VI is valid
  vi_output = d$fit_vb$runset$procs$proc_output(1)
  vi_valid = TRUE
  for (i in 1:length(vi_output)){
    line = vi_output[i]
    if(grepl(".+not guaranteed to be meaningful.+",line)){
      vi_valid = FALSE
    }
  }
  
  if (vi_valid){
  dd = d$fit_vb$summary() %>% select(variable, median) %>%  pivot_wider(names_from = 1, values_from = 2)
  dd$n = n
  dd$method = 'VI_med'
  ddd = d$fit_vb$output()
  df = bind_rows(df,dd)
  
  dd = d$fit_vb$summary() %>% select(variable, q5) %>%  pivot_wider(names_from = 1, values_from = 2)
  dd$n = n
  dd$method = 'VI_q5'
  df = bind_rows(df,dd)
  
  dd = d$fit_vb$summary() %>% select(variable, q95) %>%  pivot_wider(names_from = 1, values_from = 2)
  dd$n = n
  dd$method = 'VI_q95'
  df = bind_rows(df,dd)
  } else {
    print(paste0("Invalid MF in ", n))
  }
  
  ##### Checking MCMC
  if (check_convergence(d$res_mcmc)){
  dd = d$res_mcmc$summary %>% select(variable, median) %>%  pivot_wider(names_from = 1, values_from = 2)
  dd$n = n
  dd$method = 'MCMC_med'
  df = bind_rows(df,dd)
  
  dd = d$res_mcmc$summary %>% select(variable, q5) %>%  pivot_wider(names_from = 1, values_from = 2)
  dd$n = n
  dd$method = 'MCMC_q5'
  df = bind_rows(df,dd)
  
  dd = d$res_mcmc$summary %>% select(variable, q95) %>%  pivot_wider(names_from = 1, values_from = 2)
  dd$n = n
  dd$method = 'MCMC_q95'
  df = bind_rows(df,dd)
  } else{
    print(paste0("No MCMC convergence in ", n))
  }
}

df2 = df %>% select(n, method, starts_with('b'), starts_with('s')) %>% 
  separate(col=method, into = c('quant', 'q'), sep="_") 

dd = df2 %>% filter(quant == 'MCMC') %>% 
  #select(-starts_with('s')) %>% 
  pivot_longer(cols=-(1:3)) %>% #puts b,s,sigma row by row
  pivot_wider(names_from = q, values_from = value)

### Data Frame with true values
 

dft = data.frame(
  names = c("b[1]", "b[2]","b[3]","s[1]","s[2]","s[3]"),
  vals = c(b, 1/S[1,1],1/S[2,2], 1/S[3,3])
)

#filter(!endsWith(method, 'MEAN')) %>% 
g1 = dd %>% filter(name !='sigma') %>% 
  ggplot(aes(x=n^2)) +
    geom_pointrange(aes(y=med, ymin=q5, ymax=q95, col=name)) +
  #geom_point(size=2) +
    geom_hline(data = dft, aes(yintercept=vals,col=names))+
    coord_cartesian(ylim=c(-0.25,1.35)) +
    labs(x= 'N', y='Estimates') + 
    guides(col=guide_legend(title="Parameters")) +
    ggpubr::theme_pubr()

gi = g1 + coord_cartesian(ylim=c(0.28,0.32)) + 
  labs(y='Estimates (focus on b[3])') + 
  ggpubr::theme_pubr(legend = 'none')
gi
if (FALSE){
  ggsave(filename = 'R/larger_simulation/simulation.pdf', ggarrange(g1,gi, ncol = 1, heights = c(1.5,1)), width = 6, height = 8)
  ggsave(filename = '~/Dropbox/Apps/Overleaf/bayes_calib/figures/simulation.pdf', ggarrange(g1,gi, ncol = 1, heights = c(1.5,1)), width = 6, height = 8)
}

 dd %>% filter(name=='sigma') %>% 
ggplot() + geom_point(aes(x=n, y=med, col=name))  +
  geom_pointrange(aes(x=n, y=med, ymin=q5, ymax=q95, col=name, shape=name)) +
  geom_hline(data = dft, aes(yintercept=vals,col=names)) +
  coord_cartesian(ylim=c(0,0.05)) 
 
dfti = NULL
for (n in 1:30) {
 if (!is.null(res[[n]])){
  print(res[[n]]$res_mcmc$time$total)
 }
}

res[[20]]$fit_vb$code()
# [1] "//Fitting parameters to an ellipsoid (principal axis)"         
# [2] "//We model hat_X = S(X - b)"                                   
# [3] "//Note that S is S^-1 in the direct model"                     
# [4] "//X are the observations"                                      
# [5] "data {"                                                        
# [6] "  int<lower=0> N;"                                             
# [7] "  int<lower=0> D; //Dimensionality (tested for 2 and 3D)"      
# [8] "  matrix[N, D] X; //Data Matrix (assumed centered around zero)"
# [9] "}"                                                             
# [10] ""                                                              
# [11] "parameters {"                                                  
# [12] "  row_vector[D] b; // The centers of the circle "              
# [13] "  real<lower=0> sigma;// The noise in xy-Direction"            
# [14] "  row_vector<lower=0>[D] s; //Different Scales in x, y, and z" 
# [15] "}"                                                             
# [16] ""                                                              
# [17] "model {"                                                       
# [18] "  real mu;"                                                    
# [19] "  sigma ~ normal(0,0.2);"                                      
# [20] "  s ~ lognormal(0, 0.5);"                                      
# [21] "  b ~ normal(0,1);"                                            
# [22] "  for (i in 1:N){"                                             
# [23] "    mu = 0;"                                                   
# [24] "    for (j in 1:D){"                                           
# [25] "      mu += (s[j]*(X[i,j] - b[j]))^2;"                         
# [26] "    }"                                                         
# [27] "    target += normal_lpdf(1 | sqrt(mu), sigma);"               
# [28] "  }"                                                           
# [29] "}"                                                             
# [30] ""                                                              
# [31] ""                       