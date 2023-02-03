data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  vector[N] z;
}
parameters {
  row_vector[3] b;
  row_vector<lower=0>[3] s;
  real<lower=0> sigma;
  vector<lower=0, upper = pi()>[N-1] theta;
  vector<lower=0, upper = 2*pi()>[N-1] phi;
  //Maybe replace theta sometime
  //https://mc-stan.org/docs/2_28/stan-users-guide/unit-vectors-and-rotations.html
}
transformed parameters {
  matrix[N, 3] g;
  //Initial conditions theta_0 = phi_0 = 0
  g[1,1] = b[1] + s[1] * 0;
  g[1,2] = b[2] + s[2] * 0; 
  g[1,3] = b[3] + s[3] * 1; 
  for(i in 2:N){
    g[i,1] = b[1] + s[1] * sin(theta[i-1])*cos(phi[i-1]);
    g[i,2] = b[2] + s[2] * sin(theta[i-1])*sin(phi[i-1]);
    g[i,3] = b[3] + s[3] * cos(theta[i-1]);
  } 
}
model {
  b ~ normal(0,1);
  sigma ~ normal(0,0.2);
  s ~ lognormal(0, 0.5);
  x ~ normal(g[,1],sigma);
  y ~ normal(g[,2],sigma);
  z ~ normal(g[,3],sigma);
}
