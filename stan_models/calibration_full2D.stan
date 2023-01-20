data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  row_vector[2] b;
  row_vector<lower=0>[2] s;
  real<lower=0> sigma;// The noise in xy-Direction
  vector<lower=-pi(), upper = pi()>[N-1] theta;
  //Maybe replace theta sometime
  //https://mc-stan.org/docs/2_28/stan-users-guide/unit-vectors-and-rotations.html
}
transformed parameters {
  matrix[N, 2] A;
  A[1,1] = b[1] + s[1];
  A[1,2] = b[2]; 
  for(i in 2:N){
    A[i,1] = b[1] + s[1] * cos(theta[i-1]);
    A[i,2] = b[2] + s[2] * sin(theta[i-1]);
  } 
}
model {
  b ~ normal(0,1);
  sigma ~ normal(0,0.2);
  s ~ lognormal(0, 0.5);
  x ~ normal(A[,1],sigma);
  y ~ normal(A[,2],sigma);
}
