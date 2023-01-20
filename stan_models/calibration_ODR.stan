//Fitting parameters to an ellipsoid (principal axis)
//We model hat_X = S(X - b)
//Note that S is S^-1 in the direct model
//X are the observations
data {
  int<lower=0> N;
  int<lower=0> D; //Dimensionality (tested for 2 and 3D)
  matrix[N, D] X; //Data Matrix (assumed centered around zero)
}

parameters {
  row_vector[D] b; // The centers of the circle 
  real<lower=0> sigma;// The noise in xy-Direction
  row_vector<lower=0>[D] s; //Different Scales in x, y, and z
}

model {
  real mu;
  sigma ~ normal(0,0.2);
  s ~ lognormal(0, 0.5);
  b ~ normal(0,1);
  for (i in 1:N){
    mu = 0;
    for (j in 1:D){
      mu += (s[j]*(X[i,j] - b[j]))^2;
    }
    target += normal_lpdf(1 | sqrt(mu), sigma);
  }
}


