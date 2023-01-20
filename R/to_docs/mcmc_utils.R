get_2D_Data <- function(b = c(0,0), N = 50, S = diag(c(1,1)), sigma=1/50., random=FALSE) {
  phi = seq(0,2*pi, length.out=N+1)[1:N]
  if (random){
    phi = runif(N, 0,2*pi)[1:N]
  } 
  x = matrix(c(cos(phi), sin(phi)), nrow=2, byrow = TRUE)
  x_obs = S %*% x + b
  acc = t(x_obs)
  return(acc * matrix(rnorm(2*N,1,sigma), nrow=N))
}

get_3D_Data <- function(b = c(0,0,0), 
                        N, S = diag(c(1,1,1)), 
                        sigma=1/50., 
                        phi_max = 2*pi,
                        theta_max = pi,
                        random=FALSE) {
  phis = seq(0,phi_max, length.out=N+1)[1:N]
  thetas = seq(0.1,theta_max-0.1, length.out=N+1)[1:N]
  if (random){
    phis = runif(N, 0,2*pi)[1:N]
    thetas = runif(N, 0,pi)[1:N]
  }
  ret = matrix(NA, nrow = N*N, ncol=3)
  i = 0
  for (phi in phis){
    for (theta in thetas){
      i = i + 1
      A = matrix(c(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)), nrow=3, byrow = TRUE)
      acc = S %*% A + b
      acc = t(acc) * rnorm(3,1,sigma)
      ret[i,] = acc
    }
  }
  return(ret)
}

calibrate <- function(fit, acc) {
  b_off = spread_draws(fit, b[i]) %>% 
    mutate(i = factor(i)) %>% 
    group_by(i) %>% summarise(median(b))
  b_off = b_off$`median(b)`
  s_off = spread_draws(fit, s[i]) %>% 
    mutate(i = factor(i)) %>% 
    group_by(i) %>% summarise(median(s))
  s_off = s_off$`median(s)`
  cal_acc = data.frame(
    acc_cal_1 = s_off[1]*(acc[,1] - b_off[1]),
    acc_cal_2 = s_off[2]*(acc[,2] - b_off[2]),
    acc_cal_3 = s_off[3]*(acc[,3] - b_off[3])
  )
  return(cal_acc)
}

get_calibration_coef <- function(fit) {
  b_off = spread_draws(fit, b[i]) %>% 
    mutate(i = factor(i)) %>% 
    group_by(i) %>% summarise(median(b))
  b_off = b_off$`median(b)`
  s_off = spread_draws(fit, s[i]) %>% 
    mutate(i = factor(i)) %>% 
    group_by(i) %>% summarise(median(s))
  s_off = s_off$`median(s)`
  cal_acc = data.frame(
    b = b_off,
    s_off = s_off
  )
  return(cal_acc)
}

plot_Calibration <- function(acc, radius = 1) {
  if (TRUE){
    # Define the center and radius of the sphere
    center <- c(0, 0, 0)
    # Create a sequence of angles for the sphere
    theta <- seq(0, 2*pi, length.out = 10)
    phi <- seq(0, pi, length.out = 10)
    # Create matrices of x, y, and z coordinates for the sphere
    x <- radius * outer(cos(theta), sin(phi)) + center[1]
    y <- radius * outer(sin(theta), sin(phi)) + center[2]
    z <- radius * outer(rep(1, length(theta)), cos(phi)) + center[3]
    # Plot the sphere using the surface3d function
    surface3d(x, y, z, col = "gray", alpha = 0.1)
  }
  rgl::surface3d(x, y, z, color = "gray")
  points3d(acc[,1], acc[,2], acc[,3], col = "red", size = 5)
  rgl::text3d(1.2,0,0,"X")
  rgl::text3d(0,1.2,0,"Y")
  rgl::text3d(0,0,1.2,"Z")
  #arrow3d(c(0, 0, 0),c( 0, 0, 1), col = "red")
  # Add y-axis arrow
  #segments3d(0, 0, 0, 0, 1, 0)
  # Add z-axis arrow
  #segments3d(0, 0, 0, 0, 0, 1)
  rgl.postscript("sphere.pdf", fmt='pdf')
}

#Faster to provide num_draws
extract_info = function(fit, data, vars=NULL, num_draws=-1){
  if (is.null(vars) == FALSE){
    summary = fit$summary(vars)
  } else{
    summary = fit$summary()
  }
  res = list(
    vars = vars,
    diag = fit$diagnostic_summary(),
    summary = summary,
    code = fit$code(),
    data = data,
    time = fit$time()
  )
  if (num_draws < 0){
    num_draws = dim(fit$draws())[1]
  } else{
    num_draws = num_draws
  }
  res = c(res, num_draws=num_draws)
  return (res)
}

check_convergence = function(info){
  nd = max(info[['diag']]$num_divergent)
  rhat = max(abs(info[['summary']]$rhat))
  ess_bulk = max(abs(info[['summary']]$ess_bulk))
  ess_tail = max(abs(info[['summary']]$ess_tail))
  num = info[['num_draws']]
  print(paste(nd/num, rhat, ess_bulk/num, sep=' '))
  return (nd/num < 0.2 & rhat<1.1 & ess_bulk/num > 0.5)
}

check_coverage = function(info, var, var_name) {
  df = info$summary 
  bi = 0
  comps = 0
  in_cof = 0
  for (i in 1:nrow(df)){
    if (startsWith(df[i,1]$variable, var_name)){
      bi = bi + 1
      if (df[i,'q5'] < var[bi] && var[bi]  < df[i,'q95']){
        in_cof = in_cof + 1
      } 
      comps = comps+1
    }
  }
  return (in_cof / comps)
}