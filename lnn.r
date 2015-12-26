# returns an ExpressionSet, containing information about classification in featureData
# 'E_set' is an ExpressionSet
# 't_pi_prior' is the initial value for 't_pi'
# 'k_prior' is the intial value for 'k'
# 'b' is the concentration parameter of Dirichlet distribution
# 'plot' is a boolean indicating if plot frequency and approximated gamma distribution
# 'is_sim' indicates if lnn is called in a simulation
# 'verbose' indicates whether print intermediate result

# G is the number of genes,
# n is the number of test samples for every genes

lnn <- function(E_Set, t_pi_prior = c(0.33, 0.33, 0.34), k_prior = 1, b = c(2,2,2), plot = 0, is_sim = 0, verbose = 0)
{
  G = nrow(assayData(E_Set)[['exprs']])
  n = ncol(assayData(E_Set)[['exprs']])

  data = matrix(assayData(E_Set)[['exprs']], nrow = G, ncol = n)

  infinity = 1e100

  # 'sum_dgl_by_l' is an G * 1 matrix, the summation result of every row of 'data'
  sum_dgl_by_l = apply(data,1,sum)

  # 'sum_dgl_square_by_l' is an G * 1 matrix, the summation of every squared elements of 'data' by row
  sum_dgl_square_by_l = apply(data^2,1,sum)

  # lf123 returns log likelihood log f_1 log f_2 log f_3,
  # specific return value depends on the last parameter 'category'
  # 'psi' is the parameters, 
  # 'psi' = (\delta_1, k_1, \lambda_1, \nu_1,
  #          \delta_2, k_2, \lambda_2, \nu_2,
  #          0,        k_3, \lambda_3, \nu_3)
  # it is worthy to note that every 4 of them form a parameter group for each cluster
  # 'G_index' is a vector, it should be c(1, 2, 3, ..., G)
  lf123 <- function(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, category)
  {
    mu_0 = exp(psi[category*4-3])
    k = psi[category*4-2]
    alpha = exp(psi[category*4-1])
    beta = exp(psi[category*4])

    A = n/(2*(n*k+1))
    B_g = sum_dgl_by_l[G_index]/n
    C_g = beta + sum_dgl_square_by_l[G_index]/2 - (sum_dgl_by_l[G_index])^2/(2*n)
    D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) - n/2 * log(2*pi) - log(n*k+1)/2

    if (category == 1)
    {
      return (D - (n/2+alpha) * log (C_g + A * (mu_0 - B_g)^2))
    } 
    else if (category == 2)
    {
      return (D - (n/2+alpha) * log (C_g + A * (mu_0 + B_g)^2))
    }
    else
    {
      return (D - (n/2+alpha) * log (C_g + A * (B_g)^2))
    }
  }

  # it returns the expectation of \z_{g1}, \z_{g2}, \z_{g3}
  # this function is used in E-step
  # 'psi' is the parameters, 
  # 'psi' = (\delta_1, k_1, \lambda_1, \nu_1,
  #          \delta_2, k_2, \lambda_2, \nu_2,
  #          0,        k_3, \lambda_3, \nu_3)
  # 't_pi' = (\pi_1, \pi_2, \pi_3)
  # 'G_index' is a vector, it should be c(1, 2, 3, ..., G)
  get_tilde_z <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n)
  {
    t1 = t_pi[1] * exp(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 1))
    t2 = t_pi[2] * exp(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 2))
    t3 = t_pi[3] * exp(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 3))
    total = t1 + t2 + t3

    return (matrix(c(t1/total, t2/total, t3/total), nrow = length(G_index), ncol = 3))
  }

  # returns log likelihood that is in concern
  # this is the objective function to be maximized
  l_c <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    if (t_pi[1]<1e-6 || t_pi[2]<1e-6 || t_pi[3]<1e-6) return (-infinity)
    result = sum(tilde_z[,1] * lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 1)) + sum(tilde_z[,2] * lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 2)) + sum(tilde_z[,3] * lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 3)) + sum(tilde_z[,1] * log(t_pi[1])) + sum(tilde_z[,2] * log(t_pi[2])) + sum(tilde_z[,3] * log(t_pi[3]))

    #with dirichlet distribution
    result = result + lgamma(b[1]+b[2]+b[3]) - lgamma(b[1]) - lgamma(b[2]) - lgamma(b[3])
    result = result + sum((b-1) * log(t_pi))

    return (result)
  }

  # returns negative log likelihood in question
  # maximizing over l_c is equivalent to minimizing over negative_l_c
  negative_l_c <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    return (-l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z))
  }

  # returns gradient w.r.t. 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)
  gradient_l_c <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    G = length(G_index)
    gradient_temp = matrix(rep(0,12),4,3)
    for (category in 1:3)
    {
      mu_0 = exp(psi[category*4-3])
      k = psi[category*4-2]
      
      lambda = psi[category*4-1]
      nu = psi[category*4]

      alpha = exp(psi[category*4-1])
      beta = exp(psi[category*4])


      A = n/(2*(n*k+1))
      B_g = sum_dgl_by_l[G_index]/n
      C_g = beta + sum_dgl_square_by_l[G_index]/2 - (sum_dgl_by_l[G_index])^2/(2*n)
      D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) -n/2 * log(2*pi) - log(n*k+1)/2
      
      if (category == 1)
      {
        mu_0 = mu_0
      }
      else if (category == 2)
      {
        mu_0 = - mu_0
      }
      else
      {
        mu_0 = 0
      }

      d_delta = - sum(tilde_z[,category] * (2*A*(mu_0 - B_g))/(A*(mu_0 - B_g)^2 + C_g) * (n/2 + alpha) * mu_0)

      if (category == 3)
      {
        d_delta = 0
      }

      d_k = - (n * sum(tilde_z[,category])) / (2*(n*k+1)) + (n^2)/(2*(n*k+1)^2) * (n/2+alpha) *
              sum(tilde_z[,category]*(mu_0-B_g)^2/(A*(mu_0-B_g)^2+C_g))
      
      d_lambda = sum(tilde_z[,category] * alpha * (nu + digamma(n/2+alpha) - digamma(alpha))) - alpha * ( sum(tilde_z[,category] * log(A * (mu_0 - B_g)^2 + C_g)) )
      
      d_nu = sum(tilde_z[,category]) * alpha - beta * (n/2 + alpha) *  (sum(tilde_z[,category] / (A * (mu_0 - B_g)^2 + C_g)))
      
      gradient_temp[,category] = c(d_delta, d_k, d_lambda, d_nu)
    }

    return(c(gradient_temp))
  }

  # returns gradient of nagative log likelihood function
  gradient_negative_l_c <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    return (-gradient_l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z))
  }

  get_t_pi <- function(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  {
    #without Dirichlet
    # t1 = sum(tilde_z[,1]) / length(G_index)
    # t2 = sum(tilde_z[,2]) / length(G_index)
    # t3 = sum(tilde_z[,3]) / length(G_index)

    #with Dirichlet distribution
    denominator = length(G_index) + sum(b) - 3

    # print("===========================")
    # print(sum(tilde_z[,1]))
    # print(sum(tilde_z[,2]))
    # print(sum(tilde_z[,3]))
    
    t1 = (sum(tilde_z[,1]) + b[1] - 1) / denominator
    t2 = (sum(tilde_z[,2]) + b[2] - 1) / denominator
    t3 = (sum(tilde_z[,3]) + b[3] - 1) / denominator

    return (c(t1,t2,t3))
  }

  # limitations of parameters
  delta_min = -6
  delta_max = 15

  k_min = 0.001
  k_max = 1/k_min

  lambda_min = -6
  lambda_max = 6

  nu_min = -6
  nu_max = 6

  # initial values of 'psi' = (\delta_1, \delta_2, k, \lambda, \nu)

  median_dgl_by_l = apply(data, 1, median)

  temp = median(sort(median_dgl_by_l)[(ceiling(G * 0.95)):G])
  if (temp>0) delta_1 = log(temp)
  else delta_1 = delta_1_min

  temp = median(sort(median_dgl_by_l)[1:(trunc(G * 0.05))])
  if (temp<0) delta_2 = log(-temp)
  else delta_2 = delta_2_min

  k_1 = k_prior
  k_2 = k_prior
  k_3 = k_prior

  temp_tau = 1 / (apply(data, 1, mad)^2)
  lambda_1 = log((median(temp_tau)^2)/(mad(temp_tau)^2))
  lambda_2 = lambda_1
  lambda_3 = lambda_1
  
  nu_1 = log(median(temp_tau)/(mad(temp_tau)^2))
  nu_2 = nu_1
  nu_3 = nu_1

  if (plot)
  {
    hist(log10(temp_tau))
  }

  psi = c(delta_1, k_1, lambda_1, nu_1,
          delta_2, k_2, lambda_2, nu_2,
          0,       k_3, lambda_3, nu_3)

  t_pi = t_pi_prior
  G_index = c(1:G)

  tilde_z = get_tilde_z(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n)
  t_pi = get_t_pi(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)

  ###############################################################################
  ####-----------For Test Use Only-------------------------------------------####
  ###############################################################################

  # focus = 5
  # h = 0.0001

  # psi[focus] = psi[focus] - h
  # f1 = l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  # psi[focus] = psi[focus] + h
  # f2 = l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  # psi[focus] = psi[focus] + h
  # f3 = l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)
  # psi[focus] = psi[focus] - h

  # print(c("f1,f2,f3",f1,f2,f3))

  # derivative_approx = (f3 - f1)/(2*h)
  # derivative_accurate = gradient_l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)[focus]

  # print(c("difference:", derivative_accurate - derivative_approx))

  # return(0)

  ###############################################################################
  ####----END----For Test Use Only----END------------------------------------####
  ###############################################################################

  mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c, 
    t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, sum_dgl_square_by_l = sum_dgl_square_by_l, 
    G_index = G_index, n = n, tilde_z = tilde_z, method = 'L-BFGS-B', 
    lower=c(delta_min, k_min, lambda_min, nu_min, delta_min, k_min, lambda_min, nu_min, delta_min, k_min, lambda_min, nu_min), 
    upper=c(delta_max, k_max, lambda_max, nu_max, delta_max, k_max, lambda_max, nu_max, delta_max, k_max, lambda_max, nu_max),
    control = list(maxit = 10000))

  repeated_times = 0
  repeated_times_max = 1000

  while (repeated_times<repeated_times_max)
  {
    repeated_times = repeated_times + 1
    if (verbose)
    {
      print(c("repeated times:", repeated_times))
    }

    psi = mleinfo$par
    
    tilde_z = get_tilde_z(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n)

    t_pi = get_t_pi(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z)

    if (verbose)
    {
      print(psi)
      print(t_pi)
      # print(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 1)[17366])
      # print(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 2)[17366])
      # print(lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, 3)[17366])
      # print(tilde_z[17366,])
    #   print(negative_l_c(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n, tilde_z))
    }

    last_mleinfo = mleinfo
    mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c, 
      t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, sum_dgl_square_by_l = sum_dgl_square_by_l, 
      G_index = G_index, n = n, tilde_z = tilde_z, method = 'L-BFGS-B', 
      lower=c(delta_min, k_min, lambda_min, nu_min, delta_min, k_min, lambda_min, nu_min, delta_min, k_min, lambda_min, nu_min), 
      upper=c(delta_max, k_max, lambda_max, nu_max, delta_max, k_max, lambda_max, nu_max, delta_max, k_max, lambda_max, nu_max),
      control = list(maxit = 10000))

    if (abs(last_mleinfo$value - mleinfo$value)<1e-6) break
  }

  # if (plot)
  # {
  #   x = seq(0.01, max(temp_tau), length=1000)
  #   lines(x, dgamma(x,exp(mleinfo$par[3]),exp(mleinfo$par[4]))*G, col = 'green')
  # }

  tilde_z = get_tilde_z(psi, t_pi, sum_dgl_by_l, sum_dgl_square_by_l, G_index, n)

  # 3 means non differential expression
  # 1 means over expression
  # 2 means under expression

  fData(E_Set)[,2] = apply(tilde_z,1,which.max)
  if (is_sim)
  {
    fData(E_Set)[,3] = (fData(E_Set)[,2] == fData(E_Set)[,1])
  }
  # write.csv(tilde_z, file = "tilde_z.csv")

  return (list(E_Set, mleinfo, t_pi))

}
