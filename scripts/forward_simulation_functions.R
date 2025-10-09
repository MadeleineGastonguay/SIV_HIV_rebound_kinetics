require(lamW)
require(deSolve)

# Model definition ####

# Linear scale
ctl_nonlytic_identifiable_virus<- function(t, state, parameters){
  
  with(as.list(c(state, parameters)), {
    
    
    d_Tar <- (t >= t_go + washout_time) * (L - Tar*B*V - Tar*dTar) #Uninfected cells
    d_I <- (t >= t_go + washout_time) * (a + B * Tar * V / (1 + E / N_E) - I*dI) #Infected cells
    d_V <- (t >= t_go + washout_time) * (k * I - c *V) #Virus
    d_P <- (t >= t_go + washout_time) * (m + p*(1-f)*P*V/(N_P+V) - P*dP)  #Precurors
    d_E <- (t >= t_go + washout_time) * (P*p*f*V/(N_P+V) - E*dE)#Effectors
    
    return(list(c(d_Tar, d_I, d_V, d_P, d_E)))
  })
}

ctl_nonlytic_identifiable_virus_delay <- function(t, state, parameters){
  
  with(as.list(c(state, parameters)), {
    tauQ <- 1
    if(t < tauQ){
      Tar_lag <- Tar
      V_lag <- k/c*I
    }else{
      Tar_lag <- lagvalue(t - tauQ, 1)
      I_lag <- lagvalue(t - tauQ, 2)  
      V_lag <- k/c*I_lag
    }
    
    # Assuming V is proportional to I
    V <- k/c*I
    
    d_Tar <- (t >= t_go + washout_time) * (L - Tar*B*V - Tar*dTar) #Uninfected cells
    d_I <- (t >= t_go + washout_time) * (a + B * Tar_lag * V_lag * exp(-dTar * tauQ) / (1 + E / N_E) - I*dI) #Infected cells
    # d_V <- (t >= t_go + washout_time) * (k * I - c *V) #Virus
    d_P <- (t >= t_go + washout_time) * (m + p*(1-f)*P*V/(N_P+V) - P*dP)  #Precurors
    d_E <- (t >= t_go + washout_time) * (P*p*f*V/(N_P+V) - E*dE)#Effectors
    
    # return(list(c(d_Tar, d_I, d_V, d_P, d_E)))
    return(list(c(d_Tar, d_I, d_P, d_E)))
  })
}

ctl_nonlytic_identifiable_virus_eclipse<- function(t, state, parameters){
  
  with(as.list(c(state, parameters)), {
    
    
    d_Tar <- (t >= t_go + washout_time) * (L - Tar*B*V - Tar*dTar) #Uninfected cells
    d_I <- (t >= t_go + washout_time) * (a + q*Q - I*dI) #Infected cells
    d_V <- (t >= t_go + washout_time) * (k * I - c *V) #Virus
    d_P <- (t >= t_go + washout_time) * (m + p*(1-f)*P*V/(N_P+V) - P*dP)  #Precurors
    d_E <- (t >= t_go + washout_time) * (P*p*f*V/(N_P+V) - E*dE)#Effectors
    d_Q <- (t >= t_go + washout_time) * (B * Tar * V / (1 + E / N_E) - Q*dTar - Q*q)
    
    return(list(c(d_Tar, d_I, d_V, d_P, d_E, d_Q)))
  })
}


# Natural log scale
ctl_nonlytic_identifiable_virus_ls <- function(t, state, parameters){
  
  with(as.list(c(state, parameters)), {
    
    d_Tar <- (t >= t_go + washout_time) * (L*exp(-1*Tar) - B*exp(V) - dTar) #Uninfected cells
    d_I <- (t >= t_go + washout_time) * (exp(a_log-I) + B * exp(Tar + V - I) / (1 + (exp(E) / N_E)) - dI) #Infected cells
    d_V <- (t >= t_go + washout_time) * (k * exp(I - V) - c) #1/log(10)*( a*kdv*10^(-1*v) + B*kdv*10^x*(Nb/(Nb+10^z)) - dy ) #Virus
    d_P <- (t >= t_go + washout_time) * (m * exp(-1*P) + p*(1-f)*exp(V)/(N_P+exp(V)) - dP)  #Precurors
    d_E <- (t >= t_go + washout_time) * (exp(P-E)*p*f*exp(V)/(N_P+exp(V)) - dE)#Effectors
    
    return(list(c(d_Tar, d_I, d_V, d_P, d_E)))
  })
}

ctl_nonlytic_identifiable_virus_ls_delay <- function(t, state, parameters){
  
  with(as.list(c(state, parameters)), {
    tauQ <- 1
    if(t < tauQ){
      Tar_lag <- Tar
      V_lag <- log(k/c) + I
    }else{
      Tar_lag <- lagvalue(t - tauQ, 1)
      I_lag <- lagvalue(t - tauQ, 2)
      V_lag <- log(k/c) + I_lag
    }
    
    # Assuming V is in equilibrium with I 
    V <- log(k/c) + I
    
    d_Tar <- (t >= t_go + washout_time) * (L*exp(-1*Tar) - B*exp(V) - dTar) #Uninfected cells
    d_I <- (t >= t_go + washout_time) * (exp(a_log-I) + B * exp(Tar_lag + V_lag - dTar*tauQ - I) / (1 + (exp(E) / N_E)) - dI) #Infected cells
    # d_V <- (t >= t_go + washout_time) * (k * exp(I - V) - c) #1/log(10)*( a*kdv*10^(-1*v) + B*kdv*10^x*(Nb/(Nb+10^z)) - dy ) #Virus
    d_P <- (t >= t_go + washout_time) * (m * exp(-1*P) + p*(1-f)*exp(V)/(N_P+exp(V)) - dP)  #Precurors
    d_E <- (t >= t_go + washout_time) * (exp(P-E)*p*f*exp(V)/(N_P+exp(V)) - dE)#Effectors
    
    # return(list(c(d_Tar, d_I, d_V, d_P, d_E)))
    return(list(c(d_Tar, d_I, d_P, d_E)))
  })
}

ctl_nonlytic_identifiable_virus_ls_eclipse <- function(t, state, parameters){
  
  with(as.list(c(state, parameters)), {
    
    d_Tar <- (t >= t_go + washout_time) * (L*exp(-1*Tar) - B*exp(V) - dTar) #Uninfected cells
    d_I <- (t >= t_go + washout_time) * (exp(a_log-I) + q*exp(Q - I) - dI) #Infected cells
    d_V <- (t >= t_go + washout_time) * (k * exp(I - V) - c) #1/log(10)*( a*kdv*10^(-1*v) + B*kdv*10^x*(Nb/(Nb+10^z)) - dy ) #Virus
    d_P <- (t >= t_go + washout_time) * (m * exp(-1*P) + p*(1-f)*exp(V)/(N_P+exp(V)) - dP)  #Precurors
    d_E <- (t >= t_go + washout_time) * (exp(P-E)*p*f*exp(V)/(N_P+exp(V)) - dE)#Effectors
    d_Q <- (t >= t_go + washout_time) * (B * exp(Tar + V - Q) / (1 + (exp(E) / N_E)) - dTar - q)
    
    return(list(c(d_Tar, d_I, d_V, d_P, d_E, d_Q)))
  })
}


# Simulation function ####
sim_traj_final <- function(parameters, return_full, delay = F){
  
  
  log_rebound_threshold <- log(50)
  
  max_immune_level_P <- 1000 * N_E
  tauQ <- 1 # delay for delay model
  
  L <- parameters[1]
  B <- parameters[2]
  t_a <- parameters[3]
  m <- parameters[4]
  p <- parameters[5]
  N_P <- parameters[6]
  
  
  #ART time-averaged parameters
  T_ART <- L/dTar
  I_ART <- I1 / (t_a * dI)
  V_ART <- I1 / (t_a * dI) * k / c
  P_ART <- m / (dP - p * (1 - f) * V_ART / (V_ART + N_P))
  if(P_ART < 0){
    P_ART <- max_immune_level_P
    print("USING MAX IMMUNE LEVEL FOR P")
  }
  # if(P_ART > max_immune_level_P){
  if(dP<p*(1-f)*V_ART/(V_ART+N_P)){
    P_ART <- max_immune_level_P
    print("USING MAX IMMUNE LEVEL FOR P")
  }
  
  E_ART <- p*f*P_ART * V_ART / (dE * (V_ART + N_P))
  
  
  #functions
  R_0 <- L * B * k / ((1 + E_ART/N_E)*dTar*dI*c)
  if(delay) R_0 <- R_0*exp(-dTar*tauQ)
  if(!delay){
    r <- max(0, (-(c + dI) + sqrt((c - dI)^2 + 4 * c * dI * R_0))/2)
  }else{
    r <- 1/tauQ*lambertW0(R_0*tauQ*dI*c/(dI + c)*exp(tauQ*dI*c/(dI + c))) - dI*c/(dI + c)
  }
  if(R_0 > 1){
    p_surv <- (lambertW0(-R_0 * exp(-R_0)) + R_0) / (R_0 + rho - 1) 
    if(delay) p_surv <- p_surv*exp(dTar * tauQ)
    print(c(R_0, r))
  } else {
    p_surv <- 0
  }
  
  
  if(R_0 <= 1){
    a <- I1/t_a
  } else {
    a <- I1 * r / (p_surv * (exp(r * t_a / p_surv) - 1))
  }
  n_weight <- 3
  w_weight <- 2 ^(-1*(t_a * dI / p_surv)^n_weight)
  
  if(R_0 > 1){
    t_go <- max(0, t_a / p_surv - 1/r*log(w_weight * p_surv / (t_a * dI) * (exp(r*t_a/p_surv) - 1) + 1))
  } else {
    t_go <- 0
  }
  
  
  if(is.na(t_go)){
    # t_go <- 0
    print("t_go NAN!")
    print(c("pars", parameters))
    print(c("p surv", p_surv))
    print(c("inner log", w_weight * p_surv / (t_a * dI) * (exp(r*t_a/p_surv) - 1)))
    print(c("R0", R_0))
    if(!return_full){
      return(data.frame(peak = NA, final = NA, t_reb=NA, p_surv=p_surv, reactivation_time = t_a / p_surv))
    }else{
      return(list(data.frame(r = r, V_ART = V_ART, peak = NA, final = NA, t_reb=NA, p_surv=p_surv, reactivation_time = t_a / p_surv), 
                  data.frame(V = NA, Tar = NA,  I = NA,  P = NA, E = NA, post_activated = NA,  pre_activated = NA, time = NA)))
    }
  }
  if(t_go == -Inf){
    t_go <- 0
    print("t_go negative inf!")
    print(c("pars", parameters))
    print(c("p surv", p_surv))
    print(c("inner log", w_weight * p_surv / (t_a * dI) * (exp(r*t_a/p_surv) - 1)))
    print(c("R0", R_0))
  }
  
  
  Tar_0 <- T_ART
  if(R_0 > 1){
    I_0 <- max(w_weight*I_ART, I1/p_surv)
  } else {
    I_0 <- I_ART
  }
  V_0 <- k/c * I_0
  P_0 <- P_ART
  E_0 <- E_ART
  
  
  y_start <- c(Tar = Tar_0, I = I_0, V = V_0, P=P_0, E=E_0)
  
  #check for -Inf
  
  names(y_start) <- c("Tar", "I", "V", "P", "E") 
  
  pars_adjusted <- c(parameters, log(a), t_go)
  names(pars_adjusted) <- c("L", "B", "t_a", "m", "p", "N_P", "a_log", "t_go")
  
  
  if(sum(is.nan(y_start)) > 0){
    print("STARTING CONDITION NAN")
    print(y_start) 
    print(parameters)
  }
  
  if(!delay){
    times_seq <- seq(from=0, to=365, by=.1)
    out <- ode(log(y_start), times = times_seq, func = ctl_nonlytic_identifiable_virus_ls, parms = pars_adjusted,  atol=0)
  }else{
    times_seq <- seq(from=0, to=365, by=.1)
    y_start <- y_start[names(y_start) != "V"]
    
    out <- dede(log(y_start), times = times_seq, func = ctl_nonlytic_identifiable_virus_ls_delay, parms = pars_adjusted) 
    out <- cbind(out, V = log(k/c) + out[,"I"])
    out <- out[,c("time", "Tar", "I", "V", "P", "E")]
    
    # out <- dede(y_start, times = times_seq, func = ctl_nonlytic_identifiable_virus_delay, parms = pars_adjusted,  atol=0)
  }
  
  
  peak_VL <- max(out[,"V"])
  final_VL <- out[nrow(out),"V"]
  if(sum(is.nan(c(peak_VL, final_VL))) > 0){
    print("NaN values!")
    print(parameters)
    print(out)
  }
  if(sum(which(out[,"V"] > log_rebound_threshold)) > 0){
    t_reb <- times_seq[min(which(out[,"V"] > log_rebound_threshold))]
  } else {
    t_reb <- Inf
  }
  if(!return_full){
    return(data.frame(peak = peak_VL, final = final_VL, t_reb=t_reb, p_surv=p_surv, reactivation_time = t_a / p_surv))
  }
  else{
    #replace values before t_go with ART values
    print("V_art")
    print(V_ART)
    out[out[,"time"] < t_go+washout_time,"V"] <- log(V_ART)
    pre_activated <- rep("no", nrow(out))
    post_activated <- rep("no", nrow(out))
    pre_activated[out[,"time"] <= t_go+washout_time] <- "yes"
    post_activated[out[,"time"] >= t_go+washout_time] <- "yes"
    
    out_final <- data.frame(out, pre_activated = pre_activated, post_activated=post_activated)
    
    return(list(data.frame(r = r, V_ART = V_ART, peak = peak_VL, final = final_VL, t_reb=t_reb, p_surv=p_surv, reactivation_time = t_a / p_surv), out_final))
  }
  #return(out)
}


sim_traj_final_eclipse <- function(parameters, return_full){
  
  log_rebound_threshold <- log(50)
  
  max_immune_level_P <- 1000 * N_E
  
  L <- parameters[1]
  B <- parameters[2]
  t_a <- parameters[3]
  m <- parameters[4]
  p <- parameters[5]
  N_P <- parameters[6]
  
  q <- 1
  
  #ART time-averaged parameters
  T_ART <- L/dTar
  I_ART <- I1 / (t_a * dI)
  V_ART <- I1 / (t_a * dI) * k / c
  P_ART <- m / (dP - p * (1 - f) * V_ART / (V_ART + N_P))
  Q_ART <- 0
  if(P_ART < 0){
    P_ART <- max_immune_level_P
    print("USING MAX IMMUNE LEVEL FOR P")
  }
  if(P_ART > max_immune_level_P){
    P_ART <- max_immune_level_P
    print("USING MAX IMMUNE LEVEL FOR P")
  }
  
  E_ART <- p*f*P_ART * V_ART / (dE * (V_ART + N_P))
  
  
  #functions
  R_0 <- L * B * k / ((1 + E_ART/N_E)*dTar*dI*c) * (q/(q + dTar))
  # r <- max(0, (-(c + dI) + sqrt((c - dI)^2 + 4 * c * dI * R_0))/2)
  # Intermediates for calculating exponential growth rate
  a0 <- 1 - R_0
  a1 <- 1 / dI + 1 / c + 1 / (q + dTar)
  a2 <- (1 / dI) * (1 / c) + (1 / dI) * (1 / (q + dTar)) + (1 / c) * (1 / (q + dTar))
  a3 <- (1 / dI) * (1 / c) * (1 / (q + dTar))
  
  gamma <- -a2^3 / (27 * a3^3) + (a2 * a1) / (6 * a3^2) - a0 / (2 * a3)
  delta <- a1 / (3 * a3) - a2^2 / (9 * a3^2)
  
  sqrt_term <- sqrt(as.complex(gamma^2 + delta^3))
  r <- Re((gamma + sqrt_term)^(1/3) + (gamma - sqrt_term)^(1/3) - a2 / (3 * a3))
  
  if(R_0 > 1){
    p_surv <- (lambertW0(-R_0 * exp(-R_0)) + R_0) / (R_0 + rho - 1)*(1 + dTar/q)
    print(c(R_0, r))
  } else {
    p_surv <- 0
  }
  
  
  if(R_0 <= 1){
    a <- I1/t_a
  } else {
    a <- I1 * r / (p_surv * (exp(r * t_a / p_surv) - 1))
  }
  n_weight <- 3
  w_weight <- 2 ^(-1*(t_a * dI / p_surv)^n_weight)
  
  if(R_0 > 1){
    t_go <- max(0, t_a / p_surv - 1/r*log(w_weight * p_surv / (t_a * dI) * (exp(r*t_a/p_surv) - 1) + 1))
  } else {
    t_go <- 0
  }
  
  
  if(is.na(t_go)){
    t_go <- 0
    print("t_go NAN!")
    print(c("pars", parameters))
    print(c("p surv", p_surv))
    print(c("inner log", w_weight * p_surv / (t_a * dI) * (exp(r*t_a/p_surv) - 1)))
    print(c("R0", R_0))
  }
  if(t_go == -Inf){
    t_go <- 0
    print("t_go negative inf!")
    print(c("pars", parameters))
    print(c("p surv", p_surv))
    print(c("inner log", w_weight * p_surv / (t_a * dI) * (exp(r*t_a/p_surv) - 1)))
    print(c("R0", R_0))
  }
  
  
  
  Tar_0 <- T_ART
  if(R_0 > 1){
    I_0 <- max(w_weight*I_ART, I1/p_surv)
  } else {
    I_0 <- I_ART
  }
  V_0 <- k/c * I_0
  P_0 <- P_ART
  E_0 <- E_ART
  Q_0 <- 1e-3
  
  times_seq <- seq(from=0, to=365, by=.1)
  
  
  y_start <- c(Tar = Tar_0, I = I_0, V = V_0, P=P_0, E=E_0, Q = Q_0)
  
  #check for -Inf
  
  names(y_start) <- c("Tar", "I", "V", "P", "E", "Q") 
  print(y_start)
  
  pars_adjusted <- c(parameters, log(a), t_go)
  names(pars_adjusted) <- c("L", "B", "t_a", "m", "p", "N_P", "a_log", "t_go")
  
  if(sum(is.nan(y_start)) > 0){
    print("STARTING CONDITION NAN")
    print(y_start)
    print(parameters)
  }
  
  
  out <- ode(log(y_start), times = times_seq, func = ctl_nonlytic_identifiable_virus_ls_eclipse, parms = pars_adjusted,  atol=0)
  
  peak_VL <- max(out[,"V"])
  final_VL <- out[nrow(out),"V"]
  if(sum(is.nan(c(peak_VL, final_VL))) > 0){
    print("NaN values!")
    print(parameters)
    print(out)
  }
  if(sum(which(out[,"V"] > log_rebound_threshold)) > 0){
    t_reb <- times_seq[min(which(out[,"V"] > log_rebound_threshold))]
  } else {
    t_reb <- Inf
  }
  if(!return_full){
    return(data.frame(peak = peak_VL, final = final_VL, t_reb=t_reb, p_surv=p_surv, reactivation_time = t_a / p_surv, a= a, tgo = t_go))
  }
  else{
    #replace values before t_go with ART values
    print("V_art")
    print(V_ART)
    out[out[,"time"] < t_go,"V"] <- log(V_ART)
    pre_activated <- rep("no", nrow(out))
    post_activated <- rep("no", nrow(out))
    pre_activated[out[,"time"] <= t_go + washout_time] <- "yes"
    post_activated[out[,"time"] >= t_go + washout_time] <- "yes"
    
    out_final <- data.frame(out, pre_activated = pre_activated, post_activated=post_activated)
    return(list(data.frame(r = r, V_ART = V_ART, peak = peak_VL, final = final_VL, t_reb=t_reb, p_surv=p_surv, reactivation_time = t_a / p_surv, a= a, tgo = t_go), out_final))
  }
  #return(out)
}


# Functions to sample parameter value for simulations ####

# Mean_vec is the population parameters and varcov_mat is the FIM output from monolix
sample_fitting_results_scaled <- function(mean_vec, varcov_mat){
  sampled_model_pars <- mvrnorm(mu = mean_vec, Sigma=varcov_mat)
  sampled_model_pars['NPs_pop'] <- sampled_model_pars['NPs_pop'] + log( 10 ^5)
  sampled_model_pars['betas_pop'] <- sampled_model_pars['betas_pop'] + log( 10 ^(-7))
  return(sampled_model_pars)
}

make_covariate_matrix_list <- function(fit_results){
  #each study design is a vector with a 1 corresponding to treatment. By default whitney is 1 and barouch is 0 (vac and Ab)
  #each matrix is virus-specific
  
  #treatment order is TLR7, Vac, Study
  #Call the "study design" vector S (column)
  #To specify interaction effects, there must be a 1 in both positions. So the interaction matrix is S * S ^ T
  #The overall effect on a parameter is the sum over all elements of the element-wise product of S * S^ T with the 
  #treatment effect matrix (the diagonal elements are the bare treatment effects)
  #treatment effect matrix is upper-triangular by default
  
  param_names <- c('lambda', 'betas', 'ta', 'm', 'p', 'NPs')
  treatment_names_intx <- c('TLR7', 'VAC_1', 'STUD_1', 'AB_1')
  treatment_names_bare <- c('tlr7_1', 'vac_1', 'study_1', 'ab_1')
  effect_list <- list()
  for(parameter_idx in 1:length(param_names)){
    matrix_tmp <- matrix(0, nrow=4, ncol=4)
    rownames(matrix_tmp) <- treatment_names_intx
    colnames(matrix_tmp) <- treatment_names_intx
    for(row_idx in 1:4){
      for(col_idx in row_idx:4){
        if(row_idx == col_idx){
          #on the diagonal, just take bare effect
          effect_name <- paste(c('beta_', param_names[parameter_idx], '_', treatment_names_bare[col_idx]), sep="", collapse="")
          if(effect_name %in% names(fit_results)){
            matrix_tmp[row_idx, col_idx] <- fit_results[effect_name]
          }
        }
        else{
          #off diagonal, make intx name
          effect_name <- paste(c('beta_', param_names[parameter_idx], '_int', treatment_names_intx[row_idx], treatment_names_intx[col_idx]), sep="", collapse="")
          if(effect_name %in% names(fit_results)){
            matrix_tmp[row_idx, col_idx] <- fit_results[effect_name]
          }
        }
      }
    }
    effect_list[[param_names[parameter_idx]]] <- matrix_tmp
  }
  return(effect_list)
}

# Formatting population mean and random effects for parameters
get_mean_params <- function(fit_results){
  mean_param_names <- c('lambda_pop', 'betas_pop', 'ta_pop', 'm_pop', 'p_pop', 'NPs_pop')
  omega_param_names <- c('omega_lambda', 'omega_betas', 'omega_ta', 'omega_m', 'omega_p', 'omega_NPs')
  param_names <- c('lambda', 'betas', 'ta', 'm', 'p', 'NPs')
  return_vec_means <- 1:6 * 0
  return_vec_omega <- 1:6 * 0
  
  for(name_idx in 1:length(mean_param_names)){
    return_vec_means[name_idx] <- fit_results[paste(c(param_names[name_idx], '_pop'), sep="", collapse="")]
    omega_name <- paste(c('omega_', param_names[name_idx]), sep="", collapse="")
    if(omega_name %in% names(fit_results)){
      return_vec_omega[name_idx] <- fit_results[omega_name]
    }
  }
  df_return <- data.frame(mean = return_vec_means, omega = return_vec_omega, row.names=param_names)
  return(df_return)
}

sample_baseline_parameters_ls <- function(num_samples, mean_info_df){
  #means <- list(L=44.1, B = 1e-6 * 0.73, a = 1e-6 * 4.26, m = 7.39, p = 1.17, Np = 1e5 * 14.1) #no washout
  #means <- list(L=44.1, B = 1e-6 * 0.721, t_a = 1e-6 * 77.6, m = 7.39, p = 1.49, Np = 1e5 * 27.5) #3 day washout #FIXME
  #random_effects_logspace <- list(L=0.243, B=0.023, a=2.73, m=0, p=0.648, Np=1.67) #is this the variance or sd? #why does it say beta (fixed)?
  #random_effects_logspace <- list(L=0.222, B=0.16, a=1.97, m=0, p=0.868, Np=1.72) #washout 3 days #FIXME
  
  sampled_ls <- rep(NA, 6)
  for(i in 1:length(sampled_ls)){
    sampled_ls[i] <- rnorm(n=num_samples, mean=mean_info_df[i,'mean'], sd=mean_info_df[i, 'omega'])
  }
  return(sampled_ls)
}

add_covariate_effects <- function(parameter_start, covariates, treatment_effects){
  #covariates is a column vector of 0s and 1s specifying:
  #intTLR7STUD	intTLR7VAC	study	tlr7	vac
  
  #treatment_effects is a matrix, 6 X 5, such that treatment_effects * covariates gives the proportional adjustment of parameter_start
  parameter_end <- parameter_start
  covariate_dot_prod <- covariates %*% t(covariates)
  for(parameter_idx in 1:length(parameter_start)){
    effect_matrix <- treatment_effects[[parameter_idx]] 
    parameter_end[[parameter_idx]] <- parameter_start[[parameter_idx]] + sum(covariate_dot_prod * effect_matrix)
  }
  
  return(parameter_end)
}




# Function to format outputs from monolix ####

get_transformed_varcov_mat <- function(varcov_file, theta_estimates){
  #returns a variance covariance matrix transformed partially so that ODE parameters are in log-space but beta parameters
  #remain in regular space
  #REQUIRED that PopParam and VarCov have same order of parameters
  varcov <- as.matrix(read.table(varcov_file, sep=",", stringsAsFactors = FALSE, header=FALSE, row.names=1))
  colnames(varcov) <- rownames(varcov)
  varcov_transformed <- varcov
  for(colname in colnames(varcov)){
    for(rowname in rownames(varcov)){
      col_scale <- theta_estimates[colname, 'value']
      if(substr(colname, start=1, stop=5) == 'beta_'){
        col_scale <- Inf#1 #CHANGED TO INF TO SUPPRESS BETA VARIATION
      }
      row_scale <- theta_estimates[rowname, 'value']
      if(substr(rowname, start=1, stop=5) == 'beta_'){
        row_scale <- Inf#1 #CHANGED TO INF TO SUPPRESS BETA VARIATION
      }
      varcov_transformed[rowname, colname] <- varcov[rowname, colname] / (row_scale * col_scale)
    }
  }
  return(varcov_transformed)
}

# Function for formatting colorbars in plots ####
trim_tails <- function(range = c(-Inf, Inf)){
  trans_new("trim_tails", 
            transform = function(x) {
              force(range)
              desired_breaks <- extended_breaks(n = 5)(x[x >= range[1] & x <= range[2]])
              break_increment <- diff(desired_breaks)[1]
              x[x < range[1]] <- range[1] - break_increment
              x[x > range[2]] <- range[2] + break_increment
              x
            },
            inverse = function(x) x,
            
            breaks = function(x) {
              force(range)
              extended_breaks(n = 5)(x)
            },
            format = function(x) {
              force(range)
              x[1] <- paste("<", range[1])
              x[length(x)] <- paste(">", range[2])
              x
            })
  
} 
