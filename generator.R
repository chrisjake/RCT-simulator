# Load libraries
library(tidyverse)
library(parallel)
library(doParallel)
library(foreach)

# Define key functions
## trial_run generates a contingency matrix from a single Bernoulli simulation
trial_run = function(control_event_rate, effect_size, n){
  control_dead = rbinom(1, n, control_event_rate)            # Randomly generate number of control group events
  
  treat_event_rate = control_event_rate - effect_size
  treat_dead = rbinom(1, n, treat_event_rate)                # Randomly generate number of treatment group events
  
  trial_matrix = matrix(c(n - control_dead, control_dead,    # Generate contingency table
                          n - treat_dead, treat_dead),
                        nrow = 2,
                        ncol = 2,
                        byrow  = TRUE)
  
  result = data.frame(control_alive = trial_matrix[1,1],
                      control_dead = trial_matrix[1,2],
                      treat_alive = trial_matrix[2,1],
                      treat_dead = trial_matrix[2,2],
                      p_value = fisher.test(trial_matrix)$p.value)
  
  result                                           
}


## Calls trial_run a sim_num number of times, feeding it the set parameters
simulation = function(sim_num, control_event_rate, treat_effect_size, trial_n, prior_prob) {
  # Pre-test/prior probability is adjusted by control the proportion of simulations where the effect size is non-zero
  # Variables are perhaps named poorly: control_vs_control = no true effect of the intervention
  control_vs_control = do.call("rbind",
                               replicate(sim_num,
                                         trial_run(control_event_rate = control_event_rate,   
                                                   effect_size = 0,                #trial_run called with an effect size of 0 for control vs control
                                                   n = round(trial_n * prior_prob, 0)),         #number of c_c runs = total, * prior prob
                                         simplify = FALSE))
  
  # Treat_vs_control = effect 
  treat_vs_control = do.call("rbind",
                             replicate(sim_num,
                                       trial_run(control_event_rate = control_event_rate,
                                                 effect_size = treat_effect_size,
                                                 n = round(trial_n * (1 - prior_prob), 0)),     #number of t_c runs = opposite of c_c
                                       simplify = FALSE))
  #Merge simulations
  control_vs_control$type = "c_c"
  treat_vs_control$type = "t_c"
  
  data = dplyr::full_join(control_vs_control, treat_vs_control)
  data
}



# Generate data sets
## Goal is to calculate 18 different data sets with 500000 simulations in each
## All power adjustments will be made by varying n
## "Anaesthesia trials": Event rate 5%, effect size 1.5%; across powers 0.9 to 0.1 (by 0.1 increments), at p < 0.05
## "ICU trials": Event rate 20%, effect size 5%; across powers 0.9 to 0.1 (by 0.1 increments), at p < 0.05
## Prior probability is 0.5 for all trials

# sim_num, control_event_rate, treat_effect_size, trial_n, prior_prob

## Generate a list with the desired input parameters
var_anaes = data.frame(sim_num = 500000,
                       control_event_rate = 0.05,
                       effect_size = 0.015,
                       power = seq(0.1, 0.9, 0.1),
                       alpha = 0.05,
                       prior_probability = 0.5,
                       name = "anaes")

var_icu = data.frame(sim_num = 500000,
                     control_event_rate = 0.2,
                     effect_size = 0.05,
                     power = seq(0.1, 0.9, 0.1),
                     alpha = 0.05,
                     prior_probability = 0.5,
                     name = "icu")

var_list = full_join(var_anaes, var_icu) %>%
  mutate(name = paste(name, power, sep = "_")) %>%
  group_by(name) %>%
  # Use prop test to generate sample size, double it because it calculates size for one group. Rounded in simulation fn so no need here. 
  mutate(n = (power.prop.test(n = NULL,
                              p1 = control_event_rate,
                              p2 = control_event_rate - effect_size,
                              sig.level = alpha,
                              power = power)$n * 2))


saveRDS(var_list, file = "data/gen_settings.Rds")
rm(var_anaes, var_icu)


## Generate the output datasets with parallel methods
numCores = round(detectCores() * 0.6, 0) #Keep some cores for other things

doParallel::registerDoParallel(numCores)

foreach::foreach (i = 1:nrow(var_list)) %dopar% {
  params = var_list[i,]
  
  data = simulation(sim_num = params$sim_num,
                    control_event_rate = params$control_event_rate,
                    treat_effect_size = params$effect_size,
                    trial_n = params$n,
                    prior_prob = params$prior_probability)
  
  saveRDS(data, file = paste("data/raw/", params$name, ".Rds", sep = ""))
}