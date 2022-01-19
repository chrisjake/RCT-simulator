# Load libraries
library(tidyverse)
library(janitor)
library(mltools)

# Define Functions
## fun_extra takes a dataframe of trials (each row is a trial) and does some extra calculations for each 
fun_extra = function(df){
  df = df %>%
    mutate(result = case_when(type == "t_c" & p_value < var_thresh ~ "true_pos",            # True positive if there is a treatment effect and p < threshold
                              type == "t_c" & p_value >= var_thresh ~ "false_neg",          # False negative if there is a treatment effect and p > threshold
                              type == "c_c" & p_value < var_thresh ~ "false_pos",           # False positive if there is no treatment effect and p < threshold
                              type == "c_c" & p_value >= var_thresh ~ "true_neg"),          # True neg if there is no treatment effect and p > threshold
           sig = ifelse(p_value < var_thresh, 1, 0)) %>%
    ## This next bit bins p-values into var_nbin groups, above and below the var_thresh (via grouping by sig)
    ## This produces a string which is the upper and lower limit of the bin (as.char prevents dplyr getting sad that each group has different bin sizes)
    ## We then split the string to take the upper value, tidy up the trailing punctuation, and make it numeric
    ## This ensures we split p-values either side of var_thresh and don't have one bin that contains both significant and non-significant p-values
    ## Binning strategy is important and a few methods were tried (fixed range values, log scaled values, n-tiles, etc), this seemed to provide the best balance
    ## between bins with a reasonable minimal size of observations (important to minimise extreme results summary statistics for the bin) and a small range of included p-values
    group_by(sig) %>%
    mutate(p_bin = as.character(mltools::bin_data(p_value, bins = var_nbin, binType = "quantile"))) %>%
    ungroup() %>%
    mutate(p_bin = str_split(p_bin, pattern = ","),
           p_bin_l = as.numeric(gsub("^.", "", map_chr(p_bin, 1))),
           #Keep the upper bound as the "true" p_bin
           p_bin = as.numeric(gsub(".$", "", map_chr(p_bin, 2))),
           # Use the midway point between bounds for x-axis position
           p_bin_m = (p_bin_l + p_bin)/2)
  
  df
}

# fun_bin takes the output of fun_trial and returns a dataframe binned by p-value (each row is the number of trials at that p-value level), and calculates the PPV/NPV/FPR/FNR (as appropriate) for each bin
fun_bin = function(df) {
  df = df %>%
    group_by(p_bin, p_bin_m, sig, result) %>%
    summarise(n = n()) %>%
    #Get number of trials which are false or true positive at each p-value level
    mutate(false_pos = ifelse(result == "false_pos", n, NA),
           true_pos = ifelse(result == "true_pos", n, NA),
           false_neg = ifelse(result == "false_neg", n, NA),
           true_neg = ifelse(result == "true_neg", n, NA)) %>%
    ungroup() %>%
    select(-c(sig, result, n)) %>%
    group_by(p_bin, p_bin_m) %>%
    #A hacky solution to get the results for each p-value level on a single row
    summarise(true_pos = sum(true_pos, na.rm = TRUE),
              false_pos = sum(false_pos, na.rm = TRUE),
              true_neg = sum(true_neg, na.rm = TRUE),
              false_neg = sum(false_neg, na.rm = TRUE)) %>%
    #Calculate other key statistics
    mutate(PPV = true_pos / (true_pos + false_pos),
           FPR = false_pos / (true_pos + false_pos),
           NPV = true_neg / (true_neg + false_neg),
           FNR = false_neg / (true_neg + false_neg))
  df
}

# fun_sum takes a binned dataframe (from fun_bin) and generates a summary row

fun_sum = function(df){
  df = df %>%
    ungroup() %>%
    summarise(true_pos = sum(true_pos, na.rm = TRUE),
              false_pos = sum(false_pos, na.rm = TRUE),
              true_neg = sum(true_neg, na.rm = TRUE),
              false_neg = sum(false_neg, na.rm = TRUE)) %>%
    mutate(total = true_pos + true_neg + false_pos + false_neg,
           total_true = true_pos + false_neg,
           total_false = false_pos + true_neg,
           prob_pos_trial = (true_pos + false_pos) / total,
           prob_neg_trial = (true_neg + false_neg) / total,
           PPV = true_pos / (true_pos + false_pos),
           FPR = false_pos / (true_pos + false_pos),
           NPV = true_neg / (true_neg + false_neg),
           FNR = false_neg / (true_neg + false_neg),
           power = true_pos / total_true,
           across(where(is.numeric), round, digits = 3))
  df
}



# * * *



# Load and setup dataframes
## Load all the data and give each dataframe a power level
files = list.files(path = "./data/raw")

for (i in 1:length(files)) {
  file = paste("data/raw/", files[i], sep = "")
  name = str_match(files[i], ".+?(?=[.]Rds)")
  assign(name, readRDS(file)%>%
           mutate(power_level = str_match(file, "0[.][0-9]")))
}

## Put each dataframe into a list of dataframes
anaes = lapply(ls(pattern="anaes_0"), function(x) get(x))
icu = lapply(ls(pattern="icu_0"), function(x) get(x))


## Tidy up the environment
rm(i, file, files, name, list = ls(pattern="(anaes_)|(icu_)"))



# * * *



# Save edited dataframes
## Set Global Variables
var_thresh = 0.05
var_nbin = 100


## Raw data
anaes_raw = lapply(anaes, fun_extra)
icu_raw = lapply(icu, fun_extra)

rm(anaes, icu)


## Binned data
anaes_binned = lapply(anaes_raw, fun_bin)
icu_binned = lapply(icu_raw, fun_bin)


## Summary data
anaes_sum = lapply(anaes_binned, fun_sum)
icu_sum = lapply(icu_binned, fun_sum)


### All together now
anaes_sum = bind_rows(anaes_sum)
icu_sum = bind_rows(icu_sum)


##Write objects
dfs = ls(pattern = "(anaes)|(icu)")

lapply(dfs, function(x){
  saveRDS(get(x), file = paste("data/", x, ".Rds", sep = ""))
})