rm(list = ls())


# define constants applied
constant = data.frame("V" = 49000, "sat_to_conc" = 0.083, "c_f_O2" = 8.73,
                      "c_r_O2" = 8.29, "c_f_CO2" = 1.0, "xi" = 0.75)


# load packages
library(data.table)
library(ggplot2)
library(ggpubr)
library(readxl)


# DEFINE PREPROCESSING FUNCTIONS

# function to update column names to shorter terms
colnames_update <- function(df){
  colnames(df) <- c("t", "O2_added_1", "O2_added_2", "s_O2_center", "s_O2_edge", 
                    "c_CO2", "q_f", "q_r", "T", "n", "M")
  return(df)
} 

# function to setting daily data to all datapoints in the day
daily_obs = function(vec){
  tmp = vec[which(!is.na(vec))[1]]
  n = length(vec)
  new = rep(NA, n)
  for(i in 1:n){
    if(!is.na(vec[i])){tmp = vec[i]}
    new[i] = tmp
  }
  return(new)
}

# function to set NA values to the previous value
na_to_prev <- function(vec){
  for(i in 2:length(vec)){
    if(is.na(vec[i])){vec[i] = vec[i-1]}
  }
  return(vec)
}

# function to remove all NAs in our specific dataframe
na_to_prev_df = function(df){
  df$O2_added_1 = na_to_prev(df$O2_added_1)
  df$O2_added_2 = na_to_prev(df$O2_added_2)
  df$s_O2_center = na_to_prev(df$s_O2_center)
  df$s_O2_edge = na_to_prev(df$s_O2_edge)
  df$c_CO2 = na_to_prev(df$c_CO2)
  df$T = na_to_prev(df$T)
  df$q_f = na_to_prev(df$q_f)
  df$q_r = na_to_prev(df$q_r)
  df$n = daily_obs(df$n)
  df$M = daily_obs(df$M)
  return(df)
}

# negative outlier removal
clean_neg_outliers <- function(vec, tol){
  for(i in 2:length(vec)){
    if(vec[i] <= tol){vec[i] = vec[i-1]}
  }
  return(vec)
}

# positive outlier removal
clean_pos_outliers <- function(vec, tol){
  for(i in 2:length(vec)){
    if(vec[i] > tol){vec[i] = vec[i-1]}
  }
  return(vec)
}

clean_outliers = function(vec, ntol, ptol){
  vec = clean_neg_outliers(vec = vec, tol = ntol)
  vec = clean_pos_outliers(vec = vec, tol = ptol)
  return(vec)
}

# function for 5 min differentiation
new_O2_added_5min = function(vec, df){
  newt = seq(as.POSIXct(df$t[1]), as.POSIXct(df$t[length(df$t)]), by = (60*5))
  l = length(newt)
  n = length(vec)
  new = rep(NA, l)
  for(i in 1:l){
    tmp = 5*i
    diff = vec[tmp] - vec[tmp - 4]
    new[i] = diff
  }
  return(new)
}

# function for 5 min average
new_avg_5min = function(vec, df){
  newt = seq(as.POSIXct(df$t[1]), as.POSIXct(df$t[length(df$t)]), by = (60*5))
  l = length(newt)
  n = length(vec)
  new = rep(NA, l)
  for(i in 1:l){
    tmp = 5*i
    avg = mean(vec[(tmp - 4):tmp])
    new[i] = avg
  }
  return(new)
}

# function for 5 min sum
new_sum_5min = function(vec, df){
  newt = seq(as.POSIXct(df$t[1]), as.POSIXct(df$t[length(df$t)]), by = (60*5))
  l = length(newt)
  n = length(vec)
  new = rep(NA, l)
  for(i in 1:l){
    tmp = 5*i
    tot = sum(vec[(tmp - 4):tmp])
    new[i] = tot
  }
  return(new)
}

# function for 5 minute exact
new_exact_5min = function(vec, df){
  newt = seq(as.POSIXct(df$t[1]), as.POSIXct(df$t[length(df$t)]), by = (60*5))
  l = length(newt)
  n = length(vec)
  new = rep(NA, l)
  for(i in 1:l){
    tmp = 5*i
    val = vec[tmp]
    new[i] = val
  }
  return(new)
}


#henry coefficient
henry <- function(T){
  T = T + 273.15
  T0 = 293.15
  HR = 1700
  H0 = 4.11*1e4
  H = H0*exp(-HR*((1/T) - (1/T0)))
  return(H)
}

full_sat_to_conc <- function(T){
  T = T + 273.15
  T0 = 293.15
  HR = 1700
  H0 = 4.11*1e4
  H = H0*exp(-HR*((1/T) - (1/T0)))
  cs = 32*1e3*55.6*1*0.21/H
  return(cs)
}

# function for o2 saturation to o2 concentration
sat_to_conc <- function(vec, T){
  T = T + 273.15
  T0 = 293.15
  HR = 1700
  H0 = 4.11*1e4
  H = H0*exp(-HR*((1/T) - (1/T0)))
  cs = 32*1e3*55.6*1*0.21/H
  vec = 0.09*vec #cs*vec
  return(vec)
}

# function to change 1 min sampling to 5 min sampling for our df
new_df_5min = function(df){
  newt = seq(as.POSIXct(df$t[1]), as.POSIXct(df$t[length(df$t)]), by = (60*5))
  new = data.frame("t" = newt)
  new$O2_added_1 = new_O2_added_5min(df$O2_added_1, df)*0.2
  new$O2_added_2 = new_O2_added_5min(df$O2_added_2, df)*0.2
  new$c_O2_center = sat_to_conc(new_avg_5min(df$s_O2_center, df), df$T)
  new$c_O2_edge = sat_to_conc(new_avg_5min(df$s_O2_edge, df), df$T)
  new$c_CO2 = new_avg_5min(df$c_CO2, df)
  new$q_f = new_sum_5min(df$q_f, df)*0.2
  new$q_r = new_sum_5min(df$q_r, df)*0.2
  new$T = new_avg_5min(df$T, df)
  new$n = new_exact_5min(df$n, df)
  new$M = 1e-3*new_exact_5min(df$M, df)
  return(new)
}

# function for calculating h 
calc_h = function(df){
  h = df$O2_added_1 + df$O2_added_2 +
    constant$c_f_O2*df$q_f*1e-6 + #kg/min
    constant$c_r_O2*df$q_r*1e-6 - #kg/min
    df$c_O2_center*1e-6*(df$q_f + df$q_r)
  return(h)
}

# function for calculating g
calc_g = function(df){
  g = df$c_CO2*(df$q_f + df$q_r)*1e-6 - #kg/min
    constant$c_f_CO2*1e-6*df$q_f - #kg/min
    (1-constant$xi)*df$c_CO2*df$q_r*1e-6
  return(g)
}

#function for caalculating total O2 added, l(t)
calc_l = function(df){
  l = df$O2_added_1 + df$O2_added_2
  return(l)
}

# function to assess all data visually
plot_all_init = function(df){
  par(mfrow = c(2,5))
  plot(df$O2_added_1, type = "l")
  plot(df$O2_added_2, type = "l")
  plot(df$s_O2_center, type = "l")
  plot(df$s_O2_edge, type = "l")
  plot(df$c_CO2, type = "l")
  plot(df$q_f, type = "l")
  plot(df$q_r, type = "l")
  plot(df$T, type = "l")
  plot(df$n[which(!is.na(df$n))], type = "l")
  plot(df$M[which(!is.na(df$M))], type = "l")
}


all_data = function(){
  
  # load data
  A11 = fread("A11_2022-10-01_2022-12-01.csv")
  A13 = fread("A13_2023-01-17_2023-02-17.csv")
  A21.1 = fread("A21_2022-09-22_2022-11-22.csv")
  A21.2 = fread("A21_2022-12-08_2023-02-08.csv")

  # colnames update 
  A11 <- colnames_update(A11)
  A13 <- colnames_update(A13)
  A21.1 <- colnames_update(A21.1)
  A21.2 <- colnames_update(A21.2)
  

  
  # remove NAs
  A11 = na_to_prev_df(A11)
  A13 = na_to_prev_df(A13)
  A21.1 = na_to_prev_df(A21.1)
  A21.2 = na_to_prev_df(A21.2)
  
  # remove outliers in A11
  A11$O2_added_1 = clean_neg_outliers(A11$O2_added_1, 30000)
  A11$O2_added_2 = clean_neg_outliers(A11$O2_added_2, 30000)
  A11$s_O2_center = clean_neg_outliers(A11$s_O2_center, 50)
  A11$s_O2_center = clean_pos_outliers(A11$s_O2_center, 115)
  A11$s_O2_edge = clean_neg_outliers(A11$s_O2_edge, 50)
  A11$s_O2_edge = clean_pos_outliers(A11$s_O2_edge, 115)
  A11$c_CO2 = clean_neg_outliers(A11$c_CO2, 3)
  A11$c_CO2 = clean_pos_outliers(A11$c_CO2, 25)
  A11$q_f = clean_neg_outliers(A11$q_f, 15000)
  A11$q_f = clean_pos_outliers(A11$q_f, 35000)
  A11$q_r = clean_neg_outliers(A11$q_r, 15000)
  A11$q_r = clean_pos_outliers(A11$q_r, 30000)
  A11$T = clean_neg_outliers(A11$T, 7)
  
  # remove outliers in A13
  A13$O2_added_1 = clean_neg_outliers(A13$O2_added_1, 22000)
  A13$O2_added_2 = clean_neg_outliers(A13$O2_added_2, 20500)
  A13$s_O2_center = clean_neg_outliers(A13$s_O2_center, 76)
  A13$s_O2_center = clean_pos_outliers(A13$s_O2_center, 110)
  A13$s_O2_edge = clean_neg_outliers(A13$s_O2_edge, 88)
  A13$s_O2_edge = clean_pos_outliers(A13$s_O2_edge, 130)
  A13$c_CO2 = clean_neg_outliers(A13$c_CO2, 3)
  A13$c_CO2 = clean_pos_outliers(A13$c_CO2, 25)
  A13$q_f = clean_neg_outliers(A13$q_f, 15000)
  A13$q_f = clean_pos_outliers(A13$q_f, 35000)
  A13$q_r = clean_neg_outliers(A13$q_r, 15000)
  A13$q_r = clean_pos_outliers(A13$q_r, 30000)
  
  # remove outliers in A21.1
  A21.1$O2_added_1 = clean_neg_outliers(A21.1$O2_added_1, 12000)
  A21.1$O2_added_2 = clean_neg_outliers(A21.1$O2_added_2, 12500)
  A21.1$s_O2_center = clean_neg_outliers(A21.1$s_O2_center, 60)
  A21.1$s_O2_center = clean_pos_outliers(A21.1$s_O2_center, 100)
  A21.1$s_O2_edge = clean_neg_outliers(A21.1$s_O2_edge, 60)
  A21.1$s_O2_edge = clean_pos_outliers(A21.1$s_O2_edge, 110)
  A21.1$c_CO2 = clean_neg_outliers(A21.1$c_CO2, 5)
  A21.1$c_CO2 = clean_pos_outliers(A21.1$c_CO2, 25)
  A21.1$q_f = clean_neg_outliers(A21.1$q_f, 15000)
  A21.1$q_f = clean_pos_outliers(A21.1$q_f, 35000)
  A21.1$q_r = clean_neg_outliers(A21.1$q_r, 15000)
  A21.1$q_r = clean_pos_outliers(A21.1$q_r, 36000)
  A21.1$T = clean_neg_outliers(A21.1$T, 6)
  
  # remove outliers in A21.2
  #A21.2$O2_added_1 = clean_neg_outliers(A21.2$O2_added_1, 12000)
  #A21.2$O2_added_2 = clean_neg_outliers(A21.2$O2_added_2, 12500)
  A21.2$s_O2_center = clean_neg_outliers(A21.2$s_O2_center, 60)
  A21.2$s_O2_center = clean_pos_outliers(A21.2$s_O2_center, 110)
  A21.2$s_O2_edge = clean_neg_outliers(A21.2$s_O2_edge, 60)
  A21.2$s_O2_edge = clean_pos_outliers(A21.2$s_O2_edge, 120)
  #A21.2$c_CO2 = clean_neg_outliers(A21.2$c_CO2, 5)
  #A21.2$c_CO2 = clean_pos_outliers(A21.2$c_CO2, 25)
  A21.2$q_f = clean_neg_outliers(A21.2$q_f, 17000)
  A21.2$q_f = clean_pos_outliers(A21.2$q_f, 35000)
  A21.2$q_r = clean_neg_outliers(A21.2$q_r, 15000)
  A21.2$q_r = clean_pos_outliers(A21.2$q_r, 36000)
  A21.2$T = clean_neg_outliers(A21.2$T, 6)
  
  # 1 min sampling to 5 min sampling with concentrations instead of saturations
  A11 = new_df_5min(A11)
  A13 = new_df_5min(A13)
  A21.1 = new_df_5min(A21.1)
  A21.2 = new_df_5min(A21.2)
  
  # add data tag h
  A11["h"] = calc_h(A11)
  A13["h"] = calc_h(A13)
  A21.1["h"] = calc_h(A21.1)
  A21.2["h"] = calc_h(A21.2)
  
  # add data tag g
  A11["g"] = calc_g(A11)
  A13["g"] = calc_g(A13)
  A21.1["g"] = calc_g(A21.1)
  A21.2["g"] = calc_g(A21.2)
  
  # add data tag l
  A11["l"] = calc_l(A11)
  A13["l"] = calc_l(A13)
  A21.1["l"] = calc_l(A21.1)
  A21.2["l"] = calc_l(A21.2)
  
  # clean outliers in h and g
  
  return(list("A11" = A11, "A13" = A13, "A21.1" = A21.1, "A21.2" = A21.2))
}

all_data = all_data()

# exctract parameters of interest from dataframe
feed_preprocess <- function(feed.df){
  n = length(feed.df$Timestamp)
  time = feed.df$Timestamp[2:n]
  value = diff(as.numeric(feed.df$Average))*0.2
  df = data.frame("t" = time, "m" = value)
  return(df)
}


# combine feed data to rest of data
combine.df <- function(feed.df, all.df){
  start = which(all.df$t == feed.df$t[1])
  end = which(all.df$t == feed.df$t[length(feed.df$t)])
  new.all = all.df[c(start:end),]
  
  n = length(new.all$t)
  
  new.feed = rep(NA, n)
  
  for(i in 1:n){
    j = match(new.all$t[i], feed.df$t)
    new.feed[i] = feed.df$m[j]
  }
  new.feed[is.na(new.feed)] = 0
  new.feed[new.feed < 0] = 0
  
  new.all["m"] = new.feed
  
  return(new.all)
}


# function to return finished preprocessed data
complete_data <- function(){
  
  # extract A11, A13, A21.1 and A21.2 to not run all_data() every time
  #all_data = all_data()
  
  A11 = all_data$A11
  A13 = all_data$A13
  A21.1 = all_data$A21.1
  A21.2 = all_data$A21.2
  
  A11.221010 <- combine.df(feed_preprocess(fread("A11_221010_feeding.csv")), A11)
  A11.221011 = combine.df(feed_preprocess(fread("A11_221011_feeding.csv")), A11)
  A11.221012 = combine.df(feed_preprocess(fread("A11_221012_feeding.csv")), A11)
  A11.221013 = combine.df(feed_preprocess(fread("A11_221013_feeding.csv")), A11)
  A11.221116 = combine.df(feed_preprocess(read_excel("A11_221116_feeding.xlsx")), A11)
  
  A13.230120 = combine.df(feed_preprocess(fread("A13_230120_feeding.csv")), A13)
  A13.230121 = combine.df(feed_preprocess(fread("A13_230121_feeding.csv")), A13)
  A13.230214 = combine.df(feed_preprocess(fread("A13_230214_feeding.csv")), A13)
  A13.230215 = combine.df(feed_preprocess(fread("A13_230215_feeding.csv")), A13)
  
  A21.221004 = combine.df(feed_preprocess(fread("A21_221004_feeding.csv")), A21.1)
  A21.221006 = combine.df(feed_preprocess(fread("A21_221006_feeding.csv")), A21.1)
  A21.221117 = combine.df(feed_preprocess(fread("A21_221117_feeding.csv")), A21.1)
  A21.221118 = combine.df(feed_preprocess(fread("A21_221118_feeding.csv")), A21.1)
  
  A21.221227 = combine.df(feed_preprocess(fread("A21_221227_feeding.csv")), A21.2)
  A21.221228 = combine.df(feed_preprocess(fread("A21_221228_feeding.csv")), A21.2)
  A21.230116 = combine.df(feed_preprocess(fread("A21_230116_feeding.csv")), A21.2)
  A21.230117 = combine.df(feed_preprocess(fread("A21_230117_feeding.csv")), A21.2)
  
  return(list("A11.221010" = A11.221010, "A11.221011" = A11.221011, 
              "A11.221012" = A11.221012, "A11.221013" = A11.221013,
              "A11.221116" = A11.221116, "A13.230120" = A13.230120, 
              "A13.230121" = A13.230121, "A13.230214" = A13.230214, 
              "A13.230215" = A13.230215, "A21.221004" = A21.221004, 
              "A21.221006" = A21.221006, "A21.221117" = A21.221117, 
              "A21.221118" = A21.221118, "A21.221227" = A21.221227, 
              "A21.221228" = A21.221228, "A21.230116" = A21.230116, 
              "A21.230117" = A21.230117))
}

complete_data = complete_data()

for(i in 1:5){
  complete_data[[i]]["G1"] = rep(1, nrow(complete_data[[i]]))
  complete_data[[i]]["G2"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G3"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G4"] = rep(0, nrow(complete_data[[i]]))
}

for(i in 6:9){
  complete_data[[i]]["G1"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G2"] = rep(1, nrow(complete_data[[i]]))
  complete_data[[i]]["G3"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G4"] = rep(0, nrow(complete_data[[i]]))
}

for(i in 10:13){
  complete_data[[i]]["G1"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G2"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G3"] = rep(1, nrow(complete_data[[i]]))
  complete_data[[i]]["G4"] = rep(0, nrow(complete_data[[i]]))
}

for(i in 13:17){
  complete_data[[i]]["G1"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G2"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G3"] = rep(0, nrow(complete_data[[i]]))
  complete_data[[i]]["G4"] = rep(1, nrow(complete_data[[i]]))
}

head(complete_data[[13]])
# some manual outlier removal

# test: c(5,7,11,13,16)

plot(complete_data[[1]]$h, type = "l") # B
plot(complete_data[[2]]$h, type = "l") # D
plot(complete_data[[3]]$h, type = "l") # B
plot(complete_data[[4]]$h, type = "l") # B
# test set 5
plot(complete_data[[5]]$h, type = "l") # B
complete_data[[6]]$h = clean_neg_outliers(complete_data[[6]]$h, 0.3)
complete_data[[6]]$h = clean_pos_outliers(complete_data[[6]]$h, 0.8)
plot(complete_data[[6]]$h, type = "l") # F
complete_data[[7]]$h = clean_neg_outliers(complete_data[[7]]$h, 0.4)
complete_data[[7]]$h = clean_pos_outliers(complete_data[[7]]$h, 0.65)
#test set 7
plot(complete_data[[7]]$h, type = "l") # F
complete_data[[8]]$h = clean_neg_outliers(complete_data[[8]]$h, 0.35)
complete_data[[8]]$h = clean_pos_outliers(complete_data[[8]]$h, 0.55)
plot(complete_data[[8]]$h, type = "l") # D
complete_data[[9]]$h = clean_neg_outliers(complete_data[[9]]$h, 0.35)
complete_data[[9]]$h = clean_pos_outliers(complete_data[[9]]$h, 0.55)
plot(complete_data[[9]]$h, type = "l") # D
plot(complete_data[[10]]$h, type = "l") # E
complete_data[[10]]$h = clean_neg_outliers(complete_data[[10]]$h, 0.3)
complete_data[[10]]$h = clean_pos_outliers(complete_data[[10]]$h, 0.9)
plot(complete_data[[11]]$h, type = "l") # A
#test set 11
plot(complete_data[[12]]$h, type = "l") # A
# test set 13
plot(complete_data[[13]]$h, type = "l") # B
complete_data[[14]]$h = clean_neg_outliers(complete_data[[14]]$h, 0)
complete_data[[14]]$h = clean_pos_outliers(complete_data[[14]]$h, 0.8)
plot(complete_data[[14]]$h, type = "l") # F
complete_data[[15]]$h = clean_neg_outliers(complete_data[[15]]$h, 0)
complete_data[[15]]$h = clean_pos_outliers(complete_data[[15]]$h, 1)
plot(complete_data[[15]]$h, type = "l") # F
complete_data[[16]]$h = clean_neg_outliers(complete_data[[16]]$h, 0)
complete_data[[16]]$h = clean_pos_outliers(complete_data[[16]]$h, 1)
# test set 16
plot(complete_data[[16]]$h, type = "l") # F
complete_data[[17]]$h = clean_neg_outliers(complete_data[[17]]$h, 0.3)
complete_data[[17]]$h = clean_pos_outliers(complete_data[[17]]$h, 0.6)
plot(complete_data[[17]]$h, type = "l") # F

plot(complete_data[[1]]$g, type = "l") # A
plot(complete_data[[2]]$g, type = "l") # D
plot(complete_data[[3]]$g, type = "l") # C
plot(complete_data[[4]]$g, type = "l") # A
# test set 5
plot(complete_data[[5]]$g, type = "l") # A
plot(complete_data[[6]]$g, type = "l") # A
# test set 7
plot(complete_data[[7]]$g, type = "l") # F
plot(complete_data[[8]]$g, type = "l") # D
plot(complete_data[[9]]$g, type = "l") # D
plot(complete_data[[10]]$g, type = "l") # E
plot(complete_data[[11]]$g, type = "l") # A
plot(complete_data[[12]]$g, type = "l") # A
# test set 13
plot(complete_data[[13]]$g, type = "l") # B
plot(complete_data[[14]]$g, type = "l") # F
plot(complete_data[[15]]$g, type = "l") # F
# test set 16
plot(complete_data[[16]]$g, type = "l") # F
plot(complete_data[[17]]$g, type = "l") # F

complete_data[[1]]$T = clean_neg_outliers(complete_data[[1]]$T, 11.85)
plot(complete_data[[1]]$T, type = "l") # A
plot(complete_data[[2]]$T, type = "l") # D
plot(complete_data[[3]]$T, type = "l") # C
plot(complete_data[[4]]$T, type = "l") # A
complete_data[[5]]$T = clean_neg_outliers(complete_data[[5]]$T, 11.85)
plot(complete_data[[5]]$T, type = "l") # A
plot(complete_data[[6]]$T, type = "l") # A
complete_data[[7]]$T = clean_neg_outliers(complete_data[[7]]$T, 12.28)
plot(complete_data[[7]]$T, type = "l") # F
plot(complete_data[[8]]$T, type = "l") # D
plot(complete_data[[9]]$T, type = "l") # D
plot(complete_data[[10]]$T, type = "l") # E
complete_data[[11]]$T = clean_neg_outliers(complete_data[[11]]$T, 8.58)
plot(complete_data[[11]]$T, type = "l") # A
plot(complete_data[[12]]$T, type = "l") # A
plot(complete_data[[13]]$T, type = "l") # B
plot(complete_data[[14]]$T, type = "l") # F
complete_data[[15]]$T = clean_neg_outliers(complete_data[[15]]$T, 13.76)
plot(complete_data[[15]]$T, type = "l") # F
plot(complete_data[[16]]$T, type = "l") # F
complete_data[[17]]$T = clean_neg_outliers(complete_data[[17]]$T, 13.81)
plot(complete_data[[17]]$T, type = "l") # F




###############################################################################
# Index for each observation period, for SPDE modelling                       #
###############################################################################

all_in_1 = function(complete_data){
  
  # length of data list
  n = length(complete_data)
  
  # index instead of time
  for(i in 1:n){
    complete_data[[i]]$time.id = 1:nrow(complete_data[[i]])
    complete_data[[i]]$data.id = i
  }
  
  # rowbind all dataframes
  full_data = complete_data[[1]]
  for(i in 2:n){
    full_data = rbind(full_data, complete_data[[i]])
  }
  
  full_data = full_data[order(full_data$time.id),]
  
  return(full_data)
}

full_data = all_in_1(complete_data)

mean(complete_data[[1]]$M)
