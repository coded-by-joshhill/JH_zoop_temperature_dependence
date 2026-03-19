mutate(
  variance = case_when( # generate a variance variable...
    # where data is mean and there is data for SD or SE...
    data_type == "Mean" & rate_error_type == "Std dev" & !is.na(n_measurements_sampled) ~ (rate_error^2)/n_measurements_sampled,
    data_type == "Mean" & rate_error_type == "Std err" ~ (rate_error^2),
    TRUE ~ NA_real_),
  weight = case_when( # generate a weighting variable...
    # Where weight = 1/variance
    data_type == "Mean" & !is.na(variance) ~ 1/variance,
    TRUE ~ 1), # set all other observations (i.e., replicates) to 1
  weight = ifelse(is.na(weight) | is.infinite(weight), 1, weight),
  weight = pmax(weight, 1e-6), # apply a floor to prevent zero or negative weights
  weight = pmin(weight, 1000)) # apply a ceiling to prevent overly precise means from influencing model estimates

sum(dat$weight == 1000) # What is the sum of the observations hitting the cap
# 4 - great, not many
mean(dat$weight == 1000) * 100 # What % of my data is hitting the 1000 cap?
# ~0.4 %