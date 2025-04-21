library(randomForest)
library(ggplot2)
library(dplyr)

#Load data
data = read.csv("final_cleaned.csv")[-1,]
data$sasdate = as.Date(data$sasdate)
sum(is.na(data))
data = na.omit(data)
data = data %>%
  filter(sasdate <= as.Date("2019-12-01"))

set.seed(123) #for replication

#Define MSE function:
RMSE <- function(pred, truth){
  return(sqrt(mean((truth - pred)^2)))
} 
############################################
# Bagging and Random Forest for 1-step ahead forecast
############################################

# Set parameters for rolling window and 1 step-ahead forecast
horizon <- 1
test_window_size <- round(nrow(data)/3) #testing dataset of size 1/3 of entire dataset
train_window_size <- nrow(data)-test_window_size - horizon + 1


# Store predictions
pred_bagging_one <- rep(NA, test_window_size)
pred_rf_one <- rep(NA, test_window_size)

# Rolling window loop to get predictions and store
for (i in 1: test_window_size) {
  # Define training and testing data for the current window
  train_data <- data[i:(i + train_window_size - 1), ]
  test_data <- data[(nrow(data)-test_window_size+i):(nrow(data)), ]
  
  # Bagging model
  bagging_model_one <- randomForest(INDPRO ~ ., data = train_data, mtry = ncol(train_data) - 1, ntree = 500, nodesize=5)
  pred_bagging_one[i] <- predict(bagging_model_one, newdata = test_data)[1]
  
  # Random Forest model
  rf_model_one <- randomForest(INDPRO ~ ., data = train_data, ntree = 500, mtry = round((ncol(train_data)-1)/3), nodesize=5)
  pred_rf_one[i] <- predict(rf_model_one, newdata = test_data)[1]
}



# Calculate RMSE
test_data <- data[(nrow(data)-test_window_size+1):(nrow(data)), ]
actuals <- test_data$INDPRO
rmse_bagging_onestep <- RMSE(pred_bagging_one, actuals)
rmse_rf_onestep <- RMSE(pred_rf_one, actuals)

cat("Bagging RMSE for 1-step ahead:", rmse_bagging_onestep, "\n")
cat("Random Forest RMSE for 1-step ahead:", rmse_rf_onestep, "\n")

################################
# Visualization of 1-step ahead forecasts
################################
dates = test_data$sasdate
one_step_ahead_pred = data.frame(actuals, dates, pred_bagging_one, pred_rf_one)

#Bagging plot
Bag_one_plot = ggplot(one_step_ahead_pred) +
  geom_line(aes(x=dates, y=actuals), color = "blue") +
  geom_line(aes(x=dates, y=pred_bagging_one), color = 'red', linetype = 2) +
  labs(title = "1-step ahead forecasts using Bagging", x = "Date", y = "INDPRO predictions") +
  theme_minimal()

#RF plot
RF_one_plot = ggplot(one_step_ahead_pred) +
  geom_line(aes(x=dates, y=actuals), color = "blue") +
  geom_line(aes(x=dates, y=pred_rf_one), color = 'red', linetype = 2) +
  labs(title = "1-step ahead forecasts using RF", x = "Date", y = "INDPRO predictions") +
  theme_minimal()




###################################################
#variable importance plot for 1-step ahead forecast
###################################################

# Extract and normalize variable importance
var_importance <- importance(rf_model_one)
var_importance <- as.data.frame(var_importance)
var_importance$Variable <- rownames(var_importance)
var_importance <- var_importance %>%
  mutate(Importance = IncNodePurity / sum(IncNodePurity) * 100) %>%
  arrange(desc(Importance))

# Select the top 9 variables and combine the rest into "Other"
top_vars <- head(var_importance, 9)
other_vars <- tail(var_importance, nrow(var_importance) - 9)
other_importance <- sum(other_vars$Importance)

# Create a new dataframe with top variables and "Other" category
final_importance <- rbind(
  select(top_vars, Variable, Importance),
  data.frame(Variable = "Other", Importance = other_importance)
)


# Plot the stacked bar chart
plot_one = ggplot(final_importance, aes(x = 1, y = Importance, fill = reorder(Variable, -Importance))) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = NULL, y = "Importance (%)", title = "1-Step") +
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank())



#################################################################################################################################





############################################
# Bagging and Random Forest for 3-step ahead forecast
############################################

# Set parameters for rolling window and 1 step-ahead forecast
horizon <- 3
train_window_size <- nrow(data)-test_window_size - horizon + 1

# Store predictions
pred_bagging_three <- rep(NA, test_window_size)
pred_rf_three <- rep(NA, test_window_size)

# Rolling window loop to get predictions and store
for (i in 1: test_window_size) {
  # Define training and testing data for the current window
  train_data <- data[i:(i + train_window_size - 1), ]
  test_data <- data[(nrow(data)-test_window_size+i):(nrow(data)), ]
  
  # Bagging model
  bagging_model_three <- randomForest(INDPRO ~ ., data = train_data, mtry = ncol(train_data) - 1, ntree = 500, nodesize=5)
  pred_bagging_three[i] <- predict(bagging_model_three, newdata = test_data)[1]
  
  # Random Forest model
  rf_model_three <- randomForest(INDPRO ~ ., data = train_data, ntree = 500, mtry = round((ncol(train_data)-1)/3), nodesize=5)
  pred_rf_three[i] <- predict(rf_model_three, newdata = test_data)[1]
}


# Calculate RMSE
rmse_bagging_threestep <- RMSE(pred_bagging_three, actuals)
rmse_rf_threestep <- RMSE(pred_rf_three, actuals)

cat("Bagging RMSE for 3-step ahead:", rmse_bagging_threestep, "\n")
cat("Random Forest RMSE for 3-step ahead:", rmse_rf_threestep, "\n")

################################
# Visualization of 3-step ahead forecasts
################################
three_step_ahead_pred = data.frame(actuals, dates, pred_bagging_three, pred_rf_three)

#Bagging plot
Bag_three_plot = ggplot(three_step_ahead_pred) +
  geom_line(aes(x=dates, y=actuals), color = "blue") +
  geom_line(aes(x=dates, y=pred_bagging_three), color = 'red', linetype = 2) +
  labs(title = "3-step ahead forecasts using Bagging", x = "Date", y = "INDPRO predictions") +
  theme_minimal()

#RF plot
RF_three_plot = ggplot(three_step_ahead_pred) +
  geom_line(aes(x=dates, y=actuals), color = "blue") +
  geom_line(aes(x=dates, y=pred_rf_three), color = 'red', linetype = 2) +
  labs(title = "3-step ahead forecasts using RF", x = "Date", y = "INDPRO predictions") +
  theme_minimal()

###################################################
#variable importance plot for 3-step ahead forecast
###################################################

# Extract and normalize variable importance
var_importance <- importance(rf_model_three)
var_importance <- as.data.frame(var_importance)
var_importance$Variable <- rownames(var_importance)
var_importance <- var_importance %>%
  mutate(Importance = IncNodePurity / sum(IncNodePurity) * 100) %>%
  arrange(desc(Importance))

# Select the top 9 variables and combine the rest into "Other"
top_vars <- head(var_importance, 9)
other_vars <- tail(var_importance, nrow(var_importance) - 9)
other_importance <- sum(other_vars$Importance)

# Create a new dataframe with top variables and "Other" category
final_importance_three <- rbind(
  select(top_vars, Variable, Importance),
  data.frame(Variable = "Other", Importance = other_importance)
)

# Plot the stacked bar chart
plot_three = ggplot(final_importance_three, aes(x = 1, y = Importance, fill = reorder(Variable, -Importance))) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = NULL, y = "Importance (%)", title = "3-Step") +
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank())


#################################################################################################################################




############################################
# Bagging and Random Forest for 6-step ahead forecast
############################################

# Set parameters for rolling window and 1 step-ahead forecast
horizon <- 6
train_window_size <- nrow(data)-test_window_size - horizon + 1

# Store predictions
pred_bagging_six <- rep(NA, test_window_size)
pred_rf_six <- rep(NA, test_window_size)

# Rolling window loop to get predictions and store
for (i in 1: test_window_size) {
  # Define training and testing data for the current window
  train_data <- data[i:(i + train_window_size - 1), ]
  test_data <- data[(nrow(data)-test_window_size+i):(nrow(data)), ]
  
  # Bagging model
  bagging_model_six <- randomForest(INDPRO ~ ., data = train_data, mtry = ncol(train_data) - 1, ntree = 500, nodesize=5)
  pred_bagging_six[i] <- predict(bagging_model_six, newdata = test_data)[1]
  
  # Random Forest model
  rf_model_six <- randomForest(INDPRO ~ ., data = train_data, ntree = 500, mtry = round((ncol(train_data)-1)/3), nodesize=5)
  pred_rf_six[i] <- predict(rf_model_six, newdata = test_data)[1]
}

# Calculate RMSE
rmse_bagging_sixstep <- RMSE(pred_bagging_six, actuals)
rmse_rf_sixstep <- RMSE(pred_rf_six, actuals)

cat("Bagging RMSE for 6-step ahead:", rmse_bagging_sixstep, "\n")
cat("Random Forest RMSE for 6-step ahead:", rmse_rf_sixstep, "\n")

################################
# Visualization of 6-step ahead forecasts
################################
six_step_ahead_pred = data.frame(actuals, dates, pred_bagging_six, pred_rf_six)

#Bagging plot
Bag_six_plot = ggplot(six_step_ahead_pred) +
  geom_line(aes(x=dates, y=actuals), color = "blue") +
  geom_line(aes(x=dates, y=pred_bagging_six), color = 'red', linetype = 2) +
  labs(title = "6-step ahead forecasts using Bagging", x = "Date", y = "INDPRO predictions") +
  theme_minimal()

#RF plot
RF_six_plot = ggplot(six_step_ahead_pred) +
  geom_line(aes(x=dates, y=actuals), color = "blue") +
  geom_line(aes(x=dates, y=pred_rf_six), color = 'red', linetype = 2) +
  labs(title = "6-step ahead forecasts using RF", x = "Date", y = "INDPRO predictions") +
  theme_minimal()

###################################################
#variable importance plot for 3-step ahead forecast
###################################################

# Extract and normalize variable importance
var_importance <- importance(rf_model_six)
var_importance <- as.data.frame(var_importance)
var_importance$Variable <- rownames(var_importance)
var_importance <- var_importance %>%
  mutate(Importance = IncNodePurity / sum(IncNodePurity) * 100) %>%
  arrange(desc(Importance))

# Select the top 9 variables and combine the rest into "Other"
top_vars <- head(var_importance, 9)
other_vars <- tail(var_importance, nrow(var_importance) - 9)
other_importance <- sum(other_vars$Importance)

# Create a new dataframe with top variables and "Other" category
final_importance_six <- rbind(
  select(top_vars, Variable, Importance),
  data.frame(Variable = "Other", Importance = other_importance)
)

# Plot the stacked bar chart
plot_six = ggplot(final_importance_six, aes(x = 1, y = Importance, fill = reorder(Variable, -Importance))) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = NULL, y = "Importance (%)", title = "6-Step") +
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank())


#################################################################################################################################





############################################
# Bagging and Random Forest for 12-step ahead forecast
############################################

# Set parameters for rolling window and 1 step-ahead forecast
horizon <- 12
train_window_size <- nrow(data)-test_window_size - horizon + 1

# Store predictions
pred_bagging_twelve <- rep(NA, test_window_size)
pred_rf_twelve <- rep(NA, test_window_size)

# Rolling window loop to get predictions and store
for (i in 1: test_window_size) {
  # Define training and testing data for the current window
  train_data <- data[i:(i + train_window_size - 1), ]
  test_data <- data[(nrow(data)-test_window_size+i):(nrow(data)), ]
  
  # Bagging model
  bagging_model_twelve <- randomForest(INDPRO ~ ., data = train_data, mtry = ncol(train_data) - 1, ntree = 500, nodesize=5)
  pred_bagging_twelve[i] <- predict(bagging_model_twelve, newdata = test_data)[1]
  
  # Random Forest model
  rf_model_twelve <- randomForest(INDPRO ~ ., data = train_data, ntree = 500, mtry = round((ncol(train_data)-1)/3), nodesize=5)
  pred_rf_twelve[i] <- predict(rf_model_twelve, newdata = test_data)[1]
}

# Calculate RMSE
rmse_bagging_twelvestep <- RMSE(pred_bagging_twelve, actuals)
rmse_rf_twelvestep <- RMSE(pred_rf_twelve, actuals)

cat("Bagging RMSE for 12-step ahead:", rmse_bagging_twelvestep, "\n")
cat("Random Forest RMSE for 12-step ahead:", rmse_rf_twelvestep, "\n")

################################
# Visualization of 12-step ahead forecasts
################################
twelve_step_ahead_pred = data.frame(actuals, dates, pred_bagging_twelve, pred_rf_twelve)

#Bagging plot
Bag_twelve_plot = ggplot(twelve_step_ahead_pred) +
  geom_line(aes(x=dates, y=actuals), color = "blue") +
  geom_line(aes(x=dates, y=pred_bagging_twelve), color = 'red', linetype = 2) +
  labs(title = "12-step ahead forecasts using Bagging", x = "Date", y = "INDPRO predictions") +
  theme_minimal()

#RF plot
RF_twelve_plot = ggplot(twelve_step_ahead_pred) +
  geom_line(aes(x=dates, y=actuals), color = "blue") +
  geom_line(aes(x=dates, y=pred_rf_twelve), color = 'red', linetype = 2) +
  labs(title = "12-step ahead forecasts using RF", x = "Date", y = "INDPRO predictions") +
  theme_minimal()

###################################################
#variable importance plot for 12-step ahead forecast
###################################################

# Extract and normalize variable importance
var_importance <- importance(rf_model_twelve)
var_importance <- as.data.frame(var_importance)
var_importance$Variable <- rownames(var_importance)
var_importance <- var_importance %>%
  mutate(Importance = IncNodePurity / sum(IncNodePurity) * 100) %>%
  arrange(desc(Importance))

# Select the top 9 variables and combine the rest into "Other"
top_vars <- head(var_importance, 9)
other_vars <- tail(var_importance, nrow(var_importance) - 9)
other_importance <- sum(other_vars$Importance)

# Create a new dataframe with top variables and "Other" category
final_importance_twelve <- rbind(
  select(top_vars, Variable, Importance),
  data.frame(Variable = "Other", Importance = other_importance)
)

# Plot the stacked bar chart
plot_twelve = ggplot(final_importance_twelve, aes(x = 1, y = Importance, fill = reorder(Variable, -Importance))) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = NULL, y = "Importance (%)", title = "12-Step") +
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank())


#################################################################################################################################



##############
# all RMSEs  #
##############
RMSEs_for_bagging = c(rmse_bagging_onestep, rmse_bagging_threestep, rmse_bagging_sixstep, rmse_bagging_twelvestep)
RMSEs_for_RF = c(rmse_rf_onestep, rmse_rf_threestep, rmse_rf_sixstep, rmse_rf_twelvestep)


##############
# all plots  #
##############
library(patchwork)
plot_one + plot_three + plot_six + plot_twelve +
  plot_layout(nrow = 2) +
  plot_annotation(title = 'Variable Importance Plots (Top 10 Variables stacked to 100%)')

Bag_one_plot + Bag_three_plot + Bag_six_plot + Bag_twelve_plot +
  plot_layout(nrow = 2)

RF_one_plot + RF_three_plot + RF_six_plot + RF_twelve_plot +
  plot_layout(nrow = 2)




