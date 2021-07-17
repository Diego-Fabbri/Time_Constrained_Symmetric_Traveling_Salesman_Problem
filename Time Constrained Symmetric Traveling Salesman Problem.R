#Set your own working directory
setwd("C:/Users/diego/Documents/R/Projects/GitHub_Projects/Optimization/Time Constrained Symmetric Traveling Salesman Problem")

# Import lpSolve package
library(lpSolve)

#Import required packages (ompr)
library(dplyr)
library(ROI)
library(ROI.plugin.symphony)
library(ompr)
library(ompr.roi)
library(ggplot2)

#Set t0
t0 <- 0

#Set speed
speed <- 1

#Set matrix of costs (distances)
c <- matrix(c(0, 27.9, 54.6, 42.0, 56.5, 37.0, 30.9, 34.1,
             27.9, 0, 67.2, 25.6, 28.8, 48.4, 57.4, 21.6,
             54.6, 67.2, 0, 60.5, 95.8, 18.8, 60.4, 52.1,
             42.0, 25.6, 60.5, 0, 39.4, 43.1, 70.2, 12.2,
             56.5, 28.8, 95.8, 39.4, 0, 77.2, 84.5, 44.4,
             37.0, 48.4, 18.8, 43.1, 77.2, 0, 51.6, 34.0,
             30.9, 57.4, 60.4, 70.2, 84.5, 51.6, 0, 59.3,
             34.1, 21.6, 52.1, 12.2, 44.4, 34.0, 59.3, 0), nrow = 8, byrow = TRUE)

#Set n (depot + nodes to visit = problem size)
n <- ncol(c)

#Set matrix of travel duration
travel_time <- array(dim = c(n, n))

for (i in 1:n) {
  for (j in 1:n) {
    travel_time[i, j] <- c[i, j]/speed
  }
}

#Set Lower Bounds
LB <- c(0, 0, 0, 0, 0, 0, 0, 0)

#Set Upper Bounds
UB <- c(1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000)

#Set matrix M
M <- array(dim = c(n, n))

for (i in 1:n) {
  for (j in 1:n) {
    M[i, j] <- UB[i] - LB[i] + travel_time[i, j]
  }
}

#Build Model
Model <- MIPModel() %>%
  add_variable(x[i], i = 1:(n+1), type = "continuous") %>% #define variables
  add_variable(y[i,j], i = 1:(n+1), j = 1:(n+1), i!=j, type = "binary") %>%
  set_objective(x[n+1]-x[1], "min") %>% #define objective
  add_constraint(x[i] - x[1] >= travel_time[1,i], i = 2:n) %>% #define constraints
  add_constraint(x[n+1] - x[i] >= travel_time[i,1], i = 2:n) %>%
  add_constraint(x[i] - x[j] >= travel_time[i,j] -M[i, j]*y[i, j], i = 1:n, j = 1:n, i!=j) %>%
  add_constraint(x[j] - x[i] >= travel_time[i,j] -M[i, j] + M[i, j]*y[i, j], i = 1:n, j = 1:n, i!=j) %>%
  add_constraint(x[i] >= LB[i], i = 1:n) %>%
  add_constraint(x[i] <= UB[i], i = 1:n) %>%
  add_constraint(x[1] == t0) %>%
  solve_model(with_ROI(solver = "symphony", verbosity = 1))

#Model summary
##Status
print(paste("Model status is:", Model$status))

##Objective Function
print(paste("Objective value:", objective_value(Model)))

for (a in 1:(n+1)) {
    tmp_x <- get_solution(Model, x[i]) %>%
      filter(variable == "x", i == a) %>%
      select(value)
    
    if (tmp_x != 0) {
      print(paste("--->x[", (a-1), "] =", tmp_x))
    }
}




