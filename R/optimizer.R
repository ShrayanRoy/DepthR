

# Utilities for optimizing function in a Step-wise grid search manner
# In each step, We find the k-maximas of the function in current grid and
# consider union of nbds of those maximas and consider it as new grid
# step size and length of nbd decreases geometrically in each step


# Function to find k-maximas of the function
find_max_indices <- function(vec,k = 3) {

  sorted_indices <- order(vec, decreasing = TRUE)
  max_indices <- sorted_indices[1:k]

  return(max_indices)
}

find_neighboring_grid <- function(vec, step_size,lower = 1,upper = 9){

  neighboring_points <- list()

  for (elem in vec){
    neighbors <- seq(max(lower,elem - 2^(step_size)*step_size), min(elem + 2^(step_size)*step_size,upper), by = step_size)
    neighboring_points[[as.character(elem)]] <- neighbors
  }

  return(unique(unlist(neighboring_points)))
}

# Step wise Grid search optimization function
gridStepSearch <- function(f,lower_bound = 1, upper_bound = 10,show.prog = FALSE,
                                     initial_step_size = 1, min_step_size = 0.001,max_iter = 10) {

  step_size <- initial_step_size
  grid_points <- seq(lower_bound, upper_bound, by = step_size)

  for (i in 1:max_iter){

    f_values <- sapply(grid_points, f)
    max_indices <- find_max_indices(f_values)
    x_val <- grid_points[max_indices]
    f_val <- f_values[max_indices]

    if(show.prog){
      print(paste("Itration",i,":",x_val[which.max(f_val)]))
    }

    x_best <- x_val[which.max(f_val)]
    step_size <- max(step_size*0.5, min_step_size)
    grid_points <- find_neighboring_grid(x_val,step_size = step_size,
                              lower = lower_bound,upper = upper_bound)

  }

  f_best <- max(f_val)
  return(list(x = x_best, f_value = f_best))
}

# Example
f <- function(x) {
  return(-(sin(x) + sin(2*x) + cos(3*x) + log(x + 1) + sqrt(x)))
}

result <- gridStepSearch(f, lower_bound = 1, upper_bound = 10,
                          initial_step_size = 1, min_step_size = 0.001, max_iter = 20)
print(result)
curve(f(x),from = 1,to = 10)
abline(v = result$x,col = "red",main = "Plot of f",lty = 2,lwd = 2)



