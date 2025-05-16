# date: 10 may 2025
# author: manuel gonzález aramburu
# purpose: function to compute C index for cox model by inputting new data based 
#          on the values of the variables and the time and status of the subjects

cindex <- function(betas, new_x, new_y){
  
  # dats <data frame or tibble>: contains x, time and status for all subjects/ patients
  # varss <vector>: names (char vector) or positions (numeric vector) within dats of t
  # betas <data frame or tibble>: values of the estimated betas to be tested
  #   IMPORTANT: the data frame of tibble must have two columns named x and coefs,
  #              x with the names of the variables and coefs with the estimated betas
  # new_x <data frame or tibble or matrix>: values of the variables corresponding to betas that will be used to compute the c-index 
  # new_y <Surv object from survival library>: Surv object with the times and status of the observations  that will be used to compute the c-index
  #   IMPORTANT: new_x and new_y must have the same order so row 1 of new_x and row 1 of new_y correspond
  #              to the exact same individual and so on
  
  # transpose the new data to fit in tidydata
  
  new_x <- 
    as_tibble(new_x) |> 
    mutate(observation = row_number()) |> 
    pivot_longer(cols = -observation,
                 names_to = "x")
  
  # compute the risk scores for each new patient based on the betas
  
  new_x_risks <- 
    left_join(new_x, 
              betas,
              by = c("x")) |> 
    mutate(coefs = if_else(is.na(coefs), 0, coefs),
           x_x_coefs = value*coefs) |> 
    summarise(risk = sum(x_x_coefs),
              .by = observation)
  
  hole_grid <- 
    tibble(status = new_y[, "status"], 
           time = new_y[, "time"],
           new_x_risks)

  for(i_c in 1:length(new_y)){

    pivot <-
      hole_grid |>
      slice(i_c) |> 
      select(-observation) |> 
      rename_with(~ paste0(.x, "_p"))
    
    compare <-
      hole_grid |>
      slice(-i_c)
    
    if(i_c == 1){
      my_grid <- tibble(pivot, compare)
    }else{
      my_grid <- bind_rows(my_grid,
                           tibble(pivot, compare))
    }

  }
  
  my_grid <- 
    my_grid |>
    mutate(status_c = case_when(status_p == 1 & status == 1 ~ # both observations with the event
                                  if_else((time_p < time & risk_p > risk) | (time_p > time & risk_p < risk), 
                                          "c", "d"),
                                status_p == 1 & status == 0 ~ # one of them censored
                                  if_else(time_p > time, "n",
                                          if_else(risk_p > risk, "c", "d")),
                                status_p == 0 & status == 1 ~ # the other censored
                                  if_else(time_p < time, "n",
                                          if_else(risk_p < risk, "c", "d")),
                                status_p == 0 & status == 0 ~ "n")) # both observations censored
    
  cidx <-
    my_grid |>  
    count(status_c) |>
    bind_rows(tibble(status_c = c("d", "c"), 
                     n = c(0, 0))) |> 
    summarise(n = sum(n),
              .by = status_c) |>
    pivot_wider(values_from = n, names_from = status_c) |>
    mutate(cidx = c/(c+d))|>
    pull(cidx)
  
  return(cidx)
  
}