# date: 21 apr 2025
# author: manuel gonzález aramburu
# purpose: function to execute the random lasso method within cox model 
#          with user defined q1, q2 and lambda although there will be default 
#          values for q1 and  q2 as q1 = q2 = p/2 [Liang, L., Zhuang, Y., & Yu, P. L. H. (2024). Variable selection for high-dimensional incomplete data. Computational Statistics and Data Analysis, 192. https://doi.org/10.1016/j.csda.2023.107877]

rlasso_cox <- function(dats, varss, time, status, l, bar = T, B = 200, q1 = floor(p/2), q2 = floor(p/2), lstd = T){
  
  #--------------------------
  # input parameters
  #--------------------------
  
  # dats <data frame or tibble>: contains x, time and status for all subjects/ patients
  # varss <vector>: names (char vector) or positions (numeric vector) within dats of the dats' variables that will be the covariates of the model
  # time <char>: name of the variable in dats with survival time of the subjects/ patients
  # status <char>: name of the variable in dats with status of the subjects/ patients - 0 is censored always
  # B <numeric>: number of bootstrap samples that 
  # bar <boolean>: if F the function does not show the progress bar of the algorithm
  # q1 <numeric>: number of selected random variables for rlasso step 1
  # q2 <numeric>: number of selected random variables for rlasso step 2
  # l <numeric>: lambda parameter of lasso regresion 
  # lstd <boolean>: indicates if variables should be standardized inside the function (should be F if variables have been already standardized)
  
  #---------------------------------------------
  # preparation of the data for glmnet function
  #---------------------------------------------
  
  # dependent variable - y
  
  pre_y <-
    dats |> 
    select(sym(time), sym(status)) |> 
    rename(time = sym(time), 
           status = sym(status)) |> 
    mutate(status = if_else(toupper(status) %in% c("Y", "YES"), 
                            1, 0))
  
  y <- Surv(pre_y$time, pre_y$status)
  
  # independent variables - x
  
  x <-
    dats |> 
    select(all_of(varss)) |> 
    mutate_all(as.numeric) |> 
    as.matrix()
  
  #---------------------------------------------
  # additional needed variables
  #---------------------------------------------
  
  # total number of covariables in the analysis - p
  
  p <- length(varss)
  
  # total number of subjects in the analysis - n
  
  n <- nrow(dats)
  
  #----------------------------------------------
  # rlasso - step 1: importance of the parameters
  #----------------------------------------------
  
  # progress bar setup
  
  pb <- progress_bar$new(format = "STEP 1: [:bar] :percent [time of execution: :elapsedfull || left estimated time: :eta]",
                         total = B,
                         complete = "=",   # Caracteres de las iteraciones finalizadas
                         incomplete = "-", # Caracteres de las iteraciones no finalizadas
                         current = "+",    # Carácter actual
                         clear = FALSE,    # Si TRUE, borra la barra cuando termine
                         width = 100)      # Ancho de la barra de progreso  
  
  # parameter importance metric - I
  
  imp <- 
    tibble(vari = colnames(x)) 
  
  for(b in 1:B){
    
    if(bar){
      pb$tick()
    }
     
    b1i <- sample(1:n, n, replace = T) # bootstrap sample index
    varsi <- sample(1:p, q1, replace = F) # random variable selection index
    
    xi <- x[b1i, varsi] # subset selected variables and individuals
    yi <- y[b1i, ] # subset selected individuals
    
    fiti <- glmnet(xi, yi, 
                    lamda = l,
                    family = "cox",
                    standardize = lstd)
    
    # extract coeficients and append them to compute final importance metric
    
    coefi <- coef(fiti, s = l) 
    
    varindx <- coefi@i + 1 
    vari <- unlist(coefi@Dimnames[[1]])[varindx]
    coefs <- unlist(coefi@x)
    
    impi <- 
      tibble(vari, coefs)
    
    imp <- left_join(imp, impi, by ="vari")
    
    # print(paste("STEP 1 - ", b))
    
  }

  imp_final <- 
    imp |> 
    pivot_longer(cols = -vari) |> 
    mutate(value = if_else(is.na(value), 0, value)) |> 
    summarize(i = abs(sum(value))/B, 
              .by = vari)
  
  # note: glmnet gives coeficients not HRs
  
  #----------------------------------------
  # rlasso - step 2 - final model estimates
  #----------------------------------------
  
  # progress bar setup
  
  pb <- progress_bar$new(format = "STEP 2: [:bar] :percent [time of execution: :elapsedfull || left estimated time: :eta]",
                         total = B,
                         complete = "=",   # Caracteres de las iteraciones finalizadas
                         incomplete = "-", # Caracteres de las iteraciones no finalizadas
                         current = "+",    # Carácter actual
                         clear = FALSE,    # Si TRUE, borra la barra cuando termine
                         width = 100)      # Ancho de la barra de progreso  
  
  
  # parameter estimates
  
  params <- 
    tibble(vari = colnames(x)) 
  
  for(b in 1:B){
    
    if(bar){
      pb$tick()
    }
    
    b1i <- sample(1:n, n, replace = T) # bootstrap sample index
    varsi <- sample(1:p, q2, replace = F) # random variable selection index
    
    varsi <- 
      imp_final |> 
      filter(i != 0) |> 
      slice_sample(n = 350, weight_by = i) |> 
      pull(vari)
    
    xi <- x[b1i, varsi] # subset selected variables and individuals
    yi <- y[b1i, ] # subset selected individuals
    
    fiti <- glmnet(xi, yi, 
                   lamda = l,
                   family = "cox",
                   standardize = lstd)
    
    # extract coeficients and append them to compute final importance metric
    
    coefi <- coef(fiti, s = l) 
    
    varindx <- coefi@i + 1 
    vari <- unlist(coefi@Dimnames[[1]])[varindx]
    coefs <- unlist(coefi@x)
    
    paramsi <- 
      tibble(vari, coefs)
    
    params <- left_join(params, paramsi, by ="vari")
    
    # print(paste("STEP 2 - ", b))
    
  }
  
  params_final <- 
    params |> 
    pivot_longer(cols = -vari) |> 
    mutate(value = if_else(is.na(value), 0, value)) |> 
    summarize(coefs = sum(value)/B, 
              .by = vari) |> 
    mutate(HRs = exp(coefs)) |> 
    rename(x = vari)
  
  return(params_final)
  
}
