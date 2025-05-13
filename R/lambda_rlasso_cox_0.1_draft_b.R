# date: 10 may 2025
# author: manuel gonzález aramburu
# purpose: function to determine optimal lambda value within random lasso algorithm applied to cox model
#          with predefined q values equal to q1 = q2 = p/2 and 200 bootstrap samples
#          and 101 values of lambda raging from 0 to 1 by 0.01

lambda_rlasso_cox <- function(dats, varss, time, status, lambda = seq(0, 1, 0.01), B = 200, q1 = floor(p/2), q2 = floor(p/2), lstd = T){

  #--------------------------
  # input parameters
  #--------------------------

  # dats <data frame or tibble>: contains x, time and status for all subjects/ patients
  # varss <vector>: names (char vector) or positions (numeric vector) within dats of the dats' variables that will be the covariates of the model
  # time <char>: name of the variable in dats with survival time of the subjects/ patients
  # status <char>: name of the variable in dats with status of the subjects/ patients - 0 is censored always
  # B <numeric>: number of bootstrap samples that
  # q1 <numeric>: number of selected random variables for rlasso step 1
  # q2 <numeric>: number of selected random variables for rlasso step 2
  # lambda <numeric>: lambda parameter of lasso regresion
  # lstd <boolean>: indicates if variables should be standardized inside the function (should be F if variables have been already standardized)

  #---------------------------------------------
  # additional needed variables
  #---------------------------------------------

  # total number of covariables in the analysis - p

  p <- length(varss)

  # total number of subjects in the analysis - n

  n <- nrow(dats)

  #---------------------------------------------
  # cross validation samples generation
  #---------------------------------------------

  # shuffle the observations in the data base

  samples <- sample(n, n)

  # get the size of the samples

  size <- floor(n/5)

  # create the samples

  s1 <-
    dats |>
    slice(samples[1:size])

  s2 <-
    dats |>
    slice((size+1):(size*2))

  s3 <-
    dats |>
    slice((size*2+1):(size*3))

  s4 <-
    dats |>
    slice((size*3+1):(size*4))

  s5 <-
    dats |>
    slice((size*4+1):(size*5))

  #---------------------------------------------
  # cross validation process execution
  #---------------------------------------------

  # progress bar setup

  print("Beggining execution...")

  n_iter <- length(lambda)

  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [time of execution: :elapsedfull || left estimated time: :eta]",
                         total = n_iter,
                         complete = "=",   # Caracteres de las iteraciones finalizadas
                         incomplete = "-", # Caracteres de las iteraciones no finalizadas
                         current = "+",    # Carácter actual
                         clear = FALSE,    # Si TRUE, borra la barra cuando termine
                         width = 100)      # Ancho de la barra de progreso

  for(i_cv in 1:length(lambda)){

    # print(paste("#######################LAMBDA", i_cv, "#######################"))

    # first fold

    # training set

    cv1 <-
      bind_rows(s1, s2, s3, s4)

    # testing set

    xcv1 <-
      s5 |>
      select(all_of(varss))

    ycv1 <-
      s5 |>
      select(sym(time), sym(status)) |>
      rename(time = sym(time),
             status = sym(status)) |>
      mutate(status = if_else(toupper(status) %in% c("Y", "YES"),
                              1, 0))

    ycv1 <- Surv(ycv1$time, ycv1$status)

    # random lasso computation

    rl1 <- rlasso_cox(cv1, varss, time, status, l = lambda[i_cv], B = B, q1 = q1, q2 = q2, lstd = lstd, bar = F)

    # cindex computation

    cindex1 <- cindex(rl1 |> select(x, coefs), xcv1, ycv1)

    # second fold

    # training set

    cv2 <-
      bind_rows(s1, s2, s3, s5)

    # testing set

    xcv2 <-
      s4 |>
      select(all_of(varss))

    ycv2 <-
      s4 |>
      select(sym(time), sym(status)) |>
      rename(time = sym(time),
             status = sym(status)) |>
      mutate(status = if_else(toupper(status) %in% c("Y", "YES"),
                              1, 0))

    ycv2 <- Surv(ycv2$time, ycv2$status)

    # random lasso computation

    rl2 <- rlasso_cox(cv2, varss, time, status, l = lambda[i_cv], B = B, q1 = q1, q2 = q2, lstd = lstd, bar = F)

    # cindex computation

    cindex2 <- cindex(rl2 |> select(x, coefs), xcv2, ycv2)

    # third fold

    # training set

    cv3 <-
      bind_rows(s1, s2, s4, s5)

    # testing set

    xcv3 <-
      s3 |>
      select(all_of(varss))

    ycv3 <-
      s3 |>
      select(sym(time), sym(status)) |>
      rename(time = sym(time),
             status = sym(status)) |>
      mutate(status = if_else(toupper(status) %in% c("Y", "YES"),
                              1, 0))

    ycv3 <- Surv(ycv3$time, ycv3$status)

    # random lasso computation

    rl3 <- rlasso_cox(cv3, varss, time, status, l = lambda[i_cv], B = B, q1 = q1, q2 = q2, lstd = lstd, bar = F)

    # cindex computation

    cindex3 <- cindex(rl3 |> select(x, coefs), xcv3, ycv3)

    # fourth fold

    # training set

    cv4 <-
      bind_rows(s1, s4, s3, s5)

    # testing set

    xcv4 <-
      s2 |>
      select(all_of(varss))

    ycv4 <-
      s2 |>
      select(sym(time), sym(status)) |>
      rename(time = sym(time),
             status = sym(status)) |>
      mutate(status = if_else(toupper(status) %in% c("Y", "YES"),
                              1, 0))

    ycv4 <- Surv(ycv4$time, ycv4$status)

    # random lasso computation

    rl4 <- rlasso_cox(cv4, varss, time, status, l = lambda[i_cv], B = B, q1 = q1, q2 = q2, lstd = lstd, bar = F)

    # cindex computation

    cindex4 <- cindex(rl4 |> select(x, coefs), xcv4, ycv4)

    # fifth fold

    # training set

    cv5 <-
      bind_rows(s4, s2, s3, s5)

    # testing set

    xcv5 <-
      s1 |>
      select(all_of(varss))

    ycv5 <-
      s1 |>
      select(sym(time), sym(status)) |>
      rename(time = sym(time),
             status = sym(status)) |>
      mutate(status = if_else(toupper(status) %in% c("Y", "YES"),
                              1, 0))

    ycv5 <- Surv(ycv5$time, ycv5$status)

    # random lasso computation

    rl5 <- rlasso_cox(cv5, varss, time, status, l = lambda[i_cv], B = B, q1 = q1, q2 = q2, lstd = lstd, bar = F)

    # cindex computation

    cindex5 <- cindex(rl5 |> select(x, coefs), xcv5, ycv5)

    if(i_cv == 1){
      cindexes <- tibble(lambdas = lambda[i_cv],
                         iteration = 1:5,
                         cindexes = c(cindex1, cindex2, cindex3, cindex4, cindex5))
    }else{
      cindexes <-
        rbind(cindexes,
              tibble(lambdas = lambda[i_cv],
                     iteration = 1:5,
                     cindexes = c(cindex1, cindex2, cindex3, cindex4, cindex5)))

      cindexes_debug <<- cindexes
    }

    pb$tick()

  }

  return(cindexes)

}

# set.seed(1234)
#
# lambdas <- lambda_rlasso_cox(mibase, 7:706, "OS", "Exitus", lstd = F)
