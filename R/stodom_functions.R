################################################################################
# R functions of the stodom R package                                          #
# -----------------------------------------------------------------------------#
# 2022, Agroscope                                                              #
# developed by Sergei Schaub                                                   #
# e-mail: sergei.schaub@agroscope.admin.ch                                     #
# first version: August, 15, 2022                                              #
# last update: August, 15, 2022                                                #
################################################################################

#===============================================================================
# setting
#===============================================================================
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' @import pracma
#' @import tidyr

utils::globalVariables(c("x_axis", "n_var_a", "n_var_b", "pro_var_a", "pro_var_b",
                         "cum_var_a", "cum_var_b", "txtProgressBar", "setTxtProgressBar",
                         "statistic_observed", "larger", "dif_cum", "int_dif_cum",
                         "cum", "group", "cum_var_a_non_random"))


#===============================================================================
# first-order stochastic dominance
#===============================================================================
#' @title first-order stochastic dominance test
#'
#' @description This function tests for first-order stochastic dominance.
#' @usage fo_stodom(data_1, data_2, bins_size, n_draws, useed, variable_1, variable_2, type)
#' @param data_1 data 1.
#' @param data_2 data 2.
#' @param bins_size bin size.
#' @param n_draws number of draws to compute p values (default = 500).
#' @param useed user defined seed
#' @param variable_1 name of a (as a string); only for the output table (default = "a").
#' @param variable_2 name of b (as a string); only for the output table (default = "b").
#' @param type type of bootstrapped test, bootstrapping 1 and 2 of Barrett and Donald (2003) are available (default = "boot2").
#' @details This function computes the consistent test of first-order stochastic dominance following Barrett and Donald (2003). In detail, this function estimate their Kolmogorov-Smirnov type tests based on bootstrapping 2. The function was implemented as part of Schaub xxx
#' @return The function returns a list object containing the p-values of two dominance tests (i.e., variable 1 vs. variable 1 and variable 2 vs. variable 1).
#' @references Barrett, G. F., & Donald, S. G. (2003). Consistent tests for stochastic dominance. Econometrica, 71(1), 71-104.
#' @references Schaub, S. & El Benni, N. (2024). How do price (risk) changes influence farmers’ preference to reduce fertilizer application?
#' @examples
#'
#' # load stodom
#' require(stodom)
#'
#'  data_a <- rnorm(500, 3, 2)
#'  data_b <- rnorm(500, 1, 2)
#'
#' # estimate first-order stochastic dominance
#' fo_stodom(data_1 = data_a, data_2 = data_b, n_draws = 100, useed = 1, bins_size = 1)
#' @export


fo_stodom <- function(data_1, data_2, bins_size = 1, n_draws = 500, useed, variable_1 = "a", variable_2 = "b", type = "boot2") {

  bin_scaling = 1/bins_size # define scaling parameter

  object <- tibble(
    "FO" = "",
    pvalue = rep(NA, 2))

  object[1,1] <- paste0("FO_", variable_1, "_vs_", variable_2)
  object[2,1] <- paste0("FO_", variable_2, "_vs_", variable_1)

  # data_aux <- data_1_2
  n_draws <- n_draws

  for (z in 1:2) {

    # adjust data
    if (z == 1) {


      sample_aux_1 <- data_1
      sample_aux_2 <- data_2}

    if (z == 2) {
      sample_aux_1 <- data_2
      sample_aux_2 <- data_1}


    #-------------------------------------------------------------------------------
    # computing test statistic of observed distributions
    #-------------------------------------------------------------------------------

    min_ref <- floor(min(c(sample_aux_1, sample_aux_2))) - bins_size
    max_ref <- ceiling(max(c(sample_aux_1, sample_aux_2))) + bins_size

    data_dif <- data.frame(x_axis = seq(min_ref, max_ref, bins_size))  %>%
      mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
      distinct() %>%
      arrange(x_axis) %>%
      mutate(x_axis = as.character(x_axis)) %>%
      left_join(
        sample_aux_1 %>% as.data.frame() %>% rename(x_axis = 1) %>%
          mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
          group_by(x_axis) %>%
          summarize(n_var_a = n()) %>%
          arrange(x_axis) %>%
          mutate(pro_var_a = n_var_a/length(sample_aux_1),
                 x_axis = as.character(x_axis)),
        by = "x_axis") %>%
      left_join(
        sample_aux_2 %>% as.data.frame() %>% rename(x_axis = 1) %>%
          mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
          group_by(x_axis) %>%
          summarize(n_var_b = n()) %>%
          arrange(x_axis) %>%
          mutate(pro_var_b = n_var_b/length(sample_aux_2),
                 x_axis = as.character(x_axis)),
        by = "x_axis") %>%
      mutate(pro_var_a = ifelse(!is.na(pro_var_a), pro_var_a, 0),
             pro_var_b = ifelse(!is.na(pro_var_b), pro_var_b, 0),
             cum_var_a = cumsum(pro_var_a),
             cum_var_b = cumsum(pro_var_b),
             x_axis = as.numeric(x_axis)) %>%
      dplyr::select(x_axis, cum_var_a, cum_var_b) %>%
      mutate(dif_cum = cum_var_a - cum_var_b)



    # get maximum difference between the two ecdfs
    sub_dif <- max(data_dif$dif_cum)


    # compute test statistic
    N <- length(sample_aux_1)
    M <- length(sample_aux_2)
    statistic <- (N * M / (N + M))^0.5 * sub_dif


    #-------------------------------------------------------------------------------
    # computing test statistic of randomly sampled data
    #-------------------------------------------------------------------------------


    #////////
    # bootstrapping 2
    #////////

    if (type == "boot2") {


      stats_boots <- data.frame(
        draw = rep(NA, n_draws),
        statistic_bootstrapped = NA)


      dat_combo <- c(sample_aux_1, sample_aux_2) # combine samples

      # progress bar
      pb = txtProgressBar(min = 0, max = n_draws, initial = 0)

      for (i in 1:n_draws) {
        if (!is.na(useed)) {
          set.seed(useed + i) # set seed to make results reproducible
        }


        # randomly draw samples
        sample_a_boot <- sample(dat_combo, size = N, replace = T)
        sample_b_boot <- sample(dat_combo, size = M, replace = T)


        data_dif_boot <- data.frame(x_axis = seq(min_ref, max_ref, bins_size))  %>%
          mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
          distinct() %>%
          arrange(x_axis) %>%
          mutate(x_axis = as.character(x_axis)) %>%
          left_join(
            sample_a_boot %>% as.data.frame() %>% rename(x_axis = 1) %>%
              mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
              group_by(x_axis) %>%
              summarize(n_var_a = n()) %>%
              arrange(x_axis) %>%
              mutate(pro_var_a = n_var_a/length(sample_a_boot),
                     x_axis = as.character(x_axis)),
            by = "x_axis") %>%
          left_join(
            sample_b_boot %>% as.data.frame() %>% rename(x_axis = 1) %>%
              mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
              group_by(x_axis) %>%
              summarize(n_var_b = n()) %>%
              arrange(x_axis) %>%
              mutate(pro_var_b = n_var_b/length(sample_b_boot),
                     x_axis = as.character(x_axis)),
            by = "x_axis") %>%
          mutate(pro_var_a = ifelse(!is.na(pro_var_a), pro_var_a, 0),
                 pro_var_b = ifelse(!is.na(pro_var_b), pro_var_b, 0),
                 cum_var_a = cumsum(pro_var_a),
                 cum_var_b = cumsum(pro_var_b),
                 x_axis = as.numeric(x_axis)) %>%
          dplyr::select(x_axis, cum_var_a, cum_var_b) %>%
          mutate(dif_cum = cum_var_a - cum_var_b)




        sub_dif_sample <- max(data_dif_boot$dif_cum)

        statistic_bootstrapped <- (N * M / (N + M))^0.5 * sub_dif_sample # compute test statistic

        # store results
        stats_boots[i,1] <- i
        stats_boots[i,2] <- statistic_bootstrapped

        setTxtProgressBar(pb, i)
        close(pb)
      }

      stats_boots$statistic_observed <- statistic
      stats_boots <- stats_boots %>%
        mutate(larger = ifelse(statistic_bootstrapped > statistic_observed, 1, 0))

      sum_larger <- stats_boots %>% summarize(sum(larger)) %>% as.numeric()
      p_value <- 1/n_draws*sum_larger
      p_value

    }


    #////////
    # bootstrapping 1
    #////////

    if (type == "boot1") {

      stats_boots <- data.frame(
        draw = rep(NA, n_draws),
        statistic_bootstrapped = NA)


      dat_solo <- c(sample_aux_1) # pass dataset on

      # progress bar
      pb = txtProgressBar(min = 0, max = n_draws, initial = 0)

      for (i in 1:n_draws) {
        if (!is.na(useed)) {
          set.seed(useed + i) # set seed to make results reproducible
        }


        # randomly draw samples
        sample_a_boot <- sample(dat_solo, size = N, replace = T)


        data_dif_boot <- data.frame(x_axis = seq(min_ref, max_ref, bins_size/bin_scaling))  %>%
          distinct() %>%
          arrange(x_axis) %>%
          mutate(x_axis = as.character(x_axis)) %>%
          left_join(
            sample_a_boot %>% as.data.frame() %>% rename(x_axis = 1) %>%
              mutate(x_axis = round(x_axis,  digits = 0)) %>%
              group_by(x_axis) %>%
              summarize(n_var_a = n()) %>%
              arrange(x_axis) %>%
              mutate(pro_var_a = n_var_a/length(sample_a_boot),
                     x_axis = as.character(x_axis)),
            by = "x_axis") %>%
          left_join(
            data_dif %>% dplyr::select(x_axis, cum_var_a) %>%
              rename(cum_var_a_non_random = cum_var_a) %>%
              mutate(x_axis = round(x_axis,  digits = 0)) %>%
              arrange(x_axis) %>%
              mutate(x_axis = as.character(x_axis)),
            by = "x_axis") %>%
          mutate(pro_var_a = ifelse(!is.na(pro_var_a), pro_var_a, 0),
                 cum_var_a = cumsum(pro_var_a),
                 x_axis = as.numeric(x_axis)) %>%
          dplyr::select(x_axis, cum_var_a, cum_var_a_non_random) %>%
          mutate(dif_cum = cum_var_a - cum_var_a_non_random)


        sub_dif_sample <- max(data_dif_boot$dif_cum)

        statistic_bootstrapped <- N^0.5 * sub_dif_sample # compute test statistic

        # store results
        stats_boots[i,1] <- i
        stats_boots[i,2] <- statistic_bootstrapped

        setTxtProgressBar(pb, i)
        close(pb)
      }

      stats_boots$statistic_observed <- statistic
      stats_boots <- stats_boots %>%
        mutate(larger = ifelse(statistic_bootstrapped > statistic_observed, 1, 0))

      sum_larger <- stats_boots %>% summarize(sum(larger)) %>% as.numeric()
      p_value <- 1/n_draws*sum_larger
      p_value
    }

    #-------------------------------------------------------------------------------
    # store result
    #-------------------------------------------------------------------------------
    object[z,2] <- p_value

  }

  oject_list <- list(object, "H0: cannot reject dominance. pvalue < a, we reject that variable_1 can dominant variable_2 at the level a. See section 6 of Barrett and Donald (2003) for interpretation")

  return(oject_list)
}



#===============================================================================
# second-order stochastic dominance
#===============================================================================
#' @title second-order stochastic dominance test
#'
#' @description This function tests for second-order stochastic dominance.
#' @usage so_stodom(data_1, data_2, bins_size, n_draws, useed, variable_1, variable_2, type)
#' @param data_1 data 1.
#' @param data_2 data 2.
#' @param bins_size bin size.
#' @param n_draws number of draws to compute p values (default = 500).
#' @param useed user defined seed
#' @param variable_1 name of a (as a string); only for the output table (default = "a").
#' @param variable_2 name of b (as a string); only for the output table (default = "b").
#' @param type type of bootstrapped test, bootstrapping 1 and 2 of Barrett and Donald (2003) are available (default = "boot2").
#' @details This function computes the consistent test of second-order stochastic dominance following Barrett and Donald (2003). In detail, this function estimate their Kolmogorov-Smirnov type tests based on bootstrapping 2. The function was implemented as part of Schaub xxx
#' @return The function returns a list object containing the p-values of two dominance tests (i.e., variable 1 vs. variable 1 and variable 2 vs. variable 1).
#' @references Barrett, G. F., & Donald, S. G. (2003). Consistent tests for stochastic dominance. Econometrica, 71(1), 71-104.
#' @references Schaub, S. & El Benni, N. (2024). How do price (risk) changes influence farmers’ preference to reduce fertilizer application?
#' @examples
#'
#' # load stodom
#' require(stodom)
#'
#'  data_a <- rnorm(500, 3, 2)
#'  data_b <- rnorm(500, 1, 2)
#'
#' # estimate second-order stochastic dominance
#' so_stodom(data_1 = data_a, data_2 = data_b, n_draws = 100, useed = 1, bins_size = 1)
#' @export


so_stodom <- function(data_1, data_2 , bins_size = 1, n_draws = 500, useed, variable_1 = "a", variable_2 = "b", type = "boot2") {

  bin_scaling = 1/bins_size # define scaling parameter


  object <- tibble(
    "SO" = "",
    pvalue = rep(NA, 2))

  object[1,1] <- paste0("SO_", variable_1, "_vs_", variable_2)
  object[2,1] <- paste0("SO_", variable_2, "_vs_", variable_1)

  # data_aux <- data_1_2
  n_draws <- n_draws

  for (z in 1:2) {

    # adjust data
    if (z == 1) {

      sample_aux_1 <- data_1
      sample_aux_2 <- data_2}

    if (z == 2) {
      sample_aux_1 <- data_2
      sample_aux_2 <- data_1}

    #-------------------------------------------------------------------------------
    # computing test statistic of observed distributions
    #-------------------------------------------------------------------------------

    min_ref <- floor(min(c(sample_aux_1, sample_aux_2))) - bins_size
    max_ref <- ceiling(max(c(sample_aux_1, sample_aux_2))) + bins_size

    data_dif <- data.frame(x_axis = seq(min_ref, max_ref, bins_size))  %>%
      mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
      distinct() %>%
      arrange(x_axis) %>%
      mutate(x_axis = as.character(x_axis)) %>%
      left_join(
        sample_aux_1 %>% as.data.frame() %>% rename(x_axis = 1) %>%
          mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
          group_by(x_axis) %>%
          summarize(n_var_a = n()) %>%
          arrange(x_axis) %>%
          mutate(pro_var_a = n_var_a/length(sample_aux_1),
                 x_axis = as.character(x_axis)),
        by = "x_axis") %>%
      left_join(
        sample_aux_2 %>% as.data.frame() %>% rename(x_axis = 1) %>%
          mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
          group_by(x_axis) %>%
          summarize(n_var_b = n()) %>%
          arrange(x_axis) %>%
          mutate(pro_var_b = n_var_b/length(sample_aux_2),
                 x_axis = as.character(x_axis)),
        by = "x_axis") %>%
      mutate(pro_var_a = ifelse(!is.na(pro_var_a), pro_var_a, 0),
             pro_var_b = ifelse(!is.na(pro_var_b), pro_var_b, 0),
             cum_var_a = cumsum(pro_var_a),
             cum_var_b = cumsum(pro_var_b),
             x_axis = as.numeric(x_axis)) %>%
      dplyr::select(x_axis, cum_var_a, cum_var_b) %>%
      mutate(dif_cum = cum_var_a - cum_var_b)


    data_dif <- data_dif %>%
      mutate(int_dif_cum = cumsum(dif_cum))


    # get maximum difference between the two ecdfs
    sub_dif <- max(data_dif$int_dif_cum)


    # compute test statistic
    N <- length(sample_aux_1)
    M <- length(sample_aux_2)
    statistic <- (N * M / (N + M))^0.5 * sub_dif
    statistic


    #-------------------------------------------------------------------------------
    # computing test statistic of randomly sampled data
    #-------------------------------------------------------------------------------

    n_draws = n_draws

    #////////
    # bootstrapping 2
    #////////

    if (type == "boot2") {

      stats_boots <- data.frame(
        draw = rep(NA, n_draws),
        statistic_bootstrapped = NA)


      dat_combo <- c(sample_aux_1, sample_aux_2) # combine samples

      # progress bar
      pb = txtProgressBar(min = 0, max = n_draws, initial = 0)

      for (i in 1:n_draws) {

        set.seed(1 + i) # set seed to make results reproducible

        # randomly draw samples
        sample_a_boot <- sample(dat_combo, size = N, replace = T)
        sample_b_boot <- sample(dat_combo, size = M, replace = T)

        data_dif_boot <- data.frame(x_axis = seq(min_ref, max_ref, bins_size))  %>%
          mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
          distinct() %>%
          arrange(x_axis) %>%
          mutate(x_axis = as.character(x_axis)) %>%
          left_join(
            sample_a_boot %>% as.data.frame() %>% rename(x_axis = 1) %>%
              mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
              group_by(x_axis) %>%
              summarize(n_var_a = n()) %>%
              arrange(x_axis) %>%
              mutate(pro_var_a = n_var_a/length(sample_a_boot),
                     x_axis = as.character(x_axis)),
            by = "x_axis") %>%
          left_join(
            sample_b_boot %>% as.data.frame() %>% rename(x_axis = 1) %>%
              mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
              group_by(x_axis) %>%
              summarize(n_var_b = n()) %>%
              arrange(x_axis) %>%
              mutate(pro_var_b = n_var_b/length(sample_b_boot),
                     x_axis = as.character(x_axis)),
            by = "x_axis") %>%
          mutate(pro_var_a = ifelse(!is.na(pro_var_a), pro_var_a, 0),
                 pro_var_b = ifelse(!is.na(pro_var_b), pro_var_b, 0),
                 cum_var_a = cumsum(pro_var_a),
                 cum_var_b = cumsum(pro_var_b),
                 x_axis = as.numeric(x_axis)) %>%
          dplyr::select(x_axis, cum_var_a, cum_var_b) %>%
          mutate(dif_cum = cum_var_a - cum_var_b)

        data_dif_boot <- data_dif_boot %>%
          mutate(int_dif_cum = cumsum(dif_cum))


        # get maximum difference between the two ecdfs
        sub_dif_sample <- max(data_dif_boot$int_dif_cum)

        statistic_bootstrapped <- (N * M / (N + M))^0.5 * sub_dif_sample # compute test statistic



        # store results
        stats_boots[i,1] <- i
        stats_boots[i,2] <- statistic_bootstrapped

        setTxtProgressBar(pb, i)
        if (z == 1) {close(pb)}

      }


      stats_boots$statistic_observed <- statistic
      stats_boots <- stats_boots %>%
        mutate(larger = ifelse(statistic_bootstrapped > statistic_observed, 1, 0))

      sum_larger <- stats_boots %>% summarize(sum(larger)) %>% as.numeric()
      p_value <- 1/n_draws*sum_larger
      p_value
    }

    #////////
    # bootstrapping 1
    #////////

    if (type == "boot1") {

      stats_boots <- data.frame(
        draw = rep(NA, n_draws),
        statistic_bootstrapped = NA)


      dat_solo <- c(sample_aux_1) # pass dataset on

      # progress bar
      pb = txtProgressBar(min = 0, max = n_draws, initial = 0)

      for (i in 1:n_draws) {
        if (!is.na(useed)) {
          set.seed(useed + i) # set seed to make results reproducible
        }


        # randomly draw samples
        sample_a_boot <- sample(dat_solo, size = N, replace = T)


        data_dif_boot <- data.frame(x_axis = seq(min_ref, max_ref, bins_size/bin_scaling))  %>%
          distinct() %>%
          arrange(x_axis) %>%
          mutate(x_axis = as.character(x_axis)) %>%
          left_join(
            sample_a_boot %>% as.data.frame() %>% rename(x_axis = 1) %>%
              mutate(x_axis = round(x_axis,  digits = 0)) %>%
              group_by(x_axis) %>%
              summarize(n_var_a = n()) %>%
              arrange(x_axis) %>%
              mutate(pro_var_a = n_var_a/length(sample_a_boot),
                     x_axis = as.character(x_axis)),
            by = "x_axis") %>%
          left_join(
            data_dif %>% dplyr::select(x_axis, cum_var_a) %>%
              rename(cum_var_a_non_random = cum_var_a) %>%
              mutate(x_axis = round(x_axis,  digits = 0)) %>%
              arrange(x_axis) %>%
              mutate(x_axis = as.character(x_axis)),
            by = "x_axis") %>%
          mutate(pro_var_a = ifelse(!is.na(pro_var_a), pro_var_a, 0),
                 cum_var_a = cumsum(pro_var_a),
                 x_axis = as.numeric(x_axis)) %>%
          dplyr::select(x_axis, cum_var_a, cum_var_a_non_random) %>%
          mutate(dif_cum = cum_var_a - cum_var_a_non_random)


        data_dif_boot <- data_dif_boot %>%
          mutate(int_dif_cum = cumsum(dif_cum))


        # get maximum difference between the two ecdfs
        sub_dif_sample <- max(data_dif_boot$int_dif_cum)

        statistic_bootstrapped <- N^0.5 * sub_dif_sample # compute test statistic

        # store results
        stats_boots[i,1] <- i
        stats_boots[i,2] <- statistic_bootstrapped

        setTxtProgressBar(pb, i)
        close(pb)
      }

      stats_boots$statistic_observed <- statistic
      stats_boots <- stats_boots %>%
        mutate(larger = ifelse(statistic_bootstrapped > statistic_observed, 1, 0))

      sum_larger <- stats_boots %>% summarize(sum(larger)) %>% as.numeric()
      p_value <- 1/n_draws*sum_larger
      p_value
    }



    #-------------------------------------------------------------------------------
    # store result
    #-------------------------------------------------------------------------------
    object[z,2] <- p_value


  }

  oject_list <- list(object, "H0: cannot reject dominance. pvalue < a, we reject that variable_1 can dominant variable_2 at the level a. See section 6 of Barrett and Donald (2003) for interpretation")

  return(oject_list)
}





#===============================================================================
# get ecdf values of two variables
#===============================================================================
#' @title values of two ecdf and their cumulative difference
#'
#' @description This function computes the values of two empirical cumulative distribution function as well as their cumulative differences.
#' @usage ecdf_dat_g(data_1, data_2, bins_size)
#' @param data_1 data 1.
#' @param data_2 data 2.
#' @param bins_size bin size.
#' @details This function computes the values of two empirical cumulative distribution function as well as their cumulative differences.
#' @return The function returns a data table.
#' @examples
#'
#' # load stodom
#' require(stodom)
#'
#'  data_a <- rnorm(500, 3, 2)
#'  data_b <- rnorm(500, 1, 2)
#'
#' # compute the values of two ecdfs and their cumulative differences.
#' ecdf_dat_g(data_1 = data_a, data_2 = data_b, bins_size = 1)
#' @export


ecdf_dat_g <- function(data_1, data_2 , bins_size = 1) {

  bin_scaling = 1/bins_size # define scaling parameter

  sample_aux_1 <- data_1
  sample_aux_2 <- data_2

  min_ref <- floor(min(c(sample_aux_1, sample_aux_2))) - bins_size
  max_ref <- ceiling(max(c(sample_aux_1, sample_aux_2))) + bins_size

  data_dif <- data.frame(x_axis = seq(min_ref, max_ref, bins_size))  %>%
    mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
    distinct() %>%
    arrange(x_axis) %>%
    mutate(x_axis = as.character(x_axis)) %>%
    left_join(
      sample_aux_1 %>% as.data.frame() %>% rename(x_axis = 1) %>%
        mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
        group_by(x_axis) %>%
        summarize(n_var_a = n()) %>%
        arrange(x_axis) %>%
        mutate(pro_var_a = n_var_a/length(sample_aux_1),
               x_axis = as.character(x_axis)),
      by = "x_axis") %>%
    left_join(
      sample_aux_2 %>% as.data.frame() %>% rename(x_axis = 1) %>%
        mutate(x_axis = round(x_axis*bin_scaling,  digits = 0)/bin_scaling) %>%
        group_by(x_axis) %>%
        summarize(n_var_b = n()) %>%
        arrange(x_axis) %>%
        mutate(pro_var_b = n_var_b/length(sample_aux_2),
               x_axis = as.character(x_axis)),
      by = "x_axis") %>%
    mutate(pro_var_a = ifelse(!is.na(pro_var_a), pro_var_a, 0),
           pro_var_b = ifelse(!is.na(pro_var_b), pro_var_b, 0),
           cum_var_a = cumsum(pro_var_a),
           cum_var_b = cumsum(pro_var_b),
           x_axis = as.numeric(x_axis)) %>%
    dplyr::select(x_axis, cum_var_a, cum_var_b) %>%
    mutate(dif_cum = cum_var_a - cum_var_b)


  # adding this
  ## this is added as we know the result and rounding and putting variables in to bins can cause mall errors
  data_dif <- data_dif %>%
    mutate(int_dif_cum = cumsum(dif_cum))
  data_dif <- data_dif %>% mutate(int_dif_cum_s = int_dif_cum/abs(sum(data_dif$int_dif_cum))) # scale variable


  return(data_dif)
}




#===============================================================================
# plot ecdfs
#===============================================================================
#' @title plot ecdfs
#'
#' @description This function computes the values of two empirical cumulative distribution function and plots the values.
#' @usage ecdf_plot(data_1, data_2, bins_size)
#' @param data_1 data 1.
#' @param data_2 data 2.
#' @param bins_size bin size.
#' @details This function computes the values of two empirical cumulative distribution function and plots the values.
#' @return The function returns a plot as a ggplot2 object.
#' @examples
#'
#' # load stodom
#' require(stodom)
#'
#'  data_a <- rnorm(500, 3, 2)
#'  data_b <- rnorm(500, 1, 2)
#'
#' # plot ecdfs
#' ecdf_plot(data_1 = data_a, data_2 = data_b, bins_size = 0.1)
#' @export



ecdf_plot <- function(data_1, data_2 , bins_size = 1) {

  sample_aux_1 <- data_1
  sample_aux_2 <- data_2

  aux_1 <- ecdf_dat_g(data_1 = sample_aux_1,
                      data_2 = sample_aux_2, bins_size = bins_size)


  # prepare for plotting
  aux_2 <-
    aux_1 %>% dplyr::select(x_axis, cum_var_a, cum_var_b) %>%
    pivot_longer(!c(x_axis), names_to = "group", values_to = "cum")

  theme_set(theme_bw())

  # plot ecdfs
  fig_aux <- aux_2 %>%
    ggplot() +
    geom_hline(yintercept = 0, color = c("#525252")) +
    geom_hline(yintercept = 1, color = c("#525252")) +
    geom_line(aes(x = x_axis, y = cum,
                  linetype = group), linewidth = 0.9) +
    theme(
      axis.title = element_text(size = 13),
      axis.text  = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    scale_x_continuous(name = expression(paste("outcome (x)"))) +
    scale_y_continuous(
      breaks = c(0,1),
      name = expression(paste(widehat("F"),"(x)"[j]))) +
    scale_linetype_manual(values=c("solid", "dashed"),
                          labels = c("a", "b")) +
    guides(linetype = guide_legend(ncol = 2, title = "", order = 1))

  return(fig_aux)
}




#===============================================================================
# plot cumulative difference between two ecdfs
#===============================================================================

#' @title plot difference ecdfs
#'
#' @description This function computes the values of the cumulative difference of two empirical cumulative distribution function and plots the values.
#' @usage dif_ecdf_plot(data_1, data_2, bins_size)
#' @param data_1 data 1.
#' @param data_2 data 2.
#' @param bins_size bin size.
#' @details This function computes the values of the cumulative difference of two empirical cumulative distribution function and plots the values. This relates two showing second-order stochastic dominance.
#' @return The function returns a plot as a ggplot2 object.
#' @examples
#'
#' # load stodom
#' require(stodom)
#'
#'  data_a <- rnorm(500, 3, 2)
#'  data_b <- rnorm(500, 1, 2)
#'
#' # plot cumulative difference between two ecdfs
#' dif_ecdf_plot(data_1 = data_a, data_2 = data_b, bins_size = 0.1)
#' @export



dif_ecdf_plot <- function(data_1, data_2 , bins_size = 1) {

  sample_aux_1 <- data_1
  sample_aux_2 <- data_2

  aux_1 <- ecdf_dat_g(data_1 = sample_aux_1,
                      data_2 = sample_aux_2, bins_size = bins_size)


  theme_set(theme_bw())

  # plot cum difference of ecdfs
  fig_aux <-
    aux_1 %>%
    ggplot() +
    geom_hline(yintercept = 0, color = c("#525252")) +
    geom_line(aes(x = x_axis, y = int_dif_cum ), size = 0.9)  +
    theme(
      axis.title = element_text(size = 13),
      axis.text  = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = c(0.15, 0.85),
      legend.key.size = unit(0.5, "cm"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    scale_x_continuous(name = expression(paste("outcome (x)"))) +
    scale_y_continuous(
      breaks = 0,
      name = expression(paste(Sigma[paste("-",infinity)]^x^"*", " ", widehat("F"),"(x)"[j], " - ", widehat("F"),"(x)"[j], " ", " dx"))) +
    guides(linetype = guide_legend(ncol = 2, title = "period:", order = 1))


  fig_aux


  return(fig_aux)
}

