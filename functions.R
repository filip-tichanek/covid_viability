#' run Function
#'
#' This function loads or runs and saves (if not previously done) a model, table,
#' or any computationally intensive task. It's designed to avoid redundant computations
#' by reusing previously saved results when possible.
#'
#' @param expr An expression representing the computationally intensive task to be executed.
#' @param path A string specifying the file path where the result should be saved or loaded from.
#' @param reuse A logical flag indicating whether to attempt reusing saved results to avoid recomputation.
#' @return The result of evaluating `expr`. If `reuse` is TRUE and a saved result exists, 
#'         that result is returned; otherwise, `expr` is evaluated.
#' @examples
#' # Assuming lm_result.Rds does not exist, this will compute the linear model and save it.
#' run(lm(mpg ~ cyl, data = mtcars), path = "lm_result", reuse = TRUE)

run <- function(expr, path, reuse = TRUE) {
  fit <- NULL
  if (reuse) {
    path <- paste0(path, ".Rds")
    fit <- suppressWarnings(try(readRDS(path), silent = TRUE))
    if (inherits(fit, "try-error")) {
      fit <- NULL
    }
    }
  if (is.null(fit)) {
    fit <- eval(substitute(expr))
    if (reuse && !is.null(path) && nzchar(path)) {
      saveRDS(fit, file = path)
    }
    }
  return(fit)
  }


#' p_dir Function
#'
#' Calculates the proportion of posterior probability behind a specified
#' threshold (argument `dir` set to `>` or `<`) or probability of direction
#' (`dir` = `max`, as defined in Makowski et al., 2019). This is useful for
#' posterior analysis in Bayesian statistics.
#'
#' @param data A data frame or vector containing posterior draws (in columns or as a vector).
#' @param dir A character string specifying the direction for probability calculation:
#'            use `">"` for above threshold, `"<"` for below threshold, or `"max"` for the
#'            probability of direction as defined by Makowski et al., 2019.
#' @param tres A numeric value specifying the threshold for the calculation (e.g., 0 for
#'             no effect, 1 for odds ratio).
#' @return The calculated probability based on the direction and threshold specified.
#' @examples
#' # Assuming 'post_draws' contains posterior draws,
#' p_dir(post_draws, dir = ">", tres = 0)

p_dir <- function(data, dir, tres){
  if(dir == 'max'){
    return(
      max(length(data[data>tres]),length(data[data<tres]))/length(data)
      )
    } else if(dir == '<'){
      return(
        length(data[data<tres])/length(data)
        )
      } else if(dir == '>'){
        return(
          length(data[data>tres])/length(data)
        )
      } else {print('ERROR')}
  }


#' clean Function
#' 
#' Transliterates to ASCII, removes diacritics, punctuation, spaces, and symbols 
#' such as "+"/"-" from a string. This cleaning process makes the text more suitable 
#' for regex mining by lowering the case and removing extraneous characters.
#'
#' @param string_character A character string to be cleaned.
#' @return The cleaned text as a lower-case string without spaces, punctuation, or symbols.
#' @examples
#' clean("Example+ text - with, punctuation!")

clean <- function(string_character) {
  require(stringi)
  require(stringr)
  
  string_character %>%
    stringi::stri_trans_general("Latin-ASCII") %>% 
    stringr::str_replace_all("[[:punct:]/\\-]+", "") %>%  
    stringr::str_replace_all("\\s+", "") %>% 
    stringr::str_replace_all('\\+', "") %>% 
    stringr::str_to_lower()
}


#' Convert Excel Date Columns in a Data Frame
#'
#' Transforms specified columns in a data frame from Excel serial date numbers
#' to R date format. This function allows for multiple columns to be converted
#' simultaneously. The origin is set to account for Excel's leap year bug in 1900.
#' Columns specified are overwritten with the new date values, retaining their original names.
#'
#' @param data A data frame containing one or more columns to be converted.
#' @param dates A character vector with the names of the columns in `data` that need to be converted to date format.
#' @return A data frame with the specified columns converted to dates and overwritten.
#' @examples
#' df <- data.frame(id = 1:3, date1 = c(39338, 39339, 39340), date2 = c(39341, 39342, 39343))
#' converted_df <- excel_date(df, c("date1", "date2"))
#' @importFrom dplyr mutate
#' @importFrom base as.Date
#' @import tidyverse

excel_date <- function(data, dates) {
  require(tidyverse)
  
  data %>%
    dplyr::mutate(across(all_of(dates), ~ base::as.Date(., origin = "1899-12-30")))
}



#' Permutation Test for GLM Coefficients
#'
#' This function performs permutation tests for generalized linear models to generate 
#' empirical p-values for each model coefficient. 
#' It repeatedly permutes the response variable and recalculates the model
#' coefficients to form a distribution of coefficients under the null.
#'
#' @param model A model object returned by glm or MASS::glmn.nb
#' @param nsim The number of permutations to perform, default is 2000.
#' @return A vector of two-sided p-values for each coefficient in the model.
#' @examples
#' # Assuming 'data' is your dataframe with variables 'count', 'x1', and 'x2'
#' and you have already fitted a glm.nb model:
#' set.seed(123)
#' data <- data.frame(
#'   count = rnbinom(100, mu = 10, size = 1),
#'   x1 = rnorm(100),
#'   x2 = rnorm(100)
#' )
#' model <- MASS::glm.nb(count ~ x1 + x2, data = data)
#' p_values <- perm_model(model, nsim = 1000)
#' print(p_values)
#'
perm_model <- function(model, nsim = 2000){
  
  ## extract original data
  data <- model$model
  
  ## object where to save results
  predictor_names <- names(model$model)[-1]
  coef_null_distributions <- data.frame(matrix(ncol = length(predictor_names) + 1, nrow = nsim))
  names(coef_null_distributions) <- c("b_intercept", predictor_names)
  p_values <- vector('double', length(coef_null_distributions))
  
  ## updating models across iterations
  for (i in 1:nsim){
    data_i <- data
    data_i[, 1] <- sample(data[, 1], replace = FALSE)
    
    model_i <- update(model, data = data_i)
    
    coef_null_distributions[i, ] <- coef(model_i)
  }
  
  ## get p-values
  for(i in 1:ncol(coef_null_distributions)){
    ### add noise
    coef_null_distributions[, i] <- coef_null_distributions[, i] + rnorm(
      nsim, 
      0,
      1e-5*sd(coef_null_distributions[, i]))
    
    ### create a vector object
    vec <- coef_null_distributions[, i]
    
    ### calculate two-sided p-value
    p_values[i] <- (1 + min(
      length(vec[vec > model$coef[i]]),
      length(vec[vec < model$coef[i]])
    )) / (1 + (nsim / 2))
  }
  
  return(p_values)
}




#' Clustered Data Sampler
#'
#' The function samples clustered data by randomly selecting unique values 
#' from a specified column and replicates them with replacement. It is useful 
#' for (cluster) bootstrap
#'
#' @param data A data frame containing the data to be sampled.
#' @param id_col A string specifying the column name in `data` that contains 
#' the cluster identifier.
#' @param N An integer specifying the number of bootstrap samples to generate.
#' @param seed An integer used to set the seed for reproducibility.
#' @return A list of data frames, each containing one bootstrap sample of the 
#' clustered data.
#' @examples
#' # Assuming 'dat' is your dataframe with a column 'fam' representing clusters:
#' set.seed(16)
#' dat <- data.frame(
#'   fam = c(rep('a', 5), rep('b', 8), 'c', 'd', 'e', 'f', 'g'),
#'   outcome = c(rnorm(5, 50, 1), rnorm(8, -50, 1), rnorm(5, 0, 1))
#' )
#' new_dat <- clustdat_sampler(dat, id_col = "fam", N = 2, seed = 123)
#' print(new_dat)
#'
clustdat_sampler <- function(data, id_col, N, seed){
  set.seed(seed)
  
  reset <- list()
  
  for (i in 1:N) {
    tmp <- data.frame(id = sample(unique(data[[id_col]]), 
                                  length(unique(data[[id_col]])), 
                                  replace = TRUE))
    
    tmp <- tmp %>% 
      mutate(id_sec = factor(1:nrow(tmp))) %>% 
      left_join(data, 
                by = c("id" = id_col), 
                relationship = "many-to-many")
    
    reset[[i]] <- tmp
  }
  
  return(reset)
}



