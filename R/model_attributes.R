#' @import dplyr
#'
#'
#' @title Translates model indices to string form (e.g., "AAA")
#' @param model_ind numeric vector of model indices (e.g., 1, 4, 7)
model_numtostring <- function(model_ind) {
  model_string <- case_when(
    model_ind %% 3 == 1 ~ "A",
    model_ind %% 3 == 2 ~ "B",
    model_ind %% 3 == 0 ~ "C") %>%
    paste(., collapse = "")

  return(model_string)
}

#' @title Replacing model index with given number
#' @param model_ind numeric vector of model indices (e.g., 1, 4, 7)
#' @param new_num Replacement number to be inserted in the correct location
replace_ind <- function(model_ind, new_num) {
  model_ind[ceiling(new_num/3)] <- new_num
  return(model_ind)
}

#' @title "De"-lavaanify a parameter dataframe
#' @description
#' given a param df, reverses the lavaanify function (see [lavaan::lavaanify()]), transforming back into a string
#' @param param_df data frame detailing specified free parameters (i.e., include lhs, op, and rhs columns)
delavaanify <- function(param_df) {

  unite_df <- param_df %>%
    unite("united", lhs:rhs, sep = " ")
  delav <- "\n"
  for (row in 1:nrow(unite_df)) {
    delav <- paste0(delav, unite_df[row,], sep = "\n")
  }
  return(delav)
}

#' @title Create easily readable strings from param df
#' @description
#' given a lavaan row specification (e.g., A ~ B), transforms into a more easily readable string (e.g., B -> A)
#' @param model single row of a param df (with columns lhs, op, and rhs)
#' @returns easily readable character string
readable_lavaan <- function(model) {

  if (model$op == "~~") {
    return(paste(model$lhs, "<->", model$rhs, sep = " "))
  } else {
    return(paste(model$rhs, "->", model$lhs, sep = " "))
  }
}

#' @title Determine model validity
#' @description
#' Given a permutation dataframe and specific model_info (including antecedent/consequent information), determines if the specified model is "valid"
#' @param permute_df permutation dataframe
#' @param model_info model object, including information about model index, antecendents, and consequents
#' @returns TRUE or FALSE, indicating model validity
model_valid <- function(permute_df, model_info) {

  #generate current model
  cur_model <- permute_df[model_info$ind,]

  #first check all variables connected
  for (var in unique(c(cur_model$lhs, cur_model$rhs))) {
    var_ante_conse <- c(model_info$ante[[var]], model_info$conse[[var]])
    if (length(var_ante_conse) == 0) {
      # print(paste(var, "is not connected!"))
      return(FALSE)
    }
  }

  #next, check if there are covariances
  cov_model <- cur_model %>%
    filter(op == "~~")
  #if there are not, return true
  if (nrow(cov_model) == 0) {
    return(TRUE)
  } else {
    #otherwise, check for exo-endo covars
    for (row in 1:nrow(cov_model)) {
      #generate vars for length of lhs/rhs antes
      lhs_ante_length <- length(model_info$ante[[cov_model[row,]$lhs]])
      rhs_ante_length <- length(model_info$ante[[cov_model[row,]$rhs]])

      #if the lengths are mismatched, return false
      if ((lhs_ante_length > 0 & rhs_ante_length == 0) |
          (lhs_ante_length == 0 & rhs_ante_length > 0)) {
        # print(paste("exo-endo covar between", cov_model[row,]$lhs,
        #             "and", cov_model[row,]$rhs))
        return(FALSE)
      }
    }
  }
  #if passes all of these, return true
  return(TRUE)
}

#' @title Determine model antecedents and consequents
#' @description
#' Given a permutation dataframe and model indices to generate a specific model, loops through all variables to determine their antecedents, consequents, and direct antecedents
#' @param permute_df permutation dataframe
#' @param model_index model index represented as a numeric vector
#' @returns list containing original model index, and the model's direct antecedents, consequents, and antecedents, each their own list, sorted by variable name
ante_conse <- function(permute_df, model_index) {
  #given full permutation df and specific model index, generates named list of model antecedents and consequents (based on regression arrows)

  #extract current model
  cur_model <- permute_df[model_index, ]

  #initialize list for ante, conse, and direct ante
  dirante_list <- list()
  ante_list <- list()
  conse_list <- list()

  #for each variable...
  for (var in unique(c(cur_model$lhs, cur_model$rhs))) {
    #DETERMINE ANTECEDENTS
    #determine direct antecedents
    ante <- cur_model %>%
      filter(op == "~" & lhs == var) %>%
      .$rhs
    #set dirante_list
    dirante_list[[var]] <- ante
    #initialize past antecedent vector
    ante_past <- c()
    #while the past and original (updated) vector are not the same, continue looping
    while (!(setequal(ante, ante_past))) {
      #set past to be the same as current
      ante_past <- ante
      #update current to include 1 more step
      ante <- cur_model %>%
        filter(op == "~" & lhs %in% unique(c(ante_past, var))) %>%
        .$rhs %>%
        unique(.)

    } #will break out of loop when (updated) ante does not add new info
    #aka ante == ante_past
    #if direct ante is empty, then while loop will not run

    #once while loop is broken, set ante_list
    ante_list[[var]] <- ante

    #DETERMINE CONSEQUENTS - same logic
    #determine direct consequents
    conse <- cur_model %>%
      filter(op == "~" & rhs == var) %>%
      .$lhs
    #initialize past consequent vector
    conse_past <- c()

    #while the past and original (updated) vector are not the same, continue looping
    while (!(setequal(conse, conse_past))) {
      #set past to be the same as current
      conse_past <- conse
      #update current to include 1 more step
      conse <- cur_model %>%
        filter(op == "~" & rhs %in% unique(c(conse_past, var))) %>%
        .$lhs %>%
        unique(.)

    } #will break out of loop when (updated) conse does not add new info
    #aka conse == conse_past

    #once while loop is broken, set conse_list
    conse_list[[var]] <- conse
  }

  #returns model information as a list
  return(
    list(
      ind = model_index,
      dirante = dirante_list,
      ante = ante_list,
      conse = conse_list
    )
  )
}
