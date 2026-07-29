#' @import dplyr
#'
#' @export
#' @title Parameter Summaries of Equivalent Models
#' @description
#' `param_summary()` receives the output of [GEMini::equiv_gen()] or [GEMini::equiv_trim()] and summarizes the permutations
#' selected for each parameter across all equivalent models.
#' @param equiv_list Output of `equiv_gen()` or `equiv_trim()`, containing a list of equivalent models
#' @returns dataframe summarizing permutations selected for each parmeter across all equivalent models. Rownames are set as the original permutation specified.
param_summary <- function(equiv_list) {
  #given a list of equivalent models, creates a summary of parameters

  #extracts all equivalent model names
  all_names <- names(equiv_list$equiv_models)

  #for each name, split into it's individual characters and bind
  split_df <- data.frame()
  for (name in all_names) {
    split_string <- as.data.frame(str_split_fixed(name, "", nchar(all_names[1])))
    split_df <- rbind(split_df, split_string)
  }

  print(split_df)

  #for each column (i.e., specified parameter), count As, Bs, and Cs
  #bind to count_df with new readable_lavaan syntax
  count_df <- data.frame("configs" = c("A", "B", "C"))
  for (column in colnames(split_df)) {
    sum_a <- sum(split_df[[column]] == "A")
    sum_b <- sum(split_df[[column]] == "B")
    sum_c <- sum(split_df[[column]] == "C")

    row_spec <- equiv_list$permut_df[equiv_list$equiv_models[[1]]$ind, ][as.numeric(str_sub(column, 2)),]
    print(row_spec)

    sums_df <- as.data.frame(c(sum_a, sum_b, sum_c))
    colnames(sums_df) <- readable_lavaan(row_spec)

    count_df <- cbind(count_df, sums_df)
  }

  #clean and export
  count_clean <- count_df %>%
    column_to_rownames(var = "configs") %>%
    t() %>%
    as.data.frame()
  return(count_clean)
}

#' @export
#' @title Trimming Equivalent Models
#' @description
#' Given a list of equivalent models (i.e., the output of [GEMini::equiv_gen()]), alongside a vector specifying the "priority" for each variable
#' (such that variables with a higher priority cannot cause variables with a lower priority), `equiv_trim()` "trims" the number of equivalent models
#' to retain only those that satisfy the constraints set by the priority vector
#' @param equiv_list Output of `equiv_gen()` or `equiv_trim()`, containing a list of equivalent models
#' @param prio_vec named numeric vector indicating the priority value for each variable involved in the model.
#' @returns A truncated list of equivalent models, in the same format as the initial equiv_list argument supplied.
equiv_trim <- function(equiv_list, prio_vec) {
  #given an equiv_list, trims set of equivalent models based on priority vector and returns a copy of equiv_list (only including the trimmed models)

  #creates new column indicating whether the specified permutation is still "valid"
  trim_permut <- equiv_list$permut_df %>%
    mutate(trim_valid = ifelse(op == "~~", TRUE, prio_vec[lhs] >= prio_vec[rhs]))

  #for each implied model, creates the subset of this new trim_permut and if any of the values are false, tosses the model
  #otherwise, includes in
  trimmed_models <- list()
  for (model in equiv_list$equiv_models) {

    # print(model_numtostring(model$ind))
    permut_filter <- trim_permut[model$ind, ]

    if (sum(permut_filter$trim_valid) == nrow(permut_filter)) {
      trimmed_models[[model_numtostring(model$ind)]] <- model
    }
  }
  # View(trimmed_models)

  #create new copy (trimmed_list) and replace equiv_models with the new trimmed models
  trimmed_list <- equiv_list
  trimmed_list$equiv_models <- trimmed_models

  return(trimmed_list)
}
