#' @import dplyr
#' @import lavaan
#'
#' @export
#' @title Generate Equivalent Models
#' @description
#' `equiv_gen()` takes in a *lavaan* specification of the structural portion of an SEM and
#' generates all equivalent models according to the "replacing rule" specified by Lee and Hershberger (1990)
#' @param lav_spec *lavaan* model specification
#' @returns list containing a) a permutation dataframe, b) a list of equivalent models, and c) names of invalid models found during the generation process. For more details, see (REDACTED, in prep.)
equiv_gen <- function(lav_spec) {

  #lavaanifies to extract parameters to be estimated
  free_params <- lavaan::lavaanify(lav_spec) %>%
    filter(free != 0) %>%
    select(lhs:rhs)
  print(free_params)
  #generates permutaton df based on free_params
  permut_df <- permute_gen(free_params)
  print(permut_df)
  #generates og indexing
  og_ind <- seq(1, nrow(permut_df), 3)

  #initialize total model list
  equiv_models <- list()
  #initialize invalid models vector
  invalid_models <- c()
  #initializes stack
  model_stack <- faststack(init = 1, missing_default = NULL)

  #store og model info into equiv_models and stack
  equiv_models[[model_numtostring(og_ind)]] <- ante_conse(permut_df, og_ind)
  model_stack$push(equiv_models[[1]])

  #while the stack isn't empty...
  while (!is.null(model_stack$peek())) {

    #extract current information
    cur_info <- model_stack$pop()
    #and generate all considered models given the model info
    all_considered_models <- step_gen(permut_df, cur_info)

    #for all considered models
    for (model in names(all_considered_models)) {

      #if model already exists in equiv_models or invalid_models, skip
      if (model %in% names(equiv_models) | model %in% invalid_models) {
        next
      } else {
        #otherwise
        # print(model)
        #determine ante_conse
        model_anteconse <- ante_conse(permut_df, all_considered_models[[model]])
        #determine validity
        model_validity <- model_valid(permut_df, model_anteconse)

        #if model is valid, include model in equiv_models and add to stack
        if (model_validity) {
          equiv_models[[model_numtostring(model_anteconse$ind)]] <- model_anteconse
          model_stack$push(model_anteconse)
        } else {
          #if model is not valid, include model name in invalid_models
          invalid_models <- c(invalid_models, model_numtostring(model_anteconse$ind))
        }
      }
    }
  }

  #binds permut_df, equivalent models, and invalid models into a single list
  equiv_list <- list(permut_df, equiv_models, invalid_models)
  names(equiv_list) <- c("permut_df", "equiv_models", "invalid_models")

  return(equiv_list)
}

#' @title Create permutation dataframe
#' @description
#' Given a dataframe of free (i.e., not fixed) parameters to be estimated, generates a dataframe containing all permutations (reciprocal relationships excluded) of all parameters
#' @param param_df dataframe of parameters to be estimated, a truncated output of [lavaan::lavaanify()]
#' @returns dataframe 3 times as long (# of rows) as the input dataframe, listing all permutations of all parameters
permute_gen <- function(param_df) {
  #given free_params df, generates all permutations of estimated relations

  permut_df <- data.frame(lhs = as.character(),
                          op = as.character(),
                          rhs = as.character())
  for (row in 1:nrow(param_df)) {
    #set current row
    current_row <- param_df[row,]
    #append current row
    permut_df <- bind_rows(permut_df, current_row)

    #if current row specifies a regression...
    if (current_row$op == "~") {

      #generate flipped regression and input
      flipped_reg <- current_row %>%
        rename("rhs" = "lhs",
               "lhs" = "rhs")
      permut_df <- bind_rows(permut_df, flipped_reg)

      #generate covariance form and input
      new_covar <- current_row %>%
        mutate(op = "~~")
      permut_df <- bind_rows(permut_df, new_covar)
    }

    #if current row specifies a covariance...
    if (current_row$op == "~~") {

      #generate first regression and input
      first_reg <- current_row %>%
        mutate(op = "~")
      permut_df <- bind_rows(permut_df, first_reg)

      #generate second regression and input
      second_reg <- current_row %>%
        mutate(op = "~") %>%
        rename("rhs" = "lhs",
               "lhs" = "rhs")
      permut_df <- bind_rows(permut_df, second_reg)
    }
  }

  return(permut_df)
}

#' @title Determine limited block-recursiveness
#' @description
#' Given a model, a specific row in the model, and that model's information (i.e., regarding antecedents and consequents), determines whether the chosen focal block (defined by the specific row) obeys limited block-recursiveness
#' @param model dataframe specifying overall model (i.e., all estimated parameters)
#' @param row dataframe of 1 row length specifying focal block under investigation
#' @param model_info model info indicating information regarding variable antecedents and consequents
#' @returns If focal block fails to pass limited block-recursiveness, returns FALSE. If focal block passes limited block-recursiveness and the generated preceding block is just-identified, outputs model names of all generated names from [GEMini::pbl_jid_gen()]. Otherwise, returns TRUE.
determine_LBR <- function(model, row, model_info) {

  #determine pbl
  if (row$op == "~") {
    #if the param is a regression, then take the (unique) combined antes for both lhs and rhs, minus the direct effect from rhs -> lhs
    pbl <- union(
      model_info$ante[[row$rhs]],
      model_info$ante[[row$lhs]] %>%
        .[-which(. == row$rhs)]
    )
  } else {
    #otherwise, just take the (unique) combined antes with no modifications
    pbl <- union(
      model_info$ante[[row$rhs]],
      model_info$ante[[row$lhs]]
    )
  }
  # print("PBL:")
  # print(pbl)
  #if either the current rhs or lhs exists in PBL, toss
  if ((row$rhs %in% pbl) | (row$lhs %in% pbl)) {
    # print("PBL - FBL not block recursive!")
    return(FALSE)
  }

  #determine sbl
  if (row$op == "~") {
    #if the param is a regression, then take the (unique) combined conse for both lhs and rhs, minus the direct effect from rhs -> lhs
    sbl <- union(
      model_info$conse[[row$lhs]],
      model_info$conse[[row$rhs]] %>%
        .[-which(. == row$lhs)]
    )
  } else {
    #otherwise, just take the (unique) combined conse with no modifications
    sbl <- union(
      model_info$conse[[row$lhs]],
      model_info$conse[[row$rhs]]
    )
  }
  # print("SBL:")
  # print(sbl)
  #if either the current rhs or lhs exists in SBL, toss
  if ((row$rhs %in% sbl) | (row$lhs %in% sbl)) {
    # print("SBL - FBL not block recursive!")
    return(FALSE)
  }

  #if their exists any overlap between pbl and sbl, toss
  if (length(symdiff(pbl, sbl)) != (length(pbl) + length(sbl))) {
    # print("SBL - PBL not block recursive!")
    return(FALSE)
  }

  #if there exists any variable that is not included in pbl, fbl, or, sbl, toss
  if (!(all(union(model$rhs, model$lhs) %in% c(row$lhs, row$rhs, pbl, sbl)))) {
    # print("not all variables in PBL/FBL/SBL!")
    return(FALSE)
  }

  #if there are covariances specified:
  cov_model <- model %>%
    filter(op == "~~")
  if (nrow(cov_model) != 0) {
    #for each covariance, check if both lhs and rhs are either
    #jointly in PBL, jointly in FBL, or jointly in SBL
    #if all of these fail, the covariance must be between blocks, therefore toss
    for (cov_row in 1:nrow(cov_model)) {
      row_vars <- c(cov_model[cov_row, ]$lhs, cov_model[cov_row, ]$rhs)
      if (!(all(row_vars %in% pbl)) &
          !(all(row_vars %in% sbl)) &
          !(all(row_vars %in% c(row$lhs, row$rhs)))) {
        # print("Covariance between blocks!")
        return(FALSE)
      }
    }
  }
  # print("param passes limited-block recur!")


  #if pbl has at least 3 vars
  if (length(pbl) > 2) {
    #if estimated params involving only pbl vars show it is JID
    pbl_df <- model %>%
      mutate(pbl_param = (rhs %in% pbl) & (lhs %in% pbl))
    if (sum(pbl_df$pbl_param) == (length(pbl)*(length(pbl) - 1)/2)){
      # print("PBL is JID!")
      #run pbl_jid_gen function given pbl_df and return list of models
      pbl_gen_list <- pbl_jid_gen(pbl_df)
      return(pbl_gen_list)
    }
  }

  #if PBL is not JID, just return TRUE
  return(TRUE)
}

#' @title Generate model names from just-identified preceding block
#' @description
#' Given a model specification (dataframe of estimated parameters), generates all models involving permutations to the preceding block
#' @param pbl_df dataframe specifying overall model (i.e., all estimated parameters), with an additional column "pbl_param" indicating if the estimated parameter is "within" the preceding block
#' @returns list of model names generated by selecting all possible combinations of permutations involving variables in the preceding block
pbl_jid_gen <- function(pbl_df) {

  #initialize model_list
  model_list <- list()
  #extract rownames
  pbl_rownames <- as.numeric(rownames(pbl_df))
  #for each param in pbl_df
  for (row in 1:nrow(pbl_df)) {
    #if pbl_param is true
    if (pbl_df$pbl_param[row]) {
      #include all iterations of row number into model_list
      model_list[[paste0("param", row)]] <- c(row*3 - 2,
                                              row*3 - 1,
                                              row*3)
    } else {
      #otherwise just include the current row number
      model_list[[paste0("param", row)]] <- pbl_rownames[row]
    }
  }

  #generates all possible models
  expand_pbl <- expand.grid(model_list)
  #initialize final model list
  pbl_jid_models <- list()
  for (row in 1:nrow(expand_pbl)) {
    #for each row in df of possible models, insert into list
    pbl_jid_models[[model_numtostring(unlist(expand_pbl[row,]))]] <-
      unlist(expand_pbl[row,])
  }

  return(pbl_jid_models)
}

#' @title Apply replacing rule
#' @description
#' Given a permutation dataframe, current focal block under consideration (as a dataframe with a single row), and model info, applies the replacing rule (Lee & Hershberger, 1990)
#' @param permute_df permutation dataframe
#' @param cur_row dataframe of 1 row length specifying focal block under investigation
#' @param model_info model info indicating information regarding variable antecedents and consequents
#' @returns list of model generated, if they pass the requirements for the replacing rule
replace_rule <- function(permute_df, cur_row, model_info) {
  #initialize final list of considered models
  considered_models <- list()

  row_lhs <- cur_row$lhs
  row_op <- cur_row$op
  row_rhs <- cur_row$rhs

  #if param is a regression:
  if (row_op == "~") {
    #if all the cause dirantes are in the effect dirantes
    if (all(model_info$dirante[[row_rhs]] %in% model_info$dirante[[row_lhs]])) {
      #replace with covariance
      replace_num <- permute_df %>%
        rownames_to_column("row_ind") %>%
        filter(op == "~~") %>%
        filter(lhs %in% c(row_lhs, row_rhs) & rhs %in% c(row_lhs, row_rhs)) %>%
        .$row_ind %>%
        as.numeric(.)

      #create and insert new model_ind
      new_model_ind <- replace_ind(model_info$ind, replace_num)
      considered_models[[model_numtostring(new_model_ind)]] <- new_model_ind
    }
    #if the cause dirantes are equal to the effect dirantes (- cause)
    if (setequal(model_info$dirante[[row_rhs]],
                 model_info$dirante[[row_lhs]] %>% .[-which(. == row_rhs)])) {
      #replace with reverse reg
      replace_num <- permute_df %>%
        rownames_to_column("row_ind") %>%
        filter(op == "~") %>%
        filter(lhs == row_rhs & rhs == row_lhs) %>%
        .$row_ind %>%
        as.numeric(.)

      #create and insert new model_ind
      new_model_ind <- replace_ind(model_info$ind, replace_num)
      considered_models[[model_numtostring(new_model_ind)]] <- new_model_ind
    }
  }

  #if param is a covariance:
  if (row_op == "~~") {
    #if all rhs dirantes are in lhs dirantes, rhs is cause and lhs is effect
    if (all(model_info$dirante[[row_rhs]] %in% model_info$dirante[[row_lhs]])) {
      replace_num <- permute_df %>%
        rownames_to_column("row_ind") %>%
        filter(op == "~") %>%
        filter(lhs == row_lhs & rhs == row_rhs) %>%
        .$row_ind %>%
        as.numeric(.)

      #create and insert new model_ind
      new_model_ind <- replace_ind(model_info$ind, replace_num)
      considered_models[[model_numtostring(new_model_ind)]] <- new_model_ind
    }
    #if all lhs dirantes are in rhs dirantes, lhs is cause and rhs is effect
    if (all(model_info$dirante[[row_lhs]] %in% model_info$dirante[[row_rhs]])) {
      replace_num <- permute_df %>%
        rownames_to_column("row_ind") %>%
        filter(op == "~") %>%
        filter(lhs == row_rhs & rhs == row_lhs) %>%
        .$row_ind %>%
        as.numeric(.)

      #create and insert new model_ind
      new_model_ind <- replace_ind(model_info$ind, replace_num)
      considered_models[[model_numtostring(new_model_ind)]] <- new_model_ind
    }
  }
  #return considered_models
  return(considered_models)
}

#' @title Single step of generation
#' @description
#' Given a model (internally generated from a permutation dataframe and model info object), conducts a single "step" of generation where the focal block is iterated across all estimated parameters
#' @param permute_df permutation dataframe
#' @param model_info model info indicating information regarding variable antecedents and consequents
#' @returns list of all model names of all models potentially equivalent to the original model (does not account for validity)
step_gen <- function(permute_df, model_info) {
  #given full permutation df and specific model_index, will performe 1 step of generation and output list of all potential models

  #get current model
  cur_model <- permute_df[model_info$ind, ]

  #initialize considered models list
  considered_models <- list()

  #for each row (estimated param)
  #establish limited-block recursiveness
  for (row in 1:nrow(cur_model)) {

    cur_row <- cur_model[row, ]
    # print(paste0(cur_row$lhs, cur_row$op, cur_row$rhs))

    #run determine_LBR function
    LBR_det <- determine_LBR(cur_model, cur_row,
                             model_info)

    #if output is a list, then append to considered_models
    if (is.list(LBR_det)) {
      considered_models <- append(considered_models, LBR_det)

      #otherwise, if output is FALSE, then toss
    } else if (LBR_det == FALSE) {
      # print("skipping...")
      next
    }
    #if you get to this point, you know LBR_det is true (either list, or not false), therefore you can run the regular replacing rule
    replace_models <- replace_rule(permute_df, cur_row,
                                   model_info)
    #append to considered_models
    considered_models <- append(considered_models, replace_models)
  }
  #return final list of considered models (de-duped)
  return(considered_models[!duplicated(considered_models)])
}
