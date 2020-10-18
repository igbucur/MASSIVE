#' Function for generating the full IV model (only slab priors).
#'
#' @param J Integer number of candidate instruments.
#'
#' @return Character vector representing full IV model.
#' @export
#'
#' @examples
#' get_full_IV_model(5)
get_full_IV_model <- function(J) {
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(rep("1", J), collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}

#' Function for generating the empty IV model (only spike priors).
#'
#' @param J Integer number of candidate instruments.
#'
#' @return Character vector representing empty IV model.
#' @export
#'
#' @examples
#' get_empty_IV_model(5)
get_empty_IV_model <- function(J) {
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(rep("0", J), collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}


#' Function for generating a random IV model (combination of slab and spike priors).
#'
#' @param J Integer number of candidate instruments.
#'
#' @return Character vector representing random IV model.
#' @export
#'
#' @examples
#' get_random_IV_model(5)
get_random_IV_model <- function(J) {
  sp_or_sl <- as.character(sample(c(0, 1), J, replace = TRUE))
  paste(c(
    paste(rep("1", J), collapse = ""),
    paste(sp_or_sl, collapse = ""),
    "1", "1", "1"
  ), collapse = "|")
}

#' Function for decoding a character vector representing IV model.
#'
#' @param code Character vector encoding IV model.
#' @param sd_slab Standard deviation of slab component.
#' @param sd_spike Standard deviation of spike component.
#'
#' @return List describing prior for IV model.
#' @export
#'
#' @examples
#' decode_IV_model(code = get_random_IV_model(5), sd_slab = 1, sd_spike = 0.01)
decode_IV_model <- function(code, sd_slab, sd_spike) {
  tokens <- strsplit(code, "\\|")[[1]]
  IV_model <- stringr::str_extract_all(tokens, "[01]")
  names(IV_model) <- c('sgamma', 'salpha', 'sbeta', 'skappa_X', 'skappa_Y')
  IV_model <- lapply(IV_model, function(x) ifelse(x == "1", sd_slab, sd_spike))
  IV_model$sd_slab <- sd_slab
  IV_model$sd_spike <- sd_spike
  
  IV_model
}
