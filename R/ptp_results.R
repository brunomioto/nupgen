#' Format PTP Species Delimitation results
#'
#' @param vector
#'
#' @return
#' @export
ptp_results <- function(vector) {
  new_ptp <- vector %>%
    stringr::str_split("\\\n") %>%
    unlist() %>%
    dplyr::tibble() %>%
    dplyr::rename("label" = 1) %>%
    dplyr::filter(label != "") %>%
    dplyr::pull() %>%
    stringr::str_split("\\,") %>%
    unlist() %>%
    dplyr::tibble() %>%
    dplyr::rename("label" = 1) %>%
    dplyr::mutate(species = ifelse(stringr::str_detect(label, "\\#|Species"), stringr::str_extract(label, "[:digit:]*(?=\\s\\()") , NA),
           label = ifelse(is.na(species), stringr::str_remove_all(stringr::str_trim(label),"\\'"), NA)
    ) %>%
    tidyr::fill(species) %>%
    tidyr::drop_na()

  return(new_ptp)
}
