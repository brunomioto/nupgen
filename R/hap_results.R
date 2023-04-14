#' Format DNAsp Haplotype results
#' @description Splits a vector into a tibble with two columns, haplotype and sequence
#' @param vector A character vector
#' @return A tibble with two columns, haplotype and sequence
#' @export
#'
#' @examples
#' hap_results(vector)
hap_results <- function(vector) {
  new_hap <- vector %>%
    stringr::str_split("\\\n") %>%
    unlist() %>%
    dplyr::tibble() %>%
    dplyr::rename("label" = 1) %>%
    dplyr::filter(label != "") %>%
    dplyr::mutate(label = stringr::str_trim(label)) %>%
    tidyr::separate(col = label, into = c("haplotype", "seq"), sep = "  ") %>%
    dplyr::mutate(haplotype = stringr::str_remove_all(haplotype, "Hap_|\\:.*"),
                  seq = stringr::str_remove_all(seq, "\\[|\\]")) %>%
    tidyr::separate_rows(seq, sep = "\\s")

  return(new_hap)
}

