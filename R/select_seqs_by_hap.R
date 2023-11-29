#' Select sequences by haplotype
#'
#' This function takes a set of haplotypes and a set of sequences, and returns a subset of sequences
#' corresponding to the unique haplotypes.
#'
#' @param haps A vector containing haplotype information.
#' @param seqs A DNA alignment object.
#'
#' @return A subset of sequences corresponding to the unique sequences by haplotype.
#'
#' @export
select_seqs_by_hap <- function(haps, seqs) {

  table1 <- nupgen::hap_results(haps) %>%
    dplyr::mutate(haplotype = as.numeric(haplotype)) %>%
    dplyr::group_by(haplotype) %>%
    dplyr::slice_head(n=1)

  seq_use <- ape::labels.DNAbin(seqs) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(value %in% table1$seq) %>%
    dplyr::pull()

  alignment_final <- seqs[seq_use]

  return(alignment_final)

}
