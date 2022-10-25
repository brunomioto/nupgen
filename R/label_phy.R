#' Label phylip files
#'
#' @param seq_fasta The FASTA alignment file with full names
#' @param tree The .nwk file generated from .phy file
#'
#' @return a .nwk object
#' @export
#'

label_phy <- function(seq_fasta, tree){

  seq_names <- ape::labels.DNAbin(seq_fasta) %>%
    dplyr::as_tibble()

  tree_label <- tree$tip.label %>%
    dplyr::as_tibble()

  tabela_junto <- seq_names %>%
    dplyr::rename("seq_longo" = "value") %>%
    dplyr::mutate("novo_nome" = substr("seq_longo", start = 1, stop = 10))

  novos_nomes <- tree_label %>%
    dplyr::left_join(tabela_junto, by = c("value" = "novo_nome")) %>%
    dplyr::pull("seq_longo")

  tree$tip.label <- novos_nomes

  return(tree)
}
