#' Label phylip files
#'
#' @param seq_fasta The FASTA alignment file with full names
#' @param tree The .nwk file generated from .phy file
#'
#' @return a .nwk object
#' @export
#'
#' @examples
#' seq_normal <- ape::read.FASTA("alinhamento_cortado_2022_10_13_mioto_outgroup_3_limpo.fas")
#' tree_exe <- ape::read.tree("file.nwk")
#' new_tree <- label_phy(seq_normal, tree_exe)
#' write.tree(new_tree, "new_tree.nwk")

label_phy <- function(seq_fasta, tree){

  seq_names <- ape::labels.DNAbin(seq_fasta) %>%
    dplyr::as_tibble()

  tree_label <- tree$tip.label %>%
    dplyr::as_tibble()

  tabela_junto <- seq_names %>%
    dplyr::rename(seq_longo = value) %>%
    dplyr::mutate(novo_nome = substr(seq_longo, start = 1, stop = 10))

  novos_nomes <- tree_label %>%
    dplyr::left_join(tabela_junto, by = c("value" = "novo_nome")) %>%
    dplyr::pull(seq_longo)

  tree$tip.label <- novos_nomes

  return(tree)
}
