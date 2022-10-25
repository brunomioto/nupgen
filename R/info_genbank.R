#' Search for info from GenBank data
#'
#' @param organism The advanced search text from GenBank, like "Aylacostoma brunneum[Organism]"
#'
#' @return
#' @export
#'
info_genbank <- function(organism = "Aylacostoma brunneum[Organism]") {

  usethis::ui_done("Buscando GenBank IDs")

  lizard <- organism #We want a character vector

  lizard_search <- rentrez::entrez_search(db="nuccore",
                                          term=lizard,
                                          retmax=500)

  lizard_search$ids #gives you the NCBI ids
  l_ids <- ceiling(length(lizard_search$ids)/100)
  usethis::ui_done(glue::glue("Baixando dados de {length(lizard_search$ids)} sequencias em {l_ids} parte(s) de ~{round(length(lizard_search$ids)/l_ids)} sequencias"))

  #create tibble
  dataset <- tibble::tibble(
    #"ids" = NA,
    #"value" = NA,
    "name" = NA,
    "organism" = NA,
    "country" = NA,
    "lat" = NA,
    "lon" = NA
  )

  dataset_ids <- tibble::tibble(
    "value" = NA,
  )

  for (i in 1:l_ids) {

    usethis::ui_done("Etapa 1")
    tibble_search <- lizard_search$ids %>%
      tibble::as_tibble() %>%
      dplyr::filter(!value %in% dataset_ids$value)

    #sample
    usethis::ui_done("Etapa 2")
    ids_sample <- dplyr::slice_sample(tibble_search,
                                      n = length(lizard_search$ids)/l_ids)

    usethis::ui_done("Etapa 3")
    lizard_seqs <- rentrez::entrez_fetch(db="nuccore",
                                         id=ids_sample %>%
                                           dplyr::pull(),
                                         rettype="gb",
                                         api_key = "5df8cea3896dfea0ecfd69412d818a42c509")

    usethis::ui_done("Etapa 4")
    write(lizard_seqs, "lizards.gb") #gets sequence to a file

    usethis::ui_done("Etapa 5")
    #seq <- read.gb::read.gb("lizards.gb")

    usethis::ui_done("Etapa 6")
    seq2 <- seq %>%
      unlist(recursive = FALSE) %>%
      tibble::enframe() %>%
      dplyr::filter(stringr::str_detect(name, "FEATURES")) %>%
      tidyr::unnest(cols = "value") %>%
      tidyr::unnest(cols = "value") %>%
      dplyr::mutate("name" = stringr::str_remove(name, "\\.FEATURES")) %>%
      dplyr::filter(Location %in% c("organism", "lat_lon", "country", "gene"),
                    !stringr::str_detect(Qualifier, "\\<|\\>")) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(
        names_from = "Location",
        values_from = "Qualifier",
      ) %>%
      tidyr::separate("lat_lon", sep = "\\s",into = c("lat","N_S", "lon","E_W")) %>%
      dplyr::mutate("lat" = ifelse("N_S" == "N", as.numeric(lat), -as.numeric(lat)),
                    "lon" = ifelse("E_W" == "N", as.numeric(lon), -as.numeric(lon)),
                    #gene = ifelse(stringr::str_detect(gene, "\\COX1"), "COI", gene)
      ) %>%
      dplyr::select(-c("N_S", "E_W"))

    usethis::ui_done("Etapa 7")
    dataset_ids <- dplyr::bind_rows(dataset_ids, ids_sample)
    dataset <- dplyr::bind_rows(dataset, seq2)

    usethis::ui_done("Pausa do servidor do GenBank")
    Sys.sleep(3)

  }
  usethis::ui_done("Finalizando")
  final_dataset <- dataset %>%
    dplyr::distinct() %>%
    tidyr::drop_na("name") %>%
    dplyr::rename("code" = "name")

  return(final_dataset)

}
