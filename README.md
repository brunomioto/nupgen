
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nupgen

<!-- badges: start -->
<!-- badges: end -->

Pacote de uso interno do laboratório Nupgen DBC/Nupélia

## Instalação

Você pode instalar o pacote nupgen direto do Github usando:

``` r
# install.packages("devtools")
devtools::install_github("brunomioto/nupgen")
```

## Uso

### Renomear arquivos PHYLIP com `label_phy()`

``` r
library(nupgen)
library(ape)

#carrega arquivo FASTA (que possui os nomes corretos)
seq <- ape::read.FASTA("sequencia_fasta.fas")
#carrega arquivo .nwk (que possui os nomes cortados)
tree <- ape::read.tree("arvore.nwk")

#Função para renomear tip labels da árvore
new_tree <- label_phy(seq, tree)

#Confere se os novos nomes estão certos
new_tree$tip.label

#Salva em um novo arquivo .nwk
ape::write.tree(new_tree, "nova_tree.nwk")
```

### Obter informações de sequências do GenBank com `info_genbank()`

As informações buscadas são:

-   Código do GenBank
-   Nome do organismo
-   País de origem (pode conter mais informações)
-   Latitude
-   Longitude
-   Gene

``` r
library(nupgen)

ayla <- info_genbank(organism = "Aylacostoma brunneum[Organism]")

ayla %>% 
  head()
#> # A tibble: 6 x 6
#>   name     organism             country    lat   lon gene 
#>   <chr>    <chr>                <chr>    <dbl> <dbl> <chr>
#> 1 JQ236700 Aylacostoma brunneum Paraguay -27.4 -55.8 COI  
#> 2 JQ236703 Aylacostoma brunneum Paraguay -27.4 -55.8 COI  
#> 3 KU168373 Aylacostoma brunneum Paraguay -27.4 -55.8 <NA> 
#> 4 JQ236702 Aylacostoma brunneum Paraguay -27.4 -55.8 COI  
#> 5 KF918858 Aylacostoma brunneum Paraguay -27.4 -55.8 cytb 
#> 6 KU168374 Aylacostoma brunneum Paraguay -27.4 -55.8 <NA>
```
