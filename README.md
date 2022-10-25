
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

### Renomear arquivos PHYLIP

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
write.tree(new_tree, "nova_tree.nwk")
```