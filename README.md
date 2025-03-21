
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nupgen <a href="https://brunomioto.github.io/nupgen/"><img src="man/figures/logo.png" align="right" height="106" /></a>

<!-- badges: start -->
<!-- badges: end -->

Pacote de uso interno do laboratório Nupgen DBC/Nupélia

## Instalação

Você pode instalar o pacote nupgen direto do Github usando:

``` r
# install.packages("remotes")
remotes::install_github("brunomioto/nupgen")
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

- Código do GenBank
- Nome do organismo
- País de origem (pode conter mais informações)
- Latitude
- Longitude
- Gene

``` r
library(nupgen)

ayla <- info_genbank(organism = "Aylacostoma brunneum[Organism]")

ayla %>% 
  head()
#> # A tibble: 6 × 6
#>   code     organism             country   lat   lon gene 
#>   <chr>    <chr>                <lgl>   <dbl> <dbl> <chr>
#> 1 JQ236700 Aylacostoma brunneum NA      -27.4 -55.8 COI  
#> 2 JQ236701 Aylacostoma brunneum NA      -27.4 -55.8 COI  
#> 3 KU168374 Aylacostoma brunneum NA      -27.4 -55.8 <NA> 
#> 4 KU168375 Aylacostoma brunneum NA      -27.4 -55.8 <NA> 
#> 5 JQ236702 Aylacostoma brunneum NA      -27.4 -55.8 COI  
#> 6 JQ236704 Aylacostoma brunneum NA      -27.4 -55.8 COI
```

### Formatar resultados da Delimitação de Espécies do PTP (ou bPTP) com `ptp_results()`

Copie os resultados do PTP e salve como um vetor (Só abrir aspas e colar
todo o resultado)

``` r
vetor_ptp <- "# Max likilhood partition 
Species 1 (support = 0.329)
     MZ051026_Characidium_pellucidum_voucher_SU08-1290,MZ051108_Characidium_pellucidum_voucher_SU08-1290_4,MZ051140_Characidium_pellucidum_voucher_SU08-1290_3

Species 2 (support = 0.668)
     MG936837_Characidium_marshi_voucher_stri-3611,MG936835_Characidium_marshi_voucher_stri-11821,MG936836_Characidium_marshi_voucher_stri-4069,MG936834_Characidium_marshi_voucher_stri-6730

Species 3 (support = 0.738)
     KU288836_Characidium_rachovii_voucher_MG_ZV-P_172-2,JX111702_Characidium_rachovii_voucher_UNMDP-T_0289,KU288852_Characidium_rachovii_voucher_MG_ZV-P_172-3,KU288835_Characidium_rachovii_voucher_MG_ZV-P_172-1,KU288853_Characidium_rachovii_voucher_MG_ZV-P_172-4,KU288854_Characidium_rachovii_voucher_MG_ZV-P_172-5,KU289033_Characidium_rachovii_voucher_MG_ZV-P_309

Species 4 (support = 0.188)
     MK464095_Characidium_zebra_isolate_JCB77,MK464130_Characidium_zebra_isolate_JCB286,MK464156_Characidium_zebra_isolate_JCB354,MK464155_Characidium_zebra_isolate_JCB352,MK464153_Characidium_zebra_isolate_JCB349,MK464093_Characidium_zebra_isolate_JCB70,MK464047_Characidium_fasciatum_isolate_CIUnB884_a,MK464134_Characidium_zebra_isolate_JCB298,MK464048_Characidium_fasciatum_isolate_CIUnB884_b,MK464137_Characidium_zebra_isolate_JCB311,MK464135_Characidium_zebra_isolate_JCB299,MK464092_Characidium_zebra_isolate_JCB69,MK464125_Characidium_zebra_isolate_JCB220,MK464040_Characidium_zebra_isolate_CIUnB846_2,MK464152_Characidium_zebra_isolate_JCB348,MK464138_Characidium_zebra_isolate_JCB315,MK464131_Characidium_zebra_isolate_JCB287,MK464126_Characidium_zebra_isolate_JCB221,MK464139_Characidium_zebra_isolate_JCB316,MK464154_Characidium_zebra_isolate_JCB351,MK464094_Characidium_zebra_isolate_JCB76,MK464039_Characidium_zebra_isolate_CIUnB846_1,MK464132_Characidium_zebra_isolate_JCB294,MK464133_Characidium_zebra_isolate_JCB297,MK464140_Characidium_zebra_isolate_JCB317
"
```

Depois, utilize esse vetor dentro da função `ptp_results()`

``` r
ptp_results(vetor_ptp)
#> # A tibble: 39 × 2
#>    label                                               species
#>    <chr>                                               <chr>  
#>  1 MZ051026_Characidium_pellucidum_voucher_SU08-1290   1      
#>  2 MZ051108_Characidium_pellucidum_voucher_SU08-1290_4 1      
#>  3 MZ051140_Characidium_pellucidum_voucher_SU08-1290_3 1      
#>  4 MG936837_Characidium_marshi_voucher_stri-3611       2      
#>  5 MG936835_Characidium_marshi_voucher_stri-11821      2      
#>  6 MG936836_Characidium_marshi_voucher_stri-4069       2      
#>  7 MG936834_Characidium_marshi_voucher_stri-6730       2      
#>  8 KU288836_Characidium_rachovii_voucher_MG_ZV-P_172-2 3      
#>  9 JX111702_Characidium_rachovii_voucher_UNMDP-T_0289  3      
#> 10 KU288852_Characidium_rachovii_voucher_MG_ZV-P_172-3 3      
#> # ℹ 29 more rows
```

### Formatar resultados da análise de haplótipos do DNAsp com `hap_results()`

Copie os resultados análise de haplótipos do DNAsp e salve como um vetor
(Só abrir aspas e colar a parte dos resultados que contém os haplótipos
e nome das sequências)

``` r
vetor_hap <- "    Hap_1: 1  [ANGBF36849-19|Eigenmannia_viresc]
    Hap_2: 1  [ANGBF36850-19|Eigenmannia_viresc]
    Hap_3: 2  [ANGBF56488-19|Eigenmannia_trilin BFFDF015-19|Eigenmannia_trilinea]
    Hap_4: 2  [BSB261-10|Eigenmannia_virescens| BSB264-10|Eigenmannia_virescens|]
    Hap_5: 2  [BSB262-10|Eigenmannia_virescens| BSB260-10|Eigenmannia_virescens|]
    Hap_6: 1  [BSB263-10|Eigenmannia_virescens|]
    Hap_7: 4  [BSFFA646-07|Eigenmannia_humboldt BSFFA647-07|Eigenmannia_humboldt ANGBF36843-19|Eigenmannia_humbol ANGBF36844-19|Eigenmannia_humbol]
    Hap_8: 1  [CIUA1137-21|Eigenmannia_virescen]
"
```

Depois, utilize esse vetor dentro da função `hap_results()`

``` r
hap_results(vetor_hap)
#> # A tibble: 14 × 2
#>    haplotype seq                             
#>    <chr>     <chr>                           
#>  1 1         ANGBF36849-19|Eigenmannia_viresc
#>  2 2         ANGBF36850-19|Eigenmannia_viresc
#>  3 3         ANGBF56488-19|Eigenmannia_trilin
#>  4 3         BFFDF015-19|Eigenmannia_trilinea
#>  5 4         BSB261-10|Eigenmannia_virescens|
#>  6 4         BSB264-10|Eigenmannia_virescens|
#>  7 5         BSB262-10|Eigenmannia_virescens|
#>  8 5         BSB260-10|Eigenmannia_virescens|
#>  9 6         BSB263-10|Eigenmannia_virescens|
#> 10 7         BSFFA646-07|Eigenmannia_humboldt
#> 11 7         BSFFA647-07|Eigenmannia_humboldt
#> 12 7         ANGBF36843-19|Eigenmannia_humbol
#> 13 7         ANGBF36844-19|Eigenmannia_humbol
#> 14 8         CIUA1137-21|Eigenmannia_virescen
```

### Selecionar sequencias únicas de cada haplotipo com `select_seqs_by_hap()`

Copie os resultados análise de haplótipos do DNAsp e salve como um vetor
(Só abrir aspas e colar a parte dos resultados que contém os haplótipos
e nome das sequências).

Carregue o alinhamento (usando `ape::read.FASTA()`) e salve como um
objeto

Depois, utilize o vetor e o alinhamento dentro da função
`select_seqs_by_hap()`

``` r
selected_seqs <- select_seqs_by_hap(haps, seqs)
```

Por fim, salve o alinhamento com as sequências selecionadas utilizando
`ape::write.FASTA()`

## Criar arquivo Arlequin a partir de um arquivo FASTA com `create_arlequin()`

A função create_arlequin() cria um arquivo Arlequin (.arp) a partir de
um arquivo FASTA com sequências de DNA. Para isso, basta carregar o
arquivo FASTA, um arquivo de grupos e usar a função.

- O arquivo de grupos deve ser um data.frame com 2 colunas: “group”
  (nome dos grupos) e “name” (nome das sequências)

- Caso tenha dúvidas, o pacote tem um arquivo de exemplo, rode:
  `nupgen::groups_example`

``` r
fasta <- read.dna("arquivo_fasta.fas",format = "fasta")

grupos <- read.csv("arquivo_grupos.csv")

create_arlequin(fasta, grupos)
```
