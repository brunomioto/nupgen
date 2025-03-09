#' Create Arlequin Input File (.arp)
#'
#' This function generates an Arlequin input file (.arp) from a DNA alignment
#' and a groups data frame. The alignment should be in the DNAbin format, and
#' the groups data frame should have two columns: 'group' and 'name'. The function
#' checks the validity of the inputs and creates an Arlequin file in the specified
#' directory. If the output file already exists, the user will be prompted to delete it.
#' This script was adapted from Josh Banta (source)
#'
#' @param fasta A DNAbin object representing the DNA alignment.
#' @param groups A data.frame with two columns: 'group' and 'name', representing
#'               the sample groups.
#' @param output.dir A character string specifying the directory where the
#'                   output file will be saved. Default is the current working directory.
#' @author Josh Banta
#' @source https://sites.google.com/site/thebantalab/tutorials?authuser=0#h.39zz6drphh3u
#'
#' @return This function does not return any value. It creates a file named 'output.arp'
#'         in the specified directory.
#'
#' @importFrom cli cli_abort cli_alert_info cli_alert_success
#' @export
create_arlequin <- function(fasta,
                            groups,
                            output.dir = ".") {

  if (!inherits(fasta, "DNAbin")) {
  cli::cli_abort("The alignment must be a DNAbin object.
                                                 Use the function `ape::read.dna()`")
  }

  if (!is.data.frame(groups)) {
    cli::cli_abort("The provided groups file is not a data.frame")
  }

  if (!ncol(groups) == 2) {
    cli::cli_abort("The groups file must have two columns")
  }

  if (!identical(colnames(groups)[1:2], c("group", "name"))) {
    cli::cli_abort("The first two columns of the groups file must be 'group' and 'name'")
  }

  dat.matrix <- as.character(fasta)

  if(file.exists(paste0(output.dir, "/output.arp"))){

    cli::cli_alert_info("Before proceeding, you need to delete the existing output.arp file. Confirm?")

    resposta <- readline(prompt = "Type 'y' for yes or 'n' for no: ")

    if (tolower(resposta) == "y") {

      cli::cli_alert_success("File deleted")

      file.remove(paste0(output.dir, "/output.arp"))

    } else {

      cli::cli_abort("Delete the output.arp file and try again")

    }
  }

  cli::cli_alert_info("Creating .arp file")

  outfile <- cbind(1,rownames(dat.matrix))

  colnames(outfile) <- c("group", "name")

  dat.matrix[dat.matrix == "-"] <- "?"

  dat.matrix[dat.matrix == "N"] <- "?"

  dat.matrix[dat.matrix == "n"] <- "?"

  x2 <- dat.matrix

  colnames(x2) <- paste("locus", seq(ncol(x2)), sep = "")

  x3 <- rbind(x2, x2)
  row.names(x3) <- NULL

  x3 <- cbind(x3, 1)

  colnames(x3)[ncol(x3)] <- "pop"

  for(i in 1:nrow(x2)){

    x3[((2*i)-1), 1:(ncol(x3)-1)] <- x2[i,]
    x3[((2*i)), 1:(ncol(x3)-1)] <- x2[i, ]

    x3[((2*i)-1),"pop"] <- rownames(x2)[i]
    x3[((2*i)),"pop"] <- rownames(x2)[i]


  }

  x4 <- as.matrix(x3)

  is.odd <- function(x) x %% 2 != 0
  is.even <- function(x) x %% 2 == 0

  for(i in 1:(ncol(x4)-1)){

    for(j in 1:nrow(x4)){

      if(!is.na(x4[j,i])){

        if(x4[j,i] == "y"){

          if(is.odd(j)){

            x4[j,i] <- "?"

          }else{

            x4[j,i] <- "?"

          }

        }

        if(x4[j,i] == "r"){

          if(is.odd(j)){

            x4[j,i] <- "?"

          }else{

            x4[j,i] <- "?"

          }

        }

        if(x4[j,i] == "w"){

          if(is.odd(j)){

            x4[j,i] <- "?"

          }else{

            x4[j,i] <- "?"

          }

        }

        if(x4[j,i] == "s"){

          if(is.odd(j)){

            x4[j,i] <- "?"

          }else{

            x4[j,i] <- "?"

          }

        }

        if(x4[j,i] == "k"){

          if(is.odd(j)){

            x4[j,i] <- "?"

          }else{

            x4[j,i] <- "?"

          }

        }

        if(x4[j,i] == "m"){

          if(is.odd(j)){

            x4[j,i] <- "?"

          }else{

            x4[j,i] <- "?"

          }

        }

      }

    }

  }

  rownames(x4) <- x4[,"pop"]

  x4 <- x4[order(rownames(x4)),]

  x5 <- x4[,-which(colnames(x4) == "pop")]

  x6 <- NULL
  for(i in 1:nrow(x5)){

    x6 <- rbind(x6, noquote((paste(x5[i,], collapse = " "))))

  }

  rownames(x6) <- rownames(x5)

  x7 <- x6[which(is.odd(1:nrow(x6))),]

  cdata <- groups

  cli::cli_alert_info("Saving .arp file")

  crlf <- "\r\n"

  outfile = paste0(output.dir, "/output.arp")

  cat('[Profile]', file=outfile, append=FALSE, crlf)

  cat('', file=outfile, append=TRUE, crlf)

  cat('Title="data"', file=outfile, append=TRUE, crlf)

  write.table(paste("NBSamples=",length(table(cdata[,1])), sep = ""),
              append = TRUE, paste0(output.dir, "/output.arp"),
              row.names = FALSE, col.names = FALSE, quote = FALSE, crlf)

  cat("", file=outfile, append=TRUE, crlf)

  cat('DataType=DNA', file=outfile, append=TRUE, crlf)

  cat('GenotypicData=0', file=outfile, append=TRUE, crlf)

  cat('LocusSeparator=WHITESPACE', file=outfile, append=TRUE, crlf)

  cat("", file=outfile, append=TRUE, crlf)

  cat("[Data]", file=outfile, append=TRUE, crlf)

  cat("[[Samples]]", file=outfile, append=TRUE, sep = "")

  for(i in 1:length(names(table(cdata[,1])))){

    cat("", file=outfile, append=TRUE, crlf)

    cat("", file=outfile, append=TRUE, crlf)

    cat('SampleName=', file=outfile, append=TRUE, sep = "")

    write.table(names(table(cdata[,1]))[i], append = TRUE,
                paste0(output.dir, "/output.arp"),
                row.names = FALSE, col.names = FALSE,
                quote = TRUE, sep = "", crlf)

    write.table(paste('SampleSize=',
                      length(which(tolower(cdata[,1]) == tolower(names(table(cdata[,1])))[i])), sep = ""),
                append = TRUE, paste0(output.dir, "/output.arp"), row.names = FALSE, col.names = FALSE,
                quote = FALSE, sep = "", crlf)

    cat('SampleData={', file=outfile, append=TRUE, sep = "")

    cat("", file=outfile, append=TRUE, crlf)

    cdata.sub <- cdata[which(tolower(cdata[,1]) == tolower(names(table(cdata[,1]))[i])),]

    for(j in 1:nrow(cdata.sub)){

      cat("", file=outfile, append=TRUE, crlf)

      x6.sub <- x7[which(names(x7) == cdata.sub[j,2])]

      write.table(paste(tolower(cdata.sub[j,2]), 1, toupper(x6.sub), sep = " "),
                  append = TRUE,
                  paste0(output.dir, "/output.arp"),
                  row.names = FALSE, col.names = FALSE,
                  quote = FALSE, sep = "", crlf)
    }

    cat("", file=outfile, append=TRUE, crlf)

    cat('}', file=outfile, append=TRUE, sep = "")

  }

  cat('', file=outfile, append=TRUE, crlf)

  cat('', file=outfile, append=TRUE, crlf)

  cat('[[Structure]]', file=outfile, append=TRUE, crlf)

  cat('StructureName="Test"', file=outfile, append=TRUE, crlf)

  cat('NbGroups=1', file=outfile, append=TRUE, crlf)

  cat('', file=outfile, append=TRUE, crlf)

  cat('Group={', file=outfile, append=TRUE, crlf)

  cat('', file=outfile, append=TRUE, crlf)

  for(i in 1:length(names(table(cdata[,1])))){

    write.table(names(table(cdata[,1]))[i], append = TRUE,
                paste0(output.dir, "/output.arp"),
                row.names = FALSE, col.names = FALSE,
                quote = TRUE, sep = "", crlf)

  }

  cat('', file=outfile, append=TRUE, crlf)

  cat('}', file=outfile, append=TRUE, sep = "")

}
