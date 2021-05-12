# Normalize the scores in a bedgraph file

read_bedgraph <- function(bg.path){

    bg.path <- as.character(bg.path)
    as.data.frame(read.table(bg.path, header=F))

}

min_max_norm <- function(x) {

    (x - min(x)) / (max(x) - min(x))

  }


normalize_scores <- function(bg.df){

    big.df.norm <- big.df
    big.df$V4 <- lapply(big.df$V4, min_max_norm)
    big.df.norm

}

write_df_as_bedgraph <- function(bg.df.norm, output.path){

    output.path <- as.character(output.path)
    write.table(big.df.norm, output.path, row.names=F, col.names=F,
                quote=F, sep='\t')


}

main <- function(){

    bg.df <- read_bedgraph(snakemake$input)
    bg.df.norm <- normalize_scores(bg.df)
    write_df_as_bedgraph(big.df.norm, snakemake$output)

}


if ( !interactive() ){

    main()

}