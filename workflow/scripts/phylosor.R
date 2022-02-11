library( ape )
#library( phytools )
library( picante )
library( dplyr )
library( R.utils )
library( lubridate )

toDataFrame <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}

phylosorTable <- function( tree, metadata, locations, window ) {
    date_seq <- seq( min( metadata$date_collected ), max( metadata$date_collected ), by="days" )
    count <- 0
    for( i in (window+1):length( date_seq ) ) {
        ptm <- proc.time()
        comm <- matrix( 0, nrow=length( locations ), ncol=length( tree$tip.label ) )

        for( j in 1:length( locations ) ) {
            comm[j,] <- tree$tip.label %in% (md %>%
                                            filter( (site == locations[[j]]) & (date_collected %in% date_seq[(i-window):i]) ) %>%
                                            pull( accession_id ) )
        }
        if( max( comm ) == 0 ){
            printf( "%s being skipped. No sequences\n", date_seq[i] )
            next
        }

        colnames( comm ) <- tree$tip.label
        rownames( comm ) <- locations

        count <- count + 1

        res <- phylosor( comm, tree )
        res_df <- toDataFrame( res )
        counts <- data.frame( rowSums( comm ) )
        res_df <- merge( res_df, counts, by.x="row", by.y=0 )
        res_df <- merge( res_df, counts, by.x="col", by.y=0 )
        colnames( res_df )[4] <- "row_count"
        colnames( res_df )[5] <- "col_count"
        res_df$date <- date_seq[i-window]

        if( count == 1 ) {
            return_df <- res_df
        } else {
            return_df <- rbind( return_df, res_df )
        }

        printf( "Completed %d of %d comparisons (took %.1f seconds).\n", i, length( date_seq ), (proc.time() - ptm)["elapsed"] )
    }
    return( return_df )
}

t <- read.tree( snakemake@input[["tree"]] ) %>% root( outgroup = "EPI_ISL_402125", resolve.root = TRUE )
# Hardcode these for now, but should follow the rest of the program in no assuming them.
md <- read.csv( snakemake@input[["metadata"]] ) %>%
  select( accession_id, date_collected, site ) %>%
  filter( accession_id %in% t$tip.label ) %>%
  mutate( date_collected = as.Date( date_collected ) )
save.image( file="/gpfs/home/natem/analysis/2022.02.08_phylosor/monterey/.RData")
results <- phylosorTable( t, md, snakemake@params[["pair_list"]], as.integer( snakemake@params[["window_size"]] ) )
write.csv( results, snakemake@output[["results"]] )
