library( ape )
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

tree_loc <- "/Users/natem/Dropbox (Scripps Research)/Personal/Code/Python/monterey/test/test_phylosor/SanDiego-BajaCalifornia.actual.1.tree"
md_loc <- "/Users/natem/Dropbox (Scripps Research)/Personal/Code/Python/monterey/resources/combined_md.csv"


t <- read.tree( tree_loc ) %>% root( outgroup = "EPI_ISL_402125", resolve.root = TRUE )
# Hardcode these for now, but should follow the rest of the program in no assuming them.
md <- read.csv( md_loc ) %>%
  select( accession_id, date_collected, site ) %>%
  filter( accession_id %in% t$tip.label ) %>%
  mutate( date_collected = as.Date( date_collected ) )

locations <- c( "San Diego", "Baja California")
date_seq <- seq( as.Date( "2020-11-01" ), as.Date( "2020-12-01" ), by="days")

comm <- matrix( 0, nrow=length( locations ), ncol=length( t$tip.label ) )

for( j in 1:length( locations ) ) {
    comm[j,] <- t$tip.label %in% (md %>%
      filter( (site == locations[[j]]) & (date_collected %in% date_seq) ) %>%
      pull( accession_id ) )
}
if( max( comm ) == 0 ){
    printf( "%s being skipped. No sequences\n", date_seq[1] )
    next
}

colnames( comm ) <- t$tip.label
rownames( comm ) <- locations

print( comm )

ptm <- proc.time()
res <- phylosor( comm, t )
print( res )
res_df <- toDataFrame( res )
counts <- data.frame( rowSums( comm ) )
res_df <- merge( res_df, counts, by.x="row", by.y=0 )
res_df <- merge( res_df, counts, by.x="col", by.y=0 )
colnames( res_df )[4] <- "row_count"
colnames( res_df )[5] <- "col_count"
res_df$date <- date_seq[1]
print( res_df )
printf( "Completed comparisons (took %.1f seconds).\n", (proc.time() - ptm)["elapsed"] )


