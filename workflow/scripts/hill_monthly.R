library( ape )
library( phytools )
library( hillR )
library( dplyr )
library( R.utils )
library( lubridate )

phylosorTable <- function( tree, metadata, locations, window ) {
    date_seq <- seq( min( metadata$date_collected ), max( metadata$date_collected ), by="months" )
    count <- 0
    for( i in 2:length( date_seq ) ) {
        ptm <- proc.time()
        comm <- matrix( 0, nrow=length( locations ), ncol=length( tree$tip.label ) )

        for( j in 1:length( locations ) ) {
            comm[j,] <- tree$tip.label %in% (md %>%
                                            filter( (site == locations[[j]]) & between(date_collected, date_seq[i-1], date_seq[i]) ) %>%
                                            pull( accession_id ) )
        }
        if( max( comm ) == 0 ){
            printf( "%s being skipped. No sequences\n", date_seq[i] )
            next
        }

        colnames( comm ) <- tree$tip.label
        rownames( comm ) <- locations

        count <- count + 1

        res_df <- hill_phylo_parti_pairwise( comm, tree, q=1 )
        res_df$count1 <- rowSums( comm )[locations[1]]
        res_df$count2 <- rowSums( comm )[locations[2]]
        res_df$date <- date_seq[i]

        if( count == 1 ) {
            return_df <- res_df
        } else {
            return_df <- rbind( return_df, res_df )
        }

        printf( "Completed %d of %d comparisons (took %.1f seconds).\n", i, length( date_seq ), (proc.time() - ptm)["elapsed"] )
    }
    return( return_df )
}

args <- commandArgs(asValue=TRUE, excludeReserved=TRUE, defaults=list( shuffle=FALSE ) )[-1]

print( args )

loc_pair <- c( args$pair1, args$pair2 )

t <- read.tree( args$tree ) %>% root( outgroup = "EPI_ISL_402125", resolve.root = TRUE )
# hillR throws a fit if there are some edgelengths that are NaN so I replace them with 0 because thats probably what they are.
t$edge.length[is.na(t$edge.length)] <- 0
# Hardcode these for now, but should follow the rest of the program in no assuming them.
md <- read.csv( args$metadata ) %>%
  select( accession_id, date_collected, site ) %>%
  filter( accession_id %in% t$tip.label ) %>%
  mutate( date_collected = as.Date( date_collected ) )
md$week <- floor_date( md$date_collected, unit="week" )

if ( args$shuffle ) {
  md <- md %>%
    group_by( week ) %>%
    mutate( site_original=site, site=sample( site ) ) %>%
    ungroup()
}

results <- phylosorTable( t, md, loc_pair, as.integer( args$window_size ) )
write.csv( results, args$output )
