library( dplyr )
library( ape )
library( R.utils )
library( lubridate )

# Args in the following order:
#   [6]: --tree: locations of tree file
#   [7]: --metadata: locations of metadata csv
#   [8]: --id: column in metadata that corresponds to tip labels
#   [9]: --date: column in metadata that corresponds to date

args <- commandArgs(asValue=TRUE, excludeReserved=TRUE )[-1]

t <- read.tree( args$tree )
md <- read.csv( args$metadata )

md <- md %>%
  select( args$id, args$date ) %>%
  filter( args$id %in% t$tip.label )
md[[args$date]] <- as.Date( md[[args$date]] )
md$week <- floor_date( md[[args$date]], unit="week" )

md <- md %>%
  group_by( week ) %>%
  mutate( id_original=!!args$id, !!args$id <- sample( !!args$id ) ) %>%
  ungroup()

write.csv( md, "shuffled_md.csv" )