#!/usr/bin/env bash

##### OPTIMAL TRANSPORT ####
# wot command line interface

wot optimal_transport --matrix matrix_MM.h5ad --cell_days cells_day.txt --growth_iters 3 --lambda1 1 --lambda2 50 --epsilon 0.05 --verbose
wot trajectory --tmap tmaps --cell_set cell_sets.gmt --day 30 --embedding embedding_coord.txt
wot transition_table --tmap tmaps --cell_set cell_sets.gmt --start_time 0 --end_time 30
