import argparse
import gzip
import csv

### get the args
parser = argparse.ArgumentParser()
parser.add_argument("-v", type=str, help="input vcf file", required=True )
parser.add_argument("-s", type=str, help="sample to population file", required=True )
parser.add_argument("-g", type=int, help="boolean (0 = use sample reads| 1 = use genotypes)", default = 0 )
parser.add_argument("-r", type=float, help="uniform recombination rate per bp in morgans/bp (float)", default = 1e-8 )
parser.add_argument("-m", type=int, help="minimum distance between successive snps (int,bp)", default = 1000 )
parser.add_argument("-o", type=str, help="ploidy file for ahmm input", default = "ahmm.ploidy")
parser.add_argument("--min_total", type=int, help="minimum number of reference population samples (int)", default = 10 )
parser.add_argument("--min_diff",type=float, help="minimum allele frequency difference between a pair of populations (float)", default = 0.1 )
args = parser.parse_args()

### get sample to population mappings
sample2pop = {}
with open( args.s ) as file :
    for line in file :
        (key, val) = line.split()
        sample2pop[key] = val

### IDs
sample_fields = []

### keep track of where we are on the chromosome
last_position = -10000000
current_chrom = "NA"

### have we printed the ploidy file yet
ploidy_print = 0

## read the vcf file 
with open( args.v ) as tsv :

    ### split on tab
    for line in csv.reader(tsv, delimiter="\t") :

        ### parse the header line
        if ( line[0].startswith('#CHROM') ) :

            ### get our info
            sample_files = line
            continue

        ### skip lines that start with ##
        if ( line[0].startswith('##') ) :
            continue

        ### remove sites that are not biallelic snps
        if ( len(line[3]) != 1 or len(line[4]) != 1 ) :
            continue

        ### now read the line and decide what to do with it
        if ( current_chrom != line[0] ) :
            current_chrom = line[0]
            last_position = -1000000

        ### if position too close skip
        elif ( int( line[1] ) - last_position < args.m ) :
            continue

        ### otherwise get population allele frequencies
        alts = {}
        totals = {}
        for s in range(9, len( line ) ) :
            if ( sample_files[s] in sample2pop ) :
                if ( sample2pop[sample_files[s]] not in totals ) :
                    totals[sample2pop[sample_files[s]]] = 0
                    alts[sample2pop[sample_files[s]]] = 0
                if ( line[s].startswith('0/0') ) :
                    totals[sample2pop[sample_files[s]]] += 2
                elif ( line[s].startswith('0/1') ) :
                    alts[sample2pop[sample_files[s]]] += 1
                    totals[sample2pop[sample_files[s]]] += 2
                elif ( line[s].startswith('1/1') ) :
                    alts[sample2pop[sample_files[s]]] += 2
                    totals[sample2pop[sample_files[s]]] += 2

                ### read haploid genotypes as needed
                elif ( line[s].startswith('1') ) :
                    alts[sample2pop[sample_files[s]]] += 1
                    totals[sample2pop[sample_files[s]]] += 1
                elif ( line[2].startswith('0') ) :
                    totals[sample2pop[sample_files[s]]] += 1

        ### make sure we have enough samples at this position
        too_few = 0
        freqs = []
        for population, total in totals.items() :
            if ( population != "admixed" ) :
                freqs.append( float(alts[population])/float(total) )
                if ( total < args.min_total ) :
                    too_few += 1
        if ( too_few > 0 ) :
            continue

        ### check allele frequencies
        diff = 0
        for p1 in range(0, len(freqs)-1) :
            for p2 in range(p1+1,len(freqs)) :
               if ( abs( freqs[p1] - freqs[p2] ) > args.min_diff ) :
                   diff += 1
        if ( diff == 0 ) :
            continue

        ### if we got this far, we're rpinting the site
        print( line[0], line[1], end ="\t", sep = "\t" )
        for population, total in totals.items() :
            if ( population != "admixed" ) :
                ref = total - alts[population]
                print ( ref, alts[population], end = "\t", sep = "\t" )

        rec = float ( float( line[1] ) - last_position ) * args.r
        print ( rec, end = "" )

        ### now print the read counts for each admixed sample
        for s in range(9, len( line ) ) :

            ### only admixed samples
            if ( sample_files[s] in sample2pop and sample2pop[sample_files[s]] == "admixed" ) :

                ### print genotypes
                if ( args.g == 1 ) :
                    if ( line[s].startswith('0/0') ) :
                        print( "\t", 2, "\t", 0, end = "", sep = "" )
                    elif ( line[s].startswith('0/1') ) :
                        print( "\t", 1, "\t", 1, end = "", sep = "" )
                    elif ( line[s].startswith('1/1') ) :
                        print( "\t", 0, "\t", 2, end = "", sep = "" )
                    elif ( line[s].startswith('0') ) :
                        print( "\t", 1, "\t", 0, end = "", sep = "" )
                    elif ( line[s].startswith('1') ) :
                        print( "\t", 0, "\t", 1, end = "", sep = "" )
                    else :
                        print( "\t", 0, "\t", 0, end = "", sep = "" )
                    if ( ploidy_print == 0 ) :
                        print(sample_files[s], 2, sep = "\t", file=open(args.o, "a") )

                ### print read counts
                else :
                    fields = line[s].split( sep = ":" )
                    if ( len( fields ) > 2 ) :
                        ad = fields[1].split( sep = "," )
                        print ( "\t", ad[0], "\t", ad[1], end = "", sep = "" )
                    else :
                        print ( "\t", 0, "\t", 0,  end = "", sep = "" )
                    if ( ploidy_print == 0 ) :
                        print(sample_files[s], 2, sep = "\t", file=open(args.o, "a") )

        ## newline
        print ()

        ### update last position after printing
        last_position = int( line[1] )
        ploidy_print = 1
