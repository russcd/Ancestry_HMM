=begin
 
 This is a simple script to take VCF formatted genotype data and produce input for Ancestry_HMM
 
 Assuming:
 
 1. SNPs have been prefiltered for quality.
 2. Genotypes for reference samples are valid + reasonably high confidence
 3. A uniform recombination rate per bp.
 4. All samples (references and admixed) are included in a single VCF.
 
 This script should be considered an alpha release—it has not been thoroughly validated across a range of possible input datatypes. We are currently working on more optimal solutions to this problem including principled LD pruning, k ancestral populations, and other recommendations for best practices data pre-filtering. Nonetheless, we hope this is helpful for some users.
 
=cut

use strict ;
use warnings ; 

### parameters
my $use_reads = 0 ;              ### if using read counts from admixed samples (1) or using genotypes (0)
my $rec_rate = 1e-8 ;            ### average per bp recombination rate in morgans
my $min_f_diff = 0.1 ;           ### minimum frequency difference between reference populations
my $min_called = 10 ;            ### minimum number of *alleles* in a reference panel to consider a site
my $min_dist = 10000 ;           ### minimum distance in basepairs between adjacent emitted sites
                                 ###   This will accomplish lazy LD pruning.

### ranges of 0-based column numbers for each reference population and the admixed samples. These should be set by the user prior to using this script. s
my @ref1 = ( 9..23 ) ;
my @ref2 = ( 24..36 ) ;
my @ad = ( 37..46 ) ;

### bookkeeping stuff
my $last_pos = -1 * $min_dist ;
my $chrom = "" ;

### read through SNPs
while (<STDIN>) { 

	### drop headers
	if ( $_ =~ m/^\#/ ) { 
		next ;
	}

	## chomp and split
	chomp ; 
	my @split = split ( /\t/, $_ ) ; 

	### new chrom, reset last/chrom
	if ( $chrom ne $split[0] ) { 
		$last_pos = $min_dist * -1 ; 
		$chrom = $split[0] ; 
	}

	### sites that are above the minimum distance
	if ( $split[1] - $last_pos < $min_dist ) { 
		next ;
	}

	### only biallelic SNPs, anything else will have length > 1 in vcf format
	if ( length( $split[4] ) > 1 || length ( $split[3] ) > 1 ) { 
		next ;
	}

	### get counts of each reference panel
	my $ref1_A = 0 ;
	my $ref1_count = 0 ;
	foreach my $o ( @ref1 ) {
		if ( $split[$o] =~ m/^(\d)\/(\d)/ ) { 
			$ref1_count += 2 ;
			$ref1_A += $1 + $2 ;
		}
	}

	### get counts for the vicuna panel
	my $ref2_A = 0 ;
	my $ref2_count = 0 ;
	foreach my $t ( @ref2 ) {
		if ( $split[$t] =~ m/^(\d)\/(\d)/ ) {
            $ref2_count += 2 ;
            $ref2_A += $1 + $2 ;
        }
    }

	### too few called of references
	if ( $ref1_count < $min_called || $ref2_count < $min_called ) {
		next ;
	}

	### not enough differentiation
	if ( abs( $ref1_A/$ref1_count - $ref2_A/$ref2_count ) < $min_f_diff ) {
		next ;
	}

	### now print site + hybids
	print $split[0], "\t", $split[1], "\t" ;
	print $ref1_count - $ref1_A, "\t", $ref1_A , "\t" ;
	print $ref2_count - $ref2_A, "\t", $ref2_A, "\t" ;
	print abs($split[1]-$last_pos)*$rec_rate ;
	foreach my $o ( @ad ) {
        
        if ( $use_reads == 0 ) {
            if ( $split[$o] =~ m/^(\d+)\/(\d+)/ ) {
                my $geno = $1 + $2 ;
                print "\t", 2-$geno, "\t", $geno ;
            }
            else {
                print "\t", 0, "\t", 0 ;
            }
		}
        else {
            if ( $split[$o] =~ m/^\d+\/\d+:(\d+)\,(\d+)/ ) {
                print "\t", $1, "\t", $2 ;
            }
            else {
                print "\t", 0, "\t", 0 ;
            }
        }
	}
	print "\n" ;

	### update last pos
    $last_pos = $split[1] ;
}
