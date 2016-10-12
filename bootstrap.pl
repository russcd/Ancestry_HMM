=begin

This script will produce bootstrapped estimates of the admixture 
time for a given set of samples. 

It requires the Parallel::ForkManager module, and the number of cores
to use may be specified below. 

Usage:

perl bootstrap.pl [panel_directory] < [list of panels to include] > [time estimates]

=cut

use strict ; 
use warnings ; 
use Parallel::ForkManager ;

my $a = 0.5 ; 			### ancestry proportion
my $bootstraps = 100 ; 
my $block_size = 5000 ;   	### block size to be defined in number of snps
my @blocks ; 

my $panel_dir = $ARGV[0] ;
my $boot_dir = "${panel_dir}/BOOTSTRAPS/" ;

mkdir "${panel_dir}/BOOTSTRAPS/" ; 

my @panels ; 
while (<STDIN>) { 
	chomp $_ ; 
	my ( $file, $ploidy ) = split ( /\t/, $_ ) ; 
	
	my @out ; 
	my $line = 0 ; 

	open IN, "<${panel_dir}/$file" ;
	while (<IN>) { 
		push @out, $_ ; 
		if ( $#out > $block_size ) { 
			open BLOCK, ">${boot_dir}/${file}.$line" ;
			print BLOCK @out ; 
			close BLOCK ; 
			push @blocks, "-i $ploidy ${boot_dir}/${file}.$line /dev/null " ; 
			$#out = -1 ; 
			$line ++ ; 
		}
	}
	open BLOCK, ">${boot_dir}/${file}.$line" ;
	print BLOCK @out ;
	close BLOCK ; 
	push @blocks, "-i $ploidy ${boot_dir}/${file}.$line /dev/null " ;
}

my @cmds ; 
foreach ( 1..$bootstraps ) {
	## sample from blocks with replacement
	my $c = "./ancestry_hmm -a $a --tolerance 1 " ;
	foreach ( 1..$#blocks ) { 
		$c .= $blocks[rand()*$#blocks] ;
	}
	push @cmds, $c ; 
}

my $manager = new Parallel::ForkManager( 32 );
foreach ( @cmds ) { 
	$manager->start and next ;

	my $cmd = `$_ 2>&1 | grep -e 'estimated_admixture time:' ` ; 
	$cmd =~ s/estimated_admixture time:\s+// ;
	$cmd =~ s/\t// ; 
	print $cmd ;

	$manager->finish;
}
$manager->wait_all_children;

system ("rm -r $boot_dir") ;
