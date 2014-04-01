#!/usr/bin/env perl

# Program: SolexaQA v.2.2
# Calculates quality statistics on Illumina FASTQ sequence files
# and creates visual representations of run quality
# Murray Cox, Patrick Biggs, Daniel Peterson and Mauro Truglio
# Massey University, New Zealand
# Email contact <m.p.cox@massey.ac.nz>
# April 2013

# Version 2.0: Complete rewrite of the SolexaQA code

# Version 2.1: Heat map now output as pdf. Resolves issues writing png files using R on some UNIX systems

# Version 2.2: Increment count to keep in concert with other programs in the package

# Released under GNU General Public License version 3

use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my $usage = "
$0 input_files [-p|probcutoff 0.05] [-h|phredcutoff 13] [-v|variance] [-m|minmax] [-s|sample 10000] [-b|bwa] [-d|directory path] [-sanger -solexa -illumina]\n
-p|probcutoff   probability value (between 0 and 1) at which base-calling error is considered too high (default; p = 0.05) *or*
-h|phredcutoff  Phred quality score (between 0 and 41) at which base-calling error is considered too high
-v|variance     calculate variance statistics
-m|minmax       calculate minimum and maximum error probabilities for each read position of each tile
-s|sample       number of sequences to be sampled per tile for statistics estimates (default; s = 10000)
-b|bwa          use BWA trimming algorithm
-d|directory    path to directory where output files are saved
-sanger         Sanger format (bypasses automatic format detection)
-solexa         Solexa format (bypasses automatic format detection)
-illumina       Illumina format (bypasses automatic format detection)
\n";

# global values
my $prob_cutoff;
my $phrd_cutoff;

my $variance = 0;
my $minmax   = 0;
my $sample   = 10000;

my $automatic_detection_lines = 10000;
my $sanger = 0;
my $solexa = 0;
my $illumina = 0;
my $format;
my $user_defined;

my $bwa;

my $directory;

my $is_hiseq_32;
my $is_hiseq_48;
my %hiseq_tile_32 = (
1 => 1,
2 => 2,
3 => 3,
4 => 4,
5 => 5,
6 => 6,
7 => 7,
8 => 8,
21 => 9,
22 => 10,
23 => 11,
24 => 12,
25 => 13,
26 => 14,
27 => 15,
28 => 16,
41 => 17,
42 => 18,
43 => 19,
44 => 20,
45 => 21,
46 => 22,
47 => 23,
48 => 24,
61 => 25,
62 => 26,
63 => 27,
64 => 28,
65 => 29,
66 => 30,
67 => 31,
68 => 32,
1101 => 1,
1102 => 2,
1103 => 3,
1104 => 4,
1105 => 5,
1106 => 6,
1107 => 7,
1108 => 8,
1201 => 9,
1202 => 10,
1203 => 11,
1204 => 12,
1205 => 13,
1206 => 14,
1207 => 15,
1208 => 16,
2101 => 17,
2102 => 18,
2103 => 19,
2104 => 20,
2105 => 21,
2106 => 22,
2107 => 23,
2108 => 24,
2201 => 25,
2202 => 26,
2203 => 27,
2204 => 28,
2205 => 29,
2206 => 30,
2207 => 31,
2208 => 32
);
my %hiseq_tile_48 = (
1101 => 1,
1102 => 2,
1103 => 3,
1104 => 4,
1105 => 5,
1106 => 6,
1107 => 7,
1108 => 8,
1201 => 9,
1202 => 10,
1203 => 11,
1204 => 12,
1205 => 13,
1206 => 14,
1207 => 15,
1208 => 16,
1301 => 17,
1302 => 18,
1303 => 19,
1304 => 20,
1305 => 21,
1306 => 22,
1307 => 23,
1308 => 24,
2101 => 25,
2102 => 26,
2103 => 27,
2104 => 28,
2105 => 29,
2106 => 30,
2107 => 31,
2108 => 32,
2201 => 33,
2202 => 34,
2203 => 35,
2204 => 36,
2205 => 37,
2206 => 38,
2207 => 39,
2208 => 40,
2301 => 41,
2302 => 42,
2303 => 43,
2304 => 44,
2305 => 45,
2306 => 46,
2307 => 47,
2308 => 48
);
my @hiseq_keys_32 = keys %hiseq_tile_32;
my @hiseq_keys_48 = keys %hiseq_tile_48;

# read user options
GetOptions(
"v|variance"      => \$variance,
"m|minmax"        => \$minmax,
"s|sample=i"      => \$sample,
"p|probcutoff=f"  => \$prob_cutoff,
"h|phredcutoff=f" => \$phrd_cutoff,
"b|bwa"           => \$bwa,
"d|directory=s"   => \$directory,
"sanger"          => \$sanger,
"solexa"          => \$solexa,
"illumina"        => \$illumina
);

# get user format (if supplied)
if( ($sanger && $solexa) || ($sanger && $illumina) || ($solexa && $illumina) ){
	die "Error: Please select only one of -sanger, -solexa or -illumina\n";
}

if( $sanger || $solexa || $illumina ){
	$user_defined = 1;
}

if( $sanger ){
	$format = "sanger";
}elsif( $solexa ){
	$format = "solexa";
}elsif( $illumina ){
	$format = "illumina";
}

# get files
my @files = @ARGV;

# check for presence of at least one input file
if( !$files[0] ){ die "$usage"; }

# check for appropriate sample size
if( $sample < 10000 ){
	print STDERR "Warning: Sample sizes less than 10,000 may lead to inaccurate estimates\n";
}

if( ( $sample > 10000 ) && $variance ){
	print STDERR "Warning: Running variance method with sample sizes greater than 10,000 may take a long time\n";
}

# check for user input values and convert Phred cutoff value to probability
if( defined( $prob_cutoff ) && defined( $phrd_cutoff ) ){
	die "Error: Please select p OR h cutoff, not both.\n";
    
}elsif( !defined( $prob_cutoff ) && !defined( $phrd_cutoff ) ){
	$prob_cutoff = 0.05;
    
}elsif( defined( $phrd_cutoff ) && $phrd_cutoff < 0 ){
	die "Error: Q cutoff must be greater than or equal to 0";
    
}elsif( defined( $prob_cutoff ) && ( $prob_cutoff < 0 || $prob_cutoff > 1 ) ){
	die "Error: P cutoff must be between 0 and 1";
    
}

# check for presence of R
if( !`which R 2> err.log` ){
	print STDERR "Warning: Subsidiary program R not found. Line graphs, histogram and heatmap will not be produced.\n";
}
`rm err.log`;


# loop through input files
FILE_LOOP: foreach my $input_file ( @files ){
	
	# initialize variables to store number of tiles and number of bases per read
	my $number_of_tiles;
	my $read_length;
	
	# just get filename, not full path (as returned by @ARGV)
	my @filepath = split( /\//, $input_file );
	my $filename = $filepath[$#filepath];
    
	# open input file for reading
	open( INPUT, "<$input_file" )
    or die "Error: Failure opening $input_file for reading: $!\n";
	
	# count number of lines in file
	my $line_counter = 0;
	while( <INPUT> && ($line_counter <= 4*$automatic_detection_lines) ){ # piPipes: most like will be able to determine
		$line_counter++;
	}
    
	# number of sequences
	my $number_of_sequences = $line_counter / 4;
	
	my $last_line = $line_counter;
	
	if( $last_line <= 1 ){
		print STDERR "Warning: File $input_file is empty\n";
		next FILE_LOOP;
	}
	
	# reset file pointer and line counter
	seek (INPUT, 0, 0);
	$line_counter = 0;
	
	# determine format
	if( !$user_defined ){
		$format = "";
	}
	if( !$format ){
		
		my $number_of_lines = $automatic_detection_lines;
		if( $automatic_detection_lines > $number_of_sequences ){
			$number_of_lines = $number_of_sequences;
		}
        
		$format = &get_format(*INPUT, $number_of_lines);
		if( !$format ){
			die "Error: File format cannot be determined\n";
		}
	}
	print STDOUT $format,"\n"; # piPipes: only need the format information
    exit;
	################# Once format detected, hashes are initialized ###########
	
	my %dict_q_to_Q;
	%dict_q_to_Q=&q_to_Q();
	
	my %dict_Q_to_p;
	%dict_Q_to_p=&Q_to_p();
    
	my %dict_Q_to_q;
	%dict_Q_to_q=&Q_to_q();
	#########################################################################
	
	
	if( $format eq "sanger" ){
		if( $user_defined ){
			print STDOUT "User defined format: Sanger FASTQ format\n";
		}else{
			print STDOUT "Automatic format detection: Sanger FASTQ format\n";
		}
	}elsif( $format eq "solexa" ){
		if( $user_defined ){
			print STDOUT "User defined format: Solexa FASTQ format, Illumina pipeline 1.2 or less\n";
		}else{
			print STDOUT "Automatic format detection: Solexa FASTQ format, Illumina pipeline 1.2 or less\n";
		}
	}elsif( $format eq "illumina" ){
		if( $user_defined ){
			print STDOUT "User defined format: Illumina FASTQ format, Illumina pipeline 1.3+\n";
		}else{
			print STDOUT "Automatic format detection: Illumina FASTQ format, Illumina pipeline 1.3+\n";
		}
	}
	
	# check for phred cutoff value
	if( defined( $phrd_cutoff ) ){
		$prob_cutoff = sprintf( "%.5f", $dict_Q_to_p{$phrd_cutoff} );
		print STDOUT "Info: Phred value of $phrd_cutoff converted to a probability value of $prob_cutoff\n";
	}
    
	# determine bwa trimming threshold
	my $threshold = 0;
	if( $bwa ){
		
		if( defined( $phrd_cutoff ) ){
			$threshold = $phrd_cutoff;
		}else{
			$threshold = &p_to_Q( $prob_cutoff );
		}
	}
    
	# go to fourth to last header line in file (last sequence ID line)
	my $line;
	while ( $line_counter < ( $last_line-3 ) ){
		$line = <INPUT>;
		$line_counter++;
	}
    
	# find tile number of last file entry (largest tile number)
	$number_of_tiles = 0;
	if( $line =~ /\S+\s\S+/ ){
		# Cassava 1.8 variant
		if( $line =~ /^@[\d\w\-\._]+:[\d\w]+:[\d\w\-]+:[\d\w]+:(\d+)/ ){
			$number_of_tiles = $1 + 1;
            # Sequence Read Archive variant
		}elsif( $line =~ /^@[\d\w\-\._\s]+:[\d\w]+:(\d+)/ ){
			$number_of_tiles = $1 + 1;
		}
        # All other variants
	}elsif( $line =~ /^@[\d\w\-:\._]*:+\d*:(\d*):[\.\d]+:[\.\/\#\d\w]+$/ ){
    $number_of_tiles = $1 + 1;
}

if( !$number_of_tiles ){
    print STDERR "Error: File $input_file does not match Solexa ID format (note: this may be caused by an incomplete final entry/empty terminal lines)\n";
    next FILE_LOOP;
}

# check if HiSeq data
$is_hiseq_32 = 0;
$is_hiseq_48 = 0;
if( $number_of_tiles == 69 || $number_of_tiles == 2209 ){
    $is_hiseq_32 = 1;
    $number_of_tiles = 33;
}elsif( $number_of_tiles == 2309 ){
    $is_hiseq_48 = 1;
    $number_of_tiles = 49;
}

# go to last sequence line in file and find read length
chomp( $line = <INPUT> );
$read_length = length ($line);

# set percentage of sequences to be used in stats calculations based on input sample size
my $percentage = ( $sample * ($number_of_tiles - 1) * 4 ) / $line_counter;
if( $percentage > 1 ){
    print STDERR "Warning: Desired sample size ($sample) greater than number of sequences per tile\n";
    $percentage = 1;
}

# open output quality file
my $quality_file;
if ( $directory ){
    # remove any trailing '/'
    $directory =~ s/\/\z//;
    my $file_name = $filename . ".quality";
    $quality_file = File::Spec->catpath( undef, $directory, $file_name );
}else{
    $quality_file = $filename . ".quality";
}
if( -e $quality_file ){
    die "Error: Quality file $quality_file already exists: $!\n";
}
open( QUALITY, ">$quality_file" )
or die "Error: Failure opening $quality_file for writing: $!\n";

# open output matrix file
my $matrix_file;
if ( $directory ){
    # remove any trailing '/'
    $directory =~ s/\/\z//;
    my $file_name = $filename . ".matrix";
    $matrix_file = File::Spec->catpath( undef, $directory, $file_name );
}else{
    $matrix_file = $filename . ".matrix";
}
if( -e $matrix_file ){
    die "error: matrix file $matrix_file already exists: $!\n";
}
open( MATRIX, ">$matrix_file" )
or die "Error: Failure opening $matrix_file for writing: $!\n";

# open segment table output file
my $segment_output_file;
if ( $directory ){
    # remove any trailing '/'
    $directory =~ s/\/\z//;
    my $file_name = $filename . ".segments";
    $segment_output_file = File::Spec->catpath( undef, $directory, $file_name );
}else{
    $segment_output_file = $filename . ".segments";
}
if( -e $segment_output_file ){
    die "Error: Segment output file $segment_output_file already exists: $!\n";
}
open( SEGMENT_OUTPUT, ">$segment_output_file" )
or die "Error: Failure opening $segment_output_file for writing: $!\n";

# ---------------------------------------------------

# create hashes to store statistics
my %mean_bytiles       = ();
my %mean_count_bytiles = ();
my %min_bytiles        = ();
my %max_bytiles        = ();
my %var_sum_bytiles    = ();
my %var_count_bytiles  = ();
my %trim_histogram     = ();

# initialize hashes to store arrays of probabilities
for( my $i = 0; $i < $number_of_tiles; $i++ ){
    
    my @new_array1;
    my @new_array2;
    my @new_array3;
    my @new_array4;
    my @new_array5;
    my @new_array6;
    
    for( my $j = 0; $j < $read_length; $j++ ){
		
        push( @new_array1, 0 );
        push( @new_array2, 0 );
        push( @new_array3, 1 );
        push( @new_array4, 0 );
        push( @new_array5, 0 );
        push( @new_array6, 0 );
    }
    
    $mean_bytiles{$i}        =  \@new_array1;
    $mean_count_bytiles{$i}  =  \@new_array2;
    $min_bytiles{$i}         =  \@new_array3;
    $max_bytiles{$i}         =  \@new_array4;
    $var_sum_bytiles{$i}     =  \@new_array5;
    $var_count_bytiles{$i}   =  \@new_array6;
}

for( my $j = 0; $j <= $read_length; $j++ ){
    $trim_histogram{$j} = 0;
}

# ---------------------------------------------------

# reset file pointer on input
seek (INPUT, 0, 0);
$line_counter = 0;

# step through sequences in input fastq file
while( $line = <INPUT> ){
    
    $line_counter++;
    
    # find quality score sequence ID lines
    next if ( ( $line_counter % 4 ) != 1 );
    
    # roll simulated dice to determine which sequences will be recorded for stat estimates
    if( rand() < $percentage ){
		
        # retrieve sequence information
        my @qual;
        my $tile_number;
        
        # regular expression to find tile number from Solexa fastq file sequence IDs
        $tile_number = 0;
        if( $line =~ /\S+\s\S+/ ){
			# Cassava 1.8 variant
            if( $line =~ /^@[\d\w\-\._]+:[\d\w]+:[\d\w\-]+:[\d\w]+:(\d+)/ ){
                $tile_number = $1;
                # Sequence Read Archive variant
			}elsif( $line =~ /^@[\d\w\-\._\s]+:[\d\w]+:(\d+)/ ){
				$tile_number = $1;
			}
			# All other variants
        }elsif( $line =~ /^@[\d\w\-:\._]*:+\d*:(\d*):[\.\d]+:[\.\/\#\d\w]+$/ ){
        $tile_number = $1;
    }
    
    if( !$tile_number ){
        die "Error: Lane ID at line number $line_counter not in correct Illumina format";
    }
    
    # change tile number if hiseq
    if( $is_hiseq_32 ){
        $tile_number = $hiseq_tile_32{$tile_number};
    }elsif( $is_hiseq_48 ){
        $tile_number = $hiseq_tile_48{$tile_number};
    }
    
    # go to quality scores line
    $line = <INPUT>;
    $line_counter++;
    
    $line = <INPUT>;
    $line_counter++;
    
    chomp( $line = <INPUT> );
    $line_counter++;
    
    # store quality scores as an array
    @qual = split(//, $line);
    
    # transform quality values from ASCII into Solexa format
    for( my $i = 0; $i < scalar @qual; $i++ ){
        
        $qual[$i] = $dict_q_to_Q{$qual[$i]};
    }
    
    # calculate probability values$dict_q_to_Q{$qual[$i]}
    my @prob;
    
    for( my $i = 0; $i < scalar @qual; $i++ ){
        
        $prob[$i] = $dict_Q_to_p{$qual[$i]};
    }
    
    # assign all probabilities to holders
    for( my $i = 0; $i < scalar @prob; $i++ ){
        
        $mean_bytiles{$tile_number}->[$i] += $prob[$i];
        $mean_count_bytiles{$tile_number}->[$i]++;
        
        # evaluate each probability against current min and max, store if better
        if( $minmax && ( $prob[$i] < $min_bytiles{$tile_number}->[$i] ) ){
            $min_bytiles{$tile_number}->[$i] = $prob[$i];
        }
        if( $minmax && ( $prob[$i] > $max_bytiles{$tile_number}->[$i] ) ){
            $max_bytiles{$tile_number}->[$i] = $prob[$i];
        }
        
        # tile 0 represents global probability
        $mean_bytiles{0}->[$i] += $prob[$i];
        $mean_count_bytiles{0}->[$i]++;
        
        # find global min and max
        if( $minmax && ( $prob[$i] < $min_bytiles{0}->[$i] ) ){
            $min_bytiles{0}->[$i] = $prob[$i];
        }
        if( $minmax && ( $prob[$i] > $max_bytiles{0}->[$i] ) ){
            $max_bytiles{0}->[$i] = $prob[$i];
        }
    }
    
    # calculate longest fragment
    my $longest_segment = 0;
    my $current_start   = 0;
    my $cutoff_hit      = 0;
    
    if( $bwa ){
        $longest_segment = &bwa_trim( $threshold, \@qual );
        
    }else{
        for( my $i = 0; $i < scalar @prob; $i++ ){
            
            if( $prob[$i] >= $prob_cutoff ){
                
                $cutoff_hit = 1;
                
                my $current_segment_length = $i - $current_start;
                $current_start = $i + 1;
                
                if( $current_segment_length > $longest_segment ){
                    $longest_segment = $current_segment_length;
                }
            }
        }
        if( !$cutoff_hit ){ $longest_segment = $read_length; }
    }
    
    # add to counter
    $trim_histogram{$longest_segment}++;
}
}

# calculate means and store in hashes
for( my $j = 0; $j < $read_length; $j++ ){
    for( my $i = 0; $i < $number_of_tiles; $i++){
        if( $mean_bytiles{$i}->[0] ){
            
            my $mean_pertile = $mean_bytiles{$i}->[$j] / $mean_count_bytiles{$i}->[$j];
            
            $mean_bytiles{$i}->[$j] = $mean_pertile;
        }
    }
}

if( $variance ){ # if calculating variance, loop through file again to find differences from mean
    
    # reset file pointer on input
    seek (INPUT, 0, 0);
    $line_counter = 0;
	
    # step through sequences in input fastq file
    while( $line = <INPUT> ){
        
        $line_counter++;
        
        # find quality score sequence ID lines
        next if ( ( $line_counter % 4 ) != 1 );
        
        # roll simulated dice to determine which sequences will be recorded for stat estimates
        if( rand() < $percentage ){
			
            # retrieve sequence information
            my @qual;
            my $tile_number;
            
            # regular expression to find tile number from Solexa fastq file sequence IDs
            $tile_number = 0;
            if( $line =~ /\S+\s\S+/ ){
				# Cassava 1.8 variant
                if( $line =~ /^@[\d\w\-\._]+:[\d\w]+:[\d\w\-]+:[\d\w]+:(\d+)/ ){
                    $tile_number = $1;
					# Sequence Read Archive variant
                }elsif( $line =~ /^@[\d\w\-\._\s]+:[\d\w]+:(\d+)/ ){
                    $tile_number = $1;
                }
				# All other variants
            }elsif( $line =~ /^@[\d\w\-:\._]*:+\d*:(\d*):[\.\d]+:[\.\/\#\d\w]+$/ ){
            $tile_number = $1;
        }
        
        if( !$tile_number ){
            die "Error: Lane ID at line number $line_counter not in correct Solexa format";
        }
        
        # change tile number if hiseq
        if( $is_hiseq_32 ){
            $tile_number = $hiseq_tile_32{$tile_number};
        }elsif( $is_hiseq_48 ){
            $tile_number = $hiseq_tile_48{$tile_number};
        }
        
        # go to quality scores lines
        $line = <INPUT>;
        $line_counter++;
        
        $line = <INPUT>;
        $line_counter++;
        
        chomp( $line = <INPUT> );
        $line_counter++;
		
        # store quality scores as an array
        @qual = split(//, $line);
        
        # transform quality values from ASCII into Solexa Q format
        for( my $i = 0; $i < scalar @qual; $i++ ){
            
            $qual[$i] = $dict_q_to_Q{$qual[$i]};
        }
        
        # calculate probability values
        my @prob;
        
        for( my $i = 0; $i < scalar @qual; $i++ ){
			
            $prob[$i] = $dict_Q_to_p{$qual[$i]};
        }
        
        # calculate sum of squared differences from the mean
        for( my $i = 0; $i < scalar @prob; $i++ ){
            
            my $diff_from_mean = $mean_bytiles{$tile_number}->[$i] - $prob[$i];
            
            $var_sum_bytiles{$tile_number}->[$i] += $diff_from_mean ** 2;
            $var_count_bytiles{$tile_number}->[$i]++;
            
            # tile 0 represents global probability
            $var_sum_bytiles{0}->[$i] += $diff_from_mean ** 2;
            $var_count_bytiles{0}->[$i]++;
        }
    }
}
}

# ---------------------------------------------------
# print output to quality file

# print header line of quality file
for( my $i = 0; $i < $number_of_tiles; $i++ ){
    
    if( $mean_count_bytiles{$i}->[0] ){
		
        # assign column ID (G if global, otherwise tile #)
        my $col_id;
        
        if( $i == 0 ){
            $col_id = "global";
        }else{
            $col_id = "tile_$i";
        }
        
        print QUALITY "mean_$col_id\t";
        
        if( $variance ){ print QUALITY "var_$col_id\t"; }
        
        if( $minmax ){ print QUALITY "min_$col_id\tmax_$col_id\t" };
        
    }
}
print QUALITY "\n";

# print stats for each read position for each tile to quality file
for( my $j = 0; $j < $read_length; $j++ ){
    
    for( my $i = 0; $i < $number_of_tiles; $i++){
        
        if( $mean_count_bytiles{$i}->[0] > 0 ){
            
            # print mean, variance, minimum, and maximum
            print QUALITY sprintf( "%.5f", $mean_bytiles{$i}->[$j] ), "\t";
            
            if( $variance ){
                my $var_pertile = $var_sum_bytiles{$i}->[$j] / $var_count_bytiles{$i}->[$j];
                print QUALITY sprintf( "%.5f", $var_pertile ), "\t";
            }
            
            if( $minmax ){
                
                my $min = sprintf( "%.5f", $min_bytiles{$i}->[$j] );
                my $max = sprintf( "%.5f", $max_bytiles{$i}->[$j] );
                
                print QUALITY $min, "\t", $max, "\t";
            }
        }
    }
    print QUALITY "\n";
}

# ---------------------------------------------------
# create matrix files for heat maps

# print header line with the column headers "Tile" and each read position
print MATRIX ".\t";

for( my $i = 1; $i < ( $read_length + 1 ); $i++){
    
    print MATRIX $i, "\t";
}
print MATRIX "\n";

# print stored mean probability values for each read position of each tile
for( my $i = 1; $i < $number_of_tiles; $i++ ){
    
    # check for zeroed tiles
    my $all_zeros = 1;
    for( my $j = 0; $j < $read_length; $j++ ){
        
        if( $mean_count_bytiles{$i}->[$j] > 0 ){
            $all_zeros = 0;
        }
    }
    next if $all_zeros;
    
    # change tile names if hiseq
    if( $is_hiseq_32 ){
        my $real_tile;
        foreach my $entry ( reverse(sort @hiseq_keys_32) ){
			
            if( $hiseq_tile_32{$entry} == $i && $entry < 69 ){
                $real_tile = $entry;
            }
        }
        print MATRIX "tile_", $real_tile, "\t";
    }elsif( $is_hiseq_48 ){
        my $real_tile;
        foreach my $entry ( reverse(sort @hiseq_keys_48) ){
			
            if( $hiseq_tile_48{$entry} == $i){
                $real_tile = $entry;
            }
        }
        print MATRIX "tile_", $real_tile, "\t";
    }else{
        # print tile # at the beginning of each line
        print MATRIX "tile_", $i, "\t";
    }
    
    for( my $j = 0; $j < $read_length; $j++ ){
        
        if( $mean_count_bytiles{$i}->[$j] > 0 ){
            
            print MATRIX sprintf( "%.5f", $mean_bytiles{$i}->[$j] ), "\t";
            
        }else{
            print MATRIX "\t";
        }
    }
    
    print MATRIX "\n";
}

# ---------------------------------------------------
# create table of read segments passing cutoff

my $segment_sum = 0;

for( my $j = 0; $j <= $read_length; $j++ ){
    $segment_sum += $trim_histogram{$j};
}

print SEGMENT_OUTPUT "read_length\tproportion_of_reads\n";

for( my $j = 0; $j <= $read_length; $j++ ){
    
    my $proportion = sprintf( "%.4f", ( $trim_histogram{$j} / $segment_sum ) );
    
    print SEGMENT_OUTPUT $j, "\t", $proportion, "\n";
}

# ---------------------------------------------------
# close files

close QUALITY
or die "Error: Cannot close quality file $quality_file: $!";

close MATRIX
or die "Error: Cannot close matrix file $matrix_file: $!";

close INPUT
or die "Error: Cannot clost input file $input_file: $!";

close SEGMENT_OUTPUT
or die "Error: Cannot close segment table file $segment_output_file: $!";






# ---------------------------------------------------
# create heat map from matrix file


open (MYFILE, ">temp.R");
print MYFILE "library(package=grid)
lo = function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, treeheight_col, treeheight_row, legend, annotation, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, ...){
    # Get height of colnames and length of rownames
    if(!is.null(coln[1])){
	    longest_coln = which.max(nchar(coln))
	    gp = list(fontsize = fontsize_col, ...)
	    coln_height = unit(1.1, \"grobheight\", textGrob(coln[longest_coln], rot = 90, gp = do.call(gpar, gp)))
            }
    else{
	    coln_height = unit(5, \"bigpts\")
    }
    
    if(!is.null(rown[1])){
	    longest_rown = which.max(nchar(rown))
	    gp = list(fontsize = fontsize_row, ...)
	    rown_width = unit(1.2, \"grobwidth\", textGrob(rown[longest_rown], gp = do.call(gpar, gp)))
            }
    else{
	    rown_width = unit(5, \"bigpts\")
    }
    
    gp = list(fontsize = fontsize, ...)
    # Legend position
    if(!is.na(legend[1])){
	    longest_break = which.max(nchar(names(legend)))
	    longest_break = unit(1.1, \"grobwidth\", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
            title_length = unit(1.1, \"grobwidth\", textGrob(\"Scale\", gp = gpar(fontface = \"bold\", ...)))
            legend_width = unit(12, \"bigpts\") + longest_break * 1.2
            legend_width = max(title_length, legend_width)
            }
    else{
	    legend_width = unit(0, \"bigpts\")
    }
    
    # Set main title height
    if(is.na(main)){
	    main_height = unit(0, \"npc\")
    }
    else{
	    main_height = unit(1.5, \"grobheight\", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }
    
    # Column annotations
    if(!is.na(annotation[[1]][1])){
	    # Column annotation height
	    annot_height = unit(ncol(annotation) * (8 + 2) + 2, \"bigpts\")
	    # Width of the correponding legend
	    longest_ann = which.max(nchar(as.matrix(annotation)))
	    annot_legend_width = unit(1.2, \"grobwidth\", textGrob(as.matrix(annotation)[longest_ann], gp = gpar(...))) + unit(12, \"bigpts\")
	    if(!annotation_legend){
            annot_legend_width = unit(0, \"npc\")
	    }
    }
    else{
	    annot_height = unit(0, \"bigpts\")
	    annot_legend_width = unit(0, \"bigpts\")
    }
    
    # Tree height
    treeheight_col = unit(treeheight_col, \"bigpts\") + unit(5, \"bigpts\")
    treeheight_row = unit(treeheight_row, \"bigpts\") + unit(5, \"bigpts\")
    
    # Set cell sizes
    if(is.na(cellwidth)){
	    matwidth = unit(1, \"npc\") - rown_width - legend_width - treeheight_row - annot_legend_width
    }
    else{
	    matwidth = unit(cellwidth * ncol, \"bigpts\")
    }
    
    if(is.na(cellheight)){
	    matheight = unit(1, \"npc\") - main_height - coln_height - treeheight_col - annot_height
    }
    else{
	    matheight = unit(cellheight * nrow, \"bigpts\")
    }
    
    
    # Produce layout()
    pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 5, widths = unit.c(treeheight_row, matwidth, rown_width, legend_width, annot_legend_width), heights = unit.c(main_height, treeheight_col, annot_height, matheight, coln_height)), gp = do.call(gpar, gp)))
        
        # Get cell dimensions
        pushViewport(vplayout(4, 2))
        cellwidth = convertWidth(unit(0:1, \"npc\"), \"bigpts\", valueOnly = T)[2] / ncol
        cellheight = convertHeight(unit(0:1, \"npc\"), \"bigpts\", valueOnly = T)[2] / nrow
        upViewport()
        
        # Return minimal cell dimension in bigpts to decide if borders are drawn
        mindim = min(cellwidth, cellheight)
        return(mindim)
        }

draw_dendrogram = function(hc, horizontal = T){
    h = hc\$height / max(hc\$height) / 1.05
    m = hc\$merge
    o = hc\$order
    n = length(o)
    
    m[m > 0] = n + m[m > 0]
    m[m < 0] = abs(m[m < 0])
    
    dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c(\"x\", \"y\")))
    dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
    
    for(i in 1:nrow(m)){
	    dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
	    dist[n + i, 2] = h[i]
    }
    
    draw_connection = function(x1, x2, y1, y2, y){
	    grid.lines(x = c(x1, x1), y = c(y1, y))
	    grid.lines(x = c(x2, x2), y = c(y2, y))
	    grid.lines(x = c(x1, x2), y = c(y, y))
    }
    
    if(horizontal){
	    for(i in 1:nrow(m)){
            draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
	    }
    }
    
    else{
	    gr = rectGrob()
	    pushViewport(viewport(height = unit(1, \"grobwidth\", gr), width = unit(1, \"grobheight\", gr), angle = 90))
	    dist[, 1] = 1 - dist[, 1]
	    for(i in 1:nrow(m)){
            draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
	    }
	    upViewport()
    }
}

draw_matrix = function(matrix, border_color, fmat, fontsize_number){
    n = nrow(matrix)
    m = ncol(matrix)
    x = (1:m)/m - 1/2/m
    y = 1 - ((1:n)/n - 1/2/n)
    for(i in 1:m){
	    grid.rect(x = x[i], y = y[1:n], width = 1/m, height = 1/n, gp = gpar(fill = matrix[,i], col = border_color))
	    if(attr(fmat, \"draw\")){
            grid.text(x = x[i], y = y[1:n], label = fmat[, i], gp = gpar(col = \"grey30\", fontsize = fontsize_number))
	    }
    }
}

draw_colnames = function(coln, ...){
    m = length(coln)
    x = (1:m)/m - 1/2/m
    grid.text(coln, x = x, y = unit(0.96, \"npc\"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))
}

draw_rownames = function(rown, ...){
    n = length(rown)
    y = 1 - ((1:n)/n - 1/2/n)
    grid.text(rown, x = unit(0.04, \"npc\"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))
}

draw_legend = function(color, breaks, legend, ...){
    height = min(unit(1, \"npc\"), unit(150, \"bigpts\"))
    pushViewport(viewport(x = 0, y = unit(1, \"npc\"), just = c(0, 1), height = height))
    legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
    breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
    h = breaks[-1] - breaks[-length(breaks)]
    grid.rect(x = 0, y = breaks[-length(breaks)], width = unit(10, \"bigpts\"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = \"#FFFFFF00\"))
    grid.text(names(legend), x = unit(12, \"bigpts\"), y = legend_pos, hjust = 0, gp = gpar(...))
    upViewport()
}

convert_annotations = function(annotation, annotation_colors){
    new = annotation
    for(i in 1:ncol(annotation)){
	    a = annotation[, i]
	    b = annotation_colors[[colnames(annotation)[i]]]
	    if(is.character(a) | is.factor(a)){
            a = as.character(a)
            if(length(setdiff(a, names(b))) > 0){
                stop(sprintf(\"Factor levels on variable %s do not match with annotation_colors\", colnames(annotation)[i]))
                    }
            new[, i] = b[a]
	    }
	    else{
            a = cut(a, breaks = 100)
            new[, i] = colorRampPalette(b)(100)[a]
	    }
    }
    return(as.matrix(new))
}

draw_annotations = function(converted_annotations, border_color){
    n = ncol(converted_annotations)
    m = nrow(converted_annotations)
    x = (1:m)/m - 1/2/m
    y = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
    for(i in 1:m){
	    grid.rect(x = x[i], unit(y[1:n], \"bigpts\"), width = 1/m, height = unit(8, \"bigpts\"), gp = gpar(fill = converted_annotations[i, ], col = border_color))
    }
}

draw_annotation_legend = function(annotation, annotation_colors, border_color, ...){
    y = unit(1, \"npc\")
    text_height = unit(1, \"grobheight\", textGrob(\"FGH\", gp = gpar(...)))
    for(i in names(annotation_colors)){
	    grid.text(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = \"bold\", ...))
	    y = y - 1.5 * text_height
	    if(is.character(annotation[, i]) | is.factor(annotation[, i])){
            for(j in 1:length(annotation_colors[[i]])){
                grid.rect(x = unit(0, \"npc\"), y = y, hjust = 0, vjust = 1, height = text_height, width = text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]][j]))
                grid.text(names(annotation_colors[[i]])[j], x = text_height * 1.3, y = y, hjust = 0, vjust = 1, gp = gpar(...))
                y = y - 1.5 * text_height
            }
	    }
	    else{
            yy = y - 4 * text_height + seq(0, 1, 0.02) * 4 * text_height
            h = 4 * text_height * 0.02
            grid.rect(x = unit(0, \"npc\"), y = yy, hjust = 0, vjust = 1, height = h, width = text_height, gp = gpar(col = \"#FFFFFF00\", fill = colorRampPalette(annotation_colors[[i]])(50)))
            txt = rev(range(grid.pretty(range(annotation[, i], na.rm = TRUE))))
            yy = y - c(0, 3) * text_height
            grid.text(txt, x = text_height * 1.3, y = yy, hjust = 0, vjust = 1, gp = gpar(...))
            y = y - 4.5 * text_height
	    }
	    y = y - 1.5 * text_height
    }
}
draw_main = function(text, ...){
    grid.text(text, gp = gpar(fontface = \"bold\", ...))
}

vplayout = function(x, y){
    return(viewport(layout.pos.row = x, layout.pos.col = y))
}

heatmap_motor = function(matrix, border_color, cellwidth, cellheight, tree_col, tree_row, treeheight_col, treeheight_row, filename, width, height, breaks, color, legend, annotation, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, fmat, fontsize_number, ...){
    grid.newpage()
    
    # Set layout
    mindim = lo(coln = colnames(matrix), rown = rownames(matrix), nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col,  ...)
    
    if(!is.na(filename)){
	    pushViewport(vplayout(1:5, 1:5))
	    
	    if(is.na(height)){
            height = convertHeight(unit(0:1, \"npc\"), \"inches\", valueOnly = T)[2]
	    }
	    if(is.na(width)){
            width = convertWidth(unit(0:1, \"npc\"), \"inches\", valueOnly = T)[2]
	    }
	    
	    # Get file type
	    r = regexpr(\"\\\\.[a-zA-Z]*\$\", filename)
	    if(r == -1) stop(\"Improper filename\")
            ending = substr(filename, r + 1, r + attr(r, \"match.length\"))
            
            f = switch(ending,
            pdf = function(x, ...) pdf(x, ...),
            png = function(x, ...) png(x, units = \"in\", res = 300, ...),
            jpeg = function(x, ...) jpeg(x, units = \"in\", res = 300, ...),
            jpg = function(x, ...) jpeg(x, units = \"in\", res = 300, ...),
            tiff = function(x, ...) tiff(x, units = \"in\", res = 300, compression = \"lzw\", ...),
            bmp = function(x, ...) bmp(x, units = \"in\", res = 300, ...),
            stop(\"File type should be: pdf, png, bmp, jpg, tiff\")
            )
            
            # print(sprintf(\"height:%f width:%f\", height, width))
            f(filename, height = height, width = width)
            heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, ...)
            garbage<-dev.off()
            upViewport()
            return()
            }
    
    # Omit border color if cell size is too small
    if(mindim < 3) border_color = NA
        
        # Draw title
        if(!is.na(main)){
            pushViewport(vplayout(1, 2))
            draw_main(main, fontsize = 1.3 * fontsize, ...)
            upViewport()
        }
    
    # Draw tree for the columns
    if(!is.na(tree_col[[1]][1]) & treeheight_col != 0){
	    pushViewport(vplayout(2, 2))
	    draw_dendrogram(tree_col, horizontal = T)
	    upViewport()
    }
    
    # Draw tree for the rows
    if(!is.na(tree_row[[1]][1]) & treeheight_row != 0){
	    pushViewport(vplayout(4, 1))
	    draw_dendrogram(tree_row, horizontal = F)
	    upViewport()
    }
    
    # Draw matrix
    pushViewport(vplayout(4, 2))
    draw_matrix(matrix, border_color, fmat, fontsize_number)
    upViewport()
    
    # Draw colnames
    if(length(colnames(matrix)) != 0){
	    pushViewport(vplayout(5, 2))
	    pars = list(colnames(matrix), fontsize = fontsize_col, ...)
	    do.call(draw_colnames, pars)
            upViewport()
            }
    
    # Draw rownames
    if(length(rownames(matrix)) != 0){
	    pushViewport(vplayout(4, 3))
	    pars = list(rownames(matrix), fontsize = fontsize_row, ...)
	    do.call(draw_rownames, pars)
            upViewport()
            }
    
    # Draw annotation tracks
    if(!is.na(annotation[[1]][1])){
	    pushViewport(vplayout(3, 2))
	    converted_annotation = convert_annotations(annotation, annotation_colors)
	    draw_annotations(converted_annotation, border_color)
	    upViewport()
    }
    
    # Draw annotation legend
    if(!is.na(annotation[[1]][1]) & annotation_legend){
	    if(length(rownames(matrix)) != 0){
            pushViewport(vplayout(4:5, 5))
	    }
	    else{
            pushViewport(vplayout(3:5, 5))
	    }
	    draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)
	    upViewport()
    }
    
    # Draw legend
    if(!is.na(legend[1])){
	    length(colnames(matrix))
	    if(length(rownames(matrix)) != 0){
            pushViewport(vplayout(4:5, 4))
	    }
	    else{
            pushViewport(vplayout(3:5, 4))
	    }
	    draw_legend(color, breaks, legend, fontsize = fontsize, ...)
	    upViewport()
    }
    
    
}

generate_breaks = function(x, n){
    seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
    return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

scale_colours = function(mat, col = rainbow(10), breaks = NA){
    mat = as.matrix(mat)
    return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

cluster_mat = function(mat, distance, method){
    if(!(method %in% c(\"ward\", \"single\", \"complete\", \"average\", \"mcquitty\", \"median\", \"centroid\"))){
	    stop(\"clustering method has to one form the list: 'ward', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.\")
    }
    if(!(distance[1] %in% c(\"correlation\", \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\", \"minkowski\")) & class(distance) != \"dist\"){
	    print(!(distance[1] %in% c(\"correlation\", \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\", \"minkowski\")) | class(distance) != \"dist\")
	    stop(\"distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'\")
    }
    if(distance[1] == \"correlation\"){
	    d = dist(1 - cor(t(mat)))
    }
    else{
	    if(class(distance) == \"dist\"){
            d = distance
	    }
	    else{
            d = dist(mat, method = distance)
	    }
    }
    
    return(hclust(d, method = method))
}

scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

scale_mat = function(mat, scale){
    if(!(scale %in% c(\"none\", \"row\", \"column\"))){
	    stop(\"scale argument shoud take values: 'none', 'row' or 'column'\")
    }
    mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
    return(mat)
}

generate_annotation_colours = function(annotation, annotation_colors, drop){
    if(is.na(annotation_colors)[[1]][1]){
	    annotation_colors = list()
    }
    count = 0
    for(i in 1:ncol(annotation)){
	    if(is.character(annotation[, i]) | is.factor(annotation[, i])){
            if (is.factor(annotation[, i]) & !drop){
                count = count + length(levels(annotation[, i]))
            }
            count = count + length(unique(annotation[, i]))
	    }
    }
    
    factor_colors = hsv((seq(0, 1, length.out = count + 1)[-1] +
    0.2)%%1, 0.7, 0.95)
    
    set.seed(3453)
    
    for(i in 1:ncol(annotation)){
	    if(!(colnames(annotation)[i] %in% names(annotation_colors))){
            if(is.character(annotation[, i]) | is.factor(annotation[, i])){
                n = length(unique(annotation[, i]))
                if (is.factor(annotation[, i]) & !drop){
                    n = length(levels(annotation[, i]))
                }
                ind = sample(1:length(factor_colors), n)
                annotation_colors[[colnames(annotation)[i]]] = factor_colors[ind]
                l = levels(as.factor(annotation[, i]))
                l = l[l %in% unique(annotation[, i])]
                if (is.factor(annotation[, i]) & !drop){
                    l = levels(annotation[, i])
                }
                names(annotation_colors[[colnames(annotation)[i]]]) = l
                factor_colors = factor_colors[-ind]
            }
            else{
                r = runif(1)
                annotation_colors[[colnames(annotation)[i]]] = hsv(r, c(0.1, 1), 1)
            }
	    }
    }
    return(annotation_colors)
}

kmeans_pheatmap = function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...){
    # Filter data
    if(!is.na(sd_limit)){
	    s = apply(mat, 1, sd)
	    mat = mat[s > sd_limit, ]
    }
    
    # Cluster data
    set.seed(1245678)
    km = kmeans(mat, k, iter.max = 100)
    mat2 = km\$centers
    
    # Compose rownames
    t = table(km\$cluster)
    rownames(mat2) = sprintf(\"cl%s_size_%d\", names(t), t)
    
    # Draw heatmap
    pheatmap(mat2, ...)
}

pheatmap = function(mat, color = colorRampPalette(rev(c(\"#D73027\", \"#FC8D59\", \"#FEE090\", \"#FFFFBF\", \"#E0F3F8\", \"#91BFDB\", \"#4575B4\")))(100), kmeans_k = NA, breaks = NA, border_color = \"grey60\", cellwidth = NA, cellheight = NA, scale = \"none\", cluster_rows = TRUE, cluster_cols = TRUE, clustering_distance_rows = \"euclidean\", clustering_distance_cols = \"euclidean\", clustering_method = \"complete\",  treeheight_row = ifelse(cluster_rows, 50, 0), treeheight_col = ifelse(cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA, legend_labels = NA, annotation = NA, annotation_colors = NA, annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = \"%.2f\", fontsize_number = 0.8 * fontsize, filename = NA, width = NA, height = NA, ...){

# Preprocess matrix
mat = as.matrix(mat)
mat = scale_mat(mat, scale)

# Kmeans
if(!is.na(kmeans_k)){
    # Cluster data
    km = kmeans(mat, kmeans_k, iter.max = 100)
    mat = km\$centers
    
    # Compose rownames
    t = table(km\$cluster)
    rownames(mat) = sprintf(\"cl%s_size_%d\", names(t), t)
}
else{
    km = NA
}

# Do clustering
if(cluster_rows){
    tree_row = cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)
    mat = mat[tree_row\$order, ]
}
else{
    tree_row = NA
    treeheight_row = 0
}

if(cluster_cols){
    tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)
    mat = mat[, tree_col\$order]
}
else{
    tree_col = NA
    treeheight_col = 0
}

# Format numbers to be displayed in cells
if(display_numbers){
    fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
    attr(fmat, \"draw\") = TRUE
}
else{
    fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
    attr(fmat, \"draw\") = FALSE
}


# Colors and scales
if(!is.na(legend_breaks[1]) & !is.na(legend_labels[1])){
    if(length(legend_breaks) != length(legend_labels)){
        stop(\"Lengths of legend_breaks and legend_labels must be the same\")
    }
}


if(is.na(breaks[1])){
    breaks = generate_breaks(as.vector(mat), length(color))
}
if (legend & is.na(legend_breaks[1])) {
    legend = grid.pretty(range(as.vector(breaks)))
    names(legend) = legend
}
else if(legend & !is.na(legend_breaks[1])){
    legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
    
    if(!is.na(legend_labels[1])){
        legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
        names(legend) = legend_labels
    }
    else{
        names(legend) = legend
    }
}
else {
    legend = NA
}
mat = scale_colours(mat, col = color, breaks = breaks)

# Preparing annotation colors
if(!is.na(annotation[[1]][1])){
    annotation = annotation[colnames(mat), , drop = F]
    annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
}

if(!show_rownames){
    rownames(mat) = NULL
}

if(!show_colnames){
    colnames(mat) = NULL
}

# Draw heatmap
heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, ...)
garbage<-dev.off()
invisible(list(tree_row = tree_row, tree_col = tree_col, kmeans = km))

}







filename=\"$matrix_file\"
d <- read.table(filename, header=TRUE, sep=\"\", check.names=FALSE)
#d <- d[order(d\$.),]
row.names(d) <- d\$.
d <- d[,2:ncol(d)]
dmat<-data.matrix(d)

mycol=colorRampPalette(c(\"white\",\"yellow\",\"#FF9B0D\",\"#BD2400\",\"#631300\",\"black\"))
pheatmap(dmat, cellheight=10, cellwidth=10, cluster_rows=FALSE, filename=paste(filename, '.pdf', sep=''), cluster_cols=FALSE, legend=TRUE, main=tail(strsplit(strsplit(filename, '.fastq.matrix')[[1]][1],'/')[[1]],1), color=mycol(75), breaks=c(seq(0,0.75,0.01)))

";

close (MYFILE);

system("Rscript temp.R");
system("rm temp.R");
system("rm Rplots.pdf");




# ----------------------------------------------------
# print line graph of mean probabilities for each tile at each read position from quality file

# calculate number of columns in quality file that do not contain mean probability data
my $extra_columns = 0;

if( $variance ){
    $extra_columns += 1;
}
if( $minmax ){
    $extra_columns += 2;
}

# define R script to graph quality file data



open (MYFILE, ">temp.R");
print MYFILE "library(package=grid)

filename=\"$quality_file\"
extra_columns=$extra_columns
d <- read.table(filename, header=TRUE)

pdf( paste(filename,'.pdf', sep = '') )
par(mar=c(5,6,6,2) + 0.1, oma=c(3,0,0,0), mgp = c(4, 1, 0))

varx=1:length(d[,1])
vary=d[,1]
plot(x=varx, y=vary, col='red', ylim=c(0,0.631), type='l', lwd=3, las=1, xlab='Position along read',col.lab=rgb(0,0.5,0), ylab='Mean probability of error', col.lab=rgb(0,0.5,0))
abline(h=seq(0.1,0.6,0.1), lty=1, col='grey')

legend(0.59,c(\"Error across all tiles\", \"Individual base call error per tile\"), cex=0.7, col=c('red', 'black'), lty=c(1,3), lwd=c(3,1))
Title=tail(strsplit(strsplit(filename, '.fastq.quality')[[1]][1],'/')[[1]],1)
title(paste(\'Sample: \', Title))
mtext(\"Average base call error per nucleotide position\", 3, line=1)

for(i in c(1:length(d[1,]))){
    
    if( ( i + extra_columns ) %% ( 1 + extra_columns ) == 0 ){
        
        lines(1:length(d[,i]), d[,i], lty=3, col='black')
    }
}

lines(1:length(d[,1]), d[,1], type='l', col='red')
garbage<-dev.off()



";

close (MYFILE);

system("Rscript temp.R");
system("rm temp.R");




# ----------------------------------------------------
# print histogram of read segments passing cutoff

if ( $directory ){
    # remove any trailing '/'
    $directory =~ s/\/\z//;
    my $file_name = $filename;
    $segment_output_file = File::Spec->catpath( undef, $directory, $file_name );
}else{
    $segment_output_file = $filename;
}

# define R script to graph segment table

open (MYFILE, ">temp.R");
print MYFILE "

filename <- \"$segment_output_file.segments\"
cutoff <- $prob_cutoff
output <- \"$segment_output_file.hist\"

d <- read.table(filename, header=T)
maxup=max(d[,2])
pdf( paste(output,'.pdf', sep = ''), width = 11 )
if (maxup>0.45){
    # I want to plot the lower values up to 55, then a split to 95 for the
    # last top. This should make it clear which is the highest, without
    # drowning out the other data.
    
    # I want the split to be approx 5% of the scale,
    
    # as I am to plot the ranges 0 - 55 and 95 - 140, in total 10 decades,
    lower=c(0,0.4)
    upper=c(maxup-0.05, maxup)
    # This is 10 decades. I multiply that with 2 and add 5% and get 21 units on the outer
    # Y axis:
    y_outer=(lower[2]+upper[1]-upper[2])*100
    lowspan=c(0,(2*y_outer/3)-1)
    topspan=c(lowspan[2]+5, y_outer)
    
    
    cnvrt.coords <-function(x,y=NULL){
	    # Stolen from the teachingDemos library, simplified for this use case
	    xy <- xy.coords(x,y, recycle=TRUE)
	    cusr <- par('usr')
	    cplt <- par('plt')
	    plt <- list()
	    plt\$x <- (xy\$x-cusr[1])/(cusr[2]-cusr[1])
	    plt\$y <- (xy\$y-cusr[3])/(cusr[4]-cusr[3])
	    fig <- list()
	    fig\$x <- plt\$x*(cplt[2]-cplt[1])+cplt[1]
	    fig\$y <- plt\$y*(cplt[4]-cplt[3])+cplt[3]
	    return( list(fig=fig) )
    }
    
    subplot <- function(fun, x, y=NULL){
	    # Stolen from the teachingDemos l	ibrary, simplified for this use case
	    old.par <- par(no.readonly=TRUE)
	    on.exit(par(old.par))
	    xy <- xy.coords(x,y)
	    xy <- cnvrt.coords(xy)\$fig
	    par(plt=c(xy\$x,xy\$y), new=TRUE)
	    fun
	    tmp.par <- par(no.readonly=TRUE)
	    return(invisible(tmp.par))
    }
    
    ##############################################
    #
    #
    # The main program starts here:
    #
    #
    
    # Setting up an outer wireframe for the plots.
    par(mar=c(8,6,6,3) + 0.1, oma=c(0,0,0,0), mgp = c(3, 1, 0))
    plot(c(0,1),c(0,y_outer),type='n',axes=FALSE,xlab=paste(\'Length of longest contiguous read segments with quality higher than\',cutoff),col.lab=rgb(0,0.5,0), ylab='Proportion of reads', col.lab=rgb(0,0.5,0))
    Title=tail(strsplit(strsplit(filename, '.fastq.segments')[[1]][1],'/')[[1]],1)
    title(paste(\'Sample: \',Title))
    mtext(paste(\"p cutoff = \",cutoff), 3, line=1)
    mtext(\"Sum of the segments = 1\", 1, line=6)
    # Plotting the lower range in the lower 11/21 of the plot.
    # xpd=FALSE to clip the bars
    tmp<-subplot(barplot(d[,2],col=\'blue\',ylim=lower,xpd=FALSE,las=0, names.arg = d[,1], las=1), x=c(0,1),y=lowspan)
    op <- par(no.readonly=TRUE)
    par(tmp)
    abline(h=seq(0.0,0.4,0.1), lty=1, col=\'grey\')
    par(op)
    subplot(barplot(d[,2],col=\'blue\',ylim=lower,xpd=FALSE,las=0, names.arg = d[,1], las=1), x=c(0,1),y=lowspan)
    
    
    # Plotting the upper range in the upper 9/21 of the plot, 1/21 left to
    # the split. Again xpd=FALSE, names.arg is set up to avoid having
    # the names plotted here, must be some easier way to do this but
    # this works
	
    tmp<-subplot(barplot(d[,2],col=\'blue\',ylim=c(round(upper[1], digits = 1),round(upper[2], digits = 1)+0.1), xpd=FALSE, las=1), x=c(0,1),y=topspan)
    op <- par(no.readonly=TRUE)
    par(tmp)
    abline(h=seq(round(upper[1], digits = 1),round(upper[2], digits = 1)+0.1,0.1), lty=1, col=\'grey\')
    par(op)
    subplot(barplot(d[,2],col=\'blue\',ylim=c(round(upper[1], digits = 1),round(upper[2], digits = 1)+0.1), xpd=FALSE, las=1), x=c(0,1),y=topspan)
	
    # Legend. An annoiance is that the colors comes in the opposite
    # order than in the plot.
    
    #legend(0.05,26,c(\'Reads\'), cex=0.7, col=c(\'blue\'), pch=15)
    
    # so far so good. (Just run the upper part to see the result so far)
    # Just want to make the ends of the axes a bit nicer.
    # All the following plots are in units of the outer coordinate system
    
    lowertop=lowspan[2]+(topspan[1]-lowspan[2])/2  # Where to end the lower axis
    breakheight=1   # Height of the break
    upperbot=lowertop+breakheight#(lowspan[2]+(topspan[1]-lowspan[2])/2)+breakheight# Where to start the upper axes
    markerheight=0.5 # Heightdifference for the break markers
    markerwidth=.03  # With of the break markers
    
    # Draw the break markers:
    #lines(c(0,0),c(1,lowertop))
    lines(c(markerwidth/-2,markerwidth/2),c(lowertop-markerheight/2,lowertop+markerheight/2))
    #lines(c(0,0),c(upperbot-breakheight,14))
    #lines(c(0,0),c(upperbot,maxup))
    lines(c(markerwidth/-2,markerwidth/2),c(upperbot-markerheight/2,upperbot+markerheight/2))
	
}else{
	
    par(mar=c(8,6,6,3) + 0.1, oma=c(0,0,0,0), mgp = c(3, 1, 0))
    barplot(d[,2], names.arg = d[,1], space = 0, ylim=c(0, 0.4), col=\'blue\', las=1, xlab=paste(\'Length of longest contiguous read segments with quality higher than\',cutoff),col.lab=rgb(0,0.5,0), ylab=\'Proportion of reads\', col.lab=rgb(0,0.5,0), axis.lty = 1, cex.names = 0.9 )
    abline(h=seq(0.1,0.4,0.1), lty=1, col=\'grey\')
    barplot(d[,2], add=TRUE, names.arg = d[,1], space = 0, ylim=c(0, 0.4), col=\'blue\', las=1, xlab=paste(\'Length of longest contiguous read segments with quality higher than\',cutoff),col.lab=rgb(0,0.5,0), ylab=\'Proportion of reads\', col.lab=rgb(0,0.5,0), axis.lty = 1, cex.names = 0.9 )
	
    Title=tail(strsplit(strsplit(filename, '.fastq.segments')[[1]][1],'/')[[1]],1)
    title(paste(\'Sample: \',Title))
    mtext(paste(\"p cutoff = \",cutoff), 3, line=1)
    mtext(\"Sum of the segments = 1\", 1, line=6)
    
    #legend(0.05,0.35,c(\'Reads\'), cex=0.7, col=c(\'blue\'), pch=15)
}
garbage<-dev.off()
";

close (MYFILE);

system("Rscript temp.R");
system("rm temp.R");






open (MYFILE, ">temp.R");
print MYFILE "
filename <- \"$segment_output_file.segments\"
output <- \"$segment_output_file.cumulative\"
cutoff <- $prob_cutoff

d<-read.table(filename, header=TRUE)
d\$sumcol=d\$proportion_of_reads

#d\$sumcol=d\$sumcol+d\$proportion_of_reads
sum=0

d\$ideal=1
d\$ideal[nrow(d)]=0
for(i in 1:nrow(d)){
    sum=sum+d\$proportion_of_reads[i]
    d\$sumcol[i]=1-sum}

pdf( paste(output,'.pdf', sep = ''), width = 11 )
par(mar=c(4.5,6,6,2) + 0.1, oma=c(0,0,0,0), mgp = c(3, 1, 0))
plot(x=d\$read_length, y=d\$sumcol, panel.first = abline(h=seq(0.0,1,0.2), lty=1, col='grey'),type='l', xlab=paste(\'Length of longest contiguous read segments with quality higher than\',cutoff),col.lab=rgb(0,0.5,0), ylab='Proportion of reads',col.lab=rgb(0,0.5,0), lwd=2, col='purple', las=1)
lines(d\$read_length,d\$ideal, type='l', lty=3, col='black', lwd=2)


legend(0.77,c(\"Reads\", \"Ideal\"), lty=c(1,3), cex=0.7, col=c('purple', 'black'), lwd=2)
Title=tail(strsplit(strsplit(filename, '.fastq.segments')[[1]][1],'/')[[1]],1)
title(paste(\'Sample: \',Title))
mtext(paste('p cutoff = ',cutoff), 3, line=1)
garbage<-dev.off()
";

close (MYFILE);
system("Rscript temp.R");
system("rm temp.R");

#&print_lookup_table(\%dict_Q_to_p, \%dict_Q_to_q);

print STDOUT $input_file, " QA completed\n";
}


# terminate
exit 0 or die "Error: Program $0 ended abnormally: $!\n";

# ----------------------------------------------------



# Change ASCII character to Phred/Solexa quality score
sub q_to_Q(){
    
    if( $format eq "sanger" ){
		my $num;
		my %dict_q_to_Q=();
		for($num=28; $num<=126; $num++){
			$dict_q_to_Q{chr($num)}=$num-33;}
        return %dict_q_to_Q;
		
    }else{
		my $num;
		my %dict_q_to_Q=();
		for($num=59; $num<=126; $num++){
			$dict_q_to_Q{chr($num)}=$num-64;}
		return %dict_q_to_Q;
    }}

# Change Phred/Solexa quality score to ASCII character
sub Q_to_q($){
    
    if( $format eq "sanger" ){
		my $num;
		my %dict_Q_to_q=();
		for($num=-5; $num<=93; $num++){
			$dict_Q_to_q{$num}=chr($num+33);}
        return %dict_Q_to_q;
		
    }else{
		my $num;
		my %dict_Q_to_q=();
		for($num=-5; $num<=93; $num++){
			$dict_Q_to_q{$num}=chr($num+64);}
		return %dict_Q_to_q;
    }}

# Change Phred/Solexa quality score to probability
sub Q_to_p(){
    
	if( $format eq "solexa" ){
		my $num;
		my %dict_Q_to_p=();
		for($num=-5; $num<=93; $num++){
            $dict_Q_to_p{$num}=(10**(-$num/10)) / ((10**(-$num/10))+1);}
		return %dict_Q_to_p;
	}else{
		my $num;
		my %dict_Q_to_p=();
		for($num=-5; $num<=93; $num++){
            $dict_Q_to_p{$num}=(10**(-$num/10));}
		return %dict_Q_to_p;
	}
}

# Change probability to Phred/Solexa quality score
sub p_to_Q($){
    
	my $p = shift;
    
    if( $format && $format eq "solexa" ){
        return -10 * &log10($p/(1-$p));
    }else{
        return -10 * &log10($p);
    }
}

# log10 function
sub log10($){
    
	my $number = shift;
	return log($number)/log(10);
}

# print summary of Q, q and p values
sub print_lookup_table{
	my %dict_Q_to_p=%{+shift};
	my %dict_Q_to_q=%{+shift};
    
    
	print STDOUT "Char\tQPhred\tProb\n";
	for( my $i = -5; $i <= 93; $i++ ){
		
		my $q = $dict_Q_to_q{$i};
		my $p = $dict_Q_to_p{$i};
		
		print STDOUT $q, "\t";
		print STDOUT $i, "\t";
		print STDOUT sprintf($p), "\n";
	}
}

# automatic format detection
sub get_format(*$){
	
	# set function variables
	local *FILEHANDLE = shift;
	my $number_of_sequences = shift;
	my $format = "";
	
	# set regular expressions
	my $sanger_regexp = qr/[!"#$%&'()*+,-.\/0123456789:]/;
	my $solexa_regexp = qr/[\;<=>\?]/;
	my $solill_regexp = qr/[JKLMNOPQRSTUVWXYZ\[\]\^\_\`abcdefgh]/;
	my $all_regexp = qr/[\@ABCDEFGHI]/;
	
	# set counters
	my $sanger_counter = 0;
	my $solexa_counter = 0;
	my $solill_counter = 0;
	
	# go to file start
	seek(FILEHANDLE, 0, 0);
	
	# step through quality scores
	for( my $i = 0; $i < $number_of_sequences; $i++ ){
    
    # retrieve qualities
    <FILEHANDLE>;
    <FILEHANDLE>;
    <FILEHANDLE>;
    my $qualities = <FILEHANDLE>;
    chomp($qualities);
    
    # check qualities
    if( $qualities =~ m/$sanger_regexp/ ){
    $sanger_counter = 1;
    last;
    }
    if( $qualities =~ m/$solexa_regexp/ ){
    $solexa_counter = 1;
    }
    if( $qualities =~ m/$solill_regexp/ ){
    $solill_counter = 1;
    }
	}
	
	# determine format
	if( $sanger_counter ){
    $format = "sanger";
	}elsif( !$sanger_counter && $solexa_counter ){
    $format = "solexa";
	}elsif( !$sanger_counter && !$solexa_counter && $solill_counter ){
    $format = "illumina";
	}
	
	# go to file start
	seek(FILEHANDLE, 0, 0);
	
	# return file format
	return( $format );
    }
    
    # trim sequences using the BWA algorithm
    sub bwa_trim($$){
	
	my $threshold = shift;
	my $array_ref = shift;
	
	my @array  = @{$array_ref};
	my $length = scalar @array;
	
	# only calculate if quality fails near end
	if( $array[$#array] >= $threshold ){
    return $length;
	}
	
	# run bwa equation
	my @arg;
	for( my $i = 0; $i < $length - 1; $i++ ){
    
    my $x = $i + 1;
    for( my $j = $x; $j < $length; $j++ ){	
    $arg[$x] += $threshold - $array[$j];
    }
	}
	
	# find number of 5' bases to retain
	my $index = 0;
	my $maxval = 0;
	for ( 1 .. $#arg ){
    if ( $maxval < $arg[$_] ){
    $index = $_;
    $maxval = $arg[$_];
    }
	}
	
	# number of bases to retain
	return $index;
    }
