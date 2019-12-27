#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd ;

my $CXX='clang++' ;
my $CXXFLAGS='-std=c++11' ;
my $version='0.7' ;

GetOptions(
    "compiler=s" => \$CXX,
    "flags=s" => \$CXXFLAGS,
);

# Change folder name (needed for having portable path)
my $abs_path = getcwd;
print "adjusting folder path to $abs_path in TEnrich.cpp ...\n" ;
$abs_path =~ s/\//\\\//g ;
`sed 's~\<FOLDER_INSTALL\>~\"'$abs_path'\"~g' src/TEnrich.cpp | sed 's~\<VERSION\>~\"'$version'\"~g' > src/TEnrich_clean.cpp` ;
`sed 's~\<FOLDER_INSTALL\>~\"'$abs_path'\"~g' src/functions.cpp | sed 's~\<VERSION\>~\"'$version'\"~g' > src/functions_clean.cpp` ;

# make the file using clang ++ 
print "compiling TEnrich_clean.cpp with clang++ ...\n" ;
`$CXX $CXXFLAGS src/TEnrich_clean.cpp src/binom_pval.cpp src/fisher_pval.cpp src/functions_clean.cpp -o TEnrich` ;

# unzip te data
#my $te_data = "db/hg19_TE_repmask_LTRm_s_20140131.bed.gz" ;
#my $unzip_te_data = "db/hg19_TE_repmask_LTRm_s_20140131.bed" ;
#if ( -e $unzip_te_data ){
#    print "te database already unziped, skipping ... \n" ;
#}
#else
#{
#    print "unzipping te database ... \n" ;
#    `gunzip db/hg19_TE_repmask_LTRm_s_20140131.bed.gz` ;
#}

# clean up directories
print "cleaning folders ...\n\n" ;
`rm src/*_clean.cpp` ;

print "if make.pl was successful, you should see the help appear:\n\n" ;
sleep(1) ;
print(`./TEnrich -h`) ;
