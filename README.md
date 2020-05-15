_________________________________________________________________________________
## TEnrich ##

            A code for to compute statistical enrichment of transposable elements
            on a group of bed files. 

Code was written in C++11 , by Alexandre Coudray from the
laboratory of Virology and Genetics at the EPFL in 2019.

_________________________________________________________________________________
## Download and get TEnrich ##

### Two options : ###
    
    Go to the desired folder where you want to install TEnrich and launch 
    git clone https://github.com/alexdray86/TEnrich.git

Second option : 

    Click on the "Clone and download" (green button) above

Then, Follow the instruction below to compile the executer and launch the script 

_________________________________________________________________________________
## How to build the executer ##

#### How to make the file ####
launch the following command from the root of TEnrich folder

./make.pl

Made for clang++ on MacOS, feel free to change the compiler

Default compiler : clang++ (developed with clang++-7.0) 

Default flags    : -std=c++11

#### IMPORTANT : ####
Choose wisely where you make the script, because
you should not move it afterwards (databases in TEnrich/db/ 
folder should be accessible). If you want to move it, just
re-do a make after the move.

### How to build with LINUX (tested on SCITAS - EPFL clusters) ### 

Here an exemple working with Ubuntu 16.04 and clang6.0 :

./make.pl --compiler "clang++-6.0" --flags "-std=c++11"

On Scitas, you can load module load gcc and use g++

./make.pl --compiler "g++" --flags "-std=c++11"

### BUG with MacOSX Mojave and higher ###
You may encounter this bug after making the update towards Mojave or higher
fatal error: 'wchar.h' file not found
Apparently the names of the folders containing main dependencies were changed
To solve it this is one option (worked on Mojave 10.14.6) 
    
    xcode-select --install
    sudo xcode-select --switch /Library/Developer/CommandLineTools/
    open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

#### How to get help ####
Once compiled, launch the help with :

./TEnrich -h        (or --help)

_________________________________________________________________________________
## Description ##

        The script was made to work on large number of bed files
        at once. Therefore the input should be a path to a dire-
        ctory containing bed files. Every files with *.bed exte-
        nsion inside that folder will be used by the script. The
        genome size for the enrichment can be changed (Warning:
        don't put scientific notation, e.g. 2.7e9 won't work). 

#### Usage : ####
   _________________________________________________________________________________

-i [--bed_dir] path/to/dirWithBeds [string]       every file with *.bed ext in the folder will be us
                                                  ed

-o [--out_dir] path/to/dirOut [string]            The folder is created and results written inside (
                                                  WARNING: everything is cleaned before a new run)

--single_file path/to/bed_file [string]           if this is given, it will use a single file instead
                                                  of a group of bed file to do the enrichment (cance–
                                                  ls --bed_dir option). If no index column is given,
                                                  will perform the enrichment analysis on each single
                                                  lines. Otherwise, it will group lines per name of
                                                  the feature in the column designed by --idx_col op–
                                                  tion.

--type_hypergeom ['ala_bedtools_fisher','regular']type of hypergeometric test to do. By default, will
                                                  use a similar method that bedtools fisher uses (de–
                                                  scription: https://bedtools.readthedocs.io/en/late–
                                                  st/content/tools/fisher.html). With regular option,
                                                  uses a simple hypergeometric without weighting for
                                                  TE loci length

--comp_sense ['te_in_peak','peak_in_te','auto']   Defines the direction for the comparison, 'te_in_p–
                                                  eak' or 'peak_in_te'. In 'auto' mode, it will take
                                                  the enrichment of the smaller to the bigger [OPTIO–
                                                  NAL]. Default value : 'auto'

--stat_test_type ['greater','less']               For statistical test done (Hypergeometric and Bino–
                                                  mial), tell if we want the right tail ('greater')
                                                  or the left tail ('less'), in other words the prob–
                                                  ability of having a equal or greater / equal or lo–
                                                  wer number of hits in the intersect. [OPTIONAL]. D–
                                                  efault value : 'greater'

--padj ['true','false']                           tell if you want to print the adjusted p-val (with
                                                  the Benjamin-Hochsberg correction). [OPTIONAL]. De–
                                                  fault value : 'true'

--genome_size [integer/double]                    Genome size over which the enrichment calculation
                                                  will be done. Expect only digits. [OPTIONAL]. Defa–
                                                  ult value : 3088269832

--ref_subfam [STRING]                             subfam ref file obtained with utils/make_ref_file.–
                                                  pl. If not specified, the one for hg19 in db/ fold–
                                                  er will be used [OPTIONAL]. Default value : 'db/Su–
                                                  bfam_ref_TE.txt'

--ref_fam [STRING]                                fam ref file obtained with utils/make_ref_file.pl.
                                                  If not specified, the one for hg19 in db/ folder w–
                                                  ill be used [OPTIONAL]. Default value : 'db/Fam_re–
                                                  f_TE.txt'

--te_database [STRING]                            database of TE used to make the intersection. Shou–
                                                  ld be in bed format, as returned by utils/convert_–
                                                  repeatmasker.sh. By default, uses hg19 with LTR me–
                                                  rged by J.Duc. [OPTIONAL]. Default value : 'db/hg1–
                                                  9_TE_repmask_LTRm_s_20140131.bed

## TE database ##

The RepeatMasker 4.0.5 (Library 20140131, http://www.repeatmasker.org/species/hg.html) was used to generate an in-house repeats list where fragmented ERVs were reassembled according to the following procedure: ERVs fragments were annotated either as fragmented internal parts (ERV-int) or LTRs (Long Terminal Repeats). We then computed the frequency distribution of each LTR/ERV–int pairs and compute their respective enrichment through a Wald test. LTR to an ERV-int were merged whener the LTR type was present in at least 2% of the pairs, with a p-val<0.001. The two elements had to be in the same orientation with distance<100bp. The name of the internal part was given to the resulting LTR/ERV-int merged element (e.g. HERVH-int). Fragmented ERV-int or fragmented LTRs from the same subfamily (same name in Repbase database), with the same orientation and closer than 100bp were also merged.

_________________________________________________________________________________
## Working with another species than hg19 ##

#### Work with another Species than hg19 ####
exemple for danRer10 :
- Download TE database from RepeatMasker
- convert .fa.out file with utils/convert_repeatMasker.sh :

utils/convert_repeatMasker.sh danRer10_repMask.fa.out \
    > db/danRer10_repMask406_dfam2.bed

- remove repeats if not desired 
grep -v 'Low_complexity\|Satellite\|Simple_repeat\|Unknown' \
    db/danRer10_repMask406_dfam2.bed \
    > db/danRer10_repMask406_dfam2_TEs.bed 

- launch utils/make_ref_TE.pl to make files needed by TEnrich
INDEX is the col number of the desired feature (subfam/fam) :

#### subfam should always be in field 7 (0-based) ####
utils/make_ref_TE.pl --file db/danRer10_repMask406_dfam2_TEs.bed \
    --index 6 \
    > db/danRer10_Subfam_ref_TE.txt

#### class/families should always be in field 6 (0-based) ####
utils/make_ref_TE.pl --file db/danRer10_repMask406_dfam2_TEs.bed \
    --index 7 \
    > db/danRer10_Fam_ref_TE.txt

_________________________________________________________________________________
## Versions history ##

Updates news :
- v0.2 : -> Enrichment of TE families (prev. only subfam)
         -> adjustement of make.pl for Ubuntu   
- v0.3 : -> new options added
         -> add nonTE enrichment 
- v0.4 : -> add enrichment of TE clusters (incomplete)
- v0.5 : -> new utils to be able to check foreign species
         -> fixed bug when beds of diverse size are given
- v0.6 : -> fixed bug to allow the use of relative path 
            for --bed_dir and --out_dir parameters
- v0.7 : -> add --single_file option to make enrichment 
            from a single bed file
- v1.0 : -> Improvements on the code , make it cleaner and more robust
- v1.1 : -> Fixed bug with binomial stats
- v1.2 : -> Fixed bug - errors in counting the overlaps 

