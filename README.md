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

### BUG with MacOSX Mojave and higher ###
You may encounter this bug after making the update towards Mojave or higher
fatal error: 'wchar.h' file not found
Apparently the names of the folders containing main dependencies were changed
To solve it this is one option (worked on Mojave 10.14.6) 
    
    xcode-select --install
    sudo xcode-select --switch /Library/Developer/CommandLineTools/
    open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

### How to build with LINUX (tested on SCITAS - EPFL clusters) ### 

Here an exemple working with Ubuntu 16.04 and clang6.0 :

./make.pl --compiler "clang++-6.0" --flags "-std=c++11"

On Scitas, you can load module load gcc and use g++

./make.pl --compiler "g++" --flags "-std=c++11"

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
    ./TEnrich 
    --bed_dir         path/to/dirWithBeds                     [required]
    --out_dir         path/to/dirOut                          [required] 
    --genome_size     size_genome                             [optional] 
    --comp_sense      'te_inter_peak','peak_inter_te','auto'  [optional] 
    --stat_test_type  'greater','less'                        [optional] 
    --padj            'true','false'                          [optional]
    --ref_subfam      path/to/file                            [optional]
    --ref_fam         path/to/file                            [optional]
    --ref_cluster     path/to/file                            [optional]
    --te_database     path/to/file                            [optional]
    --help (-h)         
    --version (-v) 
  
    --bed_dir 
    path/to/dirWithBeds [string] : every file with *.bed ext in the
    folder will be used

     --out_dir path/to/dirOut [string] : The folder is created and results written
    inside (WARNING: everything is cleaned before a new run)

     --genome_size [integer/double] : Genome size over which the enrichment
    calculation will be done. Expect only digits. [OPTIONAL]. Default value :
    3088269832

     --comp_sense ['te_in_peak','peak_in_te','auto'] : Defines the direction for
    the comparison, 'te_in_peak' or 'peak_in_te'. In 'auto' mode, it will take the
    enrichment of the smaller to the bigger [OPTIONAL]. Default value : 'auto'

     --stat_test_type ['greater','less'] : For statistical test done
    (Hypergeometric and Binomial), tell if we want the right tail ('greater') or
    the left tail ('less'), in other words the probability of having a equal or
    greater / equal or lower number of hits in the intersect. [OPTIONAL]. Default
    value : 'greater'

     --padj ['true','false'] : tell if you want to print the adjusted p-val (with
    the Benjamin-Hochsberg correction). [OPTIONAL]. Default value : 'true'

    --type_hypergeom ['ala_bedtools_fisher','regular'] : type of hypergeometric te-
    st to do. By default, will use a similar method that bedtools fisher uses (de-
    scription: https://bedtools.readthedocs.io/en/latest/content/tools/fisher.html).
    With regular option, uses a simple hypergeometric without weighting for TE loci 
    length.

     --ref_subfam [STRING] : subfam ref file obtained with utils/make_ref_file.pl.
    If not specified, the one for hg19 in db/ folder will be used [OPTIONAL].
    Default value : 'db/Subfam_ref_TE.txt'

     --ref_fam [STRING] : fam ref file obtained with utils/make_ref_file.pl. If
    not specified, the one for hg19 in db/ folder will be used [OPTIONAL]. Default
    value : 'db/Fam_ref_TE.txt'

     --ref_cluster [STRING] : TE cluster file obtained with
    utils/make_ref_file.pl. If not specified, the one for hg19 in db/ folder will
    be used [OPTIONAL]. Default value : 'db/Clusters_ref_TE.txt'

     --te_database [STRING] : database of TE used to make the intersection. Should
    be in bed format, as returned by utils/convert_repeatmasker.sh. By default,
    uses hg19 with LTR merged by J.Duc. [OPTIONAL]. Default value :
    'db/hg19_TE_repmask_LTRm_s_20140131.bed'  

_________________________________________________________________________________
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


