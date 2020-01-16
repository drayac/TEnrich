#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <libgen.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <iomanip>
#include "binom_pval.hpp"
#include "fisher_pval.hpp"
#include "functions.hpp" 
#include <algorithm>
#include <limits.h>

#define constant_e (2.71828)
#define PI (3.14159265359)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int main(int argc, char* argv[])
{
    ////////////////////
    // INITIALISATION //
    ////////////////////
    
    // Some constant values hard-coded
    const int LINE_MAX_SIZE = 200;      const double FDR = 0.05 ;
    const int te_data_size = 4570939 ;  bool print_padj = 1 ;
    std::string this_dir = <FOLDER_INSTALL> ;
    std::string version  = <VERSION> ; 
    
    check_folder(this_dir) ; // check this_dir exists
    
    // Get help if needed
    get_help(std::string(argv[1]), argc, version) ;

    // Get parameters
    std::vector <std::string> input = get_parameters(argc,argv) ; 
    std::string bed_path = input[0] , out_path = input[1] , type_hypergeom = input[2], comparison_direction = input[3] , comparison_type = input[4] , print_padj_str = input[5] , ref_file = input[6] , ref_file_fam = input[7] , te_data = input[9] , single_file = input[10] ;
    long int size_hg19 = std::stol(input[12]) ; 
    int idx_col = std::stoi(input[11]) ;
    if ( print_padj_str.compare("false") == 0 ){ print_padj = 0 ; }

    // Creating output folders
    create_folder(out_path,"summary_bed") ;         create_folder(out_path,"summary_te_fam") ;
    create_folder(out_path,"summary_te_subfam") ;   system("mkdir -p temp_TEnrich") ;

    // Check TE data number of fields
    std::cout << "Checking that TE data is ok\n" ;
    int expect_n_fields = 9 ;
    check_te_database(te_data, expect_n_fields ) ;
   
    // Initialize Hash Tables
    std::unordered_map<std::string, int> peak_count , peak_count_on_te , peak_count_on_te_unique, peak_len , peak_total_bp , te_inter_peak , te_inter_peak_unique, teFam_inter_peak, teFam_inter_peak_unique ; 
    std::unordered_map<std::string, std::string> summary_peak_line ;

    /////////////////////////////////////////////
    // 1) CONCATENATE MULTIPLE FILE WITH A TAG //
    /////////////////////////////////////////////
    
    std::string concat_bed = "temp_TEnrich/concat_all.bed" ;
    std::string list_files = "temp_TEnrich/temp.list_files.txt" ;
    
    if ( bed_path.compare("empty") != 0 ){ // only if multiple files are present
        concat_bed_files(bed_path, list_files, concat_bed, peak_count, peak_len, peak_total_bp, summary_peak_line, size_hg19) ;
        idx_col = 12 ;
    } else { // if working with a single-file, consider it as our concatenated bed
        concat_bed = single_file ;
        idx_col += 9 ;
    }

    ///////////////////////////////////////////////
    // 2) INTERSECT concat_bed with te_data file //
    ///////////////////////////////////////////////
   
    // bedtools intersection  
    std::string inter_bed_path = "temp_TEnrich/temp.out.bed" ;
    bedtools_intersect( "-f 0.5 -F 0.5 -e -wa -wb", // bedtools options
                        "sort -k1,1 -k2,2n -k9,9 -k10,10n", // sorting option
                        te_data, inter_bed_path, concat_bed) ;

    // Parse intersect and get all counts
    std::unordered_map<std::string, int> save_fake_list ; // fake list used with single-file
    parse_intersect(inter_bed_path, bed_path, save_fake_list, idx_col,
                    te_inter_peak , teFam_inter_peak, 
                    te_inter_peak_unique, teFam_inter_peak_unique,
                    peak_count_on_te) ;

    // Resort by peaks to count te_in_peak 
    std::unordered_map<std::string, int> peak_inter_teFam_unique, peak_inter_te_unique ;
    std::string inter_bed_sortPeaks = "temp_TEnrich/temp.out.sortPeaks.bed" ;
    sort_peaks("sort -k9,9 -k10,10n -k13", inter_bed_path, inter_bed_sortPeaks) ;
    parse_intersect_sortPeaks( inter_bed_path,    inter_bed_sortPeaks, idx_col, 
                                peak_count_on_te, peak_count_on_te_unique,
                                peak_inter_teFam_unique, peak_inter_te_unique , 
                                list_files ) ;

    //////////////////////////////////////
    // 3) COMPUTE ENRICHMENT BY SAMPLES //
    //////////////////////////////////////
  
    std::cout << "Starting enrichment analysis for sample : \n" ;
    // initialize variables 
    bool prhead_mat_all_best = 1, prhead_mat_all_best_fam = 1, prhead_summary_te  = 1 ;
    std::string matrix_path_all_best = out_path + "/matrix_padjAlaBTFisher_Subfam.txt" ;
    std::string matrix_fam = out_path + "/matrix_padjAlaBTFisher_Fam.txt" ;
    std::string line_ref;
    // nonTE parameters 
    int total_nonTE = 4433186 ; // we estimate the number of nonTE interval ~= TE interval after merge
    int total_nonTE_bp = size_hg19 - 1434913179 ; // size_hg19 - genome span of TE    
    double nonTE_avg_size = double(total_nonTE_bp) / double(total_nonTE) ;
    double nonTE_genome_ratio = double(total_nonTE_bp) / double(size_hg19) ; 

    // Create peak summary file
    std::string summary_bed_path = out_path + "/summary_bed.txt" ;
    std::ofstream summary_bed ; 
    summary_bed.open(summary_bed_path, std::fstream::app) ;

    // print summary header 
    summary_bed << "sample_name\tpeak.number\tpeak.average.size\tpeak.total.bp\tpeak.ratio.genome\tpeak.overlap.te\tpeak.ratio.overlap.te\n" ;

    // Iterate over samples (list of files in input folder) 
    std::ifstream list_f2(list_files);
    if (this_is_empty(list_f2)){
        std::cout << "Problem while opening " << list_files << ", exiting...\n" ;
        exit (EXIT_FAILURE) ;
    }
    if (list_f2.is_open()) {
        
        // Iterate over samples to compute enrichments
        while (getline(list_f2, line_ref)) {
            
            std::string tag_name = remove_ext(base_name(line_ref)) ;
            std::cout << ", " << tag_name ;
            
            // getting peak stats
            int mean_peaklen = peak_len[tag_name] , my_peak_count = peak_count[tag_name] , my_peak_count_on_te_unique = peak_count_on_te_unique[tag_name] , peak_tot_bp = peak_total_bp[tag_name] ;

            // add info to peak summary 
            std::string peak_sum_line = summary_peak_line[tag_name] ;
            summary_bed << peak_sum_line << "\t" << my_peak_count_on_te_unique << "\n" ;

            // bool to print header of summmary files 
            bool prhead_sum_all_best = 1, prhead_sum_all_best_fam = 1 ;
          
            ///////////////////////////////////////// 
            // Enrichment analysis SUBFAMILY LEVEL //
            ///////////////////////////////////////// 
            enrichment_analysis("te_subfam",
                                ref_file, size_hg19, comparison_type, tag_name,
                                my_peak_count, my_peak_count_on_te_unique,
                                peak_tot_bp, mean_peaklen, 
                                comparison_direction, te_data_size,
                                total_nonTE, total_nonTE_bp,
                                nonTE_avg_size, nonTE_genome_ratio,
                                matrix_path_all_best, out_path,
                                print_padj, prhead_mat_all_best, 
                                prhead_mat_all_best_fam, prhead_summary_te,
                                peak_inter_te_unique, te_inter_peak_unique,
                                type_hypergeom ) ;  

            ///////////////////////////////////////// 
            // Enrichment analysis FAMILY LEVEL    //
            /////////////////////////////////////////
            enrichment_analysis("te_fam",
                                ref_file_fam, size_hg19, comparison_type, tag_name,
                                my_peak_count, my_peak_count_on_te_unique,
                                peak_tot_bp, mean_peaklen, 
                                comparison_direction, te_data_size,
                                total_nonTE, total_nonTE_bp,
                                nonTE_avg_size, nonTE_genome_ratio,
                                matrix_fam, out_path,
                                print_padj, prhead_mat_all_best, 
                                prhead_mat_all_best_fam, prhead_summary_te,
                                peak_inter_teFam_unique, teFam_inter_peak_unique,
                                type_hypergeom ) ;

             prhead_summary_te = 0 ;
        }
        list_f2.close() ; summary_bed.close() ;
    }

    std::cout << "\n" ;
    // Cleaning temp directory 
    system("rm -rf temp_TEnrich/") ;
}
