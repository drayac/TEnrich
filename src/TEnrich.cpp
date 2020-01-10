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

// MAIN 
int main(int argc, char* argv[])
{
    // Some constant values hard-coded
    const int LINE_MAX_SIZE = 200;      const double FDR = 0.05 ;
    const int te_data_size = 4570939 ;  bool print_padj = 1 ;
    std::string this_dir = <FOLDER_INSTALL> ;
    std::string version  = <VERSION> ; 
    
    check_folder(this_dir) ;
    
    // Get help
    get_help(std::string(argv[1]), argc, version) ;

    // Get and initialize parameters
    std::vector <std::string> input = get_parameters(argc,argv) ; 
    std::string bed_path = input[0] , out_path = input[1] , comparison_direction = input[3] , comparison_type = input[4] , print_padj_str = input[5] , ref_file = input[6] , ref_file_fam = input[7] , te_data = input[9] , single_file = input[10] ;
    long int size_hg19 = std::stol(input[2]) ; 
    int idx_col = std::stoi(input[11]) ;
    if ( print_padj_str.compare("false") == 0 ){ print_padj = 0 ; }

    // nonTE parameters 
    int total_nonTE = 4433186 ; // we estimate the number of nonTE interval ~= TE interval after merge
    int total_nonTE_bp = size_hg19 - 1434913179 ; // size_hg19 - genome span of TE    
    double nonTE_avg_size = double(total_nonTE_bp) / double(total_nonTE) ;
    double nonTE_genome_ratio = double(total_nonTE_bp) / double(size_hg19) ; 
   
    // Create temp directory if not existant
    system("mkdir -p temp_TEnrich") ;
 
    // Creating output folders
    create_folder(out_path,"summary_bed") ;         create_folder(out_path,"summary_te_fam") ;
    create_folder(out_path,"summary_te_subfam") ;  

    // Check TE data number of fields
    std::cout << "Checking that TE data is ok\n" ;
    int expect_n_fields = 9 ;
    check_te_database(te_data, expect_n_fields ) ;
   
    // Initialize Hash Tables
    std::unordered_map<std::string, int> peak_count , peak_count_on_te , peak_count_on_te_unique, peak_len , peak_total_bp , te_inter_peak , te_inter_peak_unique, teFam_inter_peak, teFam_inter_peak_unique ; 
    std::unordered_map<std::string, std::string> summary_peak_line ;


    //////////////////////////////////////////
    // CONCATENATE MULTIPLE FILE WITH A TAG //
    //////////////////////////////////////////
    
    std::string concat_bed = "temp_TEnrich/concat_all.bed" ;
    std::string list_files = "temp_TEnrich/temp.list_files.txt" ;
    
    if ( bed_path.compare("empty") != 0 ){ // only if multiple files are present
        concat_bed_files(bed_path, list_files, concat_bed, peak_count, peak_len, peak_total_bp, summary_peak_line, size_hg19) ;
        idx_col = 12 ;
    } else { // if working with a single-file, consider it as our concatenated bed
        concat_bed = single_file ;
        idx_col += 9 ;
    }


    ////////////////////////////////////////////
    // INTERSECT concat_bed with te_data file //
    ////////////////////////////////////////////
   
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
                                peak_inter_teFam_unique,peak_inter_te_unique , 
                                list_files ) ;

    
    ////////////////////////////////////////////
    // COMPUTE ENRICHMENT BY SAMPLES          //
    ////////////////////////////////////////////
   
    // initialize variables 
    bool prhead_mat_all_best = 1, prhead_mat_all_best_fam = 1, prhead_summary_te  = 1 ;
    std::string matrix_path_all_best = out_path + "/matrix_padjAlaBTFisher_Subfam.txt" ;
    std::string matrix_fam = out_path + "/matrix_padjAlaBTFisher_Fam.txt" ;
    std::string line_ref;

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
            std::cout << "Enrichment analysis of " << tag_name << "\n" ;
            
            // getting peak stats
            int mean_peaklen = peak_len[tag_name] , my_peak_count = peak_count[tag_name] , my_peak_count_on_te_unique = peak_count_on_te_unique[tag_name] , peak_tot_bp = peak_total_bp[tag_name] ;
            double ratio_genome_peak = double(peak_tot_bp) / double(size_hg19) ;
            double peak_total_Mbp = double(peak_tot_bp) / 1000000 ;
            double peak_ratio_on_te = double(my_peak_count_on_te_unique) / double(my_peak_count) ; 
            
            // add info to peak summary 
            std::string peak_sum_line = summary_peak_line[tag_name] ;
            summary_bed << peak_sum_line << "\t" << my_peak_count_on_te_unique << "\t" << peak_ratio_on_te << "\n" ;
            // bool to print header of summmary files 
            bool prhead_sum_all_best = 1, prhead_sum_all_best_fam = 1 ;
           
            //////////////////// 
            //  TE SUBFAM     //
            ////////////////////
            
            // Begin enrichment analysis by SUBFAM 
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
                                peak_inter_te_unique, te_inter_peak_unique ) ;  
            
            /*
            std::ifstream ref_in(ref_file);
            if (this_is_empty(ref_in)){
                std::cout << "Problem while opening " << ref_file << ", exiting...\n" ;
                exit (EXIT_FAILURE) ;
            }
            if (ref_in.is_open()) {
                
                // initialize vectors
                std::vector<std::string> subfam_names;
                std::vector<double> all_pvals , all_pvals_alafisher , all_pvals_binomial, all_pvals_rev , all_pvals_alafisher_rev , all_pvals_binomial_rev , all_pvals_best , all_pvals_alafisher_best , all_pvals_binomial_best, all_total_bp_subfam, all_subfam_genome_ratio ;
                std::vector<int> all_te_inter_peak, all_te_inter_peak_unique , all_total_subfam , all_avg_subfam_size;

                // Iterate over subfam names of TEs (or sample of bed file 2) 
                std::string line;
                while (getline(ref_in, line)) { 
                    std::istringstream iss(line) ;
                    std::vector <std::string> fields ;
                    std::string field ;
                    while(std::getline(iss, field, '\t')){ 
                        fields.push_back(field);
                    }

                    std::string subfam_name = fields[0] ;
                    std::string key = subfam_name + "_" + tag_name ;
                    int n_subfam = std::stoi(fields[1]) ;
                    int total_bp_length_subfam = std::stoi(fields[2]) ;
                    int avg_subfam_size = std::stoi(fields[3]) ;
                    double ratio_genome_subfam = double(total_bp_length_subfam)/double(size_hg19) ;
                    double total_Mbp_len_subfam = double(total_bp_length_subfam) / 1000000 ;

                    all_te_inter_peak.push_back(te_inter_peak[key]) ;
                    all_te_inter_peak_unique.push_back(te_inter_peak_unique[key]) ;
                    all_total_subfam.push_back(n_subfam);
                    all_total_bp_subfam.push_back(total_Mbp_len_subfam);
                    all_avg_subfam_size.push_back(avg_subfam_size);
                    all_subfam_genome_ratio.push_back(ratio_genome_subfam);

                    // Calculation of the p-values for each 1-1 relation
                    if (te_inter_peak_unique.find(key) == te_inter_peak.end() || te_inter_peak[key] == 0){
                        subfam_names.push_back (fields[0]) ; all_pvals.push_back (1.0) ;
                        all_pvals_alafisher.push_back (1.0) ; all_pvals_binomial.push_back (1.0) ;
                        all_pvals_rev.push_back (1.0) ; all_pvals_alafisher_rev.push_back (1.0) ;
                        all_pvals_binomial_rev.push_back (1.0) ; all_pvals_best.push_back (1.0) ;
                        all_pvals_alafisher_best.push_back (1.0) ; all_pvals_binomial_best.push_back (1.0) ;
                    }
                    else
                    {
                        ////////////////////////////////////// 
                        // PVAL ENRICHMENT PEAKs AMONG TEs  //
                        //////////////////////////////////////
                        
                        // Regular HyperGeometric
                        long long a_11 = peak_inter_te_unique[key] ; //te_inter_peak[key] ;
                        long long a_12 = MAX(0L,n_subfam - a_11) ;
                        long long a_21 = MAX(0L,my_peak_count_on_te_unique - a_11) ;
                        long long a_22 = te_data_size - a_11 - a_12 - a_21 ;
                        double pval = fisher_exact(a_11, a_12, a_21, a_22, comparison_type) ;
                        all_pvals.push_back (pval) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        int total_average = mean_peaklen + avg_subfam_size ;
                        long long b_11 = peak_inter_te_unique[key] ; //te_inter_peak[key] ;
                        long long b_12 = MAX(a_12, n_subfam - b_11 ) ;
                        long long b_21 = MAX(0L,my_peak_count - b_11) ;
                        long long b_22 = MAX( int(double(size_hg19) / double(total_average)) - b_11 - b_21 - b_12, 1 )  ;
                        double pval_alafisher = fisher_exact(b_11, b_12, b_21, b_22, comparison_type) ;
                        all_pvals_alafisher.push_back (pval_alafisher) ;

                        // Binomial exact test with genome ratio as p
                        int n_tot = my_peak_count ; // number of trials
                        int X_obs = peak_inter_te_unique[key] ; //te_inter_peak[key] ;
                        double p_obs = double( mean_peaklen*X_obs ) / double(total_bp_length_subfam) ;
                        double q_obs = 1 - p_obs ;
                        double p_exp = double(total_bp_length_subfam) / double(size_hg19) ; // Prob to touch subfam by random
                        double pval_binomial = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial.push_back(pval_binomial) ;

                        //////////////////////////////////////
                        // PVAL ENRICHMENT TEs AMONG PEAKs  //
                        //////////////////////////////////////

                        // Regular HyperGeometric
                        long long ar_11 = te_inter_peak_unique[key] ;
                        long long ar_12 = MAX(0L,n_subfam - ar_11) ;
                        long long ar_21 = MAX(0L,my_peak_count_on_te_unique - ar_11) ;
                        long long ar_22 = te_data_size - ar_11 - ar_12 - ar_21 ;
                        double pval_rev = fisher_exact(ar_11, ar_12, ar_21, ar_22, comparison_type) ;
                        all_pvals_rev.push_back (pval_rev) ;

                        add_pval_auto(  comparison_direction, mean_peaklen, avg_subfam_size, 
                                        all_pvals_best, pval, pval_rev) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        total_average = avg_subfam_size ;
                        long long br_11 = te_inter_peak_unique[key] ;
                        long long br_12 = MAX(ar_12, n_subfam - br_11 ) ;
                        long long br_21 = MAX(0L,my_peak_count - br_11 ) ;
                        long long br_22 = MAX( int(double(size_hg19) / double(total_average)) - br_11 - br_21 - br_12, 1)  ;
                        double pval_alafisher_rev = fisher_exact(br_11, br_12, br_21, b_22, comparison_type) ;
                        all_pvals_alafisher_rev.push_back (pval_alafisher_rev) ;
                        
                        add_pval_auto(  comparison_direction, mean_peaklen, avg_subfam_size, 
                                        all_pvals_alafisher_best, pval_alafisher, pval_alafisher_rev) ; 
                        
                        // Binomial exact test with genome ratio as p 
                        n_tot = n_subfam ; // number of trials
                        X_obs = te_inter_peak_unique[key] ;
                        p_exp = ratio_genome_peak ; // Prob to touch subfam by random
                        double pval_binomial_rev = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial_rev.push_back(pval_binomial_rev) ;
                        subfam_names.push_back (fields[0]) ;

                        add_pval_auto(  comparison_direction, mean_peaklen, avg_subfam_size, 
                                        all_pvals_binomial_best, pval_binomial, pval_binomial_rev) ;
                    }
                }

                // Getting pval for nonTE
                subfam_names.push_back ("nonTE") ;
                int nonTE_intersect = my_peak_count - my_peak_count_on_te_unique ;
                // reg hypergeometric test
                int nonTE_a12 = total_nonTE - nonTE_intersect ;
                int nonTE_a21 = my_peak_count - nonTE_intersect ;
                int nonTE_a22 = 2*4570939 - nonTE_intersect - nonTE_a12 - nonTE_a21 ;
                double nonTE_pval_reg = fisher_exact(nonTE_intersect, nonTE_a12, nonTE_a21, nonTE_a22, comparison_type) ;
                all_pvals_best.push_back (nonTE_pval_reg) ;

                // hypergeometric alafisher
                int peak_nonTE_total_average = 276 + mean_peaklen ;
                int nonTE_b22 = MAX( int(double(size_hg19) / double(peak_nonTE_total_average)) - nonTE_intersect - nonTE_a12 - nonTE_a21 , 1)  ;
                double nonTE_pval_alafisher = fisher_exact(nonTE_intersect, nonTE_a12, nonTE_a21, nonTE_b22, comparison_type) ;
                all_pvals_alafisher_best.push_back (nonTE_pval_alafisher) ;
                
                // binomial test
                int n_tot_nonTE = my_peak_count ; // number of trials
                int X_obs_nonTE = nonTE_intersect ;
                double p_obs_nonTE = double( mean_peaklen*X_obs_nonTE ) / double(total_nonTE_bp) ;
                double q_obs_nonTE = 1 - p_obs_nonTE ;
                double p_exp_nonTE = double(total_nonTE_bp) / double(size_hg19) ; // Prob to touch fam by random
                double nonTE_pval_binomial = pbinom(X_obs_nonTE-1, n_tot_nonTE, p_exp_nonTE, comparison_type) ;
                all_pvals_binomial_best.push_back (nonTE_pval_binomial) ;

                all_te_inter_peak.push_back(nonTE_intersect) ;
                all_te_inter_peak_unique.push_back(nonTE_intersect) ;
                all_total_subfam.push_back(total_nonTE);
                all_total_bp_subfam.push_back(total_nonTE_bp);
                all_avg_subfam_size.push_back( double(nonTE_avg_size) );
                all_subfam_genome_ratio.push_back(nonTE_genome_ratio);

                // Get ajusted p-values with Benjamin-Hochsberg method 
                std::unordered_map<std::string, double> padj_hygm_reg = calc_adj_pval(&all_pvals_best, &subfam_names, print_padj) ;
                std::unordered_map<std::string, double> padj_hygm_alaFish = calc_adj_pval(&all_pvals_alafisher_best, &subfam_names, print_padj) ;
                std::unordered_map<std::string, double> padj_hygm_binom = calc_adj_pval(&all_pvals_binomial_best, &subfam_names, print_padj) ;
                
                // Printing final results with 3 adjusted p-values
                print_all_pval_adj(matrix_path_all_best, tag_name,
                        &prhead_sum_all_best, &prhead_mat_all_best, &all_pvals_binomial_best,
                        padj_hygm_reg, padj_hygm_alaFish, padj_hygm_binom,
                        &subfam_names, &all_te_inter_peak_unique,
                        &all_total_subfam,  &all_te_inter_peak,
                        my_peak_count_on_te_unique, my_peak_count, &all_total_bp_subfam,
                        &all_avg_subfam_size, &all_subfam_genome_ratio,
                        peak_total_Mbp, peak_ratio_on_te, ratio_genome_peak,
                        "te_subfam", &prhead_summary_te, out_path ) ;
            }
            ref_in.close() ;
            */

            ///////////////// 
            /*  TE FAM     */
            /////////////////
            
            // Begin enrichment analysis by FAM 
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
                                peak_inter_te_unique, te_inter_peak_unique ) ;

            // Begin enrichment analysis for sample X 
            /*
            std::ifstream ref_in_fam(ref_file_fam);
            if (this_is_empty(ref_in_fam)){
                std::cout << "Problem while opening " << ref_file_fam << ", exiting...\n" ;
                exit (EXIT_FAILURE) ;
            }
            if (ref_in_fam.is_open()) {
                
                std::string line;
                // initialize vectors
                std::vector<std::string> fam_names;
                std::vector<double> all_pvals , all_pvals_alafisher , all_pvals_binomial, all_pvals_rev , all_pvals_alafisher_rev , all_pvals_binomial_rev , all_pvals_best , all_pvals_alafisher_best , all_pvals_binomial_best ;
                std::vector<int> all_te_inter_peak, all_te_inter_peak_unique , all_total_fam , all_avg_fam_size;
                std::vector<double> all_total_bp_fam , all_fam_genome_ratio ;

                // Iterate over subfam names of TEs (or sample of bed file 2) 
                while (getline(ref_in_fam, line)) { 
                    std::istringstream iss(line) ;
                    std::vector <std::string> fields ;
                    std::string field ;
                    while(std::getline(iss, field, '\t')){ 
                        fields.push_back(field);
                    }

                    std::string fam_name = fields[0] ;
                    std::string key = fam_name + "_" + tag_name ;
                    int n_fam = std::stoi(fields[1]) ;
                    int total_bp_length_fam = std::stoi(fields[2]) ;
                    int avg_fam_size = std::stoi(fields[3]) ;
                    double ratio_genome_fam = double(total_bp_length_fam)/double(size_hg19) ;
                    double total_Mbp_len_fam = double(total_bp_length_fam) / 1000000 ;

                    all_te_inter_peak.push_back(teFam_inter_peak[key]) ;
                    all_te_inter_peak_unique.push_back(teFam_inter_peak_unique[key]) ;
                    all_total_fam.push_back(n_fam);
                    all_total_bp_fam.push_back(total_Mbp_len_fam);
                    all_avg_fam_size.push_back(avg_fam_size);
                    all_fam_genome_ratio.push_back(ratio_genome_fam);

                    // Calculation of the p-values for each 1-1 relation
                    if (teFam_inter_peak.find(key) == teFam_inter_peak.end() || teFam_inter_peak[key] == 0){
                        fam_names.push_back (fields[0]) ; all_pvals.push_back (1.0) ;
                        all_pvals_alafisher.push_back (1.0) ; all_pvals_binomial.push_back (1.0) ;
                        all_pvals_rev.push_back (1.0) ; all_pvals_alafisher_rev.push_back (1.0) ;
                        all_pvals_binomial_rev.push_back (1.0) ; all_pvals_best.push_back (1.0) ;
                        all_pvals_alafisher_best.push_back (1.0) ; all_pvals_binomial_best.push_back (1.0) ;
                    }
                    else
                    {
                        ////////////////////////////////////// 
                        // PVAL ENRICHMENT PEAKs AMONG TEs  //
                        //////////////////////////////////////
                        
                        // Regular HyperGeometric
                        int a_11 = peak_inter_teFam_unique[key] ;//teFam_inter_peak[key] ;
                        int a_12 = MAX(0L,n_fam - a_11) ;
                        int a_21 = MAX(0L,my_peak_count_on_te_unique - a_11) ;
                        int a_22 = te_data_size - a_11 - a_12 - a_21 ;
                        double pval = fisher_exact(a_11, a_12, a_21, a_22, comparison_type) ;
                        all_pvals.push_back (pval) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        int total_average = mean_peaklen + avg_fam_size ;
                        int b_11 = peak_inter_teFam_unique[key] ; //teFam_inter_peak[key] ;
                        int b_12 = MAX(a_12, n_fam - b_11 ) ;
                        int b_21 = MAX(0L,my_peak_count - b_11) ;
                        int b_22 = MAX( int(double(size_hg19) / double(total_average)) - b_11 - b_21 - b_12, 1 )  ;
                        double pval_alafisher = fisher_exact(b_11, b_12, b_21, b_22, comparison_type) ;
                        all_pvals_alafisher.push_back (pval_alafisher) ;

                        // Binomial exact test with genome ratio as p
                        int n_tot = my_peak_count ; // number of trials
                        int X_obs = peak_inter_teFam_unique[key] ; //teFam_inter_peak[key] ;
                        double p_obs = double( mean_peaklen*X_obs ) / double(total_bp_length_fam) ;
                        double q_obs = 1 - p_obs ;
                        double p_exp = double(total_bp_length_fam) / double(size_hg19) ; // Prob to touch fam by random
                        double pval_binomial = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial.push_back(pval_binomial) ;

                        //////////////////////////////////////
                        // PVAL ENRICHMENT TEs AMONG PEAKs  //
                        //////////////////////////////////////

                        // Regular HyperGeometric
                        int ar_11 = teFam_inter_peak_unique[key] ;
                        int ar_12 = MAX(0L,n_fam - ar_11) ;
                        int ar_21 = MAX(0L,my_peak_count - ar_11) ;
                        int ar_22 = te_data_size - ar_11 - ar_12 - ar_21 ;
                        double pval_rev = fisher_exact(ar_11, ar_12, ar_21, ar_22, comparison_type) ;
                        all_pvals_rev.push_back (pval_rev) ;

                        add_pval_auto(  comparison_direction, mean_peaklen, avg_fam_size, 
                                        all_pvals_best, pval, pval_rev) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        total_average = avg_fam_size ;
                        int br_11 = teFam_inter_peak_unique[key] ;
                        int br_12 = MAX(ar_12, n_fam - br_11 ) ;
                        int br_21 = MAX(0L,my_peak_count - br_11 ) ;
                        int br_22 = MAX( int(double(size_hg19) / double(total_average)) - br_11 - br_21 - br_12, 1)  ;
                        double pval_alafisher_rev = fisher_exact(br_11, br_12, br_21, b_22, comparison_type) ;
                        all_pvals_alafisher_rev.push_back (pval_alafisher_rev) ;

                        add_pval_auto(  comparison_direction, mean_peaklen, avg_fam_size, 
                                        all_pvals_alafisher_best, pval_alafisher, pval_alafisher_rev) ; 

                        // Binomial exact test with genome ratio as p 
                        n_tot = n_fam ; // number of trials
                        X_obs = teFam_inter_peak_unique[key] ;
                        p_exp = ratio_genome_peak ; // Prob to touch peak by random
                        double pval_binomial_rev = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial_rev.push_back(pval_binomial_rev) ;

                        add_pval_auto(  comparison_direction, mean_peaklen, avg_fam_size, 
                                        all_pvals_binomial_best, pval_binomial, pval_binomial_rev) ;

                        fam_names.push_back (fields[0]) ;
                    }
                }

                // Getting pval for nonTE
                fam_names.push_back ("nonTE") ;
                int nonTE_intersect = my_peak_count - my_peak_count_on_te_unique ;

                // reg hypergeometric test
                int nonTE_a12 = total_nonTE - nonTE_intersect ;
                int nonTE_a21 = my_peak_count - nonTE_intersect ;
                int nonTE_a22 = 2*4570939 - nonTE_intersect - nonTE_a12 - nonTE_a21 ;
                double nonTE_pval_reg = fisher_exact(nonTE_intersect, nonTE_a12, nonTE_a21, nonTE_a22, comparison_type) ;
                all_pvals_best.push_back (nonTE_pval_reg) ;

                // hypergeometric alafisher 
                int peak_nonTE_total_average = 276 + mean_peaklen ;
                int nonTE_b22 = MAX( int(double(size_hg19) / double(peak_nonTE_total_average)) - nonTE_intersect - nonTE_a12 - nonTE_a21 , 1)  ;
                double nonTE_pval_alafisher = fisher_exact(nonTE_intersect, nonTE_a12, nonTE_a21, nonTE_b22, comparison_type) ;
                all_pvals_alafisher_best.push_back (nonTE_pval_alafisher) ;
                
                // binomial test
                int n_tot_nonTE = my_peak_count ; // number of trials
                int X_obs_nonTE = nonTE_intersect ;
                double p_obs_nonTE = double( mean_peaklen*X_obs_nonTE ) / double(total_nonTE_bp) ;
                double q_obs_nonTE = 1 - p_obs_nonTE ;
                double p_exp_nonTE = double(total_nonTE_bp) / double(size_hg19) ; // Prob to touch fam by random
                double nonTE_pval_binomial = pbinom(X_obs_nonTE-1, n_tot_nonTE, p_exp_nonTE, comparison_type) ;
                all_pvals_binomial_best.push_back (nonTE_pval_binomial) ;

                all_te_inter_peak.push_back(nonTE_intersect) ;
                all_te_inter_peak_unique.push_back(nonTE_intersect) ;
                all_total_fam.push_back(total_nonTE);
                all_total_bp_fam.push_back(total_nonTE_bp);
                all_avg_fam_size.push_back(nonTE_avg_size);
                all_fam_genome_ratio.push_back(nonTE_genome_ratio);

                // Get ajusted p-values with Benjamin-Hochsberg method 
                std::unordered_map<std::string, double> padj_hygm_reg_fam = calc_adj_pval(&all_pvals_best, &fam_names, print_padj) ;
                std::unordered_map<std::string, double> padj_hygm_alaFish_fam = calc_adj_pval(&all_pvals_alafisher_best, &fam_names, print_padj) ;
                std::unordered_map<std::string, double> padj_hygm_binom_fam = calc_adj_pval(&all_pvals_binomial_best, &fam_names, print_padj) ;
               
                // Printing final results with 3 adjusted p-values
                print_all_pval_adj(matrix_fam, tag_name,
                        &prhead_sum_all_best_fam, &prhead_mat_all_best_fam, &all_pvals_binomial_best,
                        padj_hygm_reg_fam, padj_hygm_alaFish_fam, padj_hygm_binom_fam,
                        &fam_names, &all_te_inter_peak_unique,
                        &all_total_fam,  &all_te_inter_peak,
                        my_peak_count_on_te_unique, my_peak_count, &all_total_bp_fam,
                        &all_avg_fam_size, &all_fam_genome_ratio,
                        peak_total_Mbp, peak_ratio_on_te, ratio_genome_peak,
                        "te_fam", &prhead_summary_te, out_path ) ;
            }
            ref_in_fam.close() ;
            */
            prhead_summary_te = 0 ;
        }
        list_f2.close() ;
        

        summary_bed.close() ;
    }

    // Cleaning temp directory 
    //system("rm -rf temp_TEnrich/") ;
}
