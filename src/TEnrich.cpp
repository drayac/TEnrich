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
    const int LINE_MAX_SIZE = 200;
    const double FDR = 0.05 ;
    const int te_data_size = 4570939 ;
    bool print_padj = 1 ;
    int bed_n_field ; 
    std::string this_dir = <FOLDER_INSTALL> ;
    
    // Get help
    if ( argc > 1 ){
        std::string firstparam = std::string(argv[1]) ;
        if ( firstparam == "-h" || firstparam == "--help" ){
            print_help() ;
            exit (EXIT_FAILURE) ;
        }
        if ( firstparam == "-v" || firstparam == "--version" ){
            std::cout << "TEnrich V0.7 – written by Alexandre Coudray at the EPFL (2019)\n\n" ;
            exit (EXIT_FAILURE) ;
        }
    }

    // Get parameters
    std::vector <std::string> input = get_parameters(argc,argv) ; 
    std::string bed_path = input[0] ; std::string out_path = input[1] ;
    long int size_hg19 = std::stol(input[2]) ; 
    std::string comparison_direction = input[3] ;
    std::string comparison_type = input[4] ;
    std::string print_padj_str = input[5] ;
    std::string ref_file = input[6] ;
    std::string ref_file_fam = input[7] ;
    std::string ref_file_clust = input[8] ;
    std::string te_data = input[9] ;
    std::string single_file = input[10] ;
    int idx_col = std::stoi(input[11]) ;

    bool comp_clust = true ;

    if ( print_padj_str.compare("false") == 0 ){ print_padj = 0 ; }
    
    //check_folder(bed_path) ; 
    check_folder(this_dir) ;

    // nonTE parameters 
    int total_nonTE = 4433186 ; // we estimate the number of nonTE interval ~= TE interval after merge
    int total_nonTE_bp = size_hg19 - 1434913179 ; // size_hg19 - genome span of TE    
    double nonTE_avg_size = double(total_nonTE_bp) / double(total_nonTE) ;
    double nonTE_genome_ratio = double(total_nonTE_bp) / double(size_hg19) ; 
    
    //Check TE data number of fields
    std::cout << "Checking that TE data is ok\n" ;
    int expect_n_fields = 9 ;
    check_te_database(te_data, expect_n_fields ) ;
    
    // Create output directory if not existant
    system("mkdir -p temp_TEnrich") ;
 
    // Creating useful folders
    create_folder(out_path,"summary_bed") ; create_folder(out_path,"summary_te_fam") ;
    create_folder(out_path,"summary_te_subfam") ; create_folder(out_path,"summary_te_clust") ;

    // Initialize Hash Tables
    std::unordered_map<std::string, int> peak_count , peak_count_on_te , peak_len , peak_total_bp , te_inter_peak , te_inter_peak_unique, teFam_inter_peak, teFam_inter_peak_unique, teClust_inter_peak, teClust_inter_peak_unique ;
    std::unordered_map<std::string, std::string> summary_peak_line ;

    std::string concat_bed = "temp_TEnrich/concat_all.bed" ;
    std::string list_files = "temp_TEnrich/temp.list_files.txt" ;
    if ( bed_path.compare("empty") != 0 ){
        // CONCATENATE MULTIPLE FILE WITH A TAG
        std::cout << "Concatening all bed files...\n" << std::endl ;
        // make a list of files in the folder given as argv[1]
        std::stringstream ls_cmd ;
        ls_cmd << "ls -1 "<< bed_path << "/* > " << list_files ;
        system(&(ls_cmd.str()[0])) ;
        std::ifstream list_f(list_files);

        if (this_is_empty(list_f)){
            std::cout << "Problem while opening " << list_files << ", exiting...\n" ;
            exit (EXIT_FAILURE) ;
        }
        if (list_f.is_open()) {
            std::string line;
            while (getline(list_f, line)) {
                // initialize peak variables
                std::string tag_name = remove_ext(base_name(line)) ;
                // local variables
                int peak_n = 0 ; long long total_length = 0 ; int prev_bed_n_field ; bool first_check = 1 ;

                std::ifstream this_bed(line) ;
                if (this_is_empty(this_bed)){
                    std::cout << "Problem while opening " << line << ", exiting...\n" ;
                    exit (EXIT_FAILURE) ;
                }
                if (this_bed.is_open()) {

                    std::string line_bed;
                    std::ofstream my_out ;
                    my_out.open (concat_bed, std::fstream::app) ;

                    while (getline(this_bed, line_bed)) {
                        // Write line in concat_bed file
                        line_bed.erase(std::remove(line_bed.begin(), line_bed.end(), '\n'), line_bed.end());
                        //my_out << line_bed << "\t" << tag_name << std::endl ;

                        // Get fields to calculate peak length
                        std::istringstream iss(line_bed) ;
                        std::vector <std::string> fields ;
                        std::string field ;
                        while(std::getline(iss, field, '\t')){ 
                            fields.push_back(field);
                        }

                        // Check bed fields numbers
                        bed_n_field = fields.size() ;
                        if ( first_check ){ prev_bed_n_field = bed_n_field ; first_check = 0 ; }
                        prev_bed_n_field = bed_n_field ;

                        // WRITE OUTPUT 3 FIELD BED
                        my_out << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << tag_name << "\n" ;

                        // Record peak length
                        int peak_length = std::stoi(fields[2]) - std::stoi(fields[1]) ;
                        total_length += peak_length  ;
                        peak_n += 1 ;
                    }
                    my_out.close() ;
                }
                // calculate peak stats
                int peak_average_length = total_length / peak_n ;
                double peak_genome_ratio = double(total_length) / double(size_hg19) ;
                double peak_total_Mbp = double(total_length) / 1e6 ;

                // add stats to map hash 
                peak_count.insert({tag_name, peak_n}) ;   
                peak_len.insert({tag_name, peak_average_length}) ;
                peak_total_bp.insert({tag_name, total_length}) ;

                // insert to summary peak
                std::string summary_line = tag_name + "\t" + std::to_string(peak_n) + "\t" + std::to_string(peak_average_length) + "\t" + std::to_string(peak_total_Mbp) + "\t" + std::to_string(peak_genome_ratio) ;
                summary_peak_line.insert({tag_name, summary_line }) ;

                std::cout << tag_name << " has " << peak_n << " peaks with an average size of " << peak_average_length << " bp. Total length = " << total_length << "\n" ;
            }
            list_f.close();
        }
        idx_col = 12 ;
    }
    else
    {
        concat_bed = single_file ;
        idx_col += 9 ;
    }

    
    // Intersect concat_bed with te_data file
    std::cout << "\nIntersect with TE database...\n" << std::endl ;
    std::stringstream bedtools_cmd ;
    std::string inter_bed_path = "temp_TEnrich/temp.out.bed" ;
    bedtools_cmd << "bedtools intersect -a "<< te_data << " -b "<< concat_bed << " -f 0.5 -F 0.5 -e -wa -wb | sort -k1,1 -k2,2n -k9,9 -k10,10n > " << inter_bed_path ;
    system(&(bedtools_cmd.str()[0])) ; 

    // open list_files to write sample names (only used if single_file option provided) 
    std::ofstream fake_list ; 
    std::unordered_map<std::string, int> save_fake_list ;
    if ( bed_path.compare("empty") == 0 ){ fake_list.open(list_files, std::fstream::app) ; }
    
    // Parse intersect and get all counts
    std::cout << "Counting peaks vs TE overlaps...\n" << std::endl ;
    std::ifstream inter_bed(inter_bed_path);
    if (this_is_empty(inter_bed)){
        std::cout << "Problem while opening " << inter_bed_path << ", exiting...\n" ;
        exit (EXIT_FAILURE) ;
    }
    if (inter_bed.is_open()) {
        std::string line, prev_te, prev_peak, prev_key , prev_key_fam, prev_key_clust , prev_tag_name ;
        while (getline(inter_bed, line)) { 
            std::istringstream iss(line) ;
            std::vector <std::string> fields ;
            std::string field ;
            while(std::getline(iss, field, '\t')){ 
                fields.push_back(field);
            }
    
            int total_field = 8 ;
            std::string tag_name = fields[idx_col] ;
            if ( bed_path.compare("empty") == 0 ){
                if ( idx_col < 0 ){
                    tag_name = fields[0] + ":" + fields[1] + "-" + fields[2] ;
                    if ( save_fake_list.find(tag_name) == save_fake_list.end() ){
                        fake_list << tag_name << "\n" ;
                        save_fake_list.insert({tag_name, 1}) ;
                    }
                }
                else
                {
                    tag_name = fields[idx_col] ;
                    if ( save_fake_list.find(tag_name) == save_fake_list.end() ){
                        fake_list << tag_name << "\n" ;
                        save_fake_list.insert({tag_name, 1}) ;
                    }
                }
            }
            else
            {
                tag_name = fields[idx_col] ;
            }

            std::string key = fields[7] + "_" + tag_name ;
            std::string key_fam = fields[6] + "_" + tag_name ;
            std::string this_te = fields[0] + fields[1] + fields[2] ;
            std::string this_peak = fields[total_field+1] + fields[total_field+2] + fields[total_field+3] ; 
            te_inter_peak[key]++ ; teFam_inter_peak[key_fam]++ ; 
            
            std::string clust = fields[8] ;
            std::string key_clust = clust + "_" + tag_name ;
            teClust_inter_peak[key_clust]++ ;

            peak_count_on_te[tag_name]++ ;

            if ( this_te.compare(prev_te) != 0 ){
                te_inter_peak_unique[prev_key]++ ; teFam_inter_peak_unique[prev_key_fam]++ ;
                teClust_inter_peak_unique[prev_key_clust]++ ;
            }
            prev_te = this_te ; prev_peak = this_peak ; 
            prev_key = key ; prev_key_fam = key_fam ; 
            prev_tag_name = tag_name ;
            prev_key_clust = key_clust ; 
        }
    }
    inter_bed.close() ; fake_list.close() ;

    // Get median / quantiles for subfam 


    // Parse files and compute pval enrichments
    // open matrix out
    bool prhead_mat_all_best = 1 ; 
    std::string matrix_path_all_best = out_path + "/matrix_padjBinomial_Subfam.txt" ;
    bool prhead_mat_all_best_fam = 1 ; 
    std::string matrix_fam = out_path + "/matrix_padjBinomial_Fam.txt" ;
    bool prhead_mat_clust = 1 ;
    std::string matrix_clust = out_path + "/matrix_padjBinomial_Clusters.txt" ;

    // Create peak summary file
    std::string summary_bed_path = out_path + "/summary_bed.txt" ;
    std::ofstream summary_bed ; 
    summary_bed.open(summary_bed_path, std::fstream::app) ;

    // Initialize bool to print summary_fam/subfam header 
    bool prhead_summary_te  = 1 ;

    // print header 
    summary_bed << "sample_name\tpeak.number\tpeak.average.size\tpeak.total.bp\tpeak.ratio.genome\tpeak.overlap.te\tpeak.ratio.overlap.te\n" ;
   
    // Iterate over samples (list of files in input folder) 
    std::ifstream list_f2(list_files);
    if (this_is_empty(list_f2)){
        std::cout << "Problem while opening " << list_files << ", exiting...\n" ;
        exit (EXIT_FAILURE) ;
    }
    if (list_f2.is_open()) {
        std::string line_ref;
        while (getline(list_f2, line_ref)) {
            
            std::string tag_name = remove_ext(base_name(line_ref)) ;
            std::cout << "Enrichment analysis of " << tag_name << "\n" ;
            
            // getting peak stats
            int mean_peaklen = peak_len[tag_name] , my_peak_count = peak_count[tag_name] , my_peak_count_on_te = peak_count_on_te[tag_name] , peak_tot_bp = peak_total_bp[tag_name] ;
            double ratio_genome_peak = double(peak_tot_bp) / double(size_hg19) ;
            double peak_total_Mbp = double(peak_tot_bp) / 1000000 ;
            double peak_ratio_on_te = double(my_peak_count_on_te) / double(my_peak_count) ; 
            
            // add info to peak summary 
            std::string peak_sum_line = summary_peak_line[tag_name] ;
            summary_bed << peak_sum_line << "\t" << my_peak_count_on_te << "\t" << peak_ratio_on_te << "\n" ;
            // bool to print header of summmary files 
            bool prhead_sum_all_best = 1, prhead_sum_all_best_fam = 1 , prhead_sum_clust = 1 ;
           
            //////////////////// 
            /*  TE SUBFAM     */
            ////////////////////
            
            // Begin enrichment analysis for sample X 
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
                    if (te_inter_peak.find(key) == te_inter_peak.end() || te_inter_peak[key] == 0){
                        subfam_names.push_back (fields[0]) ; all_pvals.push_back (1.0) ;
                        all_pvals_alafisher.push_back (1.0) ; all_pvals_binomial.push_back (1.0) ;
                        all_pvals_rev.push_back (1.0) ; all_pvals_alafisher_rev.push_back (1.0) ;
                        all_pvals_binomial_rev.push_back (1.0) ; all_pvals_best.push_back (1.0) ;
                        all_pvals_alafisher_best.push_back (1.0) ; all_pvals_binomial_best.push_back (1.0) ;
                    }
                    else
                    {
                        ////////////////////////////////////// 
                        /* PVAL ENRICHMENT PEAKs AMONG TEs  */
                        //////////////////////////////////////
                        
                        // Regular HyperGeometric
                        long long a_11 = te_inter_peak[key] ;
                        long long a_12 = MAX(0L,n_subfam - a_11) ;
                        long long a_21 = MAX(0L,my_peak_count_on_te - a_11) ;
                        long long a_22 = te_data_size - a_11 - a_12 - a_21 ;
                        double pval = fisher_exact(a_11, a_12, a_21, a_22, comparison_type) ;
                        all_pvals.push_back (pval) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        int total_average = mean_peaklen + avg_subfam_size ;
                        long long b_11 = te_inter_peak[key] ;
                        long long b_12 = MAX(a_12, n_subfam - b_11 ) ;
                        long long b_21 = MAX(0L,my_peak_count - b_11) ;
                        long long b_22 = MAX( int(double(size_hg19) / double(total_average)) - b_11 - b_21 - b_12, 1 )  ;
                        double pval_alafisher = fisher_exact(b_11, b_12, b_21, b_22, comparison_type) ;
                        all_pvals_alafisher.push_back (pval_alafisher) ;

                        // Binomial exact test with genome ratio as p
                        int n_tot = my_peak_count ; // number of trials
                        int X_obs = te_inter_peak[key] ;
                        double p_obs = double( mean_peaklen*X_obs ) / double(total_bp_length_subfam) ;
                        double q_obs = 1 - p_obs ;
                        double p_exp = double(total_bp_length_subfam) / double(size_hg19) ; // Prob to touch subfam by random
                        double pval_binomial = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial.push_back(pval_binomial) ;

                        //////////////////////////////////////
                        /* PVAL ENRICHMENT TEs AMONG PEAKs  */
                        //////////////////////////////////////

                        // Regular HyperGeometric
                        long long ar_11 = te_inter_peak_unique[key] ;
                        long long ar_12 = MAX(0L,n_subfam - ar_11) ;
                        long long ar_21 = MAX(0L,my_peak_count_on_te - ar_11) ;
                        long long ar_22 = te_data_size - ar_11 - ar_12 - ar_21 ;
                        double pval_rev = fisher_exact(ar_11, ar_12, ar_21, ar_22, comparison_type) ;
                        all_pvals_rev.push_back (pval_rev) ;
                        if ( mean_peaklen > avg_subfam_size || comparison_direction.compare("te_in_peak") == 0 ) all_pvals_best.push_back (pval_rev) ;
                        if ( mean_peaklen <= avg_subfam_size || comparison_direction.compare("peak_in_te") == 0 ) all_pvals_best.push_back (pval) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        total_average = avg_subfam_size ;
                        long long br_11 = te_inter_peak_unique[key] ;
                        long long br_12 = MAX(ar_12, n_subfam - br_11 ) ;
                        long long br_21 = MAX(0L,my_peak_count - br_11 ) ;
                        long long br_22 = MAX( int(double(size_hg19) / double(total_average)) - br_11 - br_21 - br_12, 1)  ;
                        double pval_alafisher_rev = fisher_exact(br_11, br_12, br_21, b_22, comparison_type) ;
                        all_pvals_alafisher_rev.push_back (pval_alafisher_rev) ;
                        if ( mean_peaklen > avg_subfam_size || comparison_direction.compare("te_in_peak") == 0 ){
                            all_pvals_alafisher_best.push_back (pval_alafisher_rev) ;
                        }
                        if ( mean_peaklen <= avg_subfam_size || comparison_direction.compare("peak_in_te") == 0 ){
                            all_pvals_alafisher_best.push_back (pval_alafisher) ;
                        }
                        
                        // Binomial exact test with genome ratio as p 
                        n_tot = n_subfam ; // number of trials
                        X_obs = te_inter_peak_unique[key] ;
                        p_exp = ratio_genome_peak ; // Prob to touch subfam by random
                        double pval_binomial_rev = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial_rev.push_back(pval_binomial_rev) ;
                        if ( mean_peaklen > avg_subfam_size || comparison_direction.compare("te_in_peak") == 0 ) all_pvals_binomial_best.push_back (pval_binomial_rev) ;
                        if ( mean_peaklen <= avg_subfam_size || comparison_direction.compare("peak_in_te") == 0 ) all_pvals_binomial_best.push_back (pval_binomial) ;
                        subfam_names.push_back (fields[0]) ;
                    }
                }

                // Getting pval for nonTE
                subfam_names.push_back ("nonTE") ;
                int nonTE_intersect = my_peak_count - my_peak_count_on_te ;
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
                        my_peak_count_on_te, my_peak_count, &all_total_bp_subfam,
                        &all_avg_subfam_size, &all_subfam_genome_ratio,
                        peak_total_Mbp, peak_ratio_on_te, ratio_genome_peak,
                        "te_subfam", &prhead_summary_te, out_path ) ;
            }
            ref_in.close() ;

// // / // / / /
//
//
//
//
//
//
//
//
//
            if ( comp_clust ){
            //////////////////// 
            /*  TE CLUSTERS   */
            ////////////////////
            
            // Begin enrichment analysis for sample X 
            std::ifstream ref_in_clust(ref_file_clust);
            if (this_is_empty(ref_in_clust)){
                std::cout << "Problem while opening " << ref_file_clust << ", exiting...\n" ;
                exit (EXIT_FAILURE) ;
            }
            if (ref_in_clust.is_open()) {
                
                // initialize vectors
                std::vector<std::string> clust_names;
                std::vector<double> all_pvals , all_pvals_alafisher , all_pvals_binomial, all_pvals_rev , all_pvals_alafisher_rev , all_pvals_binomial_rev , all_pvals_best , all_pvals_alafisher_best , all_pvals_binomial_best, all_total_bp_clust, all_clust_genome_ratio ;
                std::vector<int> all_teClust_inter_peak, all_teClust_inter_peak_unique , all_total_clust , all_avg_clust_size;

                // Iterate over clust names of TEs (or sample of bed file 2) 
                std::string line;
                while (getline(ref_in_clust, line)) { 
                    std::istringstream iss(line) ;
                    std::vector <std::string> fields ;
                    std::string field ;
                    while(std::getline(iss, field, '\t')){ 
                        fields.push_back(field);
                    }

                    std::string clust_name = fields[0] ;
                    std::string key = clust_name + "_" + tag_name ;
                    int n_clust = std::stoi(fields[1]) ;
                    int total_bp_length_clust = std::stoi(fields[2]) ;
                    double avg_clust_size_db = std::stod(fields[3]) ;
                    int avg_clust_size = int(avg_clust_size_db) ;
                    double ratio_genome_clust = double(total_bp_length_clust)/double(size_hg19) ;
                    double total_Mbp_len_clust = double(total_bp_length_clust) / 1000000 ;

                    all_teClust_inter_peak.push_back(teClust_inter_peak[key]) ;
                    all_teClust_inter_peak_unique.push_back(teClust_inter_peak_unique[key]) ;
                    all_total_clust.push_back(n_clust);
                    all_total_bp_clust.push_back(total_Mbp_len_clust);
                    all_avg_clust_size.push_back(avg_clust_size);
                    all_clust_genome_ratio.push_back(ratio_genome_clust);

                    // Calculation of the p-values for each 1-1 relation
                    if (teClust_inter_peak.find(key) == teClust_inter_peak.end() || teClust_inter_peak[key] == 0){
                        clust_names.push_back (fields[0]) ; all_pvals.push_back (1.0) ;
                        all_pvals_alafisher.push_back (1.0) ; all_pvals_binomial.push_back (1.0) ;
                        all_pvals_rev.push_back (1.0) ; all_pvals_alafisher_rev.push_back (1.0) ;
                        all_pvals_binomial_rev.push_back (1.0) ; all_pvals_best.push_back (1.0) ;
                        all_pvals_alafisher_best.push_back (1.0) ; all_pvals_binomial_best.push_back (1.0) ;
                    }
                    else
                    {
                        ////////////////////////////////////// 
                        /* PVAL ENRICHMENT PEAKs AMONG TEs  */
                        //////////////////////////////////////
                        
                        // Regular HyperGeometric
                        long long a_11 = teClust_inter_peak[key] ;
                        long long a_12 = MAX(0L,n_clust - a_11) ;
                        long long a_21 = MAX(0L,my_peak_count_on_te - a_11) ;
                        long long a_22 = te_data_size - a_11 - a_12 - a_21 ;
                        double pval = fisher_exact(a_11, a_12, a_21, a_22, comparison_type) ;
                        all_pvals.push_back (pval) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        int total_average = mean_peaklen + avg_clust_size ;
                        long long b_11 = teClust_inter_peak[key] ;
                        long long b_12 = MAX(a_12, n_clust - b_11 ) ;
                        long long b_21 = MAX(0L,my_peak_count - b_11) ;
                        long long b_22 = MAX( int(double(size_hg19) / double(total_average)) - b_11 - b_21 - b_12, 1 )  ;
                        double pval_alafisher = fisher_exact(b_11, b_12, b_21, b_22, comparison_type) ;
                        all_pvals_alafisher.push_back (pval_alafisher) ;

                        // Binomial exact test with genome ratio as p
                        int n_tot = my_peak_count ; // number of trials
                        int X_obs = teClust_inter_peak[key] ;
                        double p_obs = double( mean_peaklen*X_obs ) / double(total_bp_length_clust) ;
                        double q_obs = 1 - p_obs ;
                        double p_exp = double(total_bp_length_clust) / double(size_hg19) ; // Prob to touch clust by random
                        double pval_binomial = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial.push_back(pval_binomial) ;

                        //////////////////////////////////////
                        /* PVAL ENRICHMENT TEs AMONG PEAKs  */
                        //////////////////////////////////////

                        // Regular HyperGeometric
                        long long ar_11 = teClust_inter_peak_unique[key] ;
                        long long ar_12 = MAX(0L,n_clust - ar_11) ;
                        long long ar_21 = MAX(0L,my_peak_count_on_te - ar_11) ;
                        long long ar_22 = te_data_size - ar_11 - ar_12 - ar_21 ;
                        double pval_rev = fisher_exact(ar_11, ar_12, ar_21, ar_22, comparison_type) ;
                        all_pvals_rev.push_back (pval_rev) ;
                        if ( mean_peaklen > avg_clust_size || comparison_direction.compare("te_in_peak") == 0 ) all_pvals_best.push_back (pval_rev) ;
                        if ( mean_peaklen <= avg_clust_size || comparison_direction.compare("peak_in_te") == 0 ) all_pvals_best.push_back (pval) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        total_average = avg_clust_size ;
                        long long br_11 = teClust_inter_peak_unique[key] ;
                        long long br_12 = MAX(ar_12, n_clust - br_11 ) ;
                        long long br_21 = MAX(0L,my_peak_count - br_11 ) ;
                        long long br_22 = MAX( int(double(size_hg19) / double(total_average)) - br_11 - br_21 - br_12, 1)  ;
                        double pval_alafisher_rev = fisher_exact(br_11, br_12, br_21, b_22, comparison_type) ;
                        all_pvals_alafisher_rev.push_back (pval_alafisher_rev) ;
                        if ( mean_peaklen > avg_clust_size || comparison_direction.compare("te_in_peak") == 0 ){
                            all_pvals_alafisher_best.push_back (pval_alafisher_rev) ;
                        }
                        if ( mean_peaklen <= avg_clust_size || comparison_direction.compare("peak_in_te") == 0 ){
                            all_pvals_alafisher_best.push_back (pval_alafisher) ;
                        }
                        
                        // Binomial exact test with genome ratio as p 
                        n_tot = n_clust ; // number of trials
                        X_obs = teClust_inter_peak_unique[key] ;
                        p_exp = ratio_genome_peak ; // Prob to touch clust by random
                        double pval_binomial_rev = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial_rev.push_back(pval_binomial_rev) ;
                        if ( mean_peaklen > avg_clust_size || comparison_direction.compare("te_in_peak") == 0 ) all_pvals_binomial_best.push_back (pval_binomial_rev) ;
                        if ( mean_peaklen <= avg_clust_size || comparison_direction.compare("peak_in_te") == 0 ) all_pvals_binomial_best.push_back (pval_binomial) ;

                        clust_names.push_back (fields[0]) ;
                    }
                }

                // Get ajusted p-values with Benjamin-Hochsberg method 
                std::unordered_map<std::string, double> padj_hygm_reg = calc_adj_pval(&all_pvals_best, &clust_names, print_padj) ;
                std::unordered_map<std::string, double> padj_hygm_alaFish = calc_adj_pval(&all_pvals_alafisher_best, &clust_names, print_padj) ;
                std::unordered_map<std::string, double> padj_hygm_binom = calc_adj_pval(&all_pvals_binomial_best, &clust_names, print_padj) ;
                
                // Printing final results with 3 adjusted p-values
                print_all_pval_adj(matrix_clust, tag_name,
                        &prhead_sum_clust, &prhead_mat_clust, &all_pvals_binomial_best,
                        padj_hygm_reg, padj_hygm_alaFish, padj_hygm_binom,
                        &clust_names, &all_teClust_inter_peak_unique,
                        &all_total_clust,  &all_teClust_inter_peak,
                        my_peak_count_on_te, my_peak_count, &all_total_bp_clust,
                        &all_avg_clust_size, &all_clust_genome_ratio,
                        peak_total_Mbp, peak_ratio_on_te, ratio_genome_peak,
                        "te_clust", &prhead_summary_te, out_path ) ;
            }
            ref_in_clust.close() ;

            }

            ////
            ///
            //
            //
            //
            //
            //
            //
            //

            ///////////////// 
            /*  TE FAM     */
            /////////////////
            // Begin enrichment analysis for sample X 
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
                        /* PVAL ENRICHMENT PEAKs AMONG TEs  */
                        //////////////////////////////////////
                        
                        // Regular HyperGeometric
                        int a_11 = teFam_inter_peak[key] ;
                        int a_12 = MAX(0L,n_fam - a_11) ;
                        int a_21 = MAX(0L,my_peak_count_on_te - a_11) ;
                        int a_22 = te_data_size - a_11 - a_12 - a_21 ;
                        double pval = fisher_exact(a_11, a_12, a_21, a_22, comparison_type) ;
                        all_pvals.push_back (pval) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        int total_average = mean_peaklen + avg_fam_size ;
                        int b_11 = teFam_inter_peak[key] ;
                        int b_12 = MAX(a_12, n_fam - b_11 ) ;
                        int b_21 = MAX(0L,my_peak_count - b_11) ;
                        int b_22 = MAX( int(double(size_hg19) / double(total_average)) - b_11 - b_21 - b_12, 1 )  ;
                        double pval_alafisher = fisher_exact(b_11, b_12, b_21, b_22, comparison_type) ;
                        all_pvals_alafisher.push_back (pval_alafisher) ;

                        // Binomial exact test with genome ratio as p
                        int n_tot = my_peak_count ; // number of trials
                        int X_obs = teFam_inter_peak[key] ;
                        double p_obs = double( mean_peaklen*X_obs ) / double(total_bp_length_fam) ;
                        double q_obs = 1 - p_obs ;
                        double p_exp = double(total_bp_length_fam) / double(size_hg19) ; // Prob to touch fam by random
                        double pval_binomial = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial.push_back(pval_binomial) ;

                        //////////////////////////////////////
                        /* PVAL ENRICHMENT TEs AMONG PEAKs  */
                        //////////////////////////////////////

                        // Regular HyperGeometric
                        int ar_11 = teFam_inter_peak_unique[key] ;
                        int ar_12 = MAX(0L,n_fam - ar_11) ;
                        int ar_21 = MAX(0L,my_peak_count - ar_11) ;
                        int ar_22 = te_data_size - ar_11 - ar_12 - ar_21 ;
                        double pval_rev = fisher_exact(ar_11, ar_12, ar_21, ar_22, comparison_type) ;
                        all_pvals_rev.push_back (pval_rev) ;
                        if ( mean_peaklen > avg_fam_size || comparison_direction.compare("te_in_peak") == 0 ) all_pvals_best.push_back (pval_rev) ;
                        if ( mean_peaklen <= avg_fam_size || comparison_direction.compare("peak_in_te") == 0 ) all_pvals_best.push_back (pval) ;

                        // HyperGeometric with genome occupency - 'ala bedtools fisher'
                        total_average = avg_fam_size ;
                        int br_11 = teFam_inter_peak_unique[key] ;
                        int br_12 = MAX(ar_12, n_fam - br_11 ) ;
                        int br_21 = MAX(0L,my_peak_count - br_11 ) ;
                        int br_22 = MAX( int(double(size_hg19) / double(total_average)) - br_11 - br_21 - br_12, 1)  ;
                        double pval_alafisher_rev = fisher_exact(br_11, br_12, br_21, b_22, comparison_type) ;
                        all_pvals_alafisher_rev.push_back (pval_alafisher_rev) ;
                        if ( mean_peaklen > avg_fam_size || comparison_direction.compare("te_in_peak") == 0 ) all_pvals_alafisher_best.push_back (pval_alafisher_rev) ;
                        if ( mean_peaklen <= avg_fam_size || comparison_direction.compare("peak_in_te") == 0 ) all_pvals_alafisher_best.push_back (pval_alafisher) ;
                        
                        // Binomial exact test with genome ratio as p 
                        n_tot = n_fam ; // number of trials
                        X_obs = teFam_inter_peak_unique[key] ;
                        p_exp = ratio_genome_peak ; // Prob to touch peak by random
                        double pval_binomial_rev = pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
                        all_pvals_binomial_rev.push_back(pval_binomial_rev) ;
                        if ( mean_peaklen > avg_fam_size || comparison_direction.compare("te_in_peak") == 0 ) all_pvals_binomial_best.push_back (pval_binomial_rev) ;
                        if ( mean_peaklen <= avg_fam_size || comparison_direction.compare("peak_in_te") == 0 ) all_pvals_binomial_best.push_back (pval_binomial) ;
                        
                        fam_names.push_back (fields[0]) ;
                    }
                }

                // Getting pval for nonTE
                fam_names.push_back ("nonTE") ;
                //int total_nonTE = 4570939 ; // we estimate the number of nonTE interval ~= TE interval
                //int total_nonTE_bp = 1260181506 ; // size_hg19 - genome span of TE   
                int nonTE_intersect = my_peak_count - my_peak_count_on_te ;

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
                        my_peak_count_on_te, my_peak_count, &all_total_bp_fam,
                        &all_avg_fam_size, &all_fam_genome_ratio,
                        peak_total_Mbp, peak_ratio_on_te, ratio_genome_peak,
                        "te_fam", &prhead_summary_te, out_path ) ;
            }
            ref_in_fam.close() ;
            prhead_summary_te = 0 ;
        }
    }
    list_f2.close() ;
    summary_bed.close() ;

    // Cleaning temp directory 
    system("rm -rf temp_TEnrich/") ;
}
