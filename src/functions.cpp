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

#define constant_e (2.71828)
#define PI (3.14159265359)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

std::string this_dir = <FOLDER_INSTALL> ;
std::string version  = <VERSION> ;

void add_param(std::string &var, std::string name, int &i, int c_argc, char *c_argv[]){
    if (std::string(c_argv[i]) == name) {
        if (i + 1 < c_argc) { // Make sure we aren't at the end of argv!
            var = c_argv[i + 1]; 
            i++ ; // Increment 'i' so we don't get the argument as the next argv[i].
        } else { 
            std::cerr << name << " option requires one argument." << std::endl;
            exit (EXIT_FAILURE) ;
        }
    }
}

std::vector <std::string> get_parameters(int c_argc, char *c_argv[]){
    if (c_argc < 1) {
        std::cerr << "Not enough parameters were given, exiting ... " << std::endl;
        exit (EXIT_FAILURE) ;
    }
    // initialize default values
    std::vector <std::string> named_params ;
    std::string bed_dir = "empty", out_dir = "empty", genome_size = "3088269832" ;
    std::string padj = "true" , comp_sense = "auto" ; 
    std::string type_hypergeom = "ala_bedtools_fisher" , stat_test_type = "greater" ;
    std::string ref_subfam = this_dir + "/db/Subfam_ref_TE.txt" ;
    std::string ref_fam = this_dir + "/db/Fam_ref_TE.txt" ;
    std::string ref_cluster = this_dir + "/db/Clusters_ref_TE.txt" ;
    std::string te_database = this_dir + "/db/hg19_TE_repmask_LTRm_s_20140131.bed.gz" ;
    std::string single_file = "empty" , idx_col = "-1" ;
    std::string destination ; 
    int i ;
    for (i = 1; i < c_argc; ++i) {
        std::string this_param = std::string(c_argv[i]) ;
        add_param(bed_dir,      "--bed_dir",    i,  c_argc, c_argv) ;
        add_param(out_dir,      "--out_dir",    i,  c_argc, c_argv) ;
        add_param(type_hypergeom,  "--type_hypergeom",i,  c_argc, c_argv) ;
        add_param(comp_sense,   "--comp_sense", i,  c_argc, c_argv) ;
        add_param(stat_test_type,"--stat_test_type",i,c_argc,c_argv) ;
        add_param(padj,         "--padj",       i,  c_argc, c_argv) ;
        add_param(ref_subfam,   "--ref_subfam", i,  c_argc, c_argv) ;
        add_param(ref_fam,      "--ref_fam",    i,  c_argc, c_argv) ;
        add_param(ref_cluster,  "--ref_cluster",i,  c_argc, c_argv) ;
        add_param(te_database,  "--te_database",i,  c_argc, c_argv) ;
        add_param(single_file,  "--single_file",i,  c_argc, c_argv) ;
        add_param(idx_col,      "--idx_col",    i,  c_argc, c_argv) ;
        add_param(genome_size,  "--genome_size",i,  c_argc, c_argv) ;

    }

    // check formats of a few parameters
    if ( out_dir.compare("empty") == 0 ){
        std::cerr << "bed and out folder should be precised ... exiting\n" ;
        exit (EXIT_FAILURE) ;
    }
    if ( not is_digits(genome_size) ){
        std::cerr << "size genome parameters should only contain digits !! exiting\n" ;
        exit (EXIT_FAILURE) ;
    }
    if ( comp_sense.compare("peak_in_te") != 0 && comp_sense.compare("te_in_peak") != 0 && comp_sense.compare("auto") != 0 ){
        std::cerr << "option comp_sense should be either 'peak_in_te', 'te_in_peak' or 'auto'. By default 'auto'. Exiting...\n\n" ; 
        exit ( EXIT_FAILURE ) ;
    }
    if ( stat_test_type.compare("greater") != 0 && stat_test_type.compare("less") != 0){
        std::cerr << "option stat_test_type should be either 'greater', 'less'. By default 'greater'. Exiting...\n\n" ; 
        exit ( EXIT_FAILURE ) ;
    }
    if ( padj.compare("true") != 0 && padj.compare("false") != 0){
        std::cerr << "option padj should be 'true' or 'false'. By default 'true'. Exiting...\n\n" ; 
        exit ( EXIT_FAILURE ) ;
    }

    // add params to output vector 
    named_params.push_back(bed_dir) ;       named_params.push_back(out_dir) ;
    named_params.push_back(type_hypergeom) ;named_params.push_back(comp_sense) ;
    named_params.push_back(stat_test_type) ;
    named_params.push_back(padj) ;          named_params.push_back(ref_subfam) ;
    named_params.push_back(ref_fam) ;       named_params.push_back(ref_cluster) ;
    named_params.push_back(te_database);    named_params.push_back(single_file) ;
    named_params.push_back(idx_col) ;       named_params.push_back(genome_size) ;
    
    return named_params ;
}

void get_help(std::string firstparam, int argc_, std::string version ){
    if ( argc_ > 1 ){
        if ( firstparam == "-h" || firstparam == "--help" ){
            print_help(50,50,"no_detail") ;
            exit (EXIT_FAILURE) ;
        }
        if ( firstparam == "-hd" ){
            print_help(50,50,"detail") ;
            exit (EXIT_FAILURE) ;
        }
        if ( firstparam == "-v" || firstparam == "--version" ){
            std::cout << "TEnrich v" << version << "– written by Alexandre Coudray at the EPFL (2019)\n\n" ;
            exit (EXIT_FAILURE) ;
        }
    }
}

void print_help_line(std::string header, std::string long_line, int X, int Y){
    int left_char = long_line.size() , line_n = 0 ;
    int start_char = 0, length_line = Y , shift = 0 ;
    do {
        start_char = line_n * Y + shift ;
        length_line = MIN(left_char, Y) ;

        std::cout << std::left << std::setw(X) << header <<  long_line.substr(start_char,length_line) ;
        header = "" ;
        if ( start_char + length_line + 1 < long_line.size() ){
            char last_char = long_line[start_char + length_line - 1] ;
            char next_char = long_line[start_char + length_line] ; 
            char sec_next_char = long_line[start_char + length_line + 1] ;
            if ( next_char != ' ' && sec_next_char == ' ' ){
                std::cout << next_char  ;
                shift += 2 ;
            }
            if ( next_char != ' ' && sec_next_char != ' ' && left_char > 2*Y && last_char != ' ' ){
                std::cout << "–"  ;
            }
            if ( next_char == ' ' && sec_next_char != ' ' ){
                shift++ ;
            }
        }
        std::cout << std::endl ;

        line_n++ ;
        left_char = long_line.size() - start_char ;

    } while ( left_char >= Y )  ;
    std::cout << std::endl ;
}


void check_te_database( std::string te_data_ziped, int n_expect_fields){
    // gzip and head first 1k line of te_database
    std::string te_data = "temp_te_head.txt" ;
    std::stringstream gzip_head_cmd ;
    gzip_head_cmd << "gzip -cd "<< te_data_ziped << " | head -1000 > " << te_data ;
    system(&(gzip_head_cmd.str()[0])) ;

    // read and count number of fields 
    std::ifstream te_tab_check(te_data);
    int n_line_check_te = 1 ;
    if (this_is_empty(te_tab_check)){
        std::cout << "Problem while opening " << te_data << ", exiting...\n" ;
        exit (EXIT_FAILURE) ;
    }
    if (te_tab_check.is_open()) {
        std::string line ;
        while (getline(te_tab_check, line)) {
            std::istringstream iss(line) ;
            std::vector <std::string> fields ;
            std::string field ;
            while(std::getline(iss, field, '\t')){
                fields.push_back(field);
            }
            int n_fields_check = fields.size() ;
            if ( n_fields_check < n_expect_fields ){
                std::cout << "TE database has not enough fields, (should contain 9), exiting...\n" ;
                exit ( EXIT_FAILURE ) ;
            }
            else if ( n_fields_check > n_expect_fields ){
                std::cout << "TE database has too many fields, (should contain 9), exiting...\n" ;
                exit ( EXIT_FAILURE ) ;
            }
            n_line_check_te++ ;
            if ( n_line_check_te > 10 ){
                 break ;
            }
        }
    }
    te_tab_check.close() ;
    system("rm temp_te_head.txt") ;
    std::cout << "TE database seems to have to right number of fields\n" ;
}


// concatenate bed files 
void concat_bed_files(  std::string bed_path, 
                        std::string list_files, 
                        std::string concat_bed, 
                        std::unordered_map<std::string, int> &peak_count, 
                        std::unordered_map<std::string, int> &peak_len, 
                        std::unordered_map<std::string, int> &peak_total_bp, 
                        std::unordered_map<std::string, std::string> &summary_peak_line, 
                        long int size_hg19)
{
    std::cout << "Concatening all bed files..." << std::endl ;
    int bed_n_field ;
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

                    // WRITE 3 FIELD BED with TAG
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
            std::cout << "total_length : " << total_length << std::endl ;

            // add stats to map hash
            peak_count.insert({tag_name, peak_n}) ;
            peak_len.insert({tag_name, peak_average_length}) ;
            peak_total_bp.insert({tag_name, total_length}) ;

            // insert to summary peak
            std::string summary_line = tag_name + "\t" + std::to_string(peak_n) + "\t" + std::to_string(peak_average_length) + "\t" + std::to_string(peak_total_Mbp) + "\t" + std::to_string(peak_genome_ratio) ;
            summary_peak_line.insert({tag_name, summary_line }) ;
        }
        list_f.close();
	}
}


void bedtools_intersect(std::string bedtools_options, 
                        std::string sort_options, 
                        std::string te_data, 
                        std::string inter_bed_path, 
                        std::string concat_bed)
{
    std::cout << "Intersect with TE database..." << std::endl ;
    std::stringstream bedtools_cmd ;
    bedtools_cmd << "bedtools intersect -a "<< te_data << " -b "<< concat_bed << " " << bedtools_options << " | " << sort_options << " > " << inter_bed_path ;
    system(&(bedtools_cmd.str()[0])) ; 
}

void sort_peaks(std::string sort_options, 
                std::string inter_bed_path, 
                std::string sort_peaks_bed)
{
    std::cout << "Re-Sort data by Peaks..." << std::endl ;
    std::stringstream sort_cmd ;
    sort_cmd << sort_options << " " << inter_bed_path << " > " << sort_peaks_bed ;
    system(&(sort_cmd.str()[0])) ; 
}


void parse_intersect(   std::string inter_bed_path, 
                        std::string bed_path, 
                        std::unordered_map<std::string,int> &save_fake_list,
                        int idx_col, 
                        std::unordered_map<std::string, int> &te_inter_peak,  
                        std::unordered_map<std::string,int> &teFam_inter_peak,  
                        std::unordered_map<std::string, int> &te_inter_peak_unique, 
                        std::unordered_map<std::string, int> &teFam_inter_peak_unique,
                        std::unordered_map<std::string, int> &peak_count_on_te )
{
    std::cout << "Counting peaks vs TE overlaps..." << std::endl ;
    
    std::ofstream fake_list ; 
    std::ifstream inter_bed(inter_bed_path);
    if (this_is_empty(inter_bed)){
        std::cout << "Problem while opening " << inter_bed_path << ", exiting...\n" ;
        exit (EXIT_FAILURE) ;
    }
    if (inter_bed.is_open()) {
        std::string line, prev_te, prev_peak, prev_key , prev_key_fam, prev_tag_name ;
        while (getline(inter_bed, line)) { 
            std::istringstream iss(line) ;
            std::vector <std::string> fields ;
            std::string field ;
            while(std::getline(iss, field, '\t')){ 
                fields.push_back(field);
            }
    
            int total_field = 8 ;
            std::string tag_name = fields[0] + ":" + fields[1] + "-" + fields[2] ;
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
            
            if ( this_te.compare(prev_te) != 0 or ( this_te.compare(prev_te) == 0 and tag_name.compare(prev_tag_name) != 0 ) ){
                te_inter_peak_unique[prev_key]++ ; teFam_inter_peak_unique[prev_key_fam]++ ;
            }

            prev_te = this_te ; prev_peak = this_peak ; 
            prev_key = key ; prev_key_fam = key_fam ; 
            prev_tag_name = tag_name ;
        }
    }
    inter_bed.close() ; fake_list.close() ;
}

void parse_intersect_sortPeaks( std::string inter_bed_path, 
                                std::string bed_path, 
                                int idx_col, 
                                std::unordered_map<std::string, int> &peak_count_on_te,
                                std::unordered_map<std::string, int> &peak_count_on_te_unique,
                                std::unordered_map<std::string, int> &peak_inter_teFam_unique,
                                std::unordered_map<std::string, int> &peak_inter_te_unique,
                                std::string list_files  )
{
    std::cout << "Counting peaks vs TE overlaps - peak side ..." << std::endl ;
    std::ofstream fake_list ; 
    if ( bed_path.compare("empty") == 0 ){ fake_list.open(list_files, std::fstream::app) ; }
    std::ifstream inter_bed(inter_bed_path);
    if (this_is_empty(inter_bed)){
        std::cout << "Problem while opening " << inter_bed_path << ", exiting...\n" ;
        exit (EXIT_FAILURE) ;
    }
    if (inter_bed.is_open()) {
        std::string line, prev_te, prev_peak, prev_key , prev_key_fam, prev_tag_name ;
        while (getline(inter_bed, line)) { 
            std::istringstream iss(line) ;
            std::vector <std::string> fields ;
            std::string field ;
            while(std::getline(iss, field, '\t')){ 
                fields.push_back(field);
            }

            int total_field = 8 ;
            std::string tag_name = fields[0] + ":" + fields[1] + "-" + fields[2] ;
            if ( bed_path.compare("empty") == 0 ){
                if ( idx_col < 0 ){
                    tag_name = fields[0] + ":" + fields[1] + "-" + fields[2] ;
                } else {
                    tag_name = fields[idx_col] ;
                }
            } else {
                tag_name = fields[idx_col] ;
            }

            std::string this_peak = fields[total_field+1] + fields[total_field+2] + fields[total_field+3] ; 
            std::string key = fields[7] + "_" + tag_name ;
            std::string key_fam = fields[6] + "_" + tag_name ;

            peak_count_on_te[tag_name]++ ;

            if ( this_peak.compare(prev_peak) != 0 or ( this_peak.compare(prev_peak) == 0 and tag_name.compare(prev_tag_name) != 0 ) ){ // write prev_peak if this one if different OR same but different sample
                peak_count_on_te_unique[prev_tag_name]++ ;
                peak_inter_te_unique[prev_key]++ ;
                peak_inter_teFam_unique[prev_key_fam]++ ;
            }

            prev_peak = this_peak ; prev_tag_name = tag_name ;
            prev_key  = key       ; prev_key_fam  = key_fam  ;
        }
    }
    inter_bed.close() ; 
}

void add_pval_auto(std::string comparison_direction, int mean_peaklen, int avg_subfam_size, 
        std::vector<double> &all_pvals_best, 
        double pval, double pval_rev)
{
    if ( comparison_direction.compare("auto") == 0 ){
        // in auto mode, compute the smallest against the largest
        if ( mean_peaklen > avg_subfam_size ){
            all_pvals_best.push_back (pval_rev) ;
        }
        if ( mean_peaklen <= avg_subfam_size ){
            all_pvals_best.push_back (pval) ;
        }
    }
    else
    {
        if ( comparison_direction.compare("te_in_peak") == 0 ){
            all_pvals_best.push_back (pval_rev) ;
        }
        else if ( comparison_direction.compare("peak_in_te") == 0 ){
            all_pvals_best.push_back (pval) ;
        }
    }
}

// Calculate adjusted p-values with Benjamin-Hochsberg formula
std::unordered_map<std::string, double> calc_adj_pval ( std::vector<double> *all_pvals, std::vector<std::string> *subfam_names, bool padj ){

    std::unordered_map<std::string, double> all_final_pvals ;

    std::vector<int> order_pval((*all_pvals).size());
    std::size_t n(0);
    std::generate(std::begin(order_pval), std::end(order_pval), [&]{ return n++; });
    std::sort(  std::begin(order_pval), 
            std::end(order_pval),
            [&](int i1, int i2) { return (*all_pvals)[i1] < (*all_pvals)[i2]; } );


    if ( padj ){
        // Get adj pval for regular pval
        for (int i=0; i<(*all_pvals).size();i++){
            double this_pval = (*all_pvals)[order_pval[i]] ;
            double rank = double(i + 1) ;
            if ( this_pval > 0 ){
                // Benjamin-Hochsberg correction
                double pval_adj = this_pval * ( double((*all_pvals).size()) / rank ) ;
                if ( pval_adj > 1.0 ) pval_adj = 1.0 ;
                all_final_pvals.insert({ (*subfam_names)[order_pval[i]], pval_adj }) ;
            }
            else
            {
                all_final_pvals.insert({ (*subfam_names)[order_pval[i]], 0 }) ;
            }
        }
        return all_final_pvals ;
    }
    else
    {
        // Print regular pval in order 
        for (int i=0; i<(*all_pvals).size();i++){
            double this_pval = (*all_pvals)[order_pval[i]] ;
            if ( this_pval > 0 ){
                all_final_pvals.insert({ (*subfam_names)[order_pval[i]], this_pval }) ;
            }
            else
            {
                all_final_pvals.insert({ (*subfam_names)[order_pval[i]], 0 }) ;
            }
        }
        return all_final_pvals ;
    }
}

void print_all_pval_adj (std::string matrix_path, 
        std::string tag_name, bool *prhead_sum1, bool *first_it,
        std::vector<double> *pvals_ref,
        std::unordered_map<std::string, double> hyperGeom_reg_padj,
        std::unordered_map<std::string, double> binomial_padj,
        std::vector<std::string> *subfam_names, std::vector<int> *all_te_inter_peak_unique,
        std::vector<int> *all_total_subfam, std::vector<int> *all_te_inter_peak,
        int my_peak_count_on_te, int my_peak_count, std::vector<double> *all_total_bp_subfam,
        std::vector<int> *all_avg_subfam_size, std::vector<double> *all_subfam_genome_ratio,
        double peak_total_Mbp, double peak_ratio_on_te, double ratio_genome_peak, 
        std::string te_mode, bool *prhead_sum_te, std::string out_path,
        std::vector<double> *all_pbinom, std::vector<long int> *all_expect,
        std::vector<std::string> *all_dir_pbinom )
{

    // Writing results
    // open summary for this chipseq sample

    std::string summary_file = out_path + "/summary_bed/" + tag_name + "_" + te_mode + ".txt" ;
    std::ofstream summary ;
    summary.open (summary_file, std::fstream::app) ;
    if ( *prhead_sum1 ){ 
        summary << "name\tpadj.hGeom\tpadj.binom\tnTE\ttot.nTE\tp\tdir.auto\texpect\tavg.sfam.size\n" ;
        *prhead_sum1 = 0 ;
    }

    // Get pval order (to print ordered table) 
    std::vector<int> order_pval((*pvals_ref).size());
    std::size_t n(0);
    std::generate(std::begin(order_pval), std::end(order_pval), [&]{ return n++; });
    std::sort(  std::begin(order_pval), 
            std::end(order_pval),
            [&](int i1, int i2) { return (*pvals_ref)[i1] < (*pvals_ref)[i2]; } );

    // Print Matrix header  
    std::ofstream mat_out ; 
    mat_out.open(matrix_path, std::fstream::app) ;
    if ( *first_it ){
        mat_out << "names\t" << (*subfam_names)[0] ;
        for ( int i=1; i<(*subfam_names).size(); i++ ){
            mat_out << "\t" << (*subfam_names)[i] ;
        }
        mat_out << "\n" ;
        *first_it = 0 ;
    }
    mat_out << tag_name ;

    // Get adj pval for regular pval
    for (int i = 0 ; i < (*subfam_names).size() ; i++ )
    {
        std::string this_subfam = (*subfam_names)[order_pval[i]] ;
        double hyGm_reg = hyperGeom_reg_padj[this_subfam] ;
        //double hyGm_alaFish = hyperGeom_alaFish_padj[this_subfam] ;
        double binomial = binomial_padj[this_subfam] ;
        summary << (*subfam_names)[order_pval[i]] << "\t"  << hyGm_reg << "\t" << binomial << "\t" << (*all_te_inter_peak_unique)[order_pval[i]] << "\t" << (*all_total_subfam)[order_pval[i]] << "\t" << (*all_pbinom)[order_pval[i]] << "\t" << (*all_expect)[order_pval[i]] << "\t" << (*all_dir_pbinom)[order_pval[i]] << "\t" << (*all_avg_subfam_size)[order_pval[i]] << "\n" ;

        // Make matrix ( requires regular order )
        std::string this_subfam_2 = (*subfam_names)[i] ;
        double hyGm_reg_2 = hyperGeom_reg_padj[this_subfam_2] ;
        //double hyGm_alaFish_2 = hyperGeom_alaFish_padj[this_subfam_2] ;
        double binomial_2 = binomial_padj[this_subfam_2] ;

        if ( binomial_padj[this_subfam_2] == 1 )
            mat_out << "\t" << 0 ;
        else if ( binomial_padj[this_subfam_2] == 0 )
            mat_out << "\t" << 300 ;
        else if ( binomial_padj[this_subfam_2] > 0 ) 
            mat_out << "\t" << (-1)*log10(binomial_padj[this_subfam_2]) ;
        else
            mat_out << "\tNA" ; 

        // open summary for the TE subfam/fam and write results 
        std::replace( this_subfam_2.begin(), this_subfam_2.end(), '/', '_'); 
        std::string summary_subfam_file = out_path + "/summary_" + te_mode + "/" + this_subfam_2 + "_summary.txt" ;
        std::ofstream summary_subfam ;
        summary_subfam.open (summary_subfam_file, std::fstream::app) ;

        if ( *prhead_sum_te ){
            summary_subfam << "name\tpadj.hGeom\tpadj.binom\tnTE\ttot.nTE\tp\tdir.auto\texpect\tavg.sfam.size\n" ;
        }
        summary_subfam << (*subfam_names)[order_pval[i]] << "\t"  << hyGm_reg << "\t" << binomial << "\t" << (*all_te_inter_peak_unique)[order_pval[i]] << "\t" << (*all_total_subfam)[order_pval[i]] << "\t" << (*all_pbinom)[order_pval[i]] << "\t" << (*all_expect)[order_pval[i]] << "\t" << (*all_dir_pbinom)[order_pval[i]] << "\t" << (*all_avg_subfam_size)[order_pval[i]] << "\n" ;
        summary_subfam.close() ;
    }
    mat_out << "\n" ;                
    summary.close() ; mat_out.close() ;
}


double compute_hypergeom(   std::unordered_map<std::string, int> &peak_inter_te_unique,
                            std::string key,
                            int n_subfam,
                            int my_peak_count,
                            long int size_hg19,
                            int total_average, 
                            const int te_data_size,
                            std::string type_hypergeom,
                            std::string comparison_type )
{
    if ( type_hypergeom.compare("regular") == 0 ){
        long long a_11 = peak_inter_te_unique[key] ; //te_inter_peak[key] ;
        long long a_12 = MAX(0L,n_subfam - a_11) ;
        long long a_21 = MAX(0L,my_peak_count - a_11) ;
        long long a_22 = te_data_size - a_11 - a_12 - a_21 ;
        return fisher_exact(a_11, a_12, a_21, a_22, comparison_type) ;
    }
    else if ( type_hypergeom.compare("ala_bedtools_fisher") == 0 ){
        long long b_11 = peak_inter_te_unique[key] ; //te_inter_peak[key] ;
        long long a_12 = MAX(0L,n_subfam - b_11) ;
        long long b_12 = MAX(a_12, n_subfam - b_11 ) ;
        long long b_21 = MAX(0L,my_peak_count - b_11) ;
        long long b_22 = MAX( int(double(size_hg19) / double(total_average)) - b_11 - b_21 - b_12, 1 )  ;
        return fisher_exact(b_11, b_12, b_21, b_22, comparison_type) ;
    }
    else
    {
        std::cout << "ERROR : variable 'type_hypergeom' has non recognized value (should be either 'regular' or 'ala_bedtools_fisher'\n" ;
        exit(EXIT_FAILURE) ;
    }
}

double compute_binom( int n_tot, 
                      int X_obs, 
                      double p_exp, 
                      std::string comparison_type )
{
    return pbinom(X_obs-1, n_tot, p_exp, comparison_type) ;
}


void enrichment_analysis(   std::string type_analysis,
                            std::string ref_file, long int size_hg19, 
                            std::string comparison_type, std::string tag_name,
                            int my_peak_count, int my_peak_count_on_te_unique,
                            int peak_tot_bp, int mean_peaklen, 
                            std::string comparison_direction, const int te_data_size, 
                            int total_nonTE, int total_nonTE_bp,
                            double nonTE_avg_size, double nonTE_genome_ratio,
                            std::string matrix_path_all_best, std::string out_path,
                            bool print_padj, bool &prhead_mat_all_best, 
                            bool &prhead_mat_all_best_fam, bool &prhead_summary_te,
                            std::unordered_map<std::string, int> &peak_inter_te_unique, 
                            std::unordered_map<std::string, int> &te_inter_peak_unique,
                            std::string type_hypergeom )
{ 
    bool prhead_sum_all_best = 1, prhead_sum_all_best_fam = 1 ;
    double peak_total_Mbp = double(peak_tot_bp) / 1000000 ;
    double peak_ratio_on_te = double(my_peak_count_on_te_unique) / double(my_peak_count) ; 
    double ratio_genome_peak = double(peak_tot_bp) / double(size_hg19) ;

    std::ifstream ref_in(ref_file);
    if (this_is_empty(ref_in)){
        std::cout << "Problem while opening " << ref_file << ", exiting...\n" ;
        exit (EXIT_FAILURE) ;
    }
    if (ref_in.is_open()) {

        // initialize vectors
        std::vector<std::string> subfam_names, all_dir_pbinom;
        std::vector<double> all_pvals , all_pvals_binomial, all_pvals_rev , all_pvals_binomial_rev , all_pvals_best , all_pvals_binomial_best, all_total_bp_subfam, all_subfam_genome_ratio, all_pbinom;
        std::vector<int> all_te_inter_peak, all_te_inter_peak_unique , all_total_subfam , all_avg_subfam_size ;
        std::vector<long int> all_expect ;

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

            all_te_inter_peak.push_back(peak_inter_te_unique[key]) ;
            all_te_inter_peak_unique.push_back(te_inter_peak_unique[key]) ;
            all_total_subfam.push_back(n_subfam);
            all_total_bp_subfam.push_back(total_Mbp_len_subfam);
            all_avg_subfam_size.push_back(avg_subfam_size);
            all_subfam_genome_ratio.push_back(ratio_genome_subfam);

            // Calculation of the p-values for each 1-1 relation
            if (te_inter_peak_unique.find(key) == te_inter_peak_unique.end() || te_inter_peak_unique[key] == 0){
                subfam_names.push_back (fields[0]) ; all_pvals.push_back (1.0) ;
                all_pvals_binomial.push_back (1.0) ; all_pvals_rev.push_back (1.0) ; 
                all_pvals_binomial_rev.push_back (1.0) ; all_pvals_best.push_back (1.0) ;
                all_pvals_binomial_best.push_back (1.0) ; all_pbinom.push_back(0.0) ;
                all_expect.push_back(0) ; all_dir_pbinom.push_back("NA") ;
            }
            else
            {
                ////////////////////////////////////// 
                // PVAL ENRICHMENT PEAKs AMONG TEs  //
                //////////////////////////////////////

                // HyperGeometric test
                int total_average = mean_peaklen + avg_subfam_size ;
                double pval = compute_hypergeom(peak_inter_te_unique, key, n_subfam,
                                                my_peak_count_on_te_unique,
                                                size_hg19, total_average, te_data_size,
                                                type_hypergeom, comparison_type ) ;
                all_pvals.push_back (pval) ;
                //std::cout << "pval : " << pval << std::endl ;

                // Binomial exact test
                double p_exp = double(total_bp_length_subfam) / double(size_hg19) ; // Prob to touch subfam by random
                double pval_binomial = compute_binom( my_peak_count, peak_inter_te_unique[key], p_exp, comparison_type) ;
                all_pvals_binomial.push_back(pval_binomial) ;
                //std::cout << "pval_binomial : " << pval_binomial << std::endl ;

                //////////////////////////////////////
                // PVAL ENRICHMENT TEs AMONG PEAKs  //
                //////////////////////////////////////

                // Regular HyperGeometric
                int total_average_rev = avg_subfam_size ;
                double pval_rev = compute_hypergeom(te_inter_peak_unique, key, n_subfam,
                                                my_peak_count_on_te_unique,
                                                size_hg19, total_average_rev, te_data_size,
                                                type_hypergeom, comparison_type ) ;
                all_pvals_rev.push_back (pval_rev) ;
                //std::cout << "pval_rev : " << pval_rev << std::endl ;

                // add required pval for 'auto' mode 
                add_pval_auto(  comparison_direction, mean_peaklen, avg_subfam_size, 
                        all_pvals_best, pval, pval_rev) ;

                // Binomial exact test with genome ratio as p 
                double p_exp_rev = ratio_genome_peak ; // Prob to touch subfam by random
                double pval_binomial_rev = compute_binom( n_subfam, te_inter_peak_unique[key], p_exp, comparison_type) ;
                all_pvals_binomial_rev.push_back(pval_binomial_rev) ;
                //std::cout << "pval_binomial_rev : " << pval_binomial_rev << std::endl ;

                // add required pval for 'auto' mode
                add_pval_auto(  comparison_direction, mean_peaklen, avg_subfam_size, 
                                all_pvals_binomial_best, pval_binomial, pval_binomial_rev) ;

                //if ( comparison_direction.compare("auto") == 0 ){
                    if ( mean_peaklen > avg_subfam_size ){
                        all_pbinom.push_back (p_exp_rev) ;
                        long int expect = (long int)(double(p_exp_rev)*double(n_subfam)) ;
                        //std::cout << "expect : " << expect << std::endl ;

                        all_expect.push_back (expect) ;
                        std::string dir_pbinom = "te_in_peak" ;
                        all_dir_pbinom.push_back(dir_pbinom) ;
                    }
                    else if ( mean_peaklen <= avg_subfam_size ){
                        all_pbinom.push_back (p_exp) ;
                        long int expect = (long int)(double(p_exp)*double(my_peak_count)) ;
                        //std::cout << "expect : " << expect << std::endl ;

                        all_expect.push_back (expect) ;
                        std::string dir_pbinom = "peak_in_te" ;
                        all_dir_pbinom.push_back(dir_pbinom) ;
                    }
                //}
 
                subfam_names.push_back (fields[0]) ;
            }
        }

        //return  Getting pval for nonTE
        subfam_names.push_back ("nonTE") ;
        int nonTE_intersect = my_peak_count - my_peak_count_on_te_unique ;
        // reg hypergeometric test
        int nonTE_a12 = total_nonTE - nonTE_intersect ;
        int nonTE_a21 = my_peak_count - nonTE_intersect ;
        int nonTE_a22 = 2*4570939 - nonTE_intersect - nonTE_a12 - nonTE_a21 ;
        double nonTE_pval_reg = fisher_exact(nonTE_intersect, nonTE_a12, nonTE_a21, nonTE_a22, comparison_type) ;
        //all_pvals_best.push_back (nonTE_pval_reg) ;

        // hypergeometric alafisher
        int peak_nonTE_total_average = 276 + mean_peaklen ;
        int nonTE_b22 = MAX( int(double(size_hg19) / double(peak_nonTE_total_average)) - nonTE_intersect - nonTE_a12 - nonTE_a21 , 1)  ;
        double nonTE_pval_alafisher = fisher_exact(nonTE_intersect, nonTE_a12, nonTE_a21, nonTE_b22, comparison_type) ;
        all_pvals_best.push_back (nonTE_pval_alafisher) ;

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
        all_pbinom.push_back(p_exp_nonTE) ;
        int expect_nonTE = int(double(p_exp_nonTE)*double(X_obs_nonTE)) ;
        all_expect.push_back(expect_nonTE) ;
        all_dir_pbinom.push_back("NA") ;

        // Get ajusted p-values with Benjamin-Hochsberg method 
        std::unordered_map<std::string, double> padj_hygm_reg = calc_adj_pval(&all_pvals_best, &subfam_names, print_padj) ;
        //std::unordered_map<std::string, double> padj_hygm_alaFish = calc_adj_pval(&all_pvals_alafisher_best, &subfam_names, print_padj) ;
        std::unordered_map<std::string, double> padj_hygm_binom = calc_adj_pval(&all_pvals_binomial_best, &subfam_names, print_padj) ;

        // Printing final results with 3 adjusted p-values
        print_all_pval_adj(matrix_path_all_best, tag_name,
                &prhead_sum_all_best, &prhead_mat_all_best, &all_pvals_binomial_best,
                padj_hygm_reg, padj_hygm_binom,
                &subfam_names, &all_te_inter_peak_unique,
                &all_total_subfam,  &all_te_inter_peak,
                my_peak_count_on_te_unique, my_peak_count, &all_total_bp_subfam,
                &all_avg_subfam_size, &all_subfam_genome_ratio,
                peak_total_Mbp, peak_ratio_on_te, ratio_genome_peak,
                type_analysis, &prhead_summary_te, out_path,
                &all_pbinom, &all_expect, &all_dir_pbinom ) ;
    }
    ref_in.close() ;
}

void print_help(int X, int Y, std::string detail){
    std::cout << std::endl ;
    std::cout << "/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\__/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\__________________________________________________/\\\\\\" << std::endl ;
    std::cout << "\\///////\\\\\\/////__\\/\\\\\\///////////__________________________________________________\\/\\\\\\     " << "V" << version << std::endl ;
    std::cout << "       \\/\\\\\\_______\\/\\\\\\__________________________________________/\\\\\\_______________\\/\\\\\\" << std::endl ;
    std::cout << "        \\/\\\\\\_______\\/\\\\\\\\\\\\\\\\\\\\\\______/\\\\/\\\\\\\\\\\\____/\\\\/\\\\\\\\\\\\\\__\\///______/\\\\\\\\\\\\\\\\_\\/\\\\\\" << std::endl ;
    std::cout << "         \\/\\\\\\_______\\/\\\\\\///////______\\/\\\\\\////\\\\\\__\\/\\\\\\/////\\\\\\__/\\\\\\___/\\\\\\//////__\\/\\\\\\\\\\\\\\\\\\\\" << std::endl ;
    std::cout << "          \\/\\\\\\_______\\/\\\\\\_____________\\/\\\\\\__\\//\\\\\\_\\/\\\\\\___\\///__\\/\\\\\\__/\\\\\\_________\\/\\\\\\/////\\\\\\" << std::endl ; 
    std::cout << "           \\/\\\\\\_______\\/\\\\\\_____________\\/\\\\\\___\\/\\\\\\_\\/\\\\\\_________\\/\\\\\\_\\//\\\\\\________\\/\\\\\\___\\/\\\\\\" << std::endl ;
    std::cout << "            \\/\\\\\\_______\\/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_\\/\\\\\\___\\/\\\\\\_\\/\\\\\\_________\\/\\\\\\__\\///\\\\\\\\\\\\\\\\_\\/\\\\\\___\\/\\\\\\" << std::endl ;
    std::cout << "             \\///________\\///////////////__\\///____\\///__\\///__________\\///_____\\////////__\\///____\\///" << std::endl ;
    std::cout <<  std::endl << std::endl ;
    std::cout << "TEnrich for multiple beds :                      " << std::endl ;
    std::cout << "\n\tTEnrich --bed_dir path/to/dirWithBeds --out_dir path/to/dirOut\n\n" ;
    std::cout << "TEnrich for one single bed :                      " << std::endl ;
    std::cout << "\n\tTEnrich --single_file path/to/file.bed --idx_col 4 --out_dir path/to/dirOut\n\n" ;
    std::cout << "for detailed help of each paramters :                      " << std::endl ;
    std::cout << "\n\tTEnrich -hd\n\n" ;
    if ( detail.compare("detail") == 0 ){
        print_help_line("--bed_dir path/to/dirWithBeds [string]","every file with *.bed ext in the folder will be used",X,Y) ;
        print_help_line("--out_dir path/to/dirOut [string]","The folder is created and results written inside (WARNING: everything is cleaned before a new run)",X,Y) ;
        print_help_line("--single_file path/to/bed_file [string]", "if this is given, it will use a single file instead of a group of bed file to do the enrichment (cancels --bed_dir option). If no index column is given, will perform the enrichment analysis on each single lines. Otherwise, it will group lines per name of the feature in the column designed by --idx_col option.", X, Y) ;
        print_help_line("--idx_col [integer]","designate the index of the column (WARNING: 0-based index, which means index 0 is the first column, index 1 is the second, etc...) to be used in file given in --single_file to group lines.",X,Y) ;
        print_help_line("--type_hypergeom ['ala_bedtools_fisher','regular']","type of hypergeometric test to do. By default, will use a similar method that bedtools fisher uses (description: https://bedtools.readthedocs.io/en/latest/content/tools/fisher.html). With regular option, uses a simple hypergeometric without weighting for TE loci length",X,Y) ; 
        print_help_line("--comp_sense ['te_in_peak','peak_in_te','auto']","Defines the direction for the comparison, 'te_in_peak' or 'peak_in_te'. In 'auto' mode, it will take the enrichment of the smaller to the bigger [OPTIONAL]. Default value : 'auto'",X,Y) ;
        print_help_line("--stat_test_type ['greater','less']","For statistical test done (Hypergeometric and Binomial), tell if we want the right tail ('greater') or the left tail ('less'), in other words the probability of having a equal or greater / equal or lower number of hits in the intersect. [OPTIONAL]. Default value : 'greater'",X,Y) ;
        print_help_line("--padj ['true','false']","tell if you want to print the adjusted p-val (with the Benjamin-Hochsberg correction). [OPTIONAL]. Default value : 'true'",X,Y) ;
        print_help_line("--genome_size [integer/double]","Genome size over which the enrichment calculation will be done. Expect only digits. [OPTIONAL]. Default value : 3088269832",X,Y) ;
        print_help_line("--ref_subfam [STRING]","subfam ref file obtained with utils/make_ref_file.pl. If not specified, the one for hg19 in db/ folder will be used [OPTIONAL]. Default value : 'db/Subfam_ref_TE.txt'",X,Y) ;
        print_help_line("--ref_fam [STRING]","fam ref file obtained with utils/make_ref_file.pl. If not specified, the one for hg19 in db/ folder will be used [OPTIONAL]. Default value : 'db/Fam_ref_TE.txt'",X,Y) ;
        print_help_line("--te_database [STRING]","database of TE used to make the intersection. Should be in bed format, as returned by utils/convert_repeatmasker.sh. By default, uses hg19 with LTR merged by J.Duc. [OPTIONAL]. Default value : 'db/hg19_TE_repmask_LTRm_s_20140131.bed",X,Y) ; 
    }

    exit (EXIT_FAILURE) ;
}

bool is_digits(const std::string &str)
{
        return std::all_of(str.begin(), str.end(), ::isdigit); // C++11
}

bool exists (const std::string& name1) {
    struct stat buffer;   
    return (stat (name1.c_str(), &buffer) == 0 ) ; 
}

void check_folder (const std::string& dir) {
    if ( exists(dir) ){
        //std::cout << "Directory exists, beginning process..." << std::endl ;
    }
    else
    {
        std::cerr << "Directory doesn't exists. Exiting..." << std::endl ;
        exit (EXIT_FAILURE) ;
    }
}

std::string remove_ext(std::string my_str){
    size_t lastindex = my_str.find_first_of("."); 
    std::string rawname = my_str.substr(0, lastindex); 
    return rawname ;
}

std::string base_name(std::string path)
{
    path.erase(std::remove(path.begin(), path.end(), '.'), path.end());
    std::string newstring = path.substr(path.find_last_of("/\\") + 1) ;
    return newstring;
}

bool this_is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

void create_folder(std::string out_path, std::string name_folder){
    std::stringstream create_folder_cmd ;
    create_folder_cmd << "mkdir -p "<< out_path << "/" << name_folder ;
    system(&(create_folder_cmd.str()[0])) ;
}
