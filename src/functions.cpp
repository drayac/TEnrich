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
#include "functions.hpp"
#include <algorithm>

#define constant_e (2.71828)
#define PI (3.14159265359)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

std::string this_dir = <FOLDER_INSTALL> ;
std::string version  = <VERSION> ;

std::vector <std::string> get_parameters(int c_argc, char *c_argv[]){
    if (c_argc < 1) {
        std::cerr << "Not enough parameters were given, exiting ... " << std::endl;
        exit (EXIT_FAILURE) ;
    }
    std::vector <std::string> named_params ;
    std::string bed_dir = "empty", out_dir = "empty", size_genome = "3088269832" ;
    std::string padj = "true" , comp_sense = "auto", stat_test_type = "greater" ;
    std::string ref_subfam = this_dir + "/db/Subfam_ref_TE.txt" ;
    std::string ref_fam = this_dir + "/db/Fam_ref_TE.txt" ;
    std::string ref_clust = this_dir + "/db/Clusters_ref_TE.txt" ;
    std::string te_database = this_dir + "/db/hg19_TE_repmask_LTRm_s_20140131.bed.gz" ;
    std::string single_file = "empty" , idx_col = "-1" ;
    std::string destination ; 
    int i ;
    for (i = 1; i < c_argc; ++i) {
        std::string this_param = std::string(c_argv[i]) ;
        if (std::string(c_argv[i]) == "--bed_dir") {
            if (i + 1 < c_argc) { // Make sure we aren't at the end of argv!
                 bed_dir = c_argv[i + 1]; i++ ; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { 
                std::cerr << "--bed_dir option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--out_dir") {
            if (i + 1 < c_argc) { 
                out_dir = c_argv[i + 1]; i++ ;
            } else { 
                std::cerr << "--out_dir option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
       if (std::string(c_argv[i]) == "--genome_size") {
            if (i + 1 < c_argc) { 
                size_genome = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--genome_size option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--comp_sense") {
            if (i + 1 < c_argc) { 
                comp_sense = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--comp_sense option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--stat_test_type") {
            if (i + 1 < c_argc) { 
                stat_test_type = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--stat_test_type option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        } 
        if (std::string(c_argv[i]) == "--padj") {
            if (i + 1 < c_argc) { 
                padj = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--padj option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--ref_subfam") {
            if (i + 1 < c_argc) { 
                ref_subfam = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--ref_subfam option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--ref_fam") {
            if (i + 1 < c_argc) { 
                ref_fam = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--ref_fam option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--ref_cluster") {
            if (i + 1 < c_argc) { 
                ref_clust = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--ref_cluster option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--te_database") {
            if (i + 1 < c_argc) { 
                te_database = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--padj option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--single_file") {
            if (i + 1 < c_argc) { 
                single_file = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--single_file option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
        if (std::string(c_argv[i]) == "--idx_col") {
            if (i + 1 < c_argc) { 
                idx_col = c_argv[i + 1]; i++ ; 
            } else { 
                std::cerr << "--idx_col option requires one argument." << std::endl;
                exit (EXIT_FAILURE) ;
            }
        }
    }
    if ( out_dir.compare("empty") == 0 ){
        std::cerr << "bed and out folder should be precised ... exiting\n" ;
        exit (EXIT_FAILURE) ;
    }
    else
    {
        named_params.push_back(bed_dir) ; named_params.push_back(out_dir) ;
    }

    if ( is_digits(size_genome) ){
        named_params.push_back(size_genome) ;
    }
    else
    {
        std::cerr << "size genome parameters should only contain digits !! exiting\n" ;
        exit (EXIT_FAILURE) ;
    }

    if ( comp_sense.compare("peak_in_te") != 0 && comp_sense.compare("te_in_peak") != 0 && comp_sense.compare("auto") != 0 ){
        std::cerr << "option comp_sense should be either 'peak_in_te', 'te_in_peak' or 'auto'. By default 'auto'. Exiting...\n\n" ; 
        exit ( EXIT_FAILURE ) ;
    }
    named_params.push_back(comp_sense) ;
    
    if ( stat_test_type.compare("greater") != 0 && stat_test_type.compare("less") != 0){
        std::cerr << "option stat_test_type should be either 'greater', 'less'. By default 'greater'. Exiting...\n\n" ; 
        exit ( EXIT_FAILURE ) ;
    }
    named_params.push_back(stat_test_type) ;
   
    if ( padj.compare("true") != 0 && padj.compare("false") != 0){
        std::cerr << "option padj should be either 'true', 'false'. By default 'true'. Exiting...\n\n" ; 
        exit ( EXIT_FAILURE ) ;
    }

    named_params.push_back(padj) ;
    named_params.push_back(ref_subfam) ;
    named_params.push_back(ref_fam) ;
    named_params.push_back(ref_clust) ;
    named_params.push_back(te_database) ;
    named_params.push_back(single_file) ;
    named_params.push_back(idx_col) ;
    
    return named_params ;
}

void get_help(std::string firstparam, int argc_, std::string version ){
    if ( argc_ > 1 ){
        if ( firstparam == "-h" || firstparam == "--help" ){
            print_help(50,50) ;
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

void print_help(int X, int Y){

	std::cout << " __________         _     __ " << std::endl;
	std::cout << "/_  __/ __/__  ____(_)___/ / " << std::endl;
	std::cout << " / / / _// _ \\/ __/ / __/ _ \\" << std::endl;
	std::cout << "/_/ /___/_//_/_/ /_/\\__/_//_/" << std::endl;
	std::cout << "                        v." << version << std::endl << std::endl ;
	std::cout << "Required options :                      " << std::endl ;
	std::cout << "./TEnrich --bed_dir path/to/dirWithBeds \\\n" ;
	std::cout << "          --out_dir path/to/dirOut \\\n" ;
    std::cout << "\n" ;
    print_help_line("--bed_dir path/to/dirWithBeds [string]","every file with *.bed ext in the folder will be used",X,Y) ;
    print_help_line("--single_file path/to/bed_file [string]", "if this is given, it will use a single file instead of a group of bed file to do the enrichment (cancels --bed_dir option). If no index column is given, will perform the enrichment analysis on each single lines. Otherwise, it will group lines per name of the feature in the column designed by --idx_col option.", X, Y) ;
    print_help_line("--idx_col [integer]","designate the index of the column (WARNING: 0-based index, which means index 0 is the first column, index 1 is the second, etc...) to be used in file given in --single_file to group lines.",X,Y) ;
    print_help_line("--out_dir path/to/dirOut [string]","The folder is created and results written inside (WARNING: everything is cleaned before a new run)",X,Y) ;
    print_help_line("--genome_size [integer/double]","Genome size over which the enrichment calculation will be done. Expect only digits. [OPTIONAL]. Default value : 3088269832",X,Y) ;
    print_help_line("--comp_sense ['te_in_peak','peak_in_te','auto']","Defines the direction for the comparison, 'te_in_peak' or 'peak_in_te'. In 'auto' mode, it will take the enrichment of the smaller to the bigger [OPTIONAL]. Default value : 'auto'",X,Y) ;
    print_help_line("--stat_test_type ['greater','less']","For statistical test done (Hypergeometric and Binomial), tell if we want the right tail ('greater') or the left tail ('less'), in other words the probability of having a equal or greater / equal or lower number of hits in the intersect. [OPTIONAL]. Default value : 'greater'",X,Y) ;
    print_help_line("--padj ['true','false']","tell if you want to print the adjusted p-val (with the Benjamin-Hochsberg correction). [OPTIONAL]. Default value : 'true'",X,Y) ;
    print_help_line("--ref_subfam [STRING]","subfam ref file obtained with utils/make_ref_file.pl. If not specified, the one for hg19 in db/ folder will be used [OPTIONAL]. Default value : 'db/Subfam_ref_TE.txt'",X,Y) ;
    print_help_line("--ref_fam [STRING]","fam ref file obtained with utils/make_ref_file.pl. If not specified, the one for hg19 in db/ folder will be used [OPTIONAL]. Default value : 'db/Fam_ref_TE.txt'",X,Y) ;
    print_help_line("--te_database [STRING]","database of TE used to make the intersection. Should be in bed format, as returned by utils/convert_repeatmasker.sh. By default, uses hg19 with LTR merged by J.Duc. [OPTIONAL]. Default value : 'db/hg19_TE_repmask_LTRm_s_20140131.bed",X,Y) ; 

	exit (EXIT_FAILURE) ;
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
        std::cout << "Directory exists, beginning process..." << std::endl ;
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
    std::cout << "before: " << path << std::endl ;
    path.erase(std::remove(path.begin(), path.end(), '.'), path.end());
    std::string newstring = path.substr(path.find_last_of("/\\") + 1) ;
    std::cout << "after: " << newstring << std::endl ;
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
        std::unordered_map<std::string, double> hyperGeom_alaFish_padj,
        std::unordered_map<std::string, double> binomial_padj,
        std::vector<std::string> *subfam_names, std::vector<int> *all_te_inter_peak_unique,
       std::vector<int> *all_total_subfam, std::vector<int> *all_te_inter_peak,
       int my_peak_count_on_te, int my_peak_count, std::vector<double> *all_total_bp_subfam,
        std::vector<int> *all_avg_subfam_size, std::vector<double> *all_subfam_genome_ratio,
       double peak_total_Mbp, double peak_ratio_on_te, double ratio_genome_peak, 
       std::string te_mode, bool *prhead_sum_te, std::string out_path ){
    
    // Writing results
    // open summary for this chipseq sample
    
    std::string summary_file = out_path + "/summary_bed/" + tag_name + "_" + te_mode + ".txt" ;
    std::ofstream summary ;
    summary.open (summary_file, std::fstream::app) ;
    if ( *prhead_sum1 ){ 
        summary << "subfam_name\tpadj.hypergeom.reg\tpadj.hypergeom.alafisher\tpadj.binomial\tte.count.with.peak\tte.total.n\tpeak.count.on.te\ttotal.peak.ov.te\ttotal.peak.count\ttotal.Mbp.te\tavg.subfam.size\tte.genome.ratio\tpeak.total.Mbp\tpeak.ratio.on.te\tratio.genome.peak\n" ;
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
        double hyGm_alaFish = hyperGeom_alaFish_padj[this_subfam] ;
        double binomial = binomial_padj[this_subfam] ;
        summary << (*subfam_names)[order_pval[i]] << "\t"  << hyGm_reg << "\t" << hyGm_alaFish << "\t" << binomial << "\t" << (*all_te_inter_peak_unique)[order_pval[i]] << "\t" << (*all_total_subfam)[order_pval[i]] << "\t" << (*all_te_inter_peak)[order_pval[i]] << "\t" << my_peak_count_on_te << "\t" << my_peak_count  << "\t" << (*all_total_bp_subfam)[order_pval[i]] << "\t" << (*all_avg_subfam_size)[order_pval[i]] << "\t" << (*all_subfam_genome_ratio)[order_pval[i]] << "\t" << peak_total_Mbp << "\t" << peak_ratio_on_te << "\t" << ratio_genome_peak << "\n" ;

        // Make matrix ( requires regular order )
        std::string this_subfam_2 = (*subfam_names)[i] ;
        double hyGm_reg_2 = hyperGeom_reg_padj[this_subfam_2] ;
        double hyGm_alaFish_2 = hyperGeom_alaFish_padj[this_subfam_2] ;
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
            summary_subfam << "sample_name\tpadj.hypergeom.reg\tpadj.hypergeom.alafisher\tpadj.binomial\tte.count.with.peak\tte.total.n\tpeak.count.on.te\ttotal.peak.ov.te\ttotal.peak.count\ttotal.Mbp.te\tavg.subfam.size\tte.genome.ratio\tpeak.total.Mbp\tpeak.ratio.on.te\tratio.genome.peak\n" ;
        }
        summary_subfam << tag_name << "\t"  << hyGm_reg_2 << "\t" << hyGm_alaFish_2 << "\t" << binomial_2 << "\t" << (*all_te_inter_peak_unique)[i] << "\t" << (*all_total_subfam)[i] << "\t" << (*all_te_inter_peak)[i] << "\t" << my_peak_count_on_te << "\t" << my_peak_count  << "\t" << (*all_total_bp_subfam)[i] << "\t" << (*all_avg_subfam_size)[i] << "\t" << (*all_subfam_genome_ratio)[i] << "\t" << peak_total_Mbp << "\t" << peak_ratio_on_te << "\t" << ratio_genome_peak << "\n" ;
        summary_subfam.close() ;
    }
    mat_out << "\n" ;                
    summary.close() ; mat_out.close() ;
}

