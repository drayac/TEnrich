#define constant_e (2.71828)
#define PI (3.14159265359)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void add_param(std::string &var, std::string name, int &i, int c_argc, char *c_argv[]) ;

std::vector <std::string> get_parameters(int , char**) ;

void get_help(std::string firstparam, int argc_, std::string version ) ;

void check_te_database( std::string te_data, int n_expect_fields) ;

bool exists (const std::string&) ;

void check_folder (const std::string&) ;

void print_help(int, int, std::string) ;

std::string remove_ext(std::string) ;

std::string base_name(std::string) ;

bool is_digits(const std::string &str) ; 

bool this_is_empty(std::ifstream&) ;

void create_folder(std::string out_path, std::string name_folder) ;

void concat_bed_files(std::string bed_path, std::string list_files, std::string concat_bed, std::unordered_map<std::string, int> &peak_count, std::unordered_map<std::string, int> &peak_len, std::unordered_map<std::string, int> &peak_total_bp, std::unordered_map<std::string, std::string> &summary_peak_line, long int size_hg19 ) ;

void bedtools_intersect(std::string bedtools_options, std::string sort_options, std::string te_data, std::string inter_bed_path, std::string concat_bed) ;

void sort_peaks(std::string sort_options, 
                std::string inter_bed_path, 
                std::string sort_peaks_bed) ;

void parse_intersect(std::string inter_bed_path, std::string bed_path, 
                        std::unordered_map<std::string,int> &save_fake_list, 
                        int idx_col,
                        std::unordered_map<std::string, int> &te_inter_peak, 
                        std::unordered_map<std::string,int> &teFam_inter_peak, 
                        std::unordered_map<std::string, int> &te_inter_peak_unique, 
                        std::unordered_map<std::string, int> &teFam_inter_peak_unique,
                        std::unordered_map<std::string, int> &peak_count_on_te) ;

void parse_intersect_sortPeaks( std::string inter_bed_path, 
                                std::string bed_path, 
                                int idx_col, 
                                std::unordered_map<std::string, int> &peak_count_on_te,
                                std::unordered_map<std::string, int> &peak_count_on_te_unique,
                                std::unordered_map<std::string, int> &peak_inter_teFam_unique,
                                std::unordered_map<std::string, int> &peak_inter_te_unique ,
                                std::string list_files ) ;

double compute_hypergeom(   std::unordered_map<std::string, int> &peak_inter_te_unique,
                            std::string key,
                            int n_subfam,
                            int my_peak_count,
                            long int size_hg19,
                            int total_average, 
                            const int te_data_size,
                            std::string type_hypergeom,
                            std::string comparison_type ) ;
double compute_binom( int n_tot, 
                      int X_obs, 
                      double p_exp, 
                      std::string comparison_type ) ;

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
                            std::string type_hypergeom ) ;

void add_pval_auto(std::string comparison_direction, int mean_peaklen, int avg_subfam_size, 
              std::vector<double> &all_pvals_best, 
              double pval, double pval_rev) ;

std::unordered_map<std::string, double> calc_adj_pval ( std::vector<double> *all_pvals, std::vector<std::string> *subfam_names, bool ) ;

void print_all_pval_adj (
        std::string matrix_path,
        std::string tag_name, 
        bool *prhead_sum1, 
        bool *first_it,
        std::vector<double> *pvals_ref,
        std::unordered_map<std::string, double> hyperGeom_reg_padj,
        std::unordered_map<std::string, double> binomial_padj,
        std::vector<std::string> *subfam_names, 
        std::vector<int> *all_te_inter_peak_unique,
        std::vector<int> *all_total_subfam, 
        std::vector<int> *all_te_inter_peak,
        int my_peak_count_on_te, 
        int my_peak_count, 
        std::vector<double> *all_total_bp_subfam,
        std::vector<int> *all_avg_subfam_size, 
        std::vector<double> *all_subfam_genome_ratio,
        double peak_total_Mbp, 
        double peak_ratio_on_te, 
        double ratio_genome_peak,
        std::string te_mode, 
        bool *prhead_sum_te, 
        std::string out_path,
        std::vector<double> *all_pbinom, 
        std::vector<long int> *all_expect,
        std::vector<std::string> *all_dir_pbinom  ) ;

