#define constant_e (2.71828)
#define PI (3.14159265359)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

std::vector <std::string> get_parameters(int , char**) ;

void get_help(std::string firstparam, int argc_, std::string version ) ;

void check_te_database( std::string te_data, int n_expect_fields) ;

bool exists (const std::string&) ;

void check_folder (const std::string&) ;

void print_help(int, int) ;

std::string remove_ext(std::string) ;

std::string base_name(std::string) ;

bool is_digits(const std::string &str) ; 

bool this_is_empty(std::ifstream&) ;

void create_folder(std::string out_path, std::string name_folder) ;

std::unordered_map<std::string, double> calc_adj_pval ( std::vector<double> *all_pvals, std::vector<std::string> *subfam_names, bool ) ;

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
        std::string te_mode, bool *prhead_sum_te, std::string out_path ) ;

