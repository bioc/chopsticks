index_db create_name_index(const SEXP names);

void do_impute(const unsigned char *snps, const int nrow, 
	       const int *subset, int nuse, index_db snp_names, SEXP Rule, 
	       double *value_a, double *value_d);
