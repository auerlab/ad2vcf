/* ad2vcf.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int ad2vcf(const char *argv[], FILE *sam_stream);
int skip_upstream_alignments(bl_vcf_t *vcf_call, FILE *sam_stream, bl_sam_buff_t *sam_buff, FILE *vcf_out_stream, ad2vcf_stats_t *stats);
int allelic_depth(bl_vcf_t *vcf_call, FILE *sam_stream, bl_sam_buff_t *sam_buff, FILE *vcf_out_stream, ad2vcf_stats_t *stats);
void update_allele_count(bl_vcf_t *vcf_call, bl_sam_t *sam_alignment, FILE *vcf_out_stream, ad2vcf_stats_t *stats);
int uchar_cmp(unsigned char *c1, unsigned char *c2);
void ad2vcf_stats_init(ad2vcf_stats_t *stats);
