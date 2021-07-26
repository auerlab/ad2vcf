/* ad2vcf.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int ad2vcf(const char *argv[], FILE *sam_stream);
int skip_upstream_alignments(bl_vcf_t *vcf_call, FILE *sam_stream, bl_sam_buff_t *sam_buff, FILE *vcf_out_stream, vcf_stats_t *vcf_stats);
int allelic_depth(bl_vcf_t *vcf_call, FILE *sam_stream, bl_sam_buff_t *sam_buff, FILE *vcf_out_stream, vcf_stats_t *vcf_stats);
void vcf_stats_update_allele_count(vcf_stats_t *vcf_stats, bl_vcf_t *vcf_call, bl_sam_t *sam_alignment, FILE *vcf_out_stream);
int uchar_cmp(unsigned char *c1, unsigned char *c2);
void vcf_stats_init(vcf_stats_t *vcf_stats);
