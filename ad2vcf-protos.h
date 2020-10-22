/* ad2vcf.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int ad2vcf(const char *argv[], FILE *sam_stream);
_Bool skip_upstream_alignments(vcf_call_t *vcf_call, FILE *sam_stream, sam_buff_t *sam_buff, FILE *vcf_out_stream, ad2vcf_stats_t *stats);
_Bool allelic_depth(vcf_call_t *vcf_call, FILE *sam_stream, sam_buff_t *sam_buff, FILE *vcf_out_stream, ad2vcf_stats_t *stats);
_Bool call_in_alignment(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment);
_Bool alignment_upstream_of_call(vcf_call_t *vcf_call, sam_alignment_t *alignment);
void update_allele_count(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment, FILE *vcf_out_stream, ad2vcf_stats_t *stats);
void sam_buff_check_order(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment);
void sam_buff_init(sam_buff_t *sam_buff, unsigned int mapq_min);
void sam_buff_add_alignment(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment);
void sam_buff_out_of_order(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment);
void vcf_out_of_order(vcf_call_t *vcf_call, char *previous_chromosome, size_t previous_pos);
void sam_buff_free_alignment(sam_buff_t *sam_buff, size_t c);
void sam_buff_shift(sam_buff_t *sam_buff, size_t c);
int uchar_cmp(unsigned char *c1, unsigned char *c2);
void stats_update_discarded(ad2vcf_stats_t *stats, sam_alignment_t *sam_alignment);
void ad2vcf_stats_init(ad2vcf_stats_t *stats);
