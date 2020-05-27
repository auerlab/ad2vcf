/* ad2vcf.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int ad2vcf(const char *argv[], FILE *sam_stream);
_Bool skip_past_alignments(vcf_call_t *vcf_call, FILE *sam_stream, sam_buff_t *sam_buff, FILE *vcf_out_stream);
_Bool allelic_depth(vcf_call_t *vcf_call, FILE *sam_stream, sam_buff_t *sam_buff, FILE *vcf_out_stream);
_Bool call_in_alignment(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment);
_Bool alignment_behind_call(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment);
void update_allele_count(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment, FILE *vcf_out_stream);
void sam_buff_check_order(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment);
void sam_buff_init(sam_buff_t *sam_buff);
void sam_buff_add_alignment(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment);
void sam_buff_out_of_order(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment);
void vcf_out_of_order(vcf_call_t *vcf_call, char *previous_chromosome, size_t previous_pos);
void sam_alignment_copy(sam_alignment_t *dest, sam_alignment_t *src);
