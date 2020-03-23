/* ad2vcf.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int ad2vcf(const char *argv[], FILE *sam_stream);
int sam_read_alignment(const char *argv[], FILE *sam_stream, sam_alignment_t *sam_alignment);
