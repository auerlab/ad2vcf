
#ifndef _SYS_STDINT_H_
#include <stdint.h>
#endif

#define     CMD_MAX                 128
#define     SAM_BUFF_MAX_ALIGNMENTS 524288

// FIXME: These should be command line arguments
#define     MAPQ_MIN                10
#define     PHRED_MIN               20
#define     PHRED_BASE              33

// FIXME: Move this to samio when complete?
typedef struct
{
    size_t  count;
    size_t  max_count;
    size_t  previous_pos;
    char    previous_rname[SAM_RNAME_MAX_CHARS + 1];
    sam_alignment_t *alignments[SAM_BUFF_MAX_ALIGNMENTS];
}   sam_buff_t;

#define     AD2VCF_STATS_INIT   { 0, 0, 0, 0, SIZE_MAX, 0, 0, 0, 0, SIZE_MAX, 0, 0 }
typedef struct
{
    size_t  total_vcf_calls,
	    total_sam_alignments,
	    discarded_sam_alignments,
	    discarded_score_sum,
	    min_discarded_score,
	    max_discarded_score,
	    total_ref_alleles,
	    total_alt_alleles,
	    total_other_alleles,
	    min_depth,
	    max_depth,
	    // Median would require an array of all VCF calls and we can
	    // always get it later from the -ad output
	    mean_depth;
}   ad2vcf_stats_t;

#include "ad2vcf-protos.h"
