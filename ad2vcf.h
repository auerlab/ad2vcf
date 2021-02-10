
#ifndef _SYS_STDINT_H_
#include <stdint.h>
#endif

#define     CMD_MAX                 128
/*
    256k was not enough for a few of the SRA CRAMs.
    NWD976804 needed more than 512k.  Bad data?
    Of 55k samples, only about 11% reached a sam buffer of > 8k.
    Set an upper limit to prevent runaway memory use and error out
    with EX_DATAERR to tell script not to retry.

    From sam-buff-stats script:
    
    Total      54982
    >4096      13723
    >8192       5824
    >16384      2795
    >32768      1612
    >65536       138
    >131072       33
    >262144       12
    >524288        0
*/
#define     SAM_BUFF_START_SIZE     4096
#define     SAM_BUFF_MAX_SIZE       524288

// FIXME: These should be command line arguments
#define     PHRED_MIN               20
#define     PHRED_BASE              33

// Yes, we actually saw a few INFO fields over 512k in some dbGap BCFs
// Match this with vcf-split
#define VCF_INFO_MAX_CHARS          1048576
#define VCF_FORMAT_MAX_CHARS        4096
#define VCF_SAMPLE_MAX_CHARS        2048

// FIXME: Move this to samio when complete?
typedef struct
{
    size_t          buff_size;;
    sam_alignment_t **alignments;
    size_t          buffered_count;
    size_t          max_count;
    size_t          previous_pos;
    char            previous_rname[SAM_RNAME_MAX_CHARS + 1];
    
    // Use 64 bits to accommodate large sums
    size_t          mapq_min,
		    mapq_low,
		    mapq_high,
		    mapq_sum,
		    reads_used,
		    total_alignments,
		    trailing_alignments,
		    discarded_alignments,
		    discarded_score_sum,
		    discarded_trailing,
		    min_discarded_score,
		    max_discarded_score,
		    unmapped_alignments;
}   sam_buff_t;

#define SAM_BUFF_MAPQ_MIN(sam_buff) ((sam_buff)->mapq_min)
#define SAM_BUFF_MAPQ_SUM(sam_buff) ((sam_buff)->mapq_sum)
#define SAM_BUFF_MAPQ_LOW(sam_buff) ((sam_buff)->mapq_low)
#define SAM_BUFF_MAPQ_HIGH(sam_buff) ((sam_buff)->mapq_high)
#define SAM_BUFF_READS_USED(sam_buff) ((sam_buff)->reads_used)

/*
 *  FIXME: This is a foster home for a random collection of unrelated stats.
 *  Find these data a permanent home.
 */

typedef struct
{
    size_t  total_vcf_calls,
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
