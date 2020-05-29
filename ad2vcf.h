
#define     CMD_MAX                 128
#define     SAM_BUFF_MAX_ALIGNMENTS 65536

// FIXME: Move this to samio when complete?
typedef struct
{
    size_t  count;
    size_t  max_count;
    size_t  previous_pos;
    char    previous_rname[SAM_RNAME_MAX_CHARS + 1];
    sam_alignment_t *alignments[SAM_BUFF_MAX_ALIGNMENTS];
}   sam_buff_t;

#include "ad2vcf-protos.h"
