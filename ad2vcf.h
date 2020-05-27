
#define     CMD_MAX                 128
#define     MAX_BUFFERED_ALIGNMENTS 32768
#define     SAM_BUFF_INIT           { 0, 0 }

// FIXME: Move this to samio when complete?
typedef struct
{
    size_t  count;
    size_t  max_count;
    sam_alignment_t *alignments[MAX_BUFFERED_ALIGNMENTS];
}   sam_buff_t;

#include "ad2vcf-protos.h"
