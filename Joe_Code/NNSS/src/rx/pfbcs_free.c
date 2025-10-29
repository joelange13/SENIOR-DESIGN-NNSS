#include "api.h"
#include "pfb_internal.h"



/*
 * pfbcs_free(cs)
 *  Output: frees state
 */
void pfbcs_free(PFBClockSync *cs){
    if (!cs) return;
    if (cs->proto) free(cs->proto);
    if (cs->bank.b) free(cs->bank.b);
    free(cs);
}
