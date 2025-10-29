#include "api.h"
/*
 * flag_present(argc, argv, flag)
 *  Output: true if argv contains the exact flag
 */
bool flag_present(int argc, char **argv, const char *flag){
    for (int i = 1; i < argc; i++) if (strcmp(argv[i], flag) == 0) return true;
    return false;
}
