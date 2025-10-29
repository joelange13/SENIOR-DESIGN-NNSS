#include "api.h"
/*
 * get_opt(argc, argv, key, defv)
 *  Output: value after key or defv if missing
 */
const char* get_opt(int argc, char **argv, const char *key, const char *defv){
    for (int i = 1; i+1 < argc; i++) if (strcmp(argv[i], key) == 0) return argv[i+1];
    return defv;
}
