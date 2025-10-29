#include "api.h"
/*
 * OpenFile(path, is_input_file)
 *  Input : path (const char*), is_input_file=true→"rb", false→"wb"
 *  Output: FILE* (opened stream) or NULL on error (prints perror)
 */
FILE* OpenFile(const char* path, bool is_input_file){
    if(is_input_file){
        FILE* ioFile = fopen(path, "rb");
        if(!ioFile){ perror("Error Opening Input File"); return NULL; }
        printf("Opened Input File.\n");
        return ioFile;
    }else{
        FILE* ioFile = fopen(path, "wb");
        if(!ioFile){ perror("Error Opening Output File"); return NULL; }
        printf("Opened Output File.\n");
        return ioFile;
    }
}
