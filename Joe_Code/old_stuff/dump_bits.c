// dump_bits.c
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(int argc, char **argv){
    if (argc < 2) {
        fprintf(stderr,"Usage: %s <file> [--msb|--lsb] [--width N] [--group G]\n", argv[0]);
        return 1;
    }
    const char *path = argv[1];
    int msb_first = 1;
    int width = 64;   // bits per line
    int group = 8;    // space every N bits

    for (int i=2; i<argc; ++i){
        if (!strcmp(argv[i],"--lsb")) msb_first = 0;
        else if (!strcmp(argv[i],"--msb")) msb_first = 1;
        else if (!strcmp(argv[i],"--width") && i+1<argc) width = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--group") && i+1<argc) group = atoi(argv[++i]);
    }

    FILE *f = fopen(path, "rb");
    if(!f){ perror("open"); return 1; }
    fseek(f,0,SEEK_END);
    long len = ftell(f);
    fseek(f,0,SEEK_SET);
    if(len <= 0){ fclose(f); return 0; }

    uint8_t *buf = (uint8_t*)malloc((size_t)len);
    if(!buf){ fclose(f); return 1; }
    if (fread(buf,1,(size_t)len,f) != (size_t)len){ perror("read"); free(buf); fclose(f); return 1; }
    fclose(f);

    int col = 0;
    int gcol = 0;

    for (long i=0; i<len; ++i){
        uint8_t b = buf[i];
        for (int k=0; k<8; ++k){
            int bit = msb_first ? ((b >> (7-k)) & 1) : ((b >> k) & 1);
            putchar(bit ? '1' : '0');
            if (++col == width) { putchar('\n'); col = 0; gcol = 0; }
            else if (group > 0 && (++gcol % group) == 0) { putchar(' '); }
        }
    }
    if (col) putchar('\n');
    free(buf);
    return 0;
}
