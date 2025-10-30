#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

static void inject_cfo(float *xIQ, size_t nin_complex, double cps)
{
    const double w = 2.0 * M_PI * cps;
    double phi = 0.0;
    for (size_t n = 0; n < nin_complex; ++n) {
        double I = xIQ[2*n + 0];
        double Q = xIQ[2*n + 1];
        double c = cos(phi);
        double s = sin(phi);
        double Ii = I * c - Q * s;
        double Qi = I * s + Q * c;
        xIQ[2*n + 0] = (float)Ii;
        xIQ[2*n + 1] = (float)Qi;
        phi += w;
        if (phi >  M_PI) phi -= 2.0 * M_PI;
        if (phi < -M_PI) phi += 2.0 * M_PI;
    }
}

int main(void)
{
    const char *in_path  = "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\goodbye.cfile";
    const char *out_path = "C:\\Users\\joebo\\OneDrive\\Desktop\\Scalable\\SENIOR-DESIGN-NNSS\\Joe_Code\\NNSS\\data\\goodbye_offset.cfile";

    // CFO value (cycles per sample)
    const double cps = 0.001;  // Change this as needed

    // Read input file
    FILE *fin = fopen(in_path, "rb");
    if (!fin) {
        perror("open input");
        return 1;
    }

    if (fseek(fin, 0, SEEK_END) != 0) { perror("fseek"); fclose(fin); return 1; }
    long nbytes = ftell(fin);
    if (nbytes < 0) { perror("ftell"); fclose(fin); return 1; }
    rewind(fin);

    size_t nfloats = (size_t)nbytes / sizeof(float);
    if (nfloats % 2 != 0) {
        fprintf(stderr, "Error: file does not contain an even number of floats (I/Q mismatch)\n");
        fclose(fin);
        return 1;
    }

    size_t nin_c = nfloats / 2;
    float *xIQ = (float*)malloc(sizeof(float) * nfloats);
    if (!xIQ) {
        fprintf(stderr, "malloc failed\n");
        fclose(fin);
        return 1;
    }

    size_t rd = fread(xIQ, sizeof(float), nfloats, fin);
    fclose(fin);
    
    if (rd != nfloats) {
        fprintf(stderr, "Short read: got %zu / %zu floats\n", rd, nfloats);
        free(xIQ);
        return 1;
    }

    printf("Read %zu complex samples from: %s\n", nin_c, in_path);

    // Inject CFO
    inject_cfo(xIQ, nin_c, cps);

    // Write to NEW output file
    FILE *fout = fopen(out_path, "wb");
    if (!fout) {
        perror("open output");
        free(xIQ);
        return 1;
    }

    size_t wr = fwrite(xIQ, sizeof(float), nfloats, fout);
    fclose(fout);
    free(xIQ);

    if (wr != nfloats) {
        fprintf(stderr, "Short write: %zu / %zu floats\n", wr, nfloats);
        return 1;
    }

    printf("Injected CFO of %.6g cycles/sample\n", cps);
    printf("Wrote %zu complex samples to: %s\n", nin_c, out_path);
    printf("\nNow run: ./build/tx_bpsk --mode rx --in %s\n", out_path);
    
    return 0;
}