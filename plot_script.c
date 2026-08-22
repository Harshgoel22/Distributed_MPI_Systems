#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#define MAX_RUNS   5
#define N_P        4
#define N_N        2
#define MAX_LINE   256

static int   P_VALUES[N_P] = {32, 48, 64, 96};
static int   N_VALUES[N_N] = {120, 240};
static double timing[N_P][N_N][MAX_RUNS];
static int    counts[N_P][N_N];

static int p_index(int P) {
    for (int i = 0; i < N_P; i++) if (P_VALUES[i] == P) return i;
    return -1;
}
static int n_index(int N) {
    for (int i = 0; i < N_N; i++) if (N_VALUES[i] == N) return i;
    return -1;
}

static int read_timing(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char line[MAX_LINE];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return 0; }
    while (fgets(line, sizeof(line), f)) {
        int run, P, N, px, py, pz;
        double t;
        if (sscanf(line, "%d,%d,%d,%d,%d,%d,%lf", &run, &P, &N, &px, &py, &pz, &t) != 7) continue;
        int pi = p_index(P);
        int ni = n_index(N);
        if (pi < 0 || ni < 0 || counts[pi][ni] >= MAX_RUNS) continue;
        timing[pi][ni][counts[pi][ni]++] = t;
    }
    fclose(f);
    return 1;
}

static void write_dat_files(const char *outdir) {
    mkdir(outdir, 0755);
    for (int pi = 0; pi < N_P; pi++) {
        for (int ni = 0; ni < N_N; ni++) {
            char path[256];
            snprintf(path, sizeof(path), "%s/P%d_N%d.dat", outdir, P_VALUES[pi], N_VALUES[ni]);
            FILE *f = fopen(path, "w");
            if (!f) continue;
            for (int r = 0; r < counts[pi][ni]; r++) fprintf(f, "%.6f\n", timing[pi][ni][r]);
            fclose(f);
        }
    }
}

static void write_gnuplot_script(const char *gp_file, const char *dat_dir) {
    FILE *f = fopen(gp_file, "w");
    if (!f) return;

    fprintf(f,
        "set terminal pngcairo enhanced font 'Arial,12' size 1000,600\n"
        "set output 'boxplot_timing.png'\n"
        "set style data boxplot\n"
        "set style boxplot range 1.5 outliers pointtype 7\n"
        "set style fill solid 0.5 border -1\n"
        "set boxwidth 0.2\n"
        "set grid ytics lc rgb '#dddddd'\n"
        "set ylabel 'Execution Time (seconds)'\n"
        "set xlabel 'Number of Processes (P)'\n"
        "set title 'N=120 (Blue) vs N=240 (Orange)'\n"
        "set yrange [0:20]\n"
        "set ytics 2\n"
        "set xtics ('32' 1, '48' 2, '64' 3, '96' 4)\n"
        "set xrange [0.5:4.5]\n"
        "unset key\n"
        "plot \\\n"
    );

    for (int pi = 0; pi < N_P; pi++) {
        fprintf(f, " '%s/P%d_N120.dat' using (%d):(column(1)) with boxplot lt 1 lc rgb '#1f77b4', \\\n",
                dat_dir, P_VALUES[pi], pi + 1);
        fprintf(f, " '%s/P%d_N240.dat' using (%d):(column(1)) with boxplot lt 1 lc rgb '#ff7f0e'%s",
                dat_dir, P_VALUES[pi], pi + 1, (pi < N_P - 1 ? ", \\\n" : "\n"));
    }
    fclose(f);
}

int main(int argc, char *argv[]) {
    const char *timing_file = (argc > 1) ? argv[1] : "timing.txt";
    const char *dat_dir = "plot_data";
    const char *gp_file = "boxplot.gp";

    memset(counts, 0, sizeof(counts));
    if (!read_timing(timing_file)) return 1;

    write_dat_files(dat_dir);
    write_gnuplot_script(gp_file, dat_dir);
    system("gnuplot boxplot.gp");

    return 0;
}