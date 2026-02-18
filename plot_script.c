#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PROCESSES 100

int main() {
    FILE *fp = fopen("timing.txt", "r");
    if (!fp) {
        printf("Cannot open timing.txt\n");
        return 1;
    }

    FILE *data = fopen("plot_data.dat", "w");
    if (!data) {
        printf("Cannot open plot_data.dat for writing\n");
        fclose(fp);
        return 1;
    }

    char line[512];
    int run, P;
    long M;
    double D1, D2, T, seed, maxD1, maxD2, time;
    
    // Track unique process counts
    int unique_P[MAX_PROCESSES];
    int num_unique = 0;

    // Parse timing.txt and extract P, Time, M
    while (fgets(line, sizeof(line), fp)) {
        if (strlen(line) <= 1) continue;  // skip empty lines
        if (line[0] == '#') continue;     // skip comments
        if (strncmp(line, "Run", 3) == 0) continue; // skip header

        int n = sscanf(line, "%d %d %ld %lf %lf %lf %lf %lf %lf %lf",
                       &run, &P, &M, &D1, &D2, &T, &seed, &maxD1, &maxD2, &time);
        if (n != 10) continue; // skip malformed lines

        fprintf(data, "%d %lf %ld\n", P, time, M);
        
        // Track unique P values
        int found = 0;
        for (int i = 0; i < num_unique; i++) {
            if (unique_P[i] == P) {
                found = 1;
                break;
            }
        }
        if (!found && num_unique < MAX_PROCESSES) {
            unique_P[num_unique++] = P;
        }
    }

    fclose(fp);
    fclose(data);

    if (num_unique == 0) {
        printf("No valid data found\n");
        return 1;
    }

    // Sort unique_P for better ordering
    for (int i = 0; i < num_unique - 1; i++) {
        for (int j = i + 1; j < num_unique; j++) {
            if (unique_P[i] > unique_P[j]) {
                int temp = unique_P[i];
                unique_P[i] = unique_P[j];
                unique_P[j] = temp;
            }
        }
    }

    // Write gnuplot script
    FILE *gp = fopen("boxplot.gnu", "w");
    if (!gp) {
        printf("Cannot open boxplot.gnu\n");
        return 1;
    }

    fprintf(gp,
        "set terminal png size 900,600\n"
        "set output 'boxplot_time_vs_process.png'\n"
        "set title 'Execution Time vs Processes'\n"
        "set xlabel 'Processes (P)'\n"
        "set ylabel 'Time (seconds)'\n"
        "set style data boxplot\n"
        "set style boxplot outliers pointtype 7\n"
        "set grid ytics\n"
    );

    // Generate xtics
    fprintf(gp, "set xtics (");
    for (int i = 0; i < num_unique; i++) {
        fprintf(gp, "'%d' %d", unique_P[i], i + 1);
        if (i < num_unique - 1) fprintf(gp, ", ");
    }
    fprintf(gp, ")\n");

    fprintf(gp, "set key outside\n");

    // Define colors array
    const char *colors[] = {"red", "blue", "green", "orange", "purple", 
                           "brown", "pink", "cyan", "magenta", "yellow"};
    int num_colors = 10;

    // Generate plot command
    fprintf(gp, "plot ");
    for (int i = 0; i < num_unique; i++) {
        fprintf(gp, "'plot_data.dat' using (column(1)==%d ? %d:1/0):2 notitle lc rgb '%s'",
                unique_P[i], i + 1, colors[i % num_colors]);
        if (i < num_unique - 1) {
            fprintf(gp, ", \\\n     ");
        }
    }
    fprintf(gp, "\n");

    fclose(gp);

    // Run gnuplot
    system("gnuplot boxplot.gnu");

    printf("Plot saved as boxplot_time_vs_process.png\n");
    printf("Found %d unique process counts: ", num_unique);
    for (int i = 0; i < num_unique; i++) {
        printf("%d ", unique_P[i]);
    }
    printf("\n");

    return 0;
}