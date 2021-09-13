#include <stdio.h>
#include <stdlib.h>

double** spk_linker(int* k, int goal, double** obs, int N, int d);
void wam_func(double** dp, double** W, int N, int d);
double l2_norm(double* v1, double* v2, int d);
void ddg_func(double** W, double** D, int N);
void lnorm_func(double** W, double** D,double** Lnorm, int N);
void jacobi_func(double** A, double** V, double* eiganvalues, int N);
void build_A_tag(double **A, double** A_tag ,int N, double c, double s, int ip, int jp);
int* pivot(double** A, int N);
double off_sqr(double** A, int N);
void eigan_union(double** V,double* eiganvalues,double** eigan,int N);
void spk_func(double** V, double** T, int N, int k);
void order_eigan(double** V,double* eiganvalues, int N);
void bubbleSort(double arr[], double arr2[], int N);
void swap(double* a, double* b);
int eigangap_heuristic(double* ev, int N);
double** k_means(int k, int N, int d, int* first, double** T);
double calculate_difference(double **current_clusters, double* vector, int cluster, int d);
int check_min_cluster(double **current_clusters, double *vector, int k, int d);
void print_mat(double** M, int dimr, int dimc);
int num_of_lines(FILE *fp);
int num_of_columns(FILE *fp);
