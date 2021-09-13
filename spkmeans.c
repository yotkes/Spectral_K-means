#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>

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


/* Main flow function that returns the specific matrix defined by goal's value */
double** spk_linker(int* z, int goal, double** obs, int N, int d)
{
    /* Variables Decleration */
    int i;
    int k;
    double *eiganvalues;
    double **W;
    double **D;
    double **Lnorm;
    double **V;
    double **T;
    double** eigan;
        
    /* Jacobi function as a stand-alone procedure (input: real-symmetric matrix) */
    if(goal == 4) /* goal == jacobi */
    {
        /* Initializing matrices for jacobi procedeure */
        
        /* Memory allocation */
        V = (double**)calloc(N,sizeof(double *));
        eiganvalues = (double*)malloc(N*sizeof(double));
        if (V == NULL || eiganvalues == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }

        for (i = 0; i < N; i++)
        {
            V[i] = (double*)calloc(N,sizeof(double));

            if (V[i] == NULL) 
            {
                printf("%s\n","An Error Has Occured");
                exit(1);
            }
        }

        eigan = (double**)malloc((N+1)*sizeof(double *));
        if (eigan == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }

        for (i = 0; i < N+1; i++)
        {
            eigan[i] = (double*)malloc(N*sizeof(double));
            if (eigan[i] == NULL) 
            {
                printf("%s\n","An Error Has Occured");
                exit(1);
            }
        }

        jacobi_func(obs,V,eiganvalues,N); /* Diagonalizing real-symmetric matrix via jacobi procedure */
        eigan_union(V,eiganvalues,eigan,N); /* Union of eiganvalues & eiganvectors */

        /* Free Memory */
        for(i=0;i<N;i++)
        {
            free(V[i]);
        } 
        free(V);
        free(eiganvalues);

        return eigan;
    }

    /* Running specific proccess defined by user's input */
    
    /* Memory allocation */
    W = (double**)malloc(N * sizeof(double *));
    if (W == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for (i = 0; i < N; i++)
    {
        W[i] = (double*)malloc(N*sizeof(double));
        
        if (W[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    /* from this point on - WAM procedure is always in use */
    wam_func(obs,W,N,d);

    if(goal == 1) /* goal == wam */
    {
        return W;
    }

    /* Memory allocation */
    D = (double**)calloc(N,sizeof(double *));
    if (D == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for (i = 0; i < N; i++)
    {
        D[i] = (double*)calloc(N,sizeof(double));
        if (D[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }  
    }

    /* from this point on - DDG procedure is always in use */
    ddg_func(W,D,N);

    if(goal == 2) /* goal == ddg */
    {   
        for(i=0;i<N;i++)
        {
            free(W[i]);
        }        
        free(W);
        
        return D;
    }

    /* Memory allocation */
    Lnorm = (double**)calloc(N,sizeof(double *));
    if (Lnorm == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }
    
    for (i = 0; i < N; i++)
    {
        Lnorm[i] = (double*)calloc(N,sizeof(double));
        if (Lnorm[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    /* from this point on - LNORM procedure is always in use */
    lnorm_func(W,D,Lnorm,N);

    if(goal == 3) /* goal == lnorm */
    {
        for(i=0;i<N;i++)
        {
            free(W[i]);
        }
        free(W);

        for(i=0;i<N;i++)
        {
            free(D[i]);
        }
        free(D);
        return Lnorm;
    }

    /* Memory allocation */
    V = (double**)calloc(N,sizeof(double *));
    eiganvalues = (double*)malloc(N*sizeof(double));
    if (V == NULL || eiganvalues == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for (i = 0; i < N; i++)
    {
        V[i] = (double*)calloc(N,sizeof(double));
        if (V[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    /* If did not return matrix till here: goal == spk */

    jacobi_func(Lnorm,V,eiganvalues,N); /* Diagonalizing L_norm via jacobi procedure */
    order_eigan(V,eiganvalues,N);   /* Ordering of eiganvalues & eiganvectors - for spkmeans purposes */

    k = *z;  /* Deriving k's value from its address */

    /* Call Eigangap Heuristic - to indicate k's value */
    if(k == 0)
    { 
        *z = eigangap_heuristic(eiganvalues,N);
        k = *z;
    }

    /* Memory allocation */
    T = (double**)calloc(N,sizeof(double *));
    if (T == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for (i=0;i<N;i++){
        T[i] = (double*)calloc(k,sizeof(double));
        if (T[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    spk_func(V,T,N,k);

    /* Free Memory */
    for(i=0;i<N;i++)
    {
        free(W[i]);
    }
    free(W);

    for(i=0;i<N;i++)
    {
        free(D[i]);
    }
    free(D);

    for(i=0;i<N;i++)
    {
        free(Lnorm[i]);
    }
    free(Lnorm);

    for(i=0;i<N;i++)
    {
        free(V[i]);
    }
    free(V);
    free(eiganvalues);

    /* Returns T if goal == spk */
    return T;
}

/* WAM function - updating W matrix (by reference) */
void wam_func(double** dp, double** W, int N, int d)
{
    /* Variables Decleration */
    int i,j;
    double norm;

    for(i=0;i<N;i++)
    {
        j=N-1;
        while(j>=i)
        {
            if(i==j)
            { 
                W[i][j]=0;
            }
            else
            {
                norm = l2_norm(dp[i],dp[j],d);
                W[i][j] = exp(-norm/2);
                W[j][i] = W[i][j];
            }
            j--;
        }
    }
}

/* Assistive function for WAM - l2 norm (or the Euclidean norm) */
double l2_norm(double* v1, double* v2, int d)
{
    /* Variables Decleration */
    int i = 0;
    double norm = 0;

    while(i<d)
    {
        norm += pow((v1[i]-v2[i]),2);
        i++;
    }
    return sqrt(norm);
}

/* DDG function - updating D matrix (by reference) */
void ddg_func(double** W, double** D, int N)
{
    /* Variables Decleration */
    int i,j; 
    double sum;

    for(i=0;i<N;i++)
    {
        sum = 0;
        for(j=0;j<N;j++)
        {
            sum += W[i][j];
        }
        D[i][i] = sum;
    }
}

/* LNORM function - updating Lnorm matrix (by reference) */
void lnorm_func(double** W, double** D, double** Lnorm, int N)
{
    /* Variables Decleration */
    double** D_nsqrt; /* D matrix raised to the power of -0.5 */
    double** T1, **T2; /* Assistive matrices - for matrices multiplication */
    int i,j,r;

    /* Memory allocation */
    D_nsqrt = (double**)calloc(N, sizeof(double *));  
    T1 = (double**)calloc(N, sizeof(double *)); 
    T2 = (double**)calloc(N, sizeof(double *));

    if (D_nsqrt == NULL || T1 == NULL || T2 == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for(i=0;i<N;i++)
    {
        D_nsqrt[i] = (double*)calloc(N, sizeof(double));
        T1[i] = (double*)calloc(N, sizeof(double));
        T2[i] = (double*)calloc(N, sizeof(double));

        if (D_nsqrt[i] == NULL || T1[i] == NULL || T2[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    /* Calculating matrix D in the power of -0.5 */
    for(i=0;i<N;i++)
    {
        if(D[i][i] != 0)
        {
            D_nsqrt[i][i] = 1/(sqrt(D[i][i]));
        }
    }

    /* First matrices multiplication - (D^-0.5) * (W) */
    for(i=0;i<N;i++)    
    {    
        for(j=0;j<N;j++)    
        {    
            for(r=0;r<N;r++)    
            {    
                T1[i][j]+= D_nsqrt[i][r]*W[r][j];    
            }
        }    
    }    

    /* Second matrices multiplication - (D^-0.5 * W) * (D^-0.5) */
    for(i=0;i<N;i++)    
    {    
        for(j=0;j<N;j++)    
        {    
            for(r=0;r<N;r++)    
            {    
                T2[i][j] += T1[i][r]*D_nsqrt[r][j];    
            }
        }
    }  

    /* Calculating matrix L_norm  */
    for(i=0;i<N;i++)    
    {    
        for(j=0;j<N;j++)    
        {
            if(i==j)
            {
                Lnorm[i][j] = 1-T2[i][j];
            }   
            else
                Lnorm[i][j] = -T2[i][j];
        }
    }

    /* Free Memory */
    for(i=0;i<N;i++)
    {
        free(D_nsqrt[i]);
    }
    free(D_nsqrt);

    for(i=0;i<N;i++)
    {
        free(T1[i]);
    }
    free(T1);

    for(i=0;i<N;i++)
    {
        free(T2[i]);
    }
    free(T2);
}

/* JACOBI function - updating V matrix (eiganvectors matrix) & eiganvalues array - both by reference */
void jacobi_func(double** A, double** V, double* eiganvalues, int N)
{   
    /* Variables Decleration */
    int i,j,ip,jp,iter_count=0;
    int sign_theta;
    double abs_theta;
    double theta, t, c, s;
    double** A_tag, **ptr;
    double** A_org, **A_tag_org;   /* Pointers for original matrices' address */
    int* pivot_arr;
    double epsilon = pow(10, -15);  /* Threshold for algorithm's convergence */ 
    
    /* Memory Allocation */
    A_tag = (double**)calloc(N, sizeof(double *)); /* A' is the iterative matrix for this procedure */
    if (A_tag == NULL)
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }
    
    for(i = 0; i < N; i++)
    {
        A_tag[i] = (double*)calloc(N, sizeof(double));
        if (A_tag[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }
    
    A_org = A;
    A_tag_org = A_tag;

    /* Iterating for convergence */
    while(iter_count < 100) 
    {
        /* Extracting i & j (pivot indexes) */
        pivot_arr = pivot(A,N);
        ip = pivot_arr[0];
        jp = pivot_arr[1];

        /* Calculating theta, sign(theta), t, c, s */
        theta = (A[jp][jp]-A[ip][ip])/(2*A[ip][jp]);
        abs_theta = theta;

        if(theta >= 0)
        {
            sign_theta = 1;
        }
        else
        {
            sign_theta = -1;
            abs_theta = -theta;
        }

        t = (sign_theta)/(abs_theta+sqrt((pow(theta,2)+1)));
        c = 1/(sqrt(pow(t,2)+1));
        s = t*c;
        

        /* Initializing V matrix's elements (eigenvectors' matrix) */
        if(iter_count == 0)
        {
            for(i = 0;i<N;i++)
            { 
                V[i][i] = 1;
            }
            V[ip][ip] = c;
            V[jp][jp]= c;
            V[ip][jp] = s;
            V[jp][ip] = -s;
        }
        else
        {
            for(i=0;i<N;i++)
            {
                double temp1,temp2;
                temp1 = V[i][ip]*c-V[i][jp]*s;
                temp2 = V[i][ip]*s+V[i][jp]*c;

                V[i][ip] = temp1;
                V[i][jp] = temp2;
            }
        }

        build_A_tag(A, A_tag, N, c, s, ip, jp); /* Build A' */

        if((off_sqr(A,N)-off_sqr(A_tag,N)) <= epsilon) /* Convergence condition */
        {
            break;
        }
        ptr = A;
        A = A_tag;
        A_tag = ptr;

        iter_count++;
    }

    if(iter_count < 100)   /* Updating A to point to A' after lastr iteration */
    {
        ptr = A;
        A = A_tag;
        A_tag = ptr;
    }
    
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            A_org[i][j] = A[i][j];
        } 
    }

    for(i=0;i<N;i++)
    {
        eiganvalues[i] = A_org[i][i];
    }

    /* Free Memory */
    for(i=0;i<N;i++)
    {
        free(A_tag_org[i]);
    }
    free(A_tag_org);
    free(pivot_arr);
}

/* Assistive function for JACOBI procedure - updating A' matrix by pivot's indices */
void build_A_tag(double **A, double **A_tag, int N, double c, double s, int ip, int jp)
{
    /* Variables Decleration */
    int r,q;

    for(r=0;r<N;r++)
    {
        for(q=r;q<N;q++)
        {
            if(q==ip && r!=ip && r!= jp)
            {
                A_tag[r][ip] = c*A[r][ip]-s*A[r][jp];
                A_tag[ip][r] = A_tag[r][ip];
            }
            else if(r==ip && q!=ip && q!=jp)
            {
                A_tag[q][ip] = c*A[q][ip]-s*A[q][jp];
                A_tag[ip][q] = A_tag[q][ip];
            }
            else if(q==jp && r!=ip && r!=jp)
            {
                A_tag[r][jp] = c*A[r][jp]+s*A[r][ip];
                A_tag[jp][r] = A_tag[r][jp];
            }
            else if(r==jp && q!=ip && q!=jp)
            {
                A_tag[q][jp] = c*A[q][jp]+s*A[q][ip];
                A_tag[jp][q] = A_tag[q][jp];
            }
            else if(r==ip && q==ip)
            {
                A_tag[ip][ip] = (c*c)*A[ip][ip]+(s*s)*A[jp][jp]-2*s*c*A[ip][jp];
            }
            else if(r==jp && q==jp)
            {
                A_tag[jp][jp] = (s*s)*A[ip][ip]+(c*c)*A[jp][jp]+2*s*c*A[ip][jp];
            }
            else if(r==ip && q==jp)
            {
                A_tag[ip][jp] = (c*c-s*s)*A[ip][jp]+(s*c)*(A[ip][ip]-A[jp][jp]);
                A_tag[jp][ip] = A_tag[ip][jp];
            }
            else
            {
                A_tag[r][q] = A[r][q];
                A_tag[q][r] = A_tag[r][q];
            }
        }   
    }
}

/* Assistive function for JACOBI procedure - finding pivot's indices */
int* pivot(double** A, int N)
{
    /* Variables Decleration */
    int ip = 0, jp = 1; /* Initialization for pivot's indices - always off-diagonal element */
    int i,j;
    double p = 0;
    double temp;
    int* pivot_arr;

    /* Memory Allocation */
    pivot_arr = (int*)malloc(2*sizeof(int));
      if (pivot_arr == NULL)
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for(i=0;i<N;i++)
    {
        for(j=i+1;j<N;j++)
        {
            if(A[i][j]!=0)
            {
                if(A[i][j]<0)
                {
                    temp = -A[i][j];
                }
                else
                    temp = A[i][j];
                if(p < temp)
                {
                    ip = i;
                    jp = j;
                    p = temp;
                }   
            }
        }
    }
    pivot_arr[0] = ip;
    pivot_arr[1] = jp;
    return pivot_arr;
}

/* Assistive function for JACOBI procedure - OFF function squared (for Jacobi's convergence condition) */
double off_sqr(double** A, int N)
{
    /* Variables Decleration */
    int i,j;
    double offA_sqr = 0;

    for(i=0;i<N;i++)
    {
        for(j=i+1;j<N;j++)
        {
            offA_sqr += 2*A[i][j]*A[i][j];   /* A is a real-symmetric matrix => A[i][j] = A[j][i] */           
        }
    }
    return offA_sqr;
}

/* Assistive function for JACOBI function - union of eiganvalues (first row) & eiganvectors for 1 matrix (in order of printing output) */
void eigan_union(double** V,double* eiganvalues,double** eigan,int N)
{
    int i,j;
    for(i=0;i<N;i++)
    {
        /* Inserting eiganvalues */
        eigan[0][i] = eiganvalues[i];
    }
    for (i=1;i<N+1;i++)
    {
        for (j=0;j<N;j++)
        {
            /* Inserting eigenvectors as row vectors */
            eigan[i][j] = V[j][i-1];
        }
    }
}
 
/* SPK function - updating T matrix (by reference) and normalizing it */
void spk_func(double** V, double** T, int N, int k)
{
    /* Variables Decleration */
    double* v_zero;
    double norm_row;
    int i,j;

    v_zero = calloc(k,sizeof(double));
    if (v_zero == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for(i = 0; i < N; i++)
    {
        for(j = 0; j < k; j++)
        {
            T[i][j] = V[i][j];
        }
    }

    for(i = 0; i < N; i++)
    {
        norm_row = l2_norm(T[i],v_zero,k);
        for(j = 0; j < k; j++)
        {
            T[i][j] = T[i][j]/norm_row;
        }
    }
    
    /* Free Memory */
    free(v_zero);
}

/* Assistive function for JACOBI function - ordering eiganvalues & eiganvectors by ascending order */
void order_eigan(double** V, double* eiganvalues, int N)
{
    /* Variables Decleration */
    double *tmp_arr; /* Pointer for assistive array */
    double **ord_eiganvec; /* Ordered eiganvectors matrix */
    int i,j,col;

    tmp_arr = (double*)malloc(N * sizeof(double));   /* Assistive array for sorting */
    ord_eiganvec = (double**)malloc(N * sizeof(double *));  /* Final ordered eiganvectors matrix */ 

    if (tmp_arr == NULL || ord_eiganvec == NULL)
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for (i = 0; i<N ; i++)
    {
        ord_eiganvec[i] = (double*)malloc(N*sizeof(double));
        if (ord_eiganvec[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }

        tmp_arr[i] = i;     /* Assistive array for keeping the original column's index */
    }

    /* Sorting eiganvalues & assistive array (in order of keeping the original matrix's columns indices) */
    bubbleSort(eiganvalues, tmp_arr, N);

    /* Sorting eiganvectors */
    for(j=0;j<N;j++)
    {
        col = tmp_arr[j];
        for(i=0;i<N;i++)
        {
            ord_eiganvec[i][j] = V[i][col];
        }
    }

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            V[i][j] = ord_eiganvec[i][j];
        } 
    }

    /* Free Memory */
    for (i = 0; i<N ; i++)
    {
        free(ord_eiganvec[i]);
    }
    free(ord_eiganvec);
    free(tmp_arr);
}

/* Bubble Sort algorithm (Stable Algorithm) - based on an assistive function: swap */
/* Code's Source: https://www.geeksforgeeks.org/bubble-sort/ */
void bubbleSort(double arr[], double arr2[], int N) 
{ 
    /* Variables Decleration */
    int i, j; 
    for (i = 0; i < N-1; i++)     
      
    /* Last i elements are already in place */ 
    for (j = 0; j < N-i-1; j++) 
        if (arr[j] > arr[j+1]) 
        {
            swap(&arr[j], &arr[j+1]);
            swap(&arr2[j], &arr2[j+1]);
        }     
} 
  
/* Assistive function for bubbleSort - swapping 2 values */
/* Code's Source: https://www.geeksforgeeks.org/bubble-sort/ */
void swap(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

/* Eigangap Heuristic - this funciton indicates k when user inputs: k=0 */
int eigangap_heuristic(double* ev, int N)
{
    /* Variables Decleration */
    double argmax = 0;
    int j,max_ind;
    int i = 0;

    max_ind = N/2;

    for(j=0;j<max_ind;j++)
    {
        if((ev[j+1]-ev[j]) > argmax)
        {
            argmax = ev[j+1]-ev[j];
            i = j;
        }
    }
    return (i+1);  /* Increasing output value by 1 according to loop range */
}

/* K-means algorithm that handles both versions (first/calculated k initial centroids) */
double** k_means(int k, int N, int d, int* first, double** T){
    /* Variables Decleration */
    int i;
    int j;
    int counter;
    int centroid;
    int prev;
    int change;
    int cluster;
    double new_val;
    int *counts;
    double **sums;
    int *members;
    double **current_clusters;
    int MAX_ITER = 300;

    /* Memory Allocation */
    counts = calloc(k, sizeof(int));
    if (counts == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    sums = calloc(k, sizeof(double*));
    if (sums == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for (i = 0; i<k; i++)
    {
        sums[i] = calloc(d,sizeof(double));
        if (sums[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    members = (int*)calloc(N, sizeof(int));
    if (members == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    current_clusters = (double**)calloc(k,sizeof(double*));
    if (current_clusters == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for (i = 0; i<k; i++)
    {
        current_clusters[i] = (double*)calloc(d,sizeof(double));
        if (current_clusters[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    for(i = 0; i<N; i++){
        members[i] = k;
    }

    /* Initializing k centroids */
    for(i = 0; i < k; i++) {
        for(j = 0; j < d; j++) {
            sums[i][j] = T[first[i]][j];
            current_clusters[i][j] = T[first[i]][j];
        }
        counts[i] = 1;
        members[first[i]] = i;
    }

    for(counter = 0; counter < MAX_ITER; counter++){
        change = 1;
        for(i = 0; i < N; i++){
            centroid = check_min_cluster(current_clusters, T[i], k, d);
            prev = members[i];
            if(prev != centroid){
                for(j = 0; j < d; j++){
                    if (prev != k){
                        sums[prev][j] -= T[i][j];
                    }
                    sums[centroid][j] += T[i][j];
                }
                if (prev != k){
                    counts[prev] -= 1;
                }
                counts[centroid] += 1;
                members[i] = centroid;
            }
        }
        for(cluster = 0; cluster < k; cluster++){
            for(j = 0; j<d; j++){
                new_val = sums[cluster][j]/counts[cluster];
                if (current_clusters[cluster][j] != new_val){
                    change = 0;
                }
                current_clusters[cluster][j] = new_val;
            }
        }

        if(change){
            break;
        }
    }
   
    /* Free Memory */
    free(counts);
    for(i=0;i<k;i++)
    {
        free(sums[i]);
    }
    free(sums);
    free(members);

    return current_clusters;
}

/* Assistive function for K-means - calculating the distance between observation and a centroid of a cluster */
double calculate_difference(double **current_clusters, double* vector, int cluster, int d)
{
    /* Variables Decleration */
    double sum = 0;
    int j;
    for (j = 0; j<d; j++){
        sum += ((vector[j] - current_clusters[cluster][j])*(vector[j] - current_clusters[cluster][j]));
    }
    return sum;
}

/* Assistive function for K-means - checks the closest cluster for a specific observation (vector) */
int check_min_cluster(double **current_clusters, double *vector, int k, int d)
{
    /* Variables Decleration */
    int i;
    int cluster;
    double new_dist;
    double dist;
    i = 0;
    dist = calculate_difference(current_clusters, vector, i, d);
    for (cluster = 1; cluster<k; cluster++){
        new_dist = calculate_difference(current_clusters, vector, cluster, d);
        if (new_dist<dist){
            i = cluster;
            dist = new_dist;
        }
    }
    return i;
}

/* Printing Matrix function */
void print_mat(double** M, int dimr, int dimc)
{   
    /* Variables Decleration */
    int i,j;
    double threshold = 0.00005;

    /* Printing the matrix's elements as row vectors */
    for(i=0;i<dimr;i++)
    {
        if(M[i][0] < 0)
        {
            /* Matrix's element is negative and its absolute value is smaller than the threshold */
            if(-M[i][0] < threshold) 
            {
                printf("0.0000");
            }
            else
            {
                printf("%.4f", M[i][0]);
            }
        }
        else
        {
            printf("%.4f", M[i][0]);
        }

        for(j=1;j<dimc;j++)
        {
            if(M[i][j] < 0)
            {
            /* Matrix's element is negative and its absolute value is smaller than the threshold */
            if(-M[i][j] < threshold)
            {
                printf(",0.0000");
            }
            else
            {
                printf("%s%.4f", ",", M[i][j]); 
            }
            }
            else
            {
                printf("%s%.4f", ",", M[i][j]); 
            }
        }
        printf("\n");
    }
}

/* Indicates # of lines for a given file */
int num_of_lines(FILE *fp){
    /* Variables Decleration */
    int ch;
    int lines = 0;

    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n'){
            lines++;
        }
    }
return (lines+1);   /* Scanning the file until last row - no new empty line  */
}

/* Indicates # of columns for a given file */
int num_of_columns(FILE *fp){
    /* Variables Decleration */
    int ch;
    int columns = 1;

    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == ','){
            columns++;
        }
        if(ch == '\n'){
            break;
        }
    }
    return columns;
}

int main(int argc, char* argv[])
{
    /* Variables Decleration */
    int k, N, d, i, j;
    char* goal;
    int igoal;
    FILE *fp;
    long int bOfFile;
    double** obs;
    double** mat;
    double** current_clusters;
    int* first;
    double n1;
    char c;

    /* Recieving the input */
    if(argc == 4)
    {
        k = atoi(argv[1]);
        goal = argv[2];
        fp = fopen(argv[3],"r");

        if(fp == NULL)   /* Checks that the file was opened properly */
        {
            printf("%s\n","An Error Has Occured");
            return 1;
        }
    }
    else 
    {
        printf("%s\n", "Invalid Input!");
        return 1;
    }

    bOfFile = ftell(fp); /* Saves the address of the beginning of the file */
    N = num_of_lines(fp);  /* Number of lines in file */

    fseek(fp, bOfFile, SEEK_SET); /* Sets the file position back to the beginning */
    if(strcmp(goal,"jacobi") != 0) /* for JACOBI procedeure the input is symmetric matrix (# of rows = # of columns) */
    {
        d = num_of_columns(fp);  /* Number of columns in file */
        fseek(fp, bOfFile, SEEK_SET); /* Sets the file position back to the beginning */
    }
    else
        d = N;   /* In case where goal = jacobi: matrix is always symmetric => dim(columns) = dim(rows) */

    /* Checks if input is valid */
    if(k >= N || d <= 0 || k < 0 || N <= 0)
    {
        printf("%s\n", "Invalid Input!");
        return 1;
    }

    /* Converting goal's string value into an integer (if valid) */
    if(strcmp(goal,"spk") == 0)
    {
        igoal = 0;
    }
    else if(strcmp(goal,"wam") == 0)
    {
        igoal = 1;
    }
    else if(strcmp(goal,"ddg") == 0)
    {
        igoal = 2;
    }
    else if(strcmp(goal,"lnorm") == 0)
    {
        igoal = 3;
    }
    else if(strcmp(goal,"jacobi") == 0)
    {
        igoal = 4;
    }
    else
    {
        printf("%s\n", "Invalid Input!");
        return 1;  
    }

    /* Memory allocation - obseravtions matrix */
    obs = (double**)malloc(N*sizeof(double*));
    if (obs == NULL)
    {
        printf("%s\n","An Error Has Occured");
        return 1;
    }

    for(i = 0; i<N; i++)
    {
        obs[i] = malloc(d*sizeof(double));
        if (obs[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            return 1;
        }
    }
    
    /* Beginning of File Scan */
    for(i=0;i<N-1;i++)
    {
        for(j=0;j<d;j++)
        {
            fscanf(fp, "%lf%c", &n1, &c);
            obs[i][j] = n1;
        }
    }

    for(j=0;j<d-1;j++)  /* Last row scan without the last element */
    {
        fscanf(fp, "%lf%c", &n1, &c);
        obs[N-1][j] = n1;
    }

    fscanf(fp, "%lf", &n1); /* Last element scan */
    obs[N-1][d-1] = n1;
    
    fclose(fp);
    /* End of File Scan */

    /* Returns a matrix defined by goal's value */
    mat = spk_linker(&k,igoal,obs,N,d);

    if(igoal == 4) /* goal == jacobi (prints eiganvalues in first row + eiganvectors' matrix) */
    {
        print_mat(mat,N+1,N); /* Printing output for JACOBI */
    }
    else if(igoal == 0) /* goal == spk (sends to K-means algorithm with first k observations as initial centroids) */
    {
        first = (int*)malloc(k*sizeof(int)); /* Allocating memory for k initial centroids' index */
        
        for(i=0;i<k;i++)
        {
            first[i] = i;
        }
        
        current_clusters = k_means(k,N,k,first,mat); /* Performing K-means algorithm */
        print_mat(current_clusters,k,k); /* Prints k clusters - output for SPK */

        /* Free Memory */
        free(first);

        for(i=0;i<k;i++)
        {
            free(current_clusters[i]);
        }
        free(current_clusters);
    }
    else    /* Printing symmetric matrix for all the other procedures */
    {
        print_mat(mat,N,N);
    }
    
    /* Free Memory */
    if(igoal == 4)  /* Jacobi's output matrix - dim(row) = N+1 */
    {
        for(i=0;i<N+1;i++)
        {
            free(mat[i]);
        }
    }
    else            /* All the other procedures output matrix - dim(row) = N */
    {
        for(i=0;i<N;i++)
        {
            free(mat[i]);
        }
    }
    free(mat);
    
    for(i=0;i<N;i++)
    {
        free(obs[i]);
    }
    free(obs);

    return 0;
}