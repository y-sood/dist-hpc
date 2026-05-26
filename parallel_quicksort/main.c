#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <time.h>

//Comparator for sequential sort
static int comp_float(const void *a, const void *b) {
    float fa = *(const float *)a, fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

//Helpers
static inline int is_power_of_two(int x) { return (x > 0) && !(x & (x-1)); }
static inline int floor_pow2(int x) { int p=1; while((p<<1)<=x) p<<=1; return p; }

//Merges two sorted array parts
static void merge_ascending(const float *a, int na, const float *b, int nb, float *out) {
    int i=0,j=0,k=0;
    while(i<na && j<nb) out[k++] = (a[i]<=b[j]) ? a[i++] : b[j++];
    while(i<na) out[k++]=a[i++];
    while(j<nb) out[k++]=b[j++];
}

//Searching first index after pivot
static int upper_bound(const float *arr, int n, float pivot) {
    int lo=0, hi=n;
    while(lo<hi){ int m=lo+(hi-lo)/2; if(arr[m]<=pivot) lo=m+1; else hi=m; }
    return lo;
}

//Pivot selection
static void select_pivot(float *medians, int gs, int strat, float *out) {
    switch(strat) {
        case 1: *out = medians[0]; break;
        
        case 2: { double s=0; for(int i=0;i<gs;i++) s+=medians[i];
                  *out=(float)(s/gs); break; }
        
        case 3: { float tmp[gs]; memcpy(tmp,medians,gs*sizeof(float));
                  qsort(tmp,gs,sizeof(float),comp_float);
                  *out=(gs%2==1)?tmp[gs/2]:0.5f*(tmp[gs/2-1]+tmp[gs/2]);
                  break; }
        
        default: *out=0.0f;
    }
}

static void _hypersort(float **segs, int *seg_n, int T, int strat) {

    //Base case - 1 thread just return sorted segment
    if (T <= 1) return;

    //Upper/Lower boundary
    int half = T / 2;

    //Global arrays for pivot selection and communication
    float  *medians  = malloc(T * sizeof(float));
    float   pivot    = 0.0f;
    float **recv_buf = malloc(T * sizeof(float *));
    int    *recv_n   = malloc(T * sizeof(int));

    //Parallel section 
    #pragma omp parallel num_threads(T) shared(segs, seg_n, medians, pivot, recv_buf, recv_n, half, strat)
    {
        //Important information for thread
        int tid     = omp_get_thread_num();
        bool is_low = (tid < half);
        int partner = is_low ? tid + half : tid - half;

        //Get my stuff from shared memory
        float *my = segs[tid];
        int my_n  = seg_n[tid];

        //Add my median to global array
        medians[tid] = (my_n > 0) ? my[my_n / 2] : 0.0f;
        //Synchronise
        #pragma omp barrier

        //Master thread of group selects pivot from medians
        #pragma omp single
        { 
            select_pivot(medians, T, strat, &pivot); 
        }

        //Back to parallel section
        //Find first index in my segment where pivot would fit
        int split   = upper_bound(my, my_n, pivot);
        //Find if I send my low or high half to partner, and how many elements I keep/send
        float *keep = is_low ? my          : my + split;
        int  keep_n = is_low ? split        : my_n - split;
        float *send = is_low ? my + split   : my;
        int  send_n = is_low ? my_n - split : split;

        //Write how many elements I send to partner in global array, and synchronise
        recv_n[tid] = send_n;
        #pragma omp barrier
        //Allocate buffer for receiving partner's half, and synchronise
        int p_send = recv_n[partner];
        float *rb  = malloc((p_send > 0 ? p_send : 1) * sizeof(float));
        recv_buf[tid] = rb;
        #pragma omp barrier

        //Transfer
        memcpy(recv_buf[partner], send, send_n * sizeof(float));
        #pragma omp barrier

        //Merge into one sorted segment
        int new_n     = keep_n + p_send;
        float *merged = malloc((new_n > 0 ? new_n : 1) * sizeof(float));
        if (is_low)
            merge_ascending(keep, keep_n, rb, p_send, merged);
        else
            merge_ascending(rb, p_send, keep, keep_n, merged);

        free(my); free(rb);
        //Update global arrays
        segs[tid]  = merged;
        seg_n[tid] = new_n;
    }

    //Cleanup
    free(medians); free(recv_buf); free(recv_n);

    //Recurse on independent halves with half threads each
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        _hypersort(segs,        seg_n,        half, strat);
        #pragma omp section
        _hypersort(segs + half, seg_n + half, half, strat);
    }
}

void global_sort(float *arr, int total_n, int T, int strat) {

    //Base case
    if (T <= 1 || total_n <= 1) {
        qsort(arr, total_n, sizeof(float), comp_float);
        return;
    }

    //Threads are always a power of two, and we can't have more threads than elements
    if (!is_power_of_two(T)) T = floor_pow2(T);
    if (T > total_n)          T = floor_pow2(total_n);
    if (T <= 1) { qsort(arr, total_n, sizeof(float), comp_float); return; }

    //Divide array into roughly equal segments for parallel sorting
    float **segs  = malloc(T * sizeof(float *)); //Array of pointers to store sorted segments
    int    *seg_n = malloc(T * sizeof(int)); //Array to store sizes of segments

    //Divide work among threads and copy segments into separate buffers
    { int base = total_n / T, extra = total_n % T, off = 0;
      for (int i = 0; i < T; i++) {
          int c = base + (i < extra ? 1 : 0);
          segs[i]  = malloc((c > 0 ? c : 1) * sizeof(float));
          memcpy(segs[i], arr + off, c * sizeof(float));
          seg_n[i] = c;
          off += c;
      }
    }

    //Begin parallel section
    #pragma omp parallel num_threads(T)
    { 
      //Sort array segments in parallel  
      int tid = omp_get_thread_num();
      qsort(segs[tid], seg_n[tid], sizeof(float), comp_float); 
    }

      //Call recursive function
      _hypersort(segs, seg_n, T, strat);

      //Merge sorted segments back into original array
      { int off = 0;
        for (int i = 0; i < T; i++) {
            memcpy(arr + off, segs[i], seg_n[i] * sizeof(float));
            free(segs[i]);
            off += seg_n[i];
        }
      }

      //Cleanup
      free(segs); free(seg_n);
}

static void gen_rand(float *a, int n)
{ for(int i=0;i<n;i++) a[i]=((float)rand()/RAND_MAX)*100.0f; }

static void print_arr(const float *a, int n)
{ for(int i=0;i<n;i++) printf("%.2f ",a[i]); printf("\n"); }

int main(int argc, char *argv[]) {
    //Read CLI arguments
    if(argc<5){ 
        fprintf(stderr,"Usage: %s <N> <T> <toggle> <strat>\n",argv[0]); return 1; 
    }
    int N=atoi(argv[1]), T=atoi(argv[2]), toggle=atoi(argv[3]), strat=atoi(argv[4]);
    
    //Generate random array
    float *arr=malloc(N*sizeof(float));
    srand((unsigned)time(NULL));
    gen_rand(arr,N);

    //Pre-sort test => count how many elements are above a certain threshold before sorting, and print array if toggle is on
    const float flag=32.33f;
    if(toggle==1){ 
        printf("Unsorted:\n"); print_arr(arr,N);
        int c=0; for(int i=0;i<N;i++) if(arr[i]>flag) c++;
        printf("Elements > %.2f (pre-sort): %d\n",flag,c); 
    }

    //Maximum threads at one time need to be maximum possible threads
    omp_set_max_active_levels(T);

    //take timings and sort
    double t0=omp_get_wtime();
    global_sort(arr,N,T,strat);
    double t1=omp_get_wtime();

    //Check sorting
    if(toggle==1){ 
        //Print array
        printf("Sorted:\n"); 
        print_arr(arr,N);
        int ok=1;
        //Check sorting
        for(int i=0;i<N-1;i++) {
            if(arr[i]>arr[i+1]) { printf("NOT SORTED at %d\n",i); 
                  ok=0; 
                  break; }
            }
        if(ok) printf("Sort check complete\n");
        int c=0; 
        //Flag check
        for(int i=0;i<N;i++) {
            if(arr[i]>flag) c++;
        }
        printf("Elements > %.2f (post-sort): %d\n",flag,c); 
    }

    //Timing measurement
    printf("Time taken: %f seconds  (N=%d, T=%d, strat=%d)\n", t1-t0, N, T, strat);
    free(arr); return 0;
}