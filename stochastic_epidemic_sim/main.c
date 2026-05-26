//SSA for simulating a malaria epidemic
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define ROOT 0

/**
 * Compute propensities for the Malaria model.
 * @param x State vector Should be of length 7!
 * @param w Result vector (propensities). Should be of length 15!
 *
 */
static inline void prop(int *x, double *w) {
	// Birth number, humans
	const double LAMBDA_H = 20;
	// Birth number, mosquitoes
	const double LAMBDA_M = 0.5;
	// Biting rate of mosquitoes
	const double B = 0.075;
	/* Probability that a bite by an infectious mosquito results in transmission
	   of disease to human*/
	const double BETA_H = 0.3;
	/* Probability that a bite results in transmission of parasite to a
	   susceptible mosquito*/
	const double BETA_M = 0.5;
	// Human mortality rate
	const double MU_H = 0.015;
	// Mosquito mortality rate
	const double MU_M = 0.02;
	// Disease induced death rate, humans
	const double DELTA_H = 0.05;
	// Disease induced death rate, mosquitoes
	const double DELTA_M = 0.15;
	// Rate of progression from exposed to infectious state, humans
	const double ALFA_H = 0.6;
	// Rate of progression from exposed to infectious state, mosquitoes
	const double ALFA_M = 0.6;
	// Recovery rate, humans
	const double R = 0.05;
	// Loss of immunity rate, humans
	const double OMEGA = 0.02;
	/* Proportion of an antibody produced by human in response to the incidence
	   of infection caused by mosquito. */
	const double NU_H = 0.5;
	/* Proportion of an antibody produced by mosquito in response to the
	   incidence of infection caused by human. */
	const double NU_M = 0.15;

	w[0] = LAMBDA_H;
	w[1] = MU_H * x[0];
	w[2] = (B * BETA_H * x[0] * x[5]) / (1 + NU_H * x[5]);
	w[3] = LAMBDA_M;
	w[4] = MU_M * x[1];
	w[5] = (B * BETA_M * x[1]*x[4]) / (1 + NU_M * x[4]);
	w[6] = MU_H * x[2];
	w[7] = ALFA_H * x[2];
	w[8] = MU_M * x[3];
	w[9] = ALFA_M * x[3];
	w[10] = (MU_H + DELTA_H) * x[4];
	w[11] = R * x[4];
	w[12] = (MU_M + DELTA_M) * x[5];
	w[13] = OMEGA * x[6];
	w[14] = MU_H * x[6];
}

//To find sum of the propensity vector
static inline double prop_sum(double* x){
    double sum = 0;
    for(int i = 0; i < 15; i++)
    {
        sum += x[i];
    }
    return sum;
}

//To compute the reaction index
static inline int compute_index(double *w, double middle){
    double sum = 0;
    for(int i = 0; i<15; i++)
    {
        sum+= w[i];
        if(sum >= middle)
        {
            return i;
        }
    };
    return -1; //Error case
}

//To print the state vector
void print_array(int *x, int size){
    printf("State vector: ");
    for(int i = 0; i < size; i++)
    {
        printf("%d ", x[i]);
    }
    printf("\n");
}

int write_array_to_file(const char* filename, int* arr, int N) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        perror("Failed to open file for writing");
        return -1;
    }
    for (int i = 0; i < N; i++) {
        fprintf(fp, "%d\n", arr[i]);
    }
    fclose(fp);
    return 0;
}

int main(int argc, char **argv){
    //Initialization of MPI
    int my_rank, size;
    //For timing
    double start_time, total_time;
    double interval_sums[4] = {0, 0, 0, 0}; //For intervals 25s -> 50s -> 75s -> 100s
    double interval_counts[4] = {0, 0, 0, 0};
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n, n_total; //Number of experiments per process, total number of experiments

    //User inputs
    if(my_rank == ROOT)
    {
        if(argc != 2)
        {
            printf("Usage: %s <number of experiments>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        n_total = atoi(argv[1]);
        n = n_total/size;// Number of experiments per process
        //Can assume n_total is divisible by size, hence n is positive
    };
    MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    //printf("I am process %d and I will do %d experiments\n", my_rank, n);

    //Decide chunk sizes
    const int chunks = 10;
    const int chunk_size = n / chunks;

    //Send buffer as before
    //Storage for results - Only x[0] is needed
    int* xpp = malloc(n*sizeof(int)); //Using memory block to store results
    //Initialise non-blocking requests
    MPI_Request send_reqs[chunks]; //One send request per chunk per process

    //Root collection buffer
    //Create memory block to store results at ROOT
    int* xpp_root = NULL;
    MPI_Request recv_reqs[size * chunks]; //Total number of chunks = Total number of receives on root
    if(my_rank == ROOT){
        xpp_root = malloc(n_total*sizeof(int));
        //Post all Irecv calls already ->
        for(int itrecv = 0; itrecv < size; ++itrecv) //For each process
        {
            for(int itc = 0; itc < chunks; ++itc) //For each chunk
            {
                //Create a receive request
                MPI_Irecv(xpp_root + itrecv*n + itc*chunk_size, chunk_size, MPI_INT, itrecv, itc, MPI_COMM_WORLD, &recv_reqs[itrecv*chunks + itc]);
            };
        };
    };

    //SSA data
    int x0[7] = {900,900,30,330,50,270,20}; //State vector initial
    int p[15][7] = {
    {1, 0, 0, 0, 0, 0, 0},
    {-1, 0, 0, 0, 0, 0, 0},
    {-1, 0, 1, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0, 0},
    {0, -1, 0, 0, 0, 0, 0},
    {0, -1, 0, 1, 0, 0, 0},
    {0, 0, -1, 0, 0, 0, 0},
    {0, 0, -1, 0, 1, 0, 0},
    {0, 0, 0, -1, 0, 0, 0},
    {0, 0, 0, -1, 0, 1, 0},
    {0, 0, 0, 0, -1, 0, 0},
    {0, 0, 0, 0, -1, 0, 1},
    {0, 0, 0, 0, 0, -1, 0},
    {1, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, -1}
    }; //Stoichiometry matrix
    int t_end = 100; //Simulation time
    double t; //Current time

    //SSA variables
    double a0, a0_inv, u1, u2, tau;
    int index;
    srand(time(NULL) + my_rank);
    int x[7]; //Allocate memory for state vector
    double w[15]; //Allocate memory for propensity vector

    //Local computation
    start_time = MPI_Wtime(); //Start timing
    
    //Monte-Carlo loop - Now per chunk
    for(int itc = 0; itc < chunks; ++itc) //Per chunk
    {
        for(int i = 0; i < chunk_size; ++i) //Per iteration fill a chunk
            {
                //For storage into buffer
                int itr = itc*chunk_size + i;

                //Initialising SSA
                t = 0;
                memcpy(x, &x0, 7*sizeof(int)); //Copy initial state vector
    
                //Initialising interval measurements
                double last_t = 0;
                double last_wall = MPI_Wtime();
                int interval_index = 0;
                double thresholds[4] = {25, 50, 75, 100};
                double interval_start_wall = last_wall;
        
                //SSA loop
                while(t <= t_end){    
                    //Compute propensities into w 
                    prop(x, w);
                    //print_array(w, 15);
                    
                    //Sum of propensities
                    a0 = prop_sum(w);
                    a0_inv = 1.0/a0;
                    
                    //Generate random numbers
                    u1 = (double)rand() / RAND_MAX;
                    u2 = (double)rand() / RAND_MAX;
                    //printf("u1: %f, u2: %f\n", u1, u2);
                    
                    //Compute time to next reaction
                    tau = -1*log(u1)*a0_inv;
                    //printf("tau: %f\n", tau);
                    
                    //Compute reaction index
                    index = compute_index(w, a0*u2);
                    
                    //Update state vector
                    for(int i=0; i<7; i++){
                        x[i] += p[index][i];
                    };
                    
                    //Update time
                    double prev_t = t;
                    t += tau;
                
                    //For time measurements (25s,50s,75s)
                    while (interval_index < 4 && prev_t < thresholds[interval_index] && t >= thresholds[interval_index]) {
                    double now_wall = MPI_Wtime();
                    interval_sums[interval_index] += now_wall - interval_start_wall;
                    interval_counts[interval_index]++;
                    interval_start_wall = now_wall;
                    interval_index++;}
                }
                
                //For case 100s
                while (interval_index < 4) {
                double now_wall = MPI_Wtime();
                interval_sums[interval_index] += now_wall - interval_start_wall;
                interval_counts[interval_index]++;
                interval_start_wall = now_wall;
                interval_index++;
                }
        
                //Store result into buffer - Per chunk now
                //printf("One experiment done : Current iteration %d\n", itr);
                xpp[itr] = x[0];
            }
    
        //Chunk is complete - Time for non_blocking send
        MPI_Isend(xpp + itc*chunk_size, chunk_size, MPI_INT, ROOT, itc, MPI_COMM_WORLD, &send_reqs[itc]);
    } //Move onto next chunk
    
    //Total time elapsed timing
    total_time = MPI_Wtime() - start_time; //End timing
    double max_time;
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
    if (my_rank == ROOT) {
        printf("Time taken: %f seconds\n", max_time);
    }

    //Find averages of intervals
    double interval_avgs[4];
    for (int i = 0; i < 4; i++) {
        interval_avgs[i] = (interval_counts[i] > 0) ? interval_sums[i] / interval_counts[i] : 0.0;
    }
    
    //Wait for all sends to be done then cleanup and post-processing
    MPI_Waitall(chunks, send_reqs, MPI_STATUS_IGNORE);
    free(xpp);

    if(my_rank == ROOT){
        MPI_Waitall(size*chunks, recv_reqs, MPI_STATUS_IGNORE);
        write_array_to_file("results.txt", xpp_root, n_total);
        free(xpp_root);
    };

    //Interval averages
    double *all_avgs = NULL;
    MPI_Win win;
    //Buffer to store all averages
    if (my_rank == ROOT) {
        all_avgs = malloc(size * 4 * sizeof(double));
    }
    MPI_Win_create(all_avgs, (my_rank == ROOT) ? size * 4 * sizeof(double) : 0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    MPI_Win_fence(0, win);
    MPI_Put(interval_avgs, 4, MPI_DOUBLE, ROOT, my_rank * 4, 4, MPI_DOUBLE, win);
    MPI_Win_fence(0, win);

    if (my_rank == ROOT) {
        printf("Average wall clock times per interval (per process):\n");
        for (int p = 0; p < size; p++) {
            printf("Process %d: ", p);
            for (int i = 0; i < 4; i++) {
                printf("%f ", all_avgs[p*4 + i]);
            }
            printf("\n");
        }
        free(all_avgs);
    }
    MPI_Win_free(&win);
    MPI_Finalize();
    
    return 0;
}