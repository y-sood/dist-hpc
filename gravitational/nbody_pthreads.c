#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>

#define EPSILON 1e-3  // Softening parameter
#define MAX_THREADS 32

typedef struct {
    double *x;  // Position x
    double *y;  // Position y
    double *mass;  // Mass
    double *vx;  // Velocity x
    double *vy;  // Velocity y
    double *brightness;  // Brightness (not used in simulation)
} Particle;

typedef struct {
    int thread_id;
    int N;
    double G;
    double delta_t;
    Particle *particles;
    double *vx_temp;
    double *vy_temp;
    pthread_mutex_t *lock;
    int *current_index;
} ThreadData;

void read_particles(const char *filename, int N, Particle *particles) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
        fread(&particles->x[i], sizeof(double), 1, file);
        fread(&particles->y[i], sizeof(double), 1, file);
        fread(&particles->mass[i], sizeof(double), 1, file);
        fread(&particles->vx[i], sizeof(double), 1, file);
        fread(&particles->vy[i], sizeof(double), 1, file);
        fread(&particles->brightness[i], sizeof(double), 1, file);
    }
    fclose(file);
}

void write_particles(const char *filename, int N, Particle *particles) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        fprintf(stderr, "Error opening file for writing: %s\n", filename);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < N; i++) {
        fwrite(&particles->x[i], sizeof(double), 1, file);
        fwrite(&particles->y[i], sizeof(double), 1, file);
        fwrite(&particles->mass[i], sizeof(double), 1, file);
        fwrite(&particles->vx[i], sizeof(double), 1, file);
        fwrite(&particles->vy[i], sizeof(double), 1, file);
        fwrite(&particles->brightness[i], sizeof(double), 1, file);
    }
    fclose(file);
}

void *compute_forces(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    Particle *particles = data->particles;
    int N = data->N;
    double G = data->G;
    double delta_t = data->delta_t;
    double *vx_temp = data->vx_temp;
    double *vy_temp = data->vy_temp;
    pthread_mutex_t *lock = data->lock;
    int *current_index = data->current_index;

    double *vx_local = (double *)calloc(N, sizeof(double));
    double *vy_local = (double *)calloc(N, sizeof(double));

    while (1) {
        pthread_mutex_lock(lock);
        int i = (*current_index)++;
        pthread_mutex_unlock(lock);

        if (i >= N) break;

        double mi_G = particles->mass[i] * G;
        double inv_mi = 1.0 / particles->mass[i];
        double ax_i = 0.0, ay_i = 0.0;
        double x_i = particles->x[i], y_i = particles->y[i];

        for (int j = i + 1; j < N; j++) {
            double inv_mj = 1.0 / particles->mass[j];
            double dx = particles->x[j] - x_i;
            double dy = particles->y[j] - y_i;
            double rij_sq = dx * dx + dy * dy;
            double rij = sqrt(rij_sq);
            double softened_r = rij + EPSILON;
            double softened_r3 = softened_r * softened_r * softened_r;
            double force = mi_G * particles->mass[j] / softened_r3;
            double fx_component = force * dx;
            double fy_component = force * dy;

            ax_i += fx_component * inv_mi;
            ay_i += fy_component * inv_mi;

            vx_local[j] -= delta_t * (fx_component * inv_mj);
            vy_local[j] -= delta_t * (fy_component * inv_mj);
        }

        vx_local[i] += delta_t * ax_i;
        vy_local[i] += delta_t * ay_i;
    }

    pthread_mutex_lock(lock);
    for (int i = 0; i < N; i++) {
        vx_temp[i] += vx_local[i];
        vy_temp[i] += vy_local[i];
    }
    pthread_mutex_unlock(lock);

    free(vx_local);
    free(vy_local);
    return NULL;
}

void computation(Particle *particles, int N, double G, double delta_t, int n_threads) {
    double *vy_temp = (double *)calloc(N, sizeof(double));
    double *vx_temp = (double *)calloc(N, sizeof(double));
    pthread_t threads[n_threads];
    ThreadData thread_data[n_threads];
    pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    int current_index = 0;

    for (int i = 0; i < n_threads; i++) {
        thread_data[i] = (ThreadData){i, N, G, delta_t, particles, vx_temp, vy_temp, &lock, &current_index};
        pthread_create(&threads[i], NULL, compute_forces, &thread_data[i]);
    }

    for (int i = 0; i < n_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    for (int i = 0; i < N; i++) {
        particles->vx[i] += vx_temp[i];
        particles->vy[i] += vy_temp[i];
        particles->x[i] += delta_t * particles->vx[i];
        particles->y[i] += delta_t * particles->vy[i];
    }

    free(vx_temp);
    free(vy_temp);
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s N filename nsteps delta_t graphics n_threads\n", argv[0]);
        return EXIT_FAILURE;
    }

    int N = atoi(argv[1]);
    char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    int n_threads = atoi(argv[6]);
    double G = 100.0 / N;

    Particle particles;
    particles.x = (double *)calloc(N, sizeof(double));
    particles.y = (double *)calloc(N, sizeof(double));
    particles.mass = (double *)calloc(N, sizeof(double));
    particles.vx = (double *)calloc(N, sizeof(double));
    particles.vy = (double *)calloc(N, sizeof(double));
    particles.brightness = (double *)calloc(N, sizeof(double));

    read_particles(filename, N, &particles);

    for (int step = 0; step < nsteps; step++) {
        computation(&particles, N, G, delta_t, n_threads);
    }

    write_particles("result.gal", N, &particles);

    free(particles.x);
    free(particles.y);
    free(particles.mass);
    free(particles.vx);
    free(particles.vy);
    free(particles.brightness);

    return EXIT_SUCCESS;
}
