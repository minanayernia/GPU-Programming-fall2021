#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
#include <pthread.h>
using namespace std;

//
//  tuned constants
//
#define density 0.0005
#define mass 0.01
#define cutoff 0.01
#define min_r (cutoff / 100)
#define dt 0.0005
#define max_threads 4

// calculate particle's bin number
int binNum(particle_t &p, int bpr)
{
  return (floor(p.x / cutoff) + bpr * floor(p.y / cutoff));
}

typedef struct{
  int starting_p;
  int ending_p;
  particle_t *particles;
  int n;
  int bpr;
  vector<particle_t *> *bins;
} thread_args;

void clearBins(vector<particle_t *> *bins, int numbins)
{
  for (int m = 0; m < numbins; m++)
    bins[m].clear();
};

void *clearBinsParallel(void *argss)
{
  thread_args *temp_arg = ((thread_args *)argss);
  for (int m = temp_arg->starting_p; m < temp_arg->ending_p; m++)
    temp_arg->bins[m].clear();

  pthread_exit(NULL);
}

void placeParticlesInBins(int n, vector<particle_t *> *bins, particle_t *particles, int bpr)
{
  for (int i = 0; i < n; i++)
    bins[binNum(particles[i], bpr)].push_back(particles + i);
}

void computeForces(int n, vector<particle_t *> *bins, particle_t *particles, int bpr)
{
  for (int p = 0; p < n; p++)
  {
    particles[p].ax = particles[p].ay = 0;

    // find current particle's bin, handle boundaries
    int cbin = binNum(particles[p], bpr);
    int lowi = -1, highi = 1, lowj = -1, highj = 1;
    if (cbin < bpr)
      lowj = 0;
    if (cbin % bpr == 0)
      lowi = 0;
    if (cbin % bpr == (bpr - 1))
      highi = 0;
    if (cbin >= bpr * (bpr - 1))
      highj = 0;

    // apply nearby forces
    for (int i = lowi; i <= highi; i++)
      for (int j = lowj; j <= highj; j++)
      {
        int nbin = cbin + i + bpr * j;
        for (int k = 0; k < bins[nbin].size(); k++)
          apply_force(particles[p], *bins[nbin][k]);
      }
  }
}

void *computeForcesParallel(void *argss)
{
  thread_args *temp_arg = ((thread_args *)argss);

  for (int p = temp_arg->starting_p; p < temp_arg->ending_p; p++)
  {
    temp_arg->particles[p].ax = temp_arg->particles[p].ay = 0;

    // find current particle's bin, handle boundaries
    int cbin = binNum(temp_arg->particles[p], temp_arg->bpr);
    int lowi = -1, highi = 1, lowj = -1, highj = 1;
    if (cbin < temp_arg->bpr)
      lowj = 0;
    if (cbin % temp_arg->bpr == 0)
      lowi = 0;
    if (cbin % temp_arg->bpr == (temp_arg->bpr - 1))
      highi = 0;
    if (cbin >= temp_arg->bpr * (temp_arg->bpr - 1))
      highj = 0;

    // apply nearby forces
    for (int i = lowi; i <= highi; i++)
      for (int j = lowj; j <= highj; j++)
      {
        int nbin = cbin + i + temp_arg->bpr * j;
        for (int k = 0; k < temp_arg->bins[nbin].size(); k++)
          apply_force(temp_arg->particles[p], *temp_arg->bins[nbin][k]);
      }
  }
  pthread_exit(NULL);
}

void moveParticles(int n, particle_t *particles)
{

  for (int p = 0; p < n; p++)
    move(particles[p]);
}

void *moveParticlesParallel(void *argss)
{
  thread_args *temp_arg = ((thread_args *)argss);

  for (int p = temp_arg->starting_p; p < temp_arg->ending_p; p++)
    move(temp_arg->particles[p]);

  pthread_exit(NULL);
}

//
//  benchmarking program
//
int main(int argc, char **argv)
{
  if (find_option(argc, argv, "-h") >= 0)
  {
    printf("Options:\n");
    printf("-h to see this help\n");
    printf("-n <int> to set the number of particles\n");
    printf("-o <filename> to specify the output file name\n");
    return 0;
  }

  int n = read_int(argc, argv, "-n", 5000);
  int threads_num = read_int(argc, argv, "-t", 8);
  if (threads_num > max_threads)
  {
    printf("no no no");
    return 1;
  }
  char *savename = read_string(argc, argv, "-o", NULL);

  pthread_t threads[threads_num];
  thread_args args[threads_num];
  // thread_args args;

  FILE *fsave = savename ? fopen(savename, "w") : NULL;
  particle_t *particles = (particle_t *)malloc(n * sizeof(particle_t));
  set_size(n);
  init_particles(n, particles);

  // create spatial bins (of size cutoff by cutoff)
  double size = sqrt(density * n);
  int bpr = ceil(size / cutoff);
  int numbins = bpr * bpr;
  vector<particle_t *> *bins = new vector<particle_t *>[numbins];

  //
  //  simulate a number of time steps
  //
  double simulation_time = read_timer();
  for (int step = 0; step < NSTEPS; step++)
  {

    // clear bins at each time step
    clearBins(bins, numbins);

    // int partition_size = numbins / threads_num + 1;

    // for (int i = 0; i < threads_num; i++)
    // {
    //   int s = i * partition_size;
    //   if (s >= numbins)
    //     continue;
    //   args[i].starting_p = s;
    //   args[i].ending_p = fmin((i + 1) * partition_size, numbins);
    //   args[i].bins = bins;

    //   int rc = pthread_create(&threads[i], NULL, clearBinsParallel, (void *)&args[i]);
    // }

    // for (int i=0; i< threads_num; i++){
    //   int64_t status = pthread_join(threads[i], NULL);
    // }

    // place particles in bins
    placeParticlesInBins(n, bins, particles, bpr);

    //
    //  compute forces
    //
    // computeForces(n, bins, particles, bpr);
    //////////////////////////////////////////////////////////st
    int compute_partition_size = n / threads_num + 1;

    for (int i = 0; i < threads_num; i++){
      int s = i * compute_partition_size;
      if (s >= n)
        continue;
      args[i].starting_p = s;
      args[i].ending_p = fmin((i + 1) * compute_partition_size, n);
      args[i].bins = bins;
      args[i].particles = particles;
      args[i].bpr = bpr;

      int rc = pthread_create(&threads[i], NULL, computeForcesParallel, (void *)&args[i]);
    }

    for (int i=0; i< threads_num; i++){
      int64_t status = pthread_join(threads[i], NULL);
    }
    //////////////////////////////////////////////////////////////end force
    //
    //  move particles
    //
    // moveParticles(n, particles);
    /////////////////////////////////////////////////st
    // int partition_size = n / threads_num + 1;

    // for (int i = 0; i < threads_num; i++)
    // {
    //   int s = i * partition_size;
    //   if (s >= n)
    //     continue;
    //   args[i].starting_p = s;
    //   args[i].ending_p = fmin((i + 1) * partition_size, n);
    //   args[i].bins = bins;
    //   args[i].particles = particles;

    //   int rc = pthread_create(&threads[i], NULL, moveParticlesParallel, (void *)&args[i]);
    // }

    // for (int i = 0; i < threads_num; i++)
    // {
    //   int64_t status = pthread_join(threads[i], NULL);
    // }
    ///////////////////////////////////////////////////////end
// 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
    //
    //  save if necessary
    //
    if (fsave && (step % SAVEFREQ) == 0)
      save(fsave, n, particles);
  }
  simulation_time = read_timer() - simulation_time;

  printf("n = %d, simulation time = %g seconds\n", n, simulation_time);

  free(particles);
  delete[] bins;
  if (fsave)
    fclose(fsave);

  return 0;
}