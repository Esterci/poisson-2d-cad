#include "appctx.h"

void save_array(double *array, AppCtx *app);

void ZeraVector(double *uold, double *unew, double *f, int ndof)
{
  for (int i = 0; i < ndof; i++)
  {
    uold[i] = 0.0;
    unew[i] = 0.0;
    f[i] = 0.0;
  }
}

// define uma macro para acessar os elementos do vetor
#define ind(i, j) (i) * (n + 2) + (j)

int main(int argc, char **argv)
{
  AppCtx app;
  //MPI_Request reqs[4]; // required variable for non-blocking calls
  MPI_Status stats[4]; // required variable for Waitall routine
  int nreqs = 0;
  AppInit(&app, argc, argv);

  int tag1 = 1, tag2 = 2;    // tag para a comunicação
  int n = app.global_n;      // numero global de pontos
  int energy = app.energy;   // energy to be injected per iteration
  int niters = app.niters;   // number of iterations
  int ndof = app.ndof;       // number of degrees of freedom
  int local_n = app.local_n; // numero de pontos em cada processo

  // Para a implementação paralela
  // Cada processo deve ter uma cópia local de uold e unew
  double *uold = (double *)malloc(ndof * sizeof(double)); //
  double *unew = (double *)malloc(ndof * sizeof(double)); //
  double *f = (double *)malloc(ndof * sizeof(double));    //
  double *tmp;

  // Inicializa os vetores
  ZeraVector(uold, unew, f, ndof);

  // Determina a posição dos pontos fonte
  int s1x = n / 4;
  int s1y = s1x;

  int s2x = 3 * n / 4;
  int s2y = s2x;

  // determine UP and DOWN neighbors
  app.neighbors[UP] = app.rank - 1;

  if (app.rank == (app.size - 1))
    app.neighbors[DOWN] = -1;

  else
    app.neighbors[DOWN] = app.rank + 1;

  // Determine os pontos de início e fim de cada processo
  int y_i = app.rank * local_n;
  int y_f = (app.rank + 1) * local_n;

  // Injeta calor nos pontos fonte
  if (s1x >= y_i && s1x < y_f)
  {
    f[ind(s1x + 1 - app.rank * local_n, s1y + 1)] = app.energy;
  };
  if (s2x >= y_i && s2x < y_f)
  {
    f[ind(s2x + 1 - app.rank * local_n, s2y + 1)] = -app.energy;
  }

  //  Variaveris auxiliares
  double h = app.h;
  double h2 = h * h;
  double error = 0.0, total_error = 0.0;

  for (int iter = 0; iter < niters; ++iter)
  {
    nreqs = 0;
    // Cada processo deve trocar as bordas com seus vizinhos
    double local_error = 0.0;

    if (app.neighbors[DOWN] != -1)
    {
      MPI_Sendrecv(&uold[ind(local_n, 1)], n, MPI_DOUBLE, app.neighbors[DOWN], tag1,
                   &uold[ind(local_n + 1, 1)], n, MPI_DOUBLE, app.neighbors[DOWN],
                   tag2, MPI_COMM_WORLD, &stats[nreqs++]);
    }

    if (app.neighbors[UP] != -1)
    {
      MPI_Sendrecv(&uold[ind(1, 1)], n, MPI_DOUBLE, app.neighbors[UP], tag2,
               &uold[ind(0, 1)], n, MPI_DOUBLE, app.neighbors[UP], tag1, 
               MPI_COMM_WORLD, &stats[nreqs++]);
    }

    // MPI_Waitall(nreqs, reqs, stats);

    // Calculo do Stencil
    for (int i = 1; i <= local_n; ++i)
    {
      for (int j = 1; j <= n; ++j)
      {
        unew[ind(i, j)] = 0.25 * (uold[ind(i - 1, j)] + uold[ind(i + 1, j)] +
                                  uold[ind(i, j - 1)] + uold[ind(i, j + 1)] +
                                  h2 * f[ind(i, j)]);
        local_error += (unew[ind(i, j)] - uold[ind(i, j)]) * (unew[ind(i, j)] - uold[ind(i, j)]);
      }
    }

    MPI_Allreduce(&local_error, &total_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    error = sqrt(total_error);

    // if (app.rank == 0)
    //   printf("iter: %d   erro: %8.8e\n", iter, error);

    if (error < app.tol)
      break;

    tmp = unew;
    unew = uold;
    uold = tmp; // swap arrays
  }

  // if (app.rank == 0)
  //   printf("erro: %8.8e\n", total_error);

  // save_array(unew, &app);

  free(uold);
  free(unew);
  free(f);

  AppFinalize(&app);

  return 0;
}