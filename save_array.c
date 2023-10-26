
#include "appctx.h"

void save_array_ascii(double *array, AppCtx *app)
{
    char filename[256];
    FILE *fp;

    if (app->rank != 0)
        return;

    sprintf(filename, "matrizes_resultantes/array.txt");
    fp = fopen(filename, "w");
    int n = app->global_n;
    fprintf(fp, "%d %d %8.8e %8.8e\n", n + 2, n + 2, app->L, app->L);
    for (int i = 0; i < (n + 2); i++)
    {
        for (int j = 0; j < (n + 2); j++)
            fprintf(fp, "%8.8e ", array[i * (n + 2) + j]);
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void save_array_binary(double *array, AppCtx *app)
{
    char filename[256];
    FILE *fp;

    if (app->rank != 0)
        return;

    sprintf(filename, "matrizes_resultantes/array.bin");
    fp = fopen(filename, "wb");
    int n = app->global_n + 2;
    double dx = app->h;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&dx, sizeof(double), 1, fp);
    fwrite(&dx, sizeof(double), 1, fp);
    fwrite(array, sizeof(double), n * n, fp);
    fclose(fp);
}

void save_array_vti(double *array, AppCtx *app)
{
    // save a parallel vtk XML image data format file
    int n = app->global_n;
    int rank = app->rank;
    int size = app->size;
    int i, j;
    char filename[256];
    FILE *fp;

    if (rank != 0)
        return;

    int x1 = 0;
    int x2 = n + 1;
    int y1 = 0;
    int y2 = n + 1;

#ifdef HAVE_MPI
    sprintf(filename, "matrizes_resultantes/array_mpi.%d.vti", rank);
#else
#ifdef OMP
    sprintf(filename, "matrizes_resultantes/array_omp.%d.vti", rank);
#else
    sprintf(filename, "matrizes_resultantes/array_serial.%d.vti", rank);
#endif
#endif
    fp = fopen(filename, "w");
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n", x1, x2, y1, y2);
    fprintf(fp, "<Piece Extent=\"%d %d %d %d 0 0\">\n", x1, x2, y1, y2);
    fprintf(fp, "<PointData>\n");
    // fprintf(fp, "<DataArray type=\"Float64\" Name=\"array\" NumberOfComponents=\"1\" format=\"binary\">\n");
    // fwrite(array, sizeof(double),app->ndof, fp);
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"array\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int i = 0; i < (n + 2); i++)
    {
        for (int j = 0; j < (n + 2); j++)
            fprintf(fp, "\n");
    }
    fprintf(fp, "</DataArray>\n");
    fprintf(fp, "</PointData>\n");
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</ImageData>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}

#define ind(i, j) (i) * (global_n + 2) + (j)
void save_array_csv(double *array, AppCtx *app)
{
    // save array in csv file

    int local_n = app->local_n;
    int global_n = app->global_n;
    char filename[256];
    FILE *fp;

#ifdef HAVE_MPI
    sprintf(filename, "matrizes_resultantes/array_mpi_%d.csv", app->rank);
#else
#ifdef HAVE_OMP
    sprintf(filename, "matrizes_resultantes/array_omp.csv");
#else
    sprintf(filename, "matrizes_resultantes/array_serial.csv");
#endif
#endif

    fp = fopen(filename, "wb");

    for (int i = 1; i <= local_n; ++i)
    {
        for (int j = 1; j <= global_n; ++j)
        {
            fprintf(fp, "%lf ,", array[ind(i, j)]);
        };
    };

    fclose(fp);
}

void save_array(double *array, AppCtx *app)
{
    switch (app->output_type)
    {
    case OUTPUT_BINARY:
        save_array_binary(array, app);
        break;
    case OUTPUT_ASCII:
        save_array_ascii(array, app);
        break;
    case OUTPUT_VTI:
        save_array_vti(array, app);
        break;
    case OUTPUT_CSV:
        save_array_csv(array, app);
        break;
    default:
        fprintf(stderr, "Tipo de saida invalida `%d'\n", app->output_type);
        AppFinalize(app);
        exit(-1);
    }
}