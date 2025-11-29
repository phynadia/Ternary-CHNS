#include "global.h"

/*!
 ****************************************************************
 * This function writes real-space data.
 *
 * At master processor, all data being gathered from 
 * all the processors using MPI_Gather.
 *
 * This approach of gathering data helps to remove complicacy from
 * creating multiple folders for each processor.
 *
 * The computing time is higher in this approach, but the time for
 * post-processing is less. 
 *
 ****************************************************************/

void write_array_real(double *u, char filename[100])
{
    int count;
    FILE *f;

    f = fopen(filename, "wb");

    count = local_nx * (Ny+2);

    MPI_Gather(u, count, MPI_DOUBLE, u1, count, 
               MPI_DOUBLE, master, MPI_COMM_WORLD);

    if(rank == master)
    {
        for(i = 0; i < Nx; i++)
        {
          for(j = 0; j < Ny; j++)
          {
           
              Idx = IDXR(i, j);
             fwrite(&u1[Idx], sizeof(double), 1, f);
            
          }
        }
    }

   fclose(f);

}//End of the function write_array_real...........................

void dump_real_mpi_IO(double *data, char filename[100])
{
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        Idx = IDXR(i, j);
        int INDEX = j + Ny * i;
          
        work[INDEX] = data[Idx];
        
      }
    }

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_File file;
    MPI_Offset len, offset;
    int local_size;

    local_size = local_nx * Ny;
    len = local_size * sizeof(double);
    offset = rank * len;

    MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);

    MPI_File_write_at_all(file, offset, work, local_size, MPI_DOUBLE,
                          MPI_STATUS_IGNORE);

    MPI_File_close(&file);

}




