#include "global.h"

/*!
********************************************************************************* 
* Program: CHNS-2D-MPI
*
* This function reads the data from real space using single master processor and 
* scatter it to all other processors.
*
**********************************************************************************/

void read_array_real(double *u, char filename[100])
{
    int count;
    FILE *f;

    f = fopen(filename, "rb");

    count = local_nx * (Ny+2);

    if(rank == master)
    {
        for(i = 0; i < Nx; i++)
        {
          for(j = 0; j < Ny; j++)
          {          
            Idx = IDXR(i, j);
    
            fread(&u1[Idx], sizeof(double), 1, f);
          }
        }
    }//Master......

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(u1, count, MPI_DOUBLE, u, count,
                MPI_DOUBLE, master, MPI_COMM_WORLD);

   fclose(f);

}// End of the function.........
