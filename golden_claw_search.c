#include "mitm.h"

/******************************************************************************/
/* search the "golden collision" */
int golden_claw_search(int maxres, u64 k1[], u64 k2[])
{
    MPI_Init(&argc, &argv);
    int rank;
    int p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    double start = wtime();
    u64 N = 1ull << n;
    for (u64 x = 0; x < N; x++) 
    {
        u64 z = f(x);
        dict_insert(z, x);
    }
    double mid = wtime();
    printf("Fill: %.1fs\n", mid - start);
    
    int nres = 0;
    u64 ncandidates = 0;
    u64 x[256];
    for (u64 z = 0; z < N; z++) 
    {
        u64 y = g(z);   
        int nx = dict_probe(y, 256, x);
        assert(nx >= 0);
        ncandidates += nx;
        for (int i = 0; i < nx; i++)
        {

            if (is_good_pair(x[i], z)) 
            {
            	if (nres == maxres)
            		return -1;
            	k1[nres] = x[i];
            	k2[nres] = z;
            	printf("SOLUTION FOUND!\n");
            	nres += 1;
            }
        }
    }
    printf("Probe: %.1fs. %" PRId64 " candidate pairs tested\n", wtime() - mid, ncandidates);
    return nres;
}