#include "mitm.h"

/******************************************************************************/
/* search the "golden collision" */
int golden_claw_search(int maxres, u64 k1[], u64 k2[])
{
    /*MPI init*/
    MPI_Init(&argc, &argv);
    int rank;
    int P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    /*On definit le type PAIR_ZX pour l'envoi MPI*/
    MPI_Datatype MPI_PAIR_ZX;
    MPI_Type_contiguous(2, MPI_UNSIGNED_LONG_LONG, &MPI_PAIR_ZX);
    MPI_Type_commit(&MPI_PAIR_ZX);é

    /*On lance le timer*/
    if(rank == 0) double start = wtime();
    u64 N = (1ull << n);

    /**************** Phase Fill *********************/
    /*On répartit les x*/
    u64 base_range = N / P;
    u64 reste = N % P;

    u64 local_range = base_range + (rank < reste ? 1 : 0);
    u64 x_start = base_range * rank  + (rank < reste ? rank : reste);
    u64 x_end = x_start + local_range;

    /***Hash shard pour répartir les z entre les processus***/
    /*On compte cmb d'élément on envoie à chaque processus*/
    int nb_thread;
    #pragma omp parallel
    {
        #pragma omp single
        nb_thread = omp_get_num_threads();
    }

    u64 counts [nb_thread][P];
    memset(counts, 0, sizeof(counts));
    struct z_dest cache[local_range];

    #pragma omp parallel
    {
        int t = omp_get_thread_num();

        #pragma omp for 
        for (u64 x = x_start; x < x_end; ++x)
        {
            u64 z = f(x);
            int dest = murmur64(z) % P;
            cache[x - x_start] = (struct z_dest){z, dest};
            counts[t][dest]++;
        }
    }

    /*Allocation mémoire*/
    struct pair_zx *send_lists[nb_thread][P];
    for (int t = 0; t < nb_thread; ++t)
        for (int i = 0; i < P; ++i)
            send_lists[t][i] = malloc(counts[t][i] * sizeof(pair_zx));

    /*Remplissage*/
    memset(counts, 0, sizeof(counts));
    #pragma omp parallel
    {
        int t = omp_get_thread_num();

        #pragma omp for 
        for (u64 x = x_start; x < x_end; ++x)
        {
            z = cache[x - x_start].z
            dest = cache[x - x_start].dest;
            u64 indice = counts[t][dest]++;
            send_lists[t][dest][indice] = (struct pair_zx){z,x};
        }
    }

    /***Communication***/
    /*Send_counts*/
    u64 send_counts[P];
    memset(send_counts, 0, sizeof(send_counts));
    for (int t = 0; i < nb_thread; ++t)
        for(int dest = 0; dest < P; ++dest)
            send_counts[dest]+= counts[t][dest];

    /*recv_counts*/
    u64 recv_counts[P];
    MPI_Alltoall(send_counts, 1, MPI_UNSIGNED_LONG_LONG,
                recv_counts, 1, MPI_UNSIGNED_LONG_LONG,
                MPI_COMM_WORLD);

    /*recv_buffer*/
    u64 total_recv = 0;
    for (int i = 0; i < P; ++i)
        total_recv += recv_counts[i];
    struct pair_zx *recv_buffer = malloc(total_recv * sizeof(struct pair_zx));

    /*send_offset*/
    u64 send_offset[P];
    send_offset[0] = 0;
    for (int i = 0; i < P; ++i)
        send_offset[i] = send_offset[i - 1] + send_counts[i - 1];

    /*recv_offset*/
    u64 recv_offset[P];
    recv_offset[0] = 0;
    for (int i = 0; i < P; ++i)
        recv_offset[i] = recv_offset[i - 1] + recv_counts[i - 1];

    /*Global buffer*/
    u64 current_offset[P];
    memcpy(current_offset, send_offset, sizeof(current_offset));
    struct pair_zx *send_buffer = malloc(send_offset[P - 1] + send_counts[P - 1] * sizeof(pair_zx));
    for (int t = 0; t < nb_thread; ++t)
        for (int dest = 0; dest < P; ++dest)
            for (u64 k = 0; k < counts[t][dest]; ++k)
                send_buffer[current_offset[dest]++] = send_lists[t][dest][k];

    /*AlltoAllv*/
    MPI_Alltoallv(send_buffer, send_counts, send_offset, MPI_PAIR_ZX,
                recv_buffer, recv_counts, recv_offset, MPI_PAIR_ZX,
                MPI_COMM_WORLD);

    /*Remplissage du dictionnaire*/
    for (u64 i = 0; i < total_recv; ++i)
    {
        u64 z = recv_buffer[i].z;
        u64 x = recv_buffer[i].x;
        dict_insert_sharded(z, x);
    }

    if (rank == 0)
    {
        double mid = wtime();
        printf("Fill: %.1fs\n", mid - start);
    }

    /****************Phase Probe*********************/
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
    MPI_Finalize();
    return nres;
}