#include "mitm.h"

/******************************************************************************/
/* search the "golden collision" */
int golden_claw_search(int maxres, u64 k1[], u64 k2[])
{
    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    /* MPI datatype  pour envoi*/
    MPI_Datatype MPI_PAIR_ZX;
    MPI_Type_contiguous(2, MPI_UNSIGNED_LONG_LONG, &MPI_PAIR_ZX);
    MPI_Type_commit(&MPI_PAIR_ZX);


    /*On lance le timer*/
    double start = 0.0;
    double mid = 0.0;
    if(rank == 0) start = wtime();
    u64 N = (1ull << n);

    /**************** Phase Fill *********************/
    /*On répartit les x*/
    u64 base_range = N / P;
    u64 reste = N % P;

    u64 local_range = base_range + (rank < reste ? 1 : 0);
    u64 x_start = base_range * rank  + (rank < reste ? rank : reste);
    u64 x_end = x_start + local_range;
    
    /*Hash shard*/
    int send_counts[P];
    memset(send_counts, 0, sizeof(send_counts));
    struct z_dest cache[local_range];
    for (u64 x = x_start; x < x_end; ++x)
    {
        u64 z = f(x);
        int dest = murmur64(z) % P;
        cache[x - x_start] = (struct z_dest){z, dest};
        send_counts[dest]++;
    }

    /*Allocation mémoire*/
    struct pair_zx *send_list[P];
    for (int i = 0; i < P; ++i)
        send_list[i] = malloc(send_counts[i] * sizeof(struct pair_zx));

    /*Remplissage*/
    memset(send_counts, 0, sizeof(send_counts));
    for (u64 x = x_start; x < x_end; ++x)
    {
        u64 z = cache[x - x_start].z;
        int dest = cache[x - x_start].dest;
        u64 indice = send_counts[dest]++;
        send_list[dest][indice] = (struct pair_zx){z,x};
    }

    /***Communication***/
    /*recv counts*/
    int recv_counts[P];
    MPI_Alltoall(send_counts, 1, MPI_INT,
                recv_counts, 1, MPI_INT,
                MPI_COMM_WORLD);

    /*send_displs & recv_displs & total_recv*/
    int send_displs[P];
    send_displs[0] = 0;
    int recv_displs[P];
    recv_displs[0] = 0;
    int total_recv = recv_counts[0];
    for (int i = 1; i < P; ++i)
    {
        send_displs[i] = send_displs[i - 1] + send_counts[i - 1];
        recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
        total_recv += recv_counts[i];
    }
        
    /*recv_buffer*/
    struct pair_zx *recv_buffer = malloc(total_recv * sizeof(struct pair_zx));

    /*send_buffer*/
    int current_offset[P];
    memcpy(current_offset, send_displs, sizeof(current_offset));
    u64 total_send = send_displs[P - 1] + send_counts[P - 1];
    struct pair_zx *send_buffer = malloc(total_send * sizeof(struct pair_zx));
    for (int dest = 0; dest < P; ++dest)
        for (u64 k = 0; k < send_counts[dest]; ++k)
            send_buffer[current_offset[dest]++] = send_list[dest][k];

    /*AlltoAllv*/
    MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_PAIR_ZX,
                recv_buffer, recv_counts, recv_displs, MPI_PAIR_ZX,
                MPI_COMM_WORLD);

    /*Préparation du dictionnaire*/
    dict_shard_setup(2 * total_recv + 1);

    /*Remplissage du dictionnaire*/
    for (u64 i = 0; i < total_recv; ++i)
    {
        u64 z = recv_buffer[i].z;
        u64 x = recv_buffer[i].x;
        dict_insert_sharded(z, x);
    }

    /*On attend tous le monde*/
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        mid = wtime();
        printf("Fill: %.1fs\n", mid - start);
    }

    /*Nettoyage mémoire*/
    for (int i = 0; i < P; ++i)
        free(send_list[i]);
    free(send_buffer);
    free(recv_buffer);

    /****************Phase Probe*********************/
    /*On cherche les solutions*/
    int nres_local = 0;
    u64 ncandidates_local = 0;
    u64 x[256];
    u64 k1_local[maxres];
    u64 k2_local[maxres];

    for (u64 z = 0; z < N; ++z) 
    {
        u64 y = g(z);   
        int nx = dict_probe(y, 256, x);
        assert(nx >= 0);
        ncandidates_local += nx;
        
        for (int i = 0; i < nx; i++)
        {
            if (is_good_pair(x[i], z)) 
            {
                if (nres_local == maxres)
                    break;
                k1_local[nres_local] = x[i];
                k2_local[nres_local] = z;
                printf("SOLUTION FOUND!\n");
                nres_local += 1;
            }
        }
    }

    /***On envoie les solutions au rank 0***/
    int *sol_counts = NULL;
    int *displs = NULL;

    if (rank == 0)
    {
        sol_counts = malloc(P * sizeof(int));
        displs = malloc(P * sizeof(int));
    }

    MPI_Gather(&nres_local, 1, MPI_INT, sol_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int total = 0;
    if (rank == 0)
    {
        for (int i = 0; i < P; i++)
        {
            displs[i] = total;
            total += sol_counts[i];
        }
    }

    u64 *k1_all = NULL;
    u64 *k2_all = NULL;
    if (rank == 0)
    {
        k1_all = malloc(total * sizeof(u64));
        k2_all = malloc(total * sizeof(u64));
    }

    MPI_Gatherv(k1_local, nres_local, MPI_UINT64_T, k1_all, sol_counts, displs, MPI_UINT64_T,
                0, MPI_COMM_WORLD);

    MPI_Gatherv(k2_local, nres_local, MPI_UINT64_T, k2_all, sol_counts, displs, MPI_UINT64_T,
                0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        memcpy(k1, k1_all, total * sizeof(u64));
        memcpy(k2, k2_all, total * sizeof(u64));
    }

    /*Stats*/
    u64 ncandidates = 0;
    MPI_Reduce(&ncandidates_local, &ncandidates, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        printf("Probe: %.1fs. %" PRId64 " candidate pairs tested\n", wtime() - mid, ncandidates);
    }
    
    /*Nettoyage mémoire*/
    if (rank == 0)
    {
        free(k1_all);
        free(k2_all);
        free(sol_counts);
        free(displs);
    }

    MPI_Type_free(&MPI_PAIR_ZX);
    return total;
}
   
