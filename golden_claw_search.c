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

    /*On lance le timer*/
    if(rank == 0) double start = wtime();

    /**************** Phase Fill *********************/
    /* Chaque processus crée son dictionnaire local */
    u64 dict_size_local = 1.125 * ((1ull << n) / P);
    dict_shard_setup(dict_size_local);

    /* On répartit les x entre les processus */
    u64 range_size = (1ull << n) / P;
    u64 x_start = range_size * rank;
    u64 x_end   = x_start + range_size;

    /* On répartit les paires (z,x) entre les processus */
    #pragma omp parallel
    {
        int t = omp_get_thread_num();
        int nb_thread = omp_get_num_threads();

        u64 INITIAL_CAPACITY = range_size / (P * P);

        /* Allocation des listes de chaque thread pour chaque process */
        entry_list send_lists[nb_thread][P];
        for (int th = 0; th < nb_thread; ++th)
        {
            for (int i = 0; i < P; ++i)
            {
                send_lists[th][i].data = malloc(INITIAL_CAPACITY * sizeof(struct pair_zx));
                send_lists[th][i].size = 0;
                send_lists[th][i].capacity = INITIAL_CAPACITY;
            }
        }

        /* Boucle principale de génération */
        #pragma omp for
        for (u64 x = x_start; x < x_end; ++x)
        {
            u64 z = f(x);
            int dest = murmur64(z) % P;

            if (send_lists[t][dest].size == send_lists[t][dest].capacity)
            {
                send_lists[t][dest].capacity *= 2;
                send_lists[t][dest].data = realloc(send_lists[t][dest].data, send_lists[t][dest].capacity * sizeof(struct pair_zx));
            }

            send_lists[t][dest].data[send_lists[t][dest].size++] = (struct pair_zx){z, x};
        }
    }

    /* Fusion des listes de tous les threads pour chaque process */
    entry_list send_lists_final[P];
    for (int i = 0; i < P; ++i)
    {
        u64 total_size = 0;
        for (int t = 0; t < nb_thread; ++t)
        {
            total_size += send_lists[t][i].size;
        }

        send_lists_final[i].data = malloc(total_size * sizeof(struct pair_zx));
        send_lists_final[i].size = 0;
        send_lists_final[i].capacity = total_size;

        /* Copie des données de chaque thread */
        for (int t = 0; t < nb_thread; ++t)
        {
            memcpy(send_lists_final[i].data + send_lists_final[i].size, send_lists[t][i].data, send_lists[t][i].size * sizeof(struct pair_zx));
            send_lists_final[i].size += send_lists[t][i].size;
            free(send_lists[t][i].data); // libération mémoire temporaire
        }
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