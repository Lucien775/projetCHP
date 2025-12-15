#include "mitm.h"

/***************************** global variables ******************************/
u64 n = 0;
u64 mask;
u64 dict_size;
struct entry *A_local;

u32 P[2][2] = {{0,0}, {0xffffffff,0xffffffff}};
u32 C[2][2];

/************************** command-line options ****************************/

void usage(char **argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--n N                       block size [default 24]\n");
        printf("--C0 N                      1st ciphertext (in hex)\n");
        printf("--C1 N                      2nd ciphertext (in hex)\n");
        printf("\n");
        printf("All arguments are required\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[4] = {
                {"n", required_argument, NULL, 'n'},
                {"C0", required_argument, NULL, '0'},
                {"C1", required_argument, NULL, '1'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        int set = 0;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) 
        {
                switch (ch){
                case 'n':
                        n = atoi(optarg);
                        mask = (1ull << n) - 1;
                        break;
                case '0':
                        set |= 1;
                        u64 c0 = strtoull(optarg, NULL, 16);
                        C[0][0] = c0 & 0xffffffff;
                        C[0][1] = c0 >> 32;
                        break;
                case '1':
                        set |= 2;
                        u64 c1 = strtoull(optarg, NULL, 16);
                        C[1][0] = c1 & 0xffffffff;
                        C[1][1] = c1 >> 32;
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        if (n == 0 || set != 3) 
        {
        	usage(argv);
        	exit(1);
        }
}

/******************************************************************************/

int main(int argc, char **argv)
{
    /* MPI init */
    MPI_Init(&argc, &argv);

    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    /* Parse command line (only rank 0 prints) */
    process_command_line_options(argc, argv);
    if (rank == 0) {
        printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n", 
               (int)n, C[0][0], C[0][1], C[1][0], C[1][1]);
    }

    /* search */
    u64 k1[16], k2[16];
    int nkey = golden_claw_search(16, k1, k2);

    /* validation only on rank 0 */
    if (rank == 0) {
        assert(nkey > 0);
        for (int i = 0; i < nkey; i++) {
            assert(f(k1[i]) == g(k2[i]));
            assert(is_good_pair(k1[i], k2[i]));
            printf("Solution found: (%" PRIx64 ", %" PRIx64 ") [checked OK]\n",
                   k1[i], k2[i]);
        }
    }

    /* cleanup */
    MPI_Finalize();
    return 0;
}
