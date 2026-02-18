#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<stdbool.h>
#include<float.h>
#include<math.h>

#define MAX_COM 8

struct pInfo {
    bool isSend; // state whether sending or receiving
    int comm_rank; // rank of process to whom we are recieving or sending
    bool isatD1orD2; // true for D1 else D2
};

// ----------------------------------------------------- utilities ---------------------------------------------------
double find_max_data(double* data_received, const int M) {
    double max_data = data_received[0];
    for(size_t i = 0; i < M; i++) {
        max_data = fmax(max_data, data_received[i]);
    }
    return max_data;
}

void prepare_data_for_D1(double* buffer_updated_for_D1, double* data_received, const int M) {
    for(size_t i = 0; i < M; i++) {
        buffer_updated_for_D1[i] = (unsigned long long) data_received[i] % 100000;
    }
}

void prepare_data_for_D2(double* buffer_updated_for_D2, double* data_received, const int M) {
    for(size_t i = 0; i < M; i++) {
        buffer_updated_for_D2[i] = data_received[i] * 100000;
    }
}

void prepare_data_at_D1(double* data_received, const int M) {
    for(size_t i = 0; i < M; i++) {
        data_received[i] = data_received[i] * data_received[i];
    }
}

void prepare_data_at_D2(double* data_received, const int M) {
    // explicitly handling negative log issue (as some values generated are less than eq 1 which will incur problems ahead)
    const double eps = DBL_MIN, logeps = log(eps);
    for(size_t i = 0; i < M; i++) {
        data_received[i] = (data_received[i] <= 0 ? logeps : log(data_received[i]));
    }
}

// -------------------------------------------------- main function ---------------------------------------------------
int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    int rank, size;
    // getting the rank & size of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(argc != 6) {
        if(rank == 0) {
            printf("Insufficient arguments!\nUsage: ./program M D1 D2 T seed\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    const int M = atoi(argv[1]), D1 = atoi(argv[2]), D2 = atoi(argv[3]), T = atoi(argv[4]), seed = atoi(argv[5]);
    srand(seed);

    if(D1 >= D2) {
        if(rank == 0) {
            printf("Invalid D1 and D1!\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // initializing data randomly
    double *data_for_sending_to_D1 = (double*)malloc(M*sizeof(double));
    double *data_for_sending_to_D2 = (double*)malloc(M*sizeof(double));
    for(int i = 0; i < M; i++) {
        data_for_sending_to_D1[i] = (double) rand() * (rank + 1) / 10000.0;
        data_for_sending_to_D2[i] = data_for_sending_to_D1[i];
    }

    // number_of_continuous_ranks_can_send_data
    int continous_ranks = D2 - D1; 
    int jump = D2 + 1;

    struct pInfo** processes_states_flow = (struct pInfo**)malloc(size * sizeof(struct pInfo*));
    for(size_t i = 0; i < size; i++) {
        processes_states_flow[i] = (struct pInfo*)malloc(MAX_COM * sizeof(struct pInfo));
        for(size_t j = 0; j < MAX_COM; j++) {
            processes_states_flow[i][j] = (struct pInfo) {false, -1, false};
        }
    }


    // -------------------------------- construction of processes states flow ---------------------------------------------
    int* free_index = (int*)malloc(size * sizeof(int));
    int* visited = (int*) malloc(size * sizeof(int));
    for(size_t i = 0; i < size; i++) {
        free_index[i] = 0;
        visited[i] = 0;
    }
    
    // processes local variables
    int valid_sender = 0, valid_senders_count = 0;


    for(int p1 = 0; p1 < size; p1+=continous_ranks) {
        for(int p2 = p1, conrun = 0; p2 < size; conrun = (conrun + 1) % continous_ranks) {
            if(visited[p2] == 1) {
                p2 += ((conrun + 1) == continous_ranks ? jump : 1);
                continue;
            }
            visited[p2] = 1; 
            
            if(p2 + D1 < size && p2 + D2 < size) {
                // for sender
                processes_states_flow[p2][free_index[p2]++] = (struct pInfo) {true, p2 + D1, true};
                processes_states_flow[p2][free_index[p2]++] = (struct pInfo) {true, p2 + D2, false};
                processes_states_flow[p2][free_index[p2]++] = (struct pInfo) {false, p2 + D1, true};
                processes_states_flow[p2][free_index[p2]++] = (struct pInfo) {false, p2 + D2, false};

                // for receiver
                processes_states_flow[p2 + D1][free_index[p2 + D1]++] = (struct pInfo) {false, p2, true}; 
                processes_states_flow[p2 + D1][free_index[p2 + D1]++] = (struct pInfo) {true, p2, true}; 
                processes_states_flow[p2 + D2][free_index[p2 + D2]++] = (struct pInfo) {false, p2, false};
                processes_states_flow[p2 + D2][free_index[p2 + D2]++] = (struct pInfo) {true, p2, false};

                valid_sender |= (rank == p2 ? 1 : 0);
                valid_senders_count += 1;
            }
            else if(p2 + D1 < size) {
                processes_states_flow[p2][free_index[p2]++] = (struct pInfo) {true, p2 + D1, true}; 
                processes_states_flow[p2][free_index[p2]++] = (struct pInfo) {false, p2 + D1, true};
                processes_states_flow[p2 + D1][free_index[p2 + D1]++] = (struct pInfo) {false, p2, true};
                processes_states_flow[p2 + D1][free_index[p2 + D1]++] = (struct pInfo) {true, p2, true};

                valid_sender |= (rank == p2 ? 1 : 0);
                valid_senders_count += 1;
            }

            p2 += ((conrun + 1) == continous_ranks ? jump : 1);
        }
    }

    double *data_received = (double*)malloc(M*sizeof(double));
    double max_of_D1 = DBL_MIN, max_of_D2 = DBL_MIN;
    MPI_Status status;
    double stime, etime, time, maxtime;

    // ----------------------------------------------- main control flow logic ---------------------------------------------
    stime = MPI_Wtime();
    for(size_t iter = 1; iter <= T; iter++) {
        for(size_t i = 0; i < MAX_COM; i++) {
            struct pInfo p = processes_states_flow[rank][i];
            if(p.comm_rank == -1) break;

            if(p.isSend == true) {
                if(p.isatD1orD2 == true) { 
                    if(p.comm_rank > rank) { // sending to rank at D1 distance far ahead
                        MPI_Send(data_for_sending_to_D1, M, MPI_DOUBLE, p.comm_rank, p.comm_rank, MPI_COMM_WORLD);
                    }
                    else { // sending to rank at D1 distance far behind
                        prepare_data_at_D1(data_received, M);
                        MPI_Send(data_received, M, MPI_DOUBLE, p.comm_rank, p.comm_rank, MPI_COMM_WORLD);
                    }
                }
                else { 
                    if(p.comm_rank > rank) { // sending to rank at D2 distance far ahead
                        MPI_Send(data_for_sending_to_D2, M, MPI_DOUBLE, p.comm_rank, p.comm_rank, MPI_COMM_WORLD);
                    }
                    else { // sending to rank at D2 distance far behind
                        prepare_data_at_D2(data_received, M);
                        MPI_Send(data_received, M, MPI_DOUBLE, p.comm_rank, p.comm_rank, MPI_COMM_WORLD);
                    }
                }
            } 
            else {
                if(p.isatD1orD2 == true) { 
                    if(p.comm_rank < rank) { // receiving from rank at D1 distance far behind
                        MPI_Recv(data_received, M, MPI_DOUBLE, p.comm_rank, rank, MPI_COMM_WORLD, &status);
                    }
                    else { // receiving from rank at D1 distance far ahead
                        MPI_Recv(data_received, M, MPI_DOUBLE, p.comm_rank, rank, MPI_COMM_WORLD, &status);

                        // updating max_of_D1 for valid sender at the end of Tth iteration and preparing data for next iterations except the last one
                        if(iter == T) max_of_D1 = fmax(max_of_D1, find_max_data(data_received, M));
                        else prepare_data_for_D1(data_for_sending_to_D1, data_received, M);
                    }
                }
                else { 
                    if(p.comm_rank < rank) { // receiving from rank at D2 distance far behind
                        MPI_Recv(data_received, M, MPI_DOUBLE, p.comm_rank, rank, MPI_COMM_WORLD, &status);
                    }
                    else { // receiving from rank at D2 distance far ahead
                        MPI_Recv(data_received, M, MPI_DOUBLE, p.comm_rank, rank, MPI_COMM_WORLD, &status);

                        // updating max_of_D2 for valid sender at the end of Tth iteration and preparing data for next iterations except the last one
                        if(iter == T) max_of_D2 = fmax(max_of_D2, find_max_data(data_received, M));
                        else prepare_data_for_D2(data_for_sending_to_D2, data_received, M);
                    }
                }
            }
        }
    }

    // valid senders except rank 0 sending max of data_recieved_D1 and max of data_recieved_D2 to rank 0
    if(rank != 0 && valid_sender == 1) {
        double sbuf[2] = {max_of_D1, max_of_D2};
        MPI_Send(sbuf, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    // merging logic (valid senders sending computed_data_recieved to rank 0)
    if(rank == 0) {
        int curr_valid_senders_count = 1;
        double recvbuf[2];

        while(curr_valid_senders_count < valid_senders_count) {
            MPI_Recv(recvbuf, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            max_of_D1 = fmax(max_of_D1, recvbuf[0]);
            max_of_D2 = fmax(max_of_D2, recvbuf[1]);
            curr_valid_senders_count += 1;
        }
    }

    etime = MPI_Wtime();
    time = etime - stime;

    // obtaining max time among all ranks
    MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        printf("%lf %lf %lf\n", max_of_D1, max_of_D2, maxtime);
    }

    // ---------------------------------- freeing the memory ----------------------------------
    for(int i = 0; i < size; i++) {
        free(processes_states_flow[i]);
    }
    free(processes_states_flow);
    free(visited);
    free(free_index);
    free(data_for_sending_to_D1);
    free(data_for_sending_to_D2);
    free(data_received);

    MPI_Finalize();
    
    return 0;
}