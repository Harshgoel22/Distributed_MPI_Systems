#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<string.h>

// These variables dont change through out the program and hence are declared as globals

int nx, ny, nz;             // sub domain size
int nx_tot, ny_tot, nz_tot; // sub domain + halo region(to store received data from neighbour)
int halo;
int myrank;
int rx,ry,rz;   // 3d coordinates of the process in the process grid
int px,py,pz;      

#define IDX(x,y,z) ((x)*(ny_tot*nz_tot) + (y)*(nz_tot) + (z))



// Helper Functions

void create_data(double **data,int F, int seed)
{
    srand(seed);
    int n = nx_tot * ny_tot * nz_tot;

    for (int i = 0; i < F; i++) {
        data[i] = (double*) malloc(n * sizeof(double));
        double *ptr = data[i];

        for (int x = 0; x < nx_tot; x++) {
            for (int y = 0; y < ny_tot; y++) {
                for (int z = 0; z < nz_tot; z++, ptr++) {
                    // Real data starts at halo and ends at halo+nx-1
                    // Left side halo region starts at 0 and ends at halo-1
                    // Right side halo region starts at halo+nx and ends at nx+2*halo-1

                    if (x < halo || x >= nx + halo || y < halo || y >= ny + halo ||     
                        z < halo || z >= nz + halo) { 
                        *ptr = 0.0;     //point in the halo region
                    }
                    else {
                        int j = (x - halo) * ny * nz + (y - halo) * nz + (z - halo);
                        *ptr = (double)rand() * (myrank + 1) / (110426.0 + i + j);  // real data point
                    }

                }
            }
        }
    }
}

/*
    Not all data points have 6 neighbours(assuming d=7).
    We create a weight mask which determines if a neighbour data point is valid or not.
    Weight=1 means neighbour is valid, 0 means not valid. 
    @see compute_d_stencil to see how this mask is used.
*/
double* create_weight_mask() {
    int n = nx_tot * ny_tot * nz_tot;
    double* weights = (double*) malloc(n * sizeof(double));

    for(int i = 0; i < n; i++) weights[i] = 0.0;

    // Only set to 1.0 if the cell is part of our "Real" domain 
    // OR if it is a halo that will be filled by a neighbor.
    for (int x = 0; x < nx_tot; x++) {
        for (int y = 0; y < ny_tot; y++) {
            for (int z = 0; z < nz_tot; z++) {
                
                double valid = 1;
                if (rx == 0      && x < halo)      valid = 0;   // first process in the X dimension(leftmost)
                if (rx == px - 1 && x >= nx+halo)  valid = 0;   // last process in the X dimension(rightmost)
                if (ry == 0      && y < halo)      valid = 0;   // first process in the Y dimension(bottom)
                if (ry == py - 1 && y >= ny+halo)  valid = 0;   // last process in the Y dimension(top)
                if (rz == 0      && z < halo)      valid = 0;   // first process in the Z dimension(front)
                if (rz == pz - 1 && z >= nz+halo)  valid = 0;   // last process in the Z dimension(back)

                weights[IDX(x, y, z)] = valid;
            }
        }
    }
    return weights;
}
/* 
    For computing average in stencil computation we need to know the number of valid neighbours.
    Number of valid neighbours=sum of neighbours for which weight is 1.
    @see compute_d_stencil to see how divisors is used.
*/
int* compute_divisors_from_weights(double* weights, int d) {
    int n = nx_tot * ny_tot * nz_tot;
    int* divisors = (int*) malloc(n * sizeof(int));
    int base_nbr = (d % 6 != 0) ? 1 : 0;

    for (int x = halo; x < nx + halo; x++) {
        for (int y = halo; y < ny + halo; y++) {
            for (int z = halo; z < nz + halo; z++) {
                
                double count = (double)base_nbr; // Start with center if required

                for (int h = 1; h <= halo; h++) {
                    count += weights[IDX(x - h, y, z)];
                    count += weights[IDX(x + h, y, z)];
                    count += weights[IDX(x, y - h, z)];
                    count += weights[IDX(x, y + h, z)];
                    count += weights[IDX(x, y, z - h)];
                    count += weights[IDX(x, y, z + h)];
                }
                divisors[IDX(x, y, z)] = (int)count;
            }
        }
    }
    return divisors;
}


void create_type_vectors(
    MPI_Datatype *sendx_face1, MPI_Datatype *sendy_face1, MPI_Datatype *sendz_face1, MPI_Datatype *recvx_face1, 
    MPI_Datatype *recvy_face1, MPI_Datatype *recvz_face1, MPI_Datatype *sendx_face2, MPI_Datatype *sendy_face2, 
    MPI_Datatype *sendz_face2, MPI_Datatype *recvx_face2, MPI_Datatype *recvy_face2, MPI_Datatype *recvz_face2)
{
    
    // Creating Type for xface
    
    int xcount = halo * ny;         // Total number of Z-axis "strips" needed to cover a single face is ny. For "halo" number of faces -> halo*ny
    int xblockLengths[xcount];      // Number of contiguous elements in each strip i.e. nz
    int xdisplacements1[xcount];    // Starting memory offsets for the internal boundary (the data that we own and send)
    int xdisplacements2[xcount];    // Starting memory offsets for the ghost cells (where we receive the data)


    // LEFT face
    for (int x = halo, idx = 0; x < 2 * halo; x++) {    // Real data starts at "halo" index. For halo exchange we need "halo" number of faces i.e. halo to halo+halo-1 
        for (int y = halo; y < halo + ny; y++, idx++) {
            xblockLengths[idx] = nz;    // Each strip is of "nz" length
            xdisplacements1[idx] = x * ny_tot * nz_tot + y * nz_tot + halo;
            xdisplacements2[idx] = (x-halo) * ny_tot * nz_tot + y * nz_tot + halo;
        }
    }
    MPI_Type_indexed(xcount, xblockLengths, xdisplacements1, MPI_DOUBLE, sendx_face1);
    MPI_Type_commit(sendx_face1);
    MPI_Type_indexed(xcount, xblockLengths, xdisplacements2, MPI_DOUBLE, recvx_face1);
    MPI_Type_commit(recvx_face1);

    // RIGHT face
    for (int x = nx, idx = 0; x < nx + halo; x++) {     // Real data ends at nx index. For halo exchange we need "halo" number of faces i.e. nx to nx+halo-1
        for (int y = halo; y < halo + ny; y++, idx++) {
            xblockLengths[idx] = nz;
            xdisplacements1[idx] = x * ny_tot * nz_tot + y * nz_tot + halo;
            xdisplacements2[idx] = (x+halo) * ny_tot * nz_tot + y * nz_tot + halo;
        }
    }
    MPI_Type_indexed(xcount, xblockLengths, xdisplacements1, MPI_DOUBLE, sendx_face2);
    MPI_Type_commit(sendx_face2);
    MPI_Type_indexed(xcount, xblockLengths, xdisplacements2, MPI_DOUBLE, recvx_face2);
    MPI_Type_commit(recvx_face2);

    /*
        Logic used for X faces in the above code applies to Y and Z faces also. Hence comments are not repeated.
    */

    // Creating Type for yface
    int ycount = halo * nx;
    int yblockLengths[ycount];
    int ydisplacements1[ycount];
    int ydisplacements2[ycount];

    // Bottom face
    for (int y = halo, idx = 0; y < 2 * halo; y++) {
        for (int x = halo; x < halo + nx; x++, idx++) {
            yblockLengths[idx] = nz;
            ydisplacements1[idx] = x * ny_tot * nz_tot + y * nz_tot + halo;
            ydisplacements2[idx] = x * ny_tot * nz_tot + (y - halo) * nz_tot + halo;
        }
    }
    MPI_Type_indexed(ycount, yblockLengths, ydisplacements1, MPI_DOUBLE, sendy_face1);
    MPI_Type_commit(sendy_face1);
    MPI_Type_indexed(ycount, yblockLengths, ydisplacements2, MPI_DOUBLE, recvy_face1);
    MPI_Type_commit(recvy_face1);

    // Top face
    for (int y = ny, idx = 0; y < ny + halo; y++) {
        for (int x = halo; x < halo + nx; x++, idx++) {
            yblockLengths[idx] = nz;
            ydisplacements1[idx] = x * ny_tot * nz_tot + y * nz_tot + halo;
            ydisplacements2[idx] = x * ny_tot * nz_tot + (y + halo) * nz_tot + halo;
        }
    }
    MPI_Type_indexed(ycount, yblockLengths, ydisplacements1, MPI_DOUBLE, sendy_face2);
    MPI_Type_commit(sendy_face2);
    MPI_Type_indexed(ycount, yblockLengths, ydisplacements2, MPI_DOUBLE, recvy_face2);
    MPI_Type_commit(recvy_face2);

    // Creating Type for zface
    int zcount = halo * nx * ny;
    int zblockLengths[zcount];
    int zdisplacements1[zcount];
    int zdisplacements2[zcount];

    // FRONT face
    for (int z = halo, idx = 0; z < 2 * halo; z++) {
        for (int x = halo; x < halo + nx; x++) {
            for (int y = halo; y < halo + ny; y++, idx++) {
                zblockLengths[idx] = 1;
                zdisplacements1[idx] = x * ny_tot * nz_tot + y * nz_tot + z;
                zdisplacements2[idx] = x * ny_tot * nz_tot + y * nz_tot + (z - halo);
            }
        }
    }
    MPI_Type_indexed(zcount, zblockLengths, zdisplacements1, MPI_DOUBLE, sendz_face1);
    MPI_Type_commit(sendz_face1);
    MPI_Type_indexed(zcount, zblockLengths, zdisplacements2, MPI_DOUBLE, recvz_face1);
    MPI_Type_commit(recvz_face1);

    // BACK face
    for (int z = nz, idx = 0; z < nz + halo; z++) {
        for (int x = halo; x < halo + nx; x++) {
            for (int y = halo; y < halo + ny; y++, idx++) {
                zblockLengths[idx] = 1;
                zdisplacements1[idx] = x * ny_tot * nz_tot + y * nz_tot + z;
                zdisplacements2[idx] = x * ny_tot * nz_tot + y * nz_tot + (z + halo);
            }
        }
    }
    MPI_Type_indexed(zcount, zblockLengths, zdisplacements1, MPI_DOUBLE, sendz_face2);
    MPI_Type_commit(sendz_face2);
    MPI_Type_indexed(zcount, zblockLengths, zdisplacements2, MPI_DOUBLE, recvz_face2);
    MPI_Type_commit(recvz_face2);
}


MPI_Request* perform_halo_exchange(
    double *data, MPI_Datatype sendx_face1, MPI_Datatype sendy_face1, MPI_Datatype sendz_face1, 
    MPI_Datatype recvx_face1, MPI_Datatype recvy_face1, MPI_Datatype recvz_face1, MPI_Datatype sendx_face2, 
    MPI_Datatype sendy_face2, MPI_Datatype sendz_face2, MPI_Datatype recvx_face2, MPI_Datatype recvy_face2, 
    MPI_Datatype recvz_face2
) {
    

    /* 
        Identify Neighbor Ranks:
            Using MPI_PROC_NULL for boundary processes ensures that the 
            global simulation boundaries do not attempt 
            to communicate with non-existent neighbors.
     */

    int nbr_xm = (rx > 0)      ? myrank - (py * pz) : MPI_PROC_NULL; 
    int nbr_xp = (rx < px - 1) ? myrank + (py * pz) : MPI_PROC_NULL; 
    int nbr_ym = (ry > 0)      ? myrank - pz        : MPI_PROC_NULL; 
    int nbr_yp = (ry < py - 1) ? myrank + pz        : MPI_PROC_NULL; 
    int nbr_zm = (rz > 0)      ? myrank - 1         : MPI_PROC_NULL; 
    int nbr_zp = (rz < pz - 1) ? myrank + 1         : MPI_PROC_NULL; 

    MPI_Request* reqs = malloc(12 * sizeof(MPI_Request));
    int r = 0;

    // x exchanges
    MPI_Irecv(data, 1, recvx_face1, nbr_xm, 0, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Isend(data, 1, sendx_face1, nbr_xm, 1, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Irecv(data, 1, recvx_face2, nbr_xp, 1, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Isend(data, 1, sendx_face2, nbr_xp, 0, MPI_COMM_WORLD, &reqs[r++]);

    // y exchanges
    MPI_Irecv(data, 1, recvy_face1, nbr_ym, 2, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Isend(data, 1, sendy_face1, nbr_ym, 3, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Irecv(data, 1, recvy_face2, nbr_yp, 3, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Isend(data, 1, sendy_face2, nbr_yp, 2, MPI_COMM_WORLD, &reqs[r++]);

    // z exchanges
    MPI_Irecv(data, 1, recvz_face1, nbr_zm, 4, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Isend(data, 1, sendz_face1, nbr_zm, 5, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Irecv(data, 1, recvz_face2, nbr_zp, 5, MPI_COMM_WORLD, &reqs[r++]);
    MPI_Isend(data, 1, sendz_face2, nbr_zp, 4, MPI_COMM_WORLD, &reqs[r++]);

    return reqs;
}


double compute_d_stencil(double *data, int x, int y, int z, int d, double *weights, int *divisors) 
{
    // if there is 6 point stencil then center value will not be included but if it is 7, it will
    // be included, so handling this below:
    double result = (d % 6 != 0) ? data[IDX(x, y, z)] : 0;

    for(int h = 1; h <= halo; h++) 
    {
        result += data[IDX(x - h, y, z)] * weights[IDX(x - h, y, z)];   
        result += data[IDX(x + h, y, z)] * weights[IDX(x + h, y, z)];
        result += data[IDX(x, y - h, z)] * weights[IDX(x, y - h, z)];
        result += data[IDX(x, y + h, z)] * weights[IDX(x, y + h, z)];
        result += data[IDX(x, y, z - h)] * weights[IDX(x, y, z - h)];
        result += data[IDX(x, y, z + h)] * weights[IDX(x, y, z + h)];
    }

    result/=(double)divisors[IDX(x, y, z)];
    return result;
}

/**
    Computes the stencil for the outer shell (boundary regions) of the local domain.
    Points located on the edges will have their stencil computation repeated by multiple loops. 
    We do not explicitly handle these overlaps because the output is overwritten with the same value.
 */

void perform_ngbr_d_stencil(double *buffer, double *data, int d, double *weights, int *divisors)
{
    // -X face (only if neighbor exists)
    for(int x = halo; x < 2 * halo; x++) {
        for(int y = halo; y < ny + halo; y++) {
            for(int z = halo; z < nz + halo; z++) {
                buffer[IDX(x, y, z)] = compute_d_stencil(data, x, y, z, d, weights, divisors);
            }
        }
    }


    // +X face
    for(int x = nx; x < nx + halo; x++) {
        for(int y = halo; y < ny + halo; y++) {
            for(int z = halo; z < nz + halo; z++) {
                buffer[IDX(x, y, z)] = compute_d_stencil(data, x, y, z, d, weights, divisors);
            }
        }
    }

    // -Y face
    for(int x = halo; x < nx + halo; x++) {
        for(int y = halo; y < 2 * halo; y++) {
            for(int z = halo; z < nz + halo; z++) {
                buffer[IDX(x, y, z)] = compute_d_stencil(data, x, y, z, d, weights, divisors);
            }
        }
    }

    // +Y face
    for(int x = halo; x < nx + halo; x++) {
        for(int y = ny; y < ny + halo; y++) {
            for(int z = halo; z < nz + halo; z++) {
                buffer[IDX(x, y, z)] = compute_d_stencil(data, x, y, z, d, weights, divisors);
            }
        }
    }

    // -Z face
    for(int x = halo; x < nx + halo; x++) {
        for(int y = halo; y < ny + halo; y++) {
            for(int z = halo; z < 2 * halo; z++) {
                buffer[IDX(x, y, z)] = compute_d_stencil(data, x, y, z, d, weights, divisors);
            }
        }
    }
    

    // +Z face
    for(int x = halo; x < nx + halo; x++) {
        for(int y = halo; y < ny + halo; y++) {
            for(int z = nz; z < nz + halo; z++) {
                buffer[IDX(x, y, z)] = compute_d_stencil(data, x, y, z, d, weights, divisors);
            }
        }
    }
    
}



int compute_isovalues(double *data, int x, int y, int z, double isovalue)
{
    // Extract scalar values from the 8 corners of the unit cube
    double corners_of_unit_cube[8] = {
        data[IDX(x,y,z)],
        data[IDX(x+1,y,z)],
        data[IDX(x,y+1,z)],
        data[IDX(x+1,y+1,z)],
        data[IDX(x,y,z+1)],
        data[IDX(x+1,y,z+1)],
        data[IDX(x,y+1,z+1)],
        data[IDX(x+1,y+1,z+1)]
    };

    // Find the minimum and maximum scalar values among the 8 corners
    double mn = corners_of_unit_cube[0], mx = corners_of_unit_cube[0];
    for(int i=1;i<8;i++){
        if(corners_of_unit_cube[i] < mn) mn = corners_of_unit_cube[i];
        if(corners_of_unit_cube[i] > mx) mx = corners_of_unit_cube[i];
    }

    // A cell contains the isosurface if the isovalue is bounded by the min and max.
    return (mn <= isovalue && mx >= isovalue);
}


/*
    Computes isovalue intersections on the boundary faces.
    To avoid double counting at process boundaries, we use a "Right-Ordering" ownership rule:
    - A data cell spanning two processes is computed ONLY by the process on the "left" (lower coordinate).
    - Consequently, the rightmost processes in the global grid (such as rx == px-1) skip these 
    boundary loops because there is no adjacent data point to complete the unit cube.
 
    Avoiding intra-process double counting of cells at shared edges:
    1. The +X face loop is primary.
    2. The +Y face loop terminates X-iterations early to avoid the X-Y edge.
    3. The +Z face loop terminates both X and Y iterations early to avoid X-Z and Y-Z edges.
 
 */


void compute_ngbr_isovalues(double *buffer, int *isovalues, int index, double isovalue)
{
    int x, y, z;

    // Handle the +X Slab (includes its own edges)
    if (rx < px - 1) {  // The right most process in X dimension does not have a neighbour, so we skip it. Same logic for Y and Z axis also.
        x = nx + halo - 1;
        for (y = halo; y < ny + halo; y++) {
            for (z = halo; z < nz + halo; z++) {
                isovalues[index] += compute_isovalues(buffer, x, y, z, isovalue);
            }
        }
    }

    // Handle the +Y Slab (excluding what +X already took to avoid double counting)
    if (ry < py - 1) {
        y = ny + halo - 1;
        // x stops at nx + halo - 2 so we don't overlap with the +X slab
        for (x = halo; x < (rx < px - 1 ? nx + halo - 1 : nx + halo); x++) {
            for (z = halo; z < nz + halo; z++) {
                isovalues[index] += compute_isovalues(buffer, x, y, z, isovalue);
            }
        }
    }

    // Handle the +Z Slab (excluding +X and +Y overlaps)
    if (rz < pz - 1) {
        z = nz + halo - 1;
        for (x = halo; x < (rx < px - 1 ? nx + halo - 1 : nx + halo); x++) {
            for (y = halo; y < (ry < py - 1 ? ny + halo - 1 : ny + halo); y++) {
                isovalues[index] += compute_isovalues(buffer, x, y, z, isovalue);
            }
        }
    }
}


void swap(double **data, double **buffer) {
    double *tmp = *data;
    *data = *buffer;
    *buffer = tmp;
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int proc;
    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    if(argc < 13) 
    {
        if(myrank == 0) printf("Insufficient arguments\n");
        MPI_Finalize();
        return 1;
    }

    // The input to the program is d, ppn, px, py, pz, nx, ny, nz, T, random seed, F, isovalue
    int d   =  atoi(argv[1]);
    int ppn =  atoi(argv[2]);
    px  =  atoi(argv[3]);
    py  =  atoi(argv[4]);
    pz  =  atoi(argv[5]);
    nx  =  atoi(argv[6]);
    ny  =  atoi(argv[7]);
    nz  =  atoi(argv[8]);    
    int T   =  atoi(argv[9]);
    int random_seed = atoi(argv[10]);
    int F = atoi(argv[11]);
    double isovalue = atof(argv[12]);


    MPI_Request* requests;
    int reqcount = 12;

    halo = d / 6;
    nx_tot = nx + 2 * halo;     // one halo region at the left and one at the right. Same logic for Y and Z dimension.
    ny_tot = ny + 2 * halo;
    nz_tot = nz + 2 * halo;
    
    rx = myrank / (py * pz);
    ry = (myrank % (py * pz)) / pz;
    rz = (myrank % (py * pz)) % pz;
    

    // data owned by each process
    double **data = malloc(F * sizeof(double*));

    // created this for efficiently performing stencil and finading isovalues
    double **buffer = malloc(F * sizeof(double*));

    // creating data
    create_data(data,F, random_seed);
    // creating weight mask for handling boundary conditions in stencil computation
    double *weights = create_weight_mask();
    // creating divisors for normalizing the stencil result
    int *divisors = compute_divisors_from_weights(weights, d);

    // creating mpi datatypes for exchanging the halo values with ease
    MPI_Datatype sendx_face1, sendy_face1, sendz_face1, recvx_face1, recvy_face1, recvz_face1,
            sendx_face2, sendy_face2, sendz_face2, recvx_face2, recvy_face2, recvz_face2;
    
    create_type_vectors(
        &sendx_face1, &sendy_face1, &sendz_face1, &recvx_face1, &recvy_face1, &recvz_face1,
        &sendx_face2, &sendy_face2, &sendz_face2, &recvx_face2, &recvy_face2, &recvz_face2
    );

    // created 1D isovalues for storing F * T isovalues ie for each field and for T number of timesteps
    int *isovalues = malloc(F * T * sizeof(int));

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    for(int f = 0; f < F; f++) {
        buffer[f] = (double*) malloc(nx_tot * ny_tot * nz_tot * sizeof(double));

        int firstDataExchange = 1;

        // perform halo exchange for d point stencil computation
        requests = perform_halo_exchange(
            data[f], sendx_face1, sendy_face1, sendz_face1, recvx_face1,
            recvy_face1, recvz_face1, sendx_face2, sendy_face2, sendz_face2, recvx_face2,
            recvy_face2, recvz_face2
        );
        
        for(int t = 0; t < T; t++) {
            // performing d point stencil computation for center 3d block cells first which
            // are not relying on halo values for computation

            // real data starts at halo and ends at nx+halo. Since we want to avoid halo regions of our real data, we start at halo+halo and end at nx+halo-halo
            for(int x = 2 * halo; x < nx; x++) {        
                for(int y = 2 * halo; y < ny; y++) {
                    for(int z = 2 * halo; z < nz; z++) {
                        buffer[f][IDX(x, y, z)] = compute_d_stencil(data[f], x, y, z, d, weights, divisors);
                    }
                }
            }

            // Wait for all processes to complete their halo exchange
            // We will need the exchanged data for perform_ngbr_d_stencil() function

            // This wait is required only for the exchange which took place before the T loop(It is required only once)
            if (firstDataExchange) {
                MPI_Waitall(reqcount, requests, MPI_STATUSES_IGNORE);
                free(requests);
                firstDataExchange = 0;
            }
        
            // performing stencil computation for neighbour cells
            perform_ngbr_d_stencil(buffer[f], data[f], d, weights, divisors);

            /*
            At this point, stencil computation is fully complete. 
            The computed values in the halo regions must be communicated to neighbours for 2 purposes:
            1. For the next iteration of stencil computation i.e. t+1 timestep
            2. For computing isovalues for neighbour cells
            */

            // performing halo exchange for updating the neighbour cells
            requests = perform_halo_exchange(
                buffer[f], sendx_face1, sendy_face1, sendz_face1, recvx_face1,
                recvy_face1, recvz_face1, sendx_face2, sendy_face2, sendz_face2, recvx_face2,
                recvy_face2, recvz_face2
            );

            // computing isovalues for cells not depedning on halos
            int index = f * T + t;
            isovalues[index] = 0;
            for(int x = halo; x < nx + halo - 1; x++) {
                for(int y = halo; y < ny + halo - 1; y++) {
                    for(int z = halo; z < nz + halo - 1; z++) {
                        isovalues[index] += compute_isovalues(buffer[f], x, y, z,isovalue);
                    }
                }
            }
            

            // Wait for halo exchanges to be done successfully
            MPI_Waitall(reqcount, requests, MPI_STATUSES_IGNORE);
            free(requests);

            // compute isovalues for neighbour cells requiring halos
            compute_ngbr_isovalues(buffer[f], isovalues, index, isovalue);

            // swapping data and buffer
            swap(&data[f], &buffer[f]);
        }
    }

    // reduce isovalues at root 0 (MPI_SUM)
    int *total_isovalues = malloc(F * T * sizeof(int));
    MPI_Reduce(isovalues, total_isovalues, F * T, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // max time across all processes
    double local_time = MPI_Wtime() - start;
    double max_time;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // root 0 printing the desired output (ie T)
    if(myrank == 0) {
        for(int t = 0; t < T; t++) {
            for(int f = 0; f < F; f++) {
                int idx = f * T + t;
                printf("%d ", total_isovalues[idx]);
            }
            printf("\n");
        }
        printf("%lf\n", max_time);
    }

    // freeing up memory
    for(int f = 0; f < F; f++) {
        free(data[f]);
        free(buffer[f]);
    }
    free(data);
    free(buffer);
    free(isovalues);
    free(total_isovalues);
    MPI_Type_free(&sendx_face1);
    MPI_Type_free(&sendy_face1);
    MPI_Type_free(&sendz_face1);

    MPI_Type_free(&recvx_face1);
    MPI_Type_free(&recvy_face1);
    MPI_Type_free(&recvz_face1);

    MPI_Type_free(&sendx_face2);
    MPI_Type_free(&sendy_face2);
    MPI_Type_free(&sendz_face2);

    MPI_Type_free(&recvx_face2);
    MPI_Type_free(&recvy_face2);
    MPI_Type_free(&recvz_face2);

    MPI_Finalize();
    return 0;
}