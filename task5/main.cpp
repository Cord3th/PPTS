#include "mpi.h"
#include <stdint.h>
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    if (argc != 4) {
        cout << "Format: A.dat B.dat C.dat" << endl;
        return 0;
    }

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    int a_tag = 1, b_tag = 2, size_tag = 3,
        new_size = pow(size, 1.0001 / 3);

    MPI_Comm cube, square, line_x, line_y, line_z;

    int dims[3] = {new_size, new_size, new_size},
        period[3] = {0, 0, 0}, coords[3],
        remain_dims[3] = {1, 1, 0};

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, false, &cube);
    MPI_Comm_rank(cube, &rank);

    MPI_Cart_coords(cube, rank, 3, coords);
    MPI_Cart_sub(cube, remain_dims, &square);

    remain_dims[1] = 0;
    MPI_Cart_sub(cube, remain_dims, &line_x);

    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(cube, remain_dims, &line_y);

    remain_dims[1] = 0;
    remain_dims[2] = 1;
    MPI_Cart_sub(cube, remain_dims, &line_z);

    int size_x, size_y, size_z;

    MPI_Comm_size(line_x, &size_x);
    MPI_Comm_size(line_y, &size_y);
    MPI_Comm_size(line_z, &size_z);

    MPI_Datatype filetype;

    uint64_t a_m, a_n, b_m, b_n;

    double file_time = 0;
    double temp_time_1, temp_time_2;
    double real_time = 0;
    double real_temp_time_1, real_temp_time_2;

    double *A, *B;
    int part_x, part_y;
    int b_part_x, b_part_y;
    int buf_sizes[2] = {-1, -1};
    int b_buf_sizes[2] = {-1, -1};

    if (coords[2] == 0) {
        temp_time_1 = MPI_Wtime();

        MPI_File file;
        MPI_File_open(square, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

        char type;
        MPI_File_read_all(file, &type, 1, MPI_CHAR, &status);
        MPI_File_read_all(file, &a_m, 1, MPI_UNSIGNED_LONG_LONG, &status);
        MPI_File_read_all(file, &a_n, 1, MPI_UNSIGNED_LONG_LONG, &status);

        temp_time_2 = MPI_Wtime();

        file_time += temp_time_2 - temp_time_1;

        if (coords[0] == 0 && coords[1] == 0) {
            cout << a_m << ' ' << new_size << endl;
        }

        real_temp_time_1 = MPI_Wtime();

        int n = a_n;
        int m = a_m;

        if (coords[0] != size_x - 1) {
            part_x = a_m / size_x;
        } else {
            part_x = a_m / size_x + a_m % size_x;
        }

        if (coords[1] != size_y - 1) {
            part_y = a_n / size_y;
        } else {
            part_y = a_n / size_y + a_n % size_y;
        }

        int start_y = coords[1] * (a_n / size_y);
        MPI_Type_create_subarray(1, &n, &part_y, &start_y, MPI_ORDER_C, MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);

        int offset = 2 * sizeof(uint64_t) + sizeof (char)
                     + sizeof(double) * coords[0] * a_n * (a_m / size_x);

        real_temp_time_2 = MPI_Wtime();
        real_time += real_temp_time_2 - real_temp_time_1;

        A = new double[part_x * part_y];

        MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

        temp_time_1 = MPI_Wtime();
        MPI_File_read_all(file, A, part_x * part_y, MPI_DOUBLE, &status);

        MPI_File_close(&file);
        temp_time_2 = MPI_Wtime();
        file_time += temp_time_2 - temp_time_1;

        real_temp_time_1 = MPI_Wtime();

        int send_to_rank;
        int send_to_coord[3] = {coords[0], coords[1], coords[1]};

        MPI_Cart_rank(cube, send_to_coord, &send_to_rank);

        buf_sizes[0] = part_x;
        buf_sizes[1] = part_y;

        if (coords[1] != 0) {
            MPI_Send(buf_sizes, 2, MPI_INT, send_to_rank, size_tag, cube);
            MPI_Send(A, part_x * part_y, MPI_DOUBLE, send_to_rank, a_tag, cube);
        }

        real_temp_time_2 = MPI_Wtime();
        real_time = real_temp_time_2 - real_temp_time_1;
    }

    real_temp_time_1 = MPI_Wtime();

    if (coords[1] == coords[2]) {
        int send_from;
        int send_from_coords[3] = {coords[0], coords[1], 0};

        MPI_Cart_rank(cube, send_from_coords, &send_from);

        if (coords[2] != 0) {
            MPI_Recv(&buf_sizes, 2, MPI_INT, send_from, size_tag, cube, &status);
            part_x = buf_sizes[0];
            part_y = buf_sizes[1];
        }
    }

    int bcast_root;
    int bcast_coords[1] = {coords[2]};

    MPI_Cart_rank(line_y, bcast_coords, &bcast_root);
    MPI_Bcast(buf_sizes, 2, MPI_INT, bcast_root, line_y);

    part_x = buf_sizes[0];
    part_y = buf_sizes[1];

    if (coords[1] == 0 && coords[2] == 0) {
    } else {
        A = new double[part_x * part_y];
    }

    if (coords[1] == coords[2]) {
        int send_from;
        int send_from_coords[3] = {coords[0], coords[1], 0};

        MPI_Cart_rank(cube, send_from_coords, &send_from);

        if (coords[2] != 0) {
            MPI_Recv(A, part_x * part_y, MPI_DOUBLE, send_from, a_tag, cube, &status);
        }

    }

    MPI_Bcast(A, part_x * part_y, MPI_DOUBLE, bcast_root, line_y);

    //matr A ready

    real_temp_time_2 = MPI_Wtime();
    real_time = real_temp_time_2 - real_temp_time_1;

    if (coords[2] == 0) {
        temp_time_1 = MPI_Wtime();

        MPI_File file;
        char type;

        MPI_File_open(square, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
        MPI_File_read_all(file, &type, 1, MPI_CHAR, &status);
        MPI_File_read_all(file, &b_m, 1, MPI_UNSIGNED_LONG_LONG, &status);
        MPI_File_read_all(file, &b_n, 1, MPI_UNSIGNED_LONG_LONG, &status);

        temp_time_2 = MPI_Wtime();
        file_time += temp_time_2 - temp_time_1;

        real_temp_time_1 = MPI_Wtime();

        int n = b_n;
        int m = b_m;

        if (coords[0] != size_x - 1) {
            b_part_x = b_m / size_x;
        } else {
            b_part_x = b_m / size_x + b_m % size_x;
        }

        if (coords[1] != size_y - 1) {
            b_part_y = b_n / size_y;
        } else {
            b_part_y = b_n / size_y + b_n % size_y;
        }

        int start_y = coords[1] * (b_n / size_y);

        MPI_Type_create_subarray(1, &n, &b_part_y, &start_y, MPI_ORDER_C, MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);

        int offset = 2 * sizeof(uint64_t) + sizeof (char)
                     + sizeof(double) * n * (b_m / size_x) * coords[0];

        real_temp_time_2 = MPI_Wtime();
        real_time = real_temp_time_2 - real_temp_time_1;

        MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

        B = new double[b_part_x * b_part_y];
        temp_time_1 = MPI_Wtime();

        MPI_File_read_all(file, B, b_part_x * b_part_y, MPI_DOUBLE, &status);

        MPI_File_close(&file);

        temp_time_2 = MPI_Wtime();
        file_time = temp_time_2 - temp_time_1;

        int send_to_rank;
        int send_to_coord[3] = {coords[0], coords[1], coords[0]};

        MPI_Cart_rank(cube, send_to_coord, &send_to_rank);

        b_buf_sizes[0] = b_part_x;
        b_buf_sizes[1] = b_part_y;

        if (coords[0] != 0) {
            MPI_Send(b_buf_sizes, 2, MPI_INT, send_to_rank, size_tag, cube);
            MPI_Send(B, b_part_x * b_part_y, MPI_DOUBLE, send_to_rank, b_tag, cube);
        }

    }

    real_temp_time_1 = MPI_Wtime();

    if (coords[0] == coords[2]) {
        int send_from;
        int send_from_coords[3] = {coords[0], coords[1], 0};

        MPI_Cart_rank(cube, send_from_coords, &send_from);

        if (coords[2] != 0) {
            MPI_Recv(&b_buf_sizes, 2, MPI_INT, send_from, size_tag, cube, &status);

            b_part_x = b_buf_sizes[0];
            b_part_y = b_buf_sizes[1];
        }
    }

    MPI_Cart_rank(line_x, bcast_coords, &bcast_root);
    MPI_Bcast(b_buf_sizes, 2, MPI_INT, bcast_root, line_x);

    b_part_x = b_buf_sizes[0];
    b_part_y = b_buf_sizes[1];

    if (coords[0] == 0 && coords[2] == 0) {
    } else {
        B = new double[b_part_x * b_part_y];
    }

    if (coords[0] == coords[2]) {
        int send_from;
        int send_from_coords[3] = {coords[0], coords[1], 0};

        MPI_Cart_rank(cube, send_from_coords, &send_from);
        if(coords[2] != 0) {
            MPI_Recv(B, b_part_x * b_part_y, MPI_DOUBLE, send_from, b_tag, cube, &status);
        }
    }

    MPI_Bcast(B, b_part_x * b_part_y, MPI_DOUBLE, bcast_root, line_x);

    // matr b reade

    double *C_local = new double[part_x * b_part_y];

    if (b_part_x != part_y) {
        cout << "Unable to multiply " << b_part_x << ' ' << part_y
             << ' ' << coords[0] << ' ' << coords[1] << ' ' << coords[2] << endl;
    } else {
        for (int i = 0; i < part_x; i++) {
            for (int j = 0; j < b_part_y; j++) {
                C_local[i * b_part_y + j] = 0;
            }
        }

        for (int i = 0; i < part_x; i++) {
            for (int k = 0; k < b_part_x; k++) {
                for (int j = 0; j < b_part_y; j++) {
                    C_local[i * b_part_y + j] += A[i * part_y + k]
                                                * B[k * b_part_y + j];
                }
            }
        }
    }

    double *C = new double[part_x * b_part_y];

    int C_rank;
    int root_coord[1] = {0};

    MPI_Cart_rank(line_z, root_coord, &C_rank);

    MPI_Reduce(C_local, C, part_x * b_part_y, MPI_DOUBLE, MPI_SUM, C_rank, line_z);

    MPI_Comm_rank(cube, &rank);

    real_temp_time_2 = MPI_Wtime();
    real_time = real_temp_time_2 - real_temp_time_1;

    if (rank == 0) {
        temp_time_1 = MPI_Wtime();

        MPI_File file;

        MPI_File_open(MPI_COMM_SELF, argv[3], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
        char type = 'd';
        MPI_File_write(file, &type, 1, MPI_CHAR, &status);
        MPI_File_write(file, &a_m, 1, MPI_UNSIGNED_LONG_LONG, &status);
        MPI_File_write(file, &b_n, 1, MPI_UNSIGNED_LONG_LONG, &status);
        MPI_File_close(&file);

        temp_time_2 = MPI_Wtime();
        file_time = temp_time_2 - temp_time_1;
    }

    if (coords[2] == 0) {
        MPI_File file;

        MPI_File_open(square, argv[3], MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

        if (coords[0] != size_x - 1) {
            part_x = a_m / size_x;
        } else {
            part_x = a_m / size_x + a_m % size_x;
        }

        if (coords[1] != size_y - 1) {
            part_y = b_n / size_y;
        } else {
            part_y = b_n / size_y + b_n % size_y;
        }

        int n = b_n;
        int start_y = coords[1] * (b_n / size_y);

        MPI_Type_create_subarray(1, &n, &part_y, &start_y, MPI_ORDER_C, MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);

        int offset = 2 * sizeof(uint64_t) + sizeof (char)
                     + sizeof(double) * n * (a_m / size_x) * coords[0];

        MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
        temp_time_1 = MPI_Wtime();

        MPI_File_write(file, C, part_x * b_part_y, MPI_DOUBLE, &status);

        MPI_File_close(&file);
        temp_time_2 = MPI_Wtime();
        file_time = temp_time_2 - temp_time_1;
    }

    double max_file_time;
    double max_real_time;

    MPI_Reduce(&file_time, &max_file_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&real_time, &max_real_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0 ) {
        cout << "File max time = " << max_file_time << endl
             << "Real max time = " << max_real_time << endl;
    }

    delete []C;
    delete []C_local;
    delete []A;
    delete []B;

    MPI_Comm_free(&cube);
    MPI_Comm_free(&square);
    MPI_Comm_free(&line_z);
    MPI_Comm_free(&line_x);
    MPI_Comm_free(&line_y);
    MPI_Finalize();
    return 0;
}
