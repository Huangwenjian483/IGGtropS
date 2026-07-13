#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cassert>
#include <chrono> // Replaced <windows.h> for high-precision timing

#define PI           3.1415926535897932  /* pi */
#define NH 26
#define NPARA 5
#define DLON 2.5
#define DLAT 2.5
#define DH 1000.0

const int nLon = (int)(360 / DLON);
const int nLat = (int)(180 / DLAT) + 1;

int mod(int a, int b) {
    int r = a % b;
    return (r < 0) ? r + b : r;
}

void multiple(const int n, const double* Psi, const double* theta, double* X)
{
    int i = 0, j = 0;
    for (i = 0; i < n; i++) {
        X[i] = 0.0;
        for (j = 0; j < n; j++)
            X[i] += Psi[j + i * n] * theta[j];
    }
}

// Linear interpolation function (for altitude)
double interp1(const double* x, const double* y, int size, double xi) {
    int i = 0;
    double t = 0.0;
    if (xi <= x[0]) return y[0];
    if (xi >= x[size - 1]) return y[size - 1];
    for (i = 0; i < size - 1; i++) {
        if (xi >= x[i] && xi <= x[i + 1]) {
            t = (xi - x[i]) / (x[i + 1] - x[i]);
            return y[i] * (1 - t) + y[i + 1] * t;
        }
    }
    return 0.0; // Should not reach here
}

int getIGS(const char* fileInName, const int numCoeffs, double* blh) {
    int row = 0;
    char line[2048];
    char* context = NULL;
    char* token = NULL;

    // Linux standard fopen
    FILE* fpIn = fopen(fileInName, "r");
    if (!fpIn) {
        printf("Failed to open file: %s!\n", fileInName);
        return 1;
    }
    else {
        printf("Successfully opened file: %s!\n", fileInName);
    }

    if (!blh) {
        printf("Memory allocation failed!\n");
        fclose(fpIn);
        return 1;
    }

    while (fgets(line, sizeof(line), fpIn)) {
        // Remove trailing newline characters
        line[strcspn(line, "\r\n")] = 0;

        context = NULL;
        // Use POSIX standard strtok_r instead of strtok_s
        token = strtok_r(line, ",", &context);
        if (!token) continue;

        for (int i = 0; i < 3; ++i) {
            token = strtok_r(NULL, ",", &context);
            if (token == NULL) {
                fprintf(stderr, "Insufficient columns at line %d!\n", row + 1);
                break;
            }
            blh[row * 3 + i] = atof(token);
        }
        row++;
    }

    fclose(fpIn);
    return 0;
}

/**
 * @brief Generates an N-by-N Discrete Cosine Transform (DCT) matrix, stored in a column-major 1D array.
 * @param n Matrix order (Input)
 * @return double* Pointer to the dynamically allocated 1D array (size: n * n)
 */
double* dctmtx(int n) {
    if (n <= 0) {
        return nullptr;
    }

    // Dynamically allocate 1D array memory in C++
    double* c = (double*)std::malloc(n * n * sizeof(double));
    if (c == nullptr) {
        std::cerr << "Error: Memory allocation failed." << std::endl;
        return nullptr;
    }

    // Precompute constant terms using functions under std:: namespace
    double c1 = sqrt(2.0 / n);
    double c2 = 1.0 / sqrt(2.0);

    // Nested loops for computation
    // r is row index (0 ~ n-1), c_idx is column index (0 ~ n-1)
    for (int r = 0; r < n; ++r) {
        for (int c_idx = 0; c_idx < n; ++c_idx) {
            double angle = M_PI * (2.0 * c_idx + 1.0) * r / (2.0 * n);
            
            // 【Column-major indexing】: Each column is arranged continuously, index is c_idx * n + r
            int index = c_idx * n + r;

            if (r == 0) {
                // Special handling for the first row
                c[index] = c1 * cos(angle) * c2;
            } else {
                // Other rows
                c[index] = c1 * cos(angle);
            }
        }
    }

    return c;
}

void IGGtropSERA_raw(const double doy, const double* blh, const double* coeff, double* ZTD) {
    int i = 0, idoy = 0, h = 0, igrid = 0;
    int idx5[5]{}, a, b;
    double lon[nLon]{}, lat[nLat]{}, normH[NH]{};
    double siteLon, siteLat, sitehgt, ix, iy, p, q, t, Y_grid[4] = { 0 };
    double Y_h[NH]{}, logY_h[NH]{};

    for (i = 0; i < nLon; i++) lon[i] = DLON * i;
    for (i = 0; i < nLat; i++) lat[i] = 90.0 - DLAT * i;
    for (i = 0; i < NH; i++) normH[i] = DH * i;

    siteLat = blh[0];
    siteLon = blh[1];
    sitehgt = blh[2];
    ix = (siteLon - lon[0]) / (lon[1] - lon[0]) + 1;
    iy = (siteLat - lat[0]) / (lat[1] - lat[0]) + 1;

    p = ix - (int)ix;
    q = iy - (int)iy;

    t = 2 * PI / 365.25 * (doy);

    for (igrid = 0; igrid < 4; igrid++) {
        a = mod((int)(ix + (igrid & 1)) - 1, nLon);
        b = (int)iy + (igrid / 2) - 1;
        if (b >= nLat) b = nLat - 1;

        for (h = 0; h < NH; h++) {
            for (i = 0; i < 5; i++) {
                idx5[i] = a + b * nLon + h * nLon * nLat + i * nLon * nLat * NH;
            }
            Y_h[h] = coeff[idx5[0]]
                + coeff[idx5[1]] * cos(t)
                + coeff[idx5[2]] * sin(t)
                + coeff[idx5[3]] * cos(2 * t)
                + coeff[idx5[4]] * sin(2 * t);
            logY_h[h] = log(Y_h[h]);
        }
        Y_grid[igrid] = exp(interp1(normH, logY_h, NH, sitehgt));
    }

    ZTD[idoy] = (1 - p) * (1 - q) * Y_grid[0]
        + p * (1 - q) * Y_grid[1]
        + (1 - p) * q * Y_grid[2]
        + p * q * Y_grid[3];
}

void IGGtropSERA_new(const double doy, const double* blh, const double* lat, const double* coeffSpa, const double* Psi, const int* c, double* ZTD) {
    int i = 0, idoy = 0, h = 0, igrid = 0;
    int idx5[5]{}, a, b;
    int num = 20;
    const int nc = 2;
    const int nLat_local = 2;
    double lon[nLon]{}, normH0[NH]{}, normH[nc]{};
    double siteLon, siteLat, sitehgt, ix, iy, p, q, t, Y_grid[4] = { 0 };
    double Y_h[NH]{}, logY_h[NH]{};
    double* coeffRec = nullptr;

    for (i = 0; i < nLon; i++) lon[i] = DLON * i;
    for (i = 0; i < NH; i++) normH0[i] = DH * i;
    for (i = 0; i < nc; i++) normH[i] = normH0[c[i]];

    coeffRec = (double*)malloc(sizeof(double) * (nLon * num));
    assert(coeffRec);
    for (i = 0; i < num; i++) {
        multiple(nLon, Psi, coeffSpa + i * nLon, coeffRec + i * nLon);
    }

    siteLat = blh[0];
    siteLon = blh[1];
    sitehgt = blh[2];
    ix = (siteLon - lon[0]) / (lon[1] - lon[0]) + 1;
    iy = (siteLat - lat[0]) / (lat[1] - lat[0]) + 1;

    p = ix - (int)ix;
    q = iy - (int)iy;

    for (idoy = 0; idoy < 1; idoy++) {
        t = 2 * PI / 365.25 * (doy);

        for (igrid = 0; igrid < 4; igrid++) {
            a = mod((int)(ix + (igrid & 1)) - 1, nLon);
            b = (int)(igrid / 2);
            if (b >= nLat_local) b = nLat_local - 1;

            for (h = 0; h < nc; h++) {
                for (i = 0; i < 5; i++) {
                    idx5[i] = a + b * nLon + h * nLon * nLat_local + i * nLon * nLat_local * nc;
                }
                Y_h[h] = coeffRec[idx5[0]]
                    + coeffRec[idx5[1]] * cos(t)
                    + coeffRec[idx5[2]] * sin(t)
                    + coeffRec[idx5[3]] * cos(2 * t)
                    + coeffRec[idx5[4]] * sin(2 * t);
                logY_h[h] = log(Y_h[h]);
            }

            Y_grid[igrid] = exp(interp1(normH, logY_h, nc, sitehgt));
        }

        ZTD[idoy] = (1 - p) * (1 - q) * Y_grid[0]
            + p * (1 - q) * Y_grid[1]
            + (1 - p) * q * Y_grid[2]
            + p * q * Y_grid[3];
    }
    free(coeffRec);
}

int count_set_bits(const uint8_t* buffer, long startBit, long bitCount) {
    int count = 0;
    for (long i = 0; i < bitCount; i++) {
        long currentBit = startBit + i;
        if ((buffer[currentBit / 8] >> (currentBit % 8)) & 1) {
            count++;
        }
    }
    return count;
}

void extract4SparseColumnsFromFiles(const char* pathB, const char* pathV, int i_mat, int j_mat, double* partial_S) {
    const int n_row = 144;
    const int n_col = 73;
    const long page_size = (long)n_row * n_col;

    int i = i_mat;
    int j = j_mat;

    FILE* fidB = fopen(pathB, "rb");
    if (!fidB) {
        perror("Failed to open B file");
        return;
    }

    long start_idx1 = (long)i * page_size + (long)j * n_row;
    long start_idx2 = (long)(i + 1) * page_size + (long)j * n_row;
    int num_bits = 2 * n_row; 

    long max_bit_idx = start_idx2 + num_bits;
    long bytes_to_read = (max_bit_idx + 7) / 8;

    uint8_t* bufferB = (uint8_t*)malloc(bytes_to_read);
    fread(bufferB, 1, bytes_to_read, fidB);
    fclose(fidB);

    int sum1 = count_set_bits(bufferB, 0, start_idx1);
    int count1 = count_set_bits(bufferB, start_idx1, num_bits);
    int sum2 = sum1 + count1 + count_set_bits(bufferB, start_idx1 + num_bits, start_idx2 - (start_idx1 + num_bits));
    int count2 = count_set_bits(bufferB, start_idx2, num_bits);

    FILE* fidV = fopen(pathV, "rb");
    if (!fidV) {
        free(bufferB);
        perror("Failed to open V file");
        return;
    }

    float* V1 = (float*)malloc(count1 * sizeof(float));
    float* V2 = (float*)malloc(count2 * sizeof(float));

    fseek(fidV, sum1 * sizeof(float), SEEK_SET);
    fread(V1, sizeof(float), count1, fidV);

    fseek(fidV, sum2 * sizeof(float), SEEK_SET);
    fread(V2, sizeof(float), count2, fidV);
    fclose(fidV);

    memset(partial_S, 0, 144 * 4 * sizeof(double));

    int v_ptr = 0;
    for (int k = 0; k < num_bits; k++) {
        if ((bufferB[(start_idx1 + k) / 8] >> ((start_idx1 + k) % 8)) & 1) {
            partial_S[k] = V1[v_ptr++];
        }
    }

    v_ptr = 0;
    for (int k = 0; k < num_bits; k++) {
        if ((bufferB[(start_idx2 + k) / 8] >> ((start_idx2 + k) % 8)) & 1) {
            partial_S[288 + k] = V2[v_ptr++];
        }
    }

    free(bufferB);
    free(V1);
    free(V2);
}

void getData_Segment(const char* pathB, const char* pathV, int m, int n, double* out_ptr) {
    long total_elements = (long)m * n;

    FILE* fidB = fopen(pathB, "rb");
    if (!fidB) return;

    long bytes_to_read = (total_elements + 7) / 8;
    uint8_t* bufferB = (uint8_t*)malloc(bytes_to_read);
    if (!bufferB) {
        fclose(fidB);
        return;
    }
    fread(bufferB, 1, bytes_to_read, fidB);
    fclose(fidB);

    int countV = count_set_bits(bufferB, 0, total_elements);
    float* bufferV = NULL;

    if (countV > 0) {
        FILE* fidV = fopen(pathV, "rb");
        if (!fidV) {
            free(bufferB);
            return;
        }
        bufferV = (float*)malloc(countV * sizeof(float));
        if (bufferV) {
            fread(bufferV, sizeof(float), countV, fidV);
        }
        fclose(fidV);
    }

    if (countV > 0 && bufferV != NULL) {
        int v_idx = 0;
        for (long i = 0; i < total_elements; i++) {
            if ((bufferB[i / 8] >> (i % 8)) & 1) {
                out_ptr[i] = (double)bufferV[v_idx++];
            } else {
                out_ptr[i] = 0.0;
            }
        }
        free(bufferV);
    } else {
        for (long i = 0; i < total_elements; i++) {
            out_ptr[i] = 0.0;
        }
    }
    free(bufferB);
}

void testPlanA(const int doy, const double* blh, const double* Psi, int nsta, double* ZTD) {
    const int nLon = 144;
    const int nLat = 73;
    const int nH = 26;
    const int nPara = 5;
    const int total_n = nLat * nH * nPara;
    const long total_elements = (long)nLon * total_n;
    int i = 0, j = 0;

    // Use Linux-style forward slash paths
    const char* dir = "./IGGtropS_BV";

    double* coeffSpa = (double*)calloc(total_elements, sizeof(double));
    double* coeffRec = (double*)calloc(total_elements, sizeof(double));

    for (i = 0; i < nsta; i++) {
        for (j = 1; j <= 5; j++) {
            char pathB[256], pathV[256];
            // Standard snprintf replaces sprintf_s
            snprintf(pathB, sizeof(pathB), "%s/IGGtropS_B%d_DCT.bin", dir, j);
            snprintf(pathV, sizeof(pathV), "%s/IGGtropS_V%d_DCT.bin", dir, j);

            long offset = (long)(j - 1) * nLon * nLat * nH;
            getData_Segment(pathB, pathV, nLon, nLat * nH, &coeffSpa[offset]);
        }

        for (j = 0; j < total_n; j++) {
            multiple(nLon, Psi, coeffSpa + j * nLon, coeffRec + j * nLon);
        }
        IGGtropSERA_raw(doy, blh + i * 3, coeffRec, ZTD + i);
    }
    free(coeffSpa);
    free(coeffRec);
}

void testPlanB(const int doy, const double* blh, const double* Psi, int nsta, double* ZTD) {
    const int nH = 26;
    const int nPara = 5;
    const double dlat = 2.5;
    const int nLon = 144;
    const int nLat = 73;
    const double dH = 1000.0;
    const int num = 20;
    const char* dir = "./IGGtropS_BV";

    double lat0 = 90.0;
    double normH0 = 0.0;
    double lat_grid[nLat] = {}, lat[2] = {};
    double* coeffSpa = nullptr;
    int i = 0, j = 0;

    for (i = 0; i < nLat; i++) lat_grid[i] = 90.0 - DLAT * i;
    coeffSpa = (double*)malloc(sizeof(double) * (nLon * num));

    for (i = 0; i < nsta; i++) {
        double iy = (blh[i * 3 + 0] - lat0) / (-dlat) + 1.0;
        double iz = (blh[i * 3 + 2] - normH0) / dH + 1.0;

        int c[2];
        c[0] = (int)floor(iz);
        c[1] = (int)ceil(iz);

        while (c[0] < 1) {
            int shift = 1 - c[0];
            c[0] += shift;
            c[1] += shift;
        }

        int b[2];
        b[0] = (int)iy; 
        b[1] = b[0] + 1;
        if (b[1] > nLat) b[1] = nLat;

        for (j = 0; j < 2; j++) {
            b[j] -= 1;
            c[j] -= 1;
            lat[j] = lat_grid[b[j]];
        }

        for (j = 1; j <= 5; j++) {
            char pathB[256], pathV[256];
            snprintf(pathB, sizeof(pathB), "%s/IGGtropS_B%d_DCT.bin", dir, j);
            snprintf(pathV, sizeof(pathV), "%s/IGGtropS_V%d_DCT.bin", dir, j);
            extract4SparseColumnsFromFiles(pathB, pathV, c[0], b[0], coeffSpa + (j - 1) * 4 * nLon);
        }

        IGGtropSERA_new(doy, blh + i * 3, lat, coeffSpa, Psi, c, ZTD + i);
    }
    free(coeffSpa);
}

void testPlanC(const int doy, const double* blh, const double* Psi, int nsta, double* ZTD) {
    const int nLon = 144;
    const int nLat = 73;
    const int nH = 26;
    const int nPara = 5;
    const int total_n = nLat * nH * nPara;
    const long total_elements = (long)nLon * total_n;
    int i = 0, j = 0;

    const char* dir = "./IGGtropS_BV";

    double* coeffSpa = (double*)calloc(total_elements, sizeof(double));
    double* coeffRec = (double*)calloc(total_elements, sizeof(double));

    for (j = 1; j <= 5; j++) {
        char pathB[256], pathV[256];
        snprintf(pathB, sizeof(pathB), "%s/IGGtropS_B%d_DCT.bin", dir, j);
        snprintf(pathV, sizeof(pathV), "%s/IGGtropS_V%d_DCT.bin", dir, j);

        long offset = (long)(j - 1) * nLon * nLat * nH;
        getData_Segment(pathB, pathV, nLon, nLat * nH, &coeffSpa[offset]);
    }

    for (j = 0; j < total_n; j++) {
        multiple(nLon, Psi, coeffSpa + j * nLon, coeffRec + j * nLon);
    }

    for (i = 0; i < nsta; i++) {
        IGGtropSERA_raw(doy, blh + i * 3, coeffRec, ZTD + i);
    }

    free(coeffSpa);
    free(coeffRec);
}

void runningTime() {
    double elapsed_time[3] = {};

    const int numCoeffs = nLon * nLat * NH * NPARA;
    const int strLenth = 256;
    const int numSta = 252;
    const int numDoy = 1;
    int doy = rand() % 365 + 1;

    // Allocated dynamically to prevent stack overflow and maintain low-spec device compatibility
    double* blh = (double*)calloc(numSta * 3, sizeof(double)); 

    double* coeffSpa = nullptr;       
    double* coeffRec = nullptr;       
    double* Dict = nullptr;           
    double* ZTD_raw = nullptr;

    const char* dir = "./IGGtropS_BV";

    char fileIGS[strLenth];  
    snprintf(fileIGS, sizeof(fileIGS), "%s/IGS_BLH.csv", dir);

    coeffSpa = (double*)malloc(sizeof(double) * numCoeffs);
    assert(coeffSpa);

    coeffRec = (double*)malloc(sizeof(double) * numCoeffs);
    assert(coeffRec);
    ZTD_raw = (double*)malloc(sizeof(double) * numSta * numDoy);
    assert(ZTD_raw);

    Dict = dctmtx(nLon);

    printf("\n--------------------------- Reading IGS Station Coordinates ---------------------------\n");
    getIGS(fileIGS, numSta, blh);

    printf("\n--------------------------- Decompressing and Computing ZTD ---------------------------\n");

    doy = 23;

    // C++11 standard high-resolution clock with excellent cross-platform compatibility
    auto start = std::chrono::high_resolution_clock::now();
    testPlanA(doy, blh, Dict, numSta, ZTD_raw);  
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    elapsed_time[0] = diff.count();
    printf("Plan A Running Time: %.6f seconds\n", elapsed_time[0]);

    start = std::chrono::high_resolution_clock::now();
    testPlanB(doy, blh, Dict, numSta, ZTD_raw); 
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    elapsed_time[1] = diff.count();
    printf("Plan B Running Time: %.6f seconds\n", elapsed_time[1]);

    start = std::chrono::high_resolution_clock::now();
    testPlanC(doy, blh, Dict, numSta, ZTD_raw); 
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    elapsed_time[2] = diff.count();
    printf("Plan C Running Time: %.6f seconds\n", elapsed_time[2]);

    // Free memory
    free(blh);
    free(coeffSpa);
    free(Dict);
    free(coeffRec);
    free(ZTD_raw);
}

int main()
{
    // Initialize random seed
    srand((unsigned int)time(NULL));
    runningTime();
    return 0;
}

// g++ -o3 -std=c++11 main.cpp -o main
// ./main