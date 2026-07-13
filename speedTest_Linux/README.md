# IGGtropS Client-side Benchmark

This repository provides a C++ implementation for evaluating the computational efficiency of the client-side IGGtropS algorithm. The program compares different strategies for reconstructing compressed IGGtropS model coefficients and computing the Zenith Total Delay (ZTD) on embedded GNSS terminals.

The benchmark was developed to evaluate the client-side algorithm proposed in our study and compare it with conventional coefficient reconstruction methods.

## Overview

Three computational schemes are implemented:

### Scheme I

- Reads the sparse IGGtropS coefficients for each station.
- Reconstructs the complete set of model coefficients.
- Computes the ZTD.

This represents the straightforward implementation without optimization and is suitable for ordinary standalone terminals.

### Scheme II

- Reads the sparse IGGtropS coefficients for each station.
- Reconstructs only the model coefficients required for the target location.
- Computes the ZTD.

This is the client-side strategy proposed in our study and significantly reduces the computational burden.

### Scheme III

- Reads the sparse coefficients only once.
- Reconstructs the complete IGGtropS coefficient set once.
- Computes the ZTD sequentially for all stations.

This strategy is suitable for large-scale network processing where many stations share the same model coefficients.

## Test Data

The benchmark uses

- 252 IGS station coordinates
- IGGtropS sparse coefficients

All required data files should be placed in

```
IGGtropS_BV/
```

including

```
IGGtropS_B1_DCT.bin
IGGtropS_B2_DCT.bin
IGGtropS_B3_DCT.bin
IGGtropS_B4_DCT.bin
IGGtropS_B5_DCT.bin

IGGtropS_V1_DCT.bin
IGGtropS_V2_DCT.bin
IGGtropS_V3_DCT.bin
IGGtropS_V4_DCT.bin
IGGtropS_V5_DCT.bin

IGS_BLH.csv
```

## Requirements

- Linux
- GCC/G++ (C++11 or later)

## Compilation

Compile using

```bash
g++ -O3 -std=c++11 main.cpp -o main
```

## Run

```bash
./main
```

The program reports the execution time of

- Scheme A (full reconstruction for every station)
- Scheme B (partial reconstruction, proposed method)
- Scheme C (batch reconstruction)

## Hardware

The benchmark in the paper was performed on

- Raspberry Pi 5
- Quad-core ARM Cortex-A76 @ 2.4 GHz
- 8 GB RAM

representing a typical low-cost embedded GNSS terminal.

## Output

Typical output includes

- Loading station coordinates
- Loading DCT dictionary
- Execution time of Scheme A
- Execution time of Scheme B
- Execution time of Scheme C

## Citation
Huang, W., Ou, J., Huo, X., Li W, Yuan, Y., Xiao, G. (2026). A New Sparse Representation Method for Tropospheric Grid Models. Satellite Navigation,

If you use this code in your research, please cite the corresponding publication describing the IGGtropS client-side algorithm.
