#!/bin/bash

# Compilación de C++ (CPU) con OpenMP y optimización Ofast
g++ -I ../include/ -I /usr/include/eigen3 -fopenmp -Ofast compute_billiard_time.cpp -o compute_billiard_cpp.o
for i in $(seq 40 40 400); do
    ./compute_billiard_cpp.o 3 2 0.408 0 1 $i >> tiempos_cpp.txt
done

# Compilación de CUDA (GPU) con optimización agresiva
nvcc -arch=sm_50 -O3 -Xcompiler -Ofast --use_fast_math -I ../include/ -o compute_billiard_cuda.o compute_billiard_time.cu -lcublas -lcusolver
for i in $(seq 40 40 400); do
    ./compute_billiard_cuda.o 3 2 0.408 0 $i >> tiempos_cuda.txt
done

# Generación de gráficos
python3 GrafTiempos.py con_opt