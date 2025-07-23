#!/bin/bash

g++ -std=c++17  -Ofast -fopenmp -I ../include -I /usr/include/eigen3 compute_billiard_video.cpp -o video.o
parallel --progress "./video.o 3 2 {1} {2} 1" ::: $(seq 0.005 0.005 4) ::: 0 90
video_k.py