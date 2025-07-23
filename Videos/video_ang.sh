#!/bin/bash

g++ -std=c++17 -Ofast -fopenmp -I ../include -I /usr/include/eigen3 compute_billiard_video.cpp -o video.o
parallel --progress "./video.o 3 2 {1} {2} 1" ::: 0.408 0.818 ::: $(0 0.05 90)
video_ang.py