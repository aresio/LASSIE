#!/bin/bash

nvcc src/c++/simulator.cu -O3 -o lassie --use_fast_math -restrict -lcusolver -lcublas --cudart static \
        -gencode arch=compute_20,code=compute_20 \
        -gencode arch=compute_30,code=compute_30 \
        -gencode arch=compute_35,code=compute_35 \
        -gencode arch=compute_50,code=compute_50
