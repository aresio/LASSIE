# About the CPU-based LSODA

In order to assess the performances of LASSIE, we compared the method with a CPU-based simulator based on the LSODA algorithm. 
This directory contains the script used for the tests, which relies on the SciPy implementation of LSODA. 

To run a simulation use the following command:

`LSODAsim.py input_dir output_dir`

where `input_dir` is the directory that contains the model to be simulated. 
The output dynamics will be stored in the directory `output_dir`.
