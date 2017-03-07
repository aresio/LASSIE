#include "Reader.h"
#include "kernels.cuh"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>     // string, stod
#include <iostream>   // cout
#include <time.h>
#include <cmath>
#include <math.h>       /* ceil */
#include <limits>
#include <cstring>
#include <thrust/transform_reduce.h>
#include <thrust/count.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/extrema.h>
#include "cublas_v2.h"
#include "cusolverDn.h"


#include <typeinfo>

using namespace std;

#ifdef _WIN32
	const string slash = "\\";
#elif _WIN64
	const string slash = "\\";
#elif __linux__
	const string slash = "/";
#else
	const string slash = "/";
#endif

#define CUDA_CHECK_RETURN(value) {											\
		cudaError_t _m_cudaStat = value;										\
		if (_m_cudaStat != cudaSuccess) {										\
			fprintf(stderr, "Error %s at line %d in file %s\n",					\
					cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
					exit(1);															\
		} }

/******************************************************************************/
/*  square  */

template<typename T>
class Simulation
{
public:
	Simulation()
{

}

	struct square
	{
		__host__ __device__ float operator()(const T & x) const
		{
			return x * x;
		}
	};

	struct checkEI
	{
		__host__ __device__ int operator()(const T & x) const
		{
			if(isnan(x) || isinf(x) || x < 0)
				return 1;
			else
				return 0;
		}
	};

	float run(string A, string B, string c_vector1, string cs_vector1,
			string t_vector1, string MX_0, string M_feed1, string modelkind1,
			string folder, int verbose, string atol_vector1, string be_step,
			string newton_iter, string newton_tol, string rkf_step, 
			string stiffness_tol, string volume1)
	{

		CUDA_CHECK_RETURN(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

		/* ******************** Reading input model files to be simulated ******************** */
		cusolverDnHandle_t handle = NULL;
		cublasHandle_t cublasHandle = NULL; 
		cudaStream_t stream = NULL;
		cusolverDnCreate(&handle);
		cublasCreate(&cublasHandle);
		cudaStreamCreate(&stream);

		cusolverDnSetStream(handle, stream);
		cublasSetStream(cublasHandle, stream);

		Reader<T>* reader = new Reader<T>;
		vector<vector<short> > left_matrix = reader -> readMatrix(A);
		vector<vector<short> > right_matrix = reader -> readMatrix(B);
		
		vector<T> c_vector = reader -> readByLine(c_vector1);
		vector<T> t_vector = reader -> readByLine(t_vector1);

		vector<T> M_0 = reader -> readByTab(MX_0);

		//vector<int> cs_vector = reader -> readByLineshort(cs_vector1);
		//vector<T> M_feed = reader -> readByTab(M_feed1);
		//int modelkind = reader -> readType(modelkind1);

		T dt_ei, dt_rkf, tollN, isStiff;
		vector<T> atol_vector;
		vector<T> M_feed;
		vector<int> cs_vector;
		int modelkind;

		int iterN;

		//reading optional files
		if(strcmp(be_step.c_str(), "NA") != 0)
		{
			dt_ei = reader -> readSingleValue(be_step);
		}
		else
		{
			dt_ei = 0.1;
		}

		if(strcmp(rkf_step.c_str(), "NA") != 0)
		{
			dt_rkf = reader -> readSingleValue(rkf_step);
		}
		else
		{
			dt_rkf = 1e-3;
		}

		if(strcmp(newton_iter.c_str(), "NA") != 0)
		{
			iterN = (int) reader -> readSingleValue(newton_iter);
		}
		else
		{
			iterN = 10000;
		}

		if(strcmp(newton_tol.c_str(), "NA") != 0)
		{
			tollN =  reader -> readSingleValue(newton_tol);
		}
		else
		{
			tollN = 1e-6;
		}

		if(strcmp(stiffness_tol.c_str(), "NA") != 0)
		{
			isStiff =  reader -> readSingleValue(stiffness_tol);
		}
		else
		{
			isStiff = 1e-6;
		}

		if(strcmp(atol_vector1.c_str(), "NA") != 0)
		{
			atol_vector = reader -> readByLine(atol_vector1);
		}
		else
		{
			for(int i = 0; i < M_0.size(); i++)
				atol_vector.push_back(1e-12);
		}

		if(strcmp(cs_vector1.c_str(), "NA") != 0)
		{
			cs_vector = reader -> readByLineshort(cs_vector1);
		}
		else
		{
			for(int i = 0; i < M_0.size(); i++)
				cs_vector.push_back(i);
		}

		if(strcmp(M_feed1.c_str(), "NA") != 0)
		{
			M_feed = reader -> readByTab(M_feed1);
		}
		else
		{
			for(int i = 0; i < M_0.size(); i++)
				M_feed.push_back(0.0);
		}

		if(strcmp(modelkind1.c_str(), "NA") != 0)
		{
			modelkind = reader -> readType(modelkind1);
		}
		else
		{
			modelkind = 0;
		}

		int bdf = 1;

		if(verbose)
		{
			cout << "dt Runge-Kutta-Fehlberg = " << dt_rkf << "\n";
			cout << "dt Backward Euler = " << dt_ei << "\n";
			cout << "Newton's iterations = " << iterN << "\n";
			cout << "Newton's tollerance = " << tollN << "\n";
			cout << "Stiffness tollerance = " << isStiff << "\n";
			cout << "Runge-Kutta-Fehlberg's tollerances = " << atol_vector[0] << "\n";

		}

		if(modelkind == 2)
		{
			cout << "\nError: the modelkind can be: deterministic/concentration or stochastic/number\n";
			cout << "***************************************************************\n\n";
			exit(-1);
		}

		T volume = 1;
		if(modelkind)
		{
			volume = reader -> readSingleValue(volume1);
			if(volume == 0)
			{
				cout << "\nError: the modelkind is stochastic/number, but there is not volume file\n";
				cout << "***************************************************************\n\n";
				exit(-1);
			}
			conversionToDeterministic(M_0, M_feed, c_vector, left_matrix, volume, verbose);
		}

		/* ********************************************************************************************************************************************* */

		/* ******************** creating structure data to save differential equations ******************** */
		// create (B-A)'
		vector<vector<short> > H = createTransposte(left_matrix, right_matrix);

		// dimensions for malloc memory
		int Nb_reaction = left_matrix.size();
		int Nb_species = left_matrix[0].size();
		int Nb_speciesSaving = cs_vector.size();
		int rowStructPoli = count_nonzero(H);
		int rowStructMoni = count_nonzero(left_matrix);
		int rowOffsetMoni = count_nonzeroRows(left_matrix);
		int Nb_times = t_vector.size();

		int lenJac = Nb_species * Nb_species;

		//polynomials
		short4* ode_poli = (short4*) malloc(rowStructPoli * sizeof(short4));
		short2* offset_poli = (short2*) malloc(Nb_species * sizeof(short2));
		create_sparse(ode_poli, H, offset_poli, 1);

		//monomial
		short4* ode_moni = (short4*) malloc(rowStructMoni * sizeof(short4));
		short2* offset_moni = (short2*) malloc(rowOffsetMoni * sizeof(short2));
		create_sparse(ode_moni, left_matrix, offset_moni, 0);
		/* ********************************************************************************************************************************************* */

		/* ******************** creation pointers from vectors ******************** */
		int* cs_vector_point = &cs_vector[0];
		T* c_vector_point = &c_vector[0];
		T* atol_vector_point = &atol_vector[0];
		T* M_0_point = &M_0[0];
		T* M_feed_point = &M_feed[0];

		T* t_vector_point = &t_vector[0];

		FILE* fd;
		string folder1 = folder + slash + "output" + slash + "SpeciesSampling";
		fd=fopen(folder1.c_str(), "w");
		if( fd==NULL )
		{
			perror("Error open output file SpeciesSampling");
			exit(-1);
		}
		fprintf(fd, "Saving species' dynamics: \n");
		for(int j=0; j < Nb_speciesSaving; j++)
		{
			fprintf(fd, "%d; ", cs_vector_point[j]);
		}
		fprintf(fd, "\n");
		fclose(fd);

		/* ********************************************************************************************************************************************* */
		/* ******************** allocation pointers on gpu ******************** */

		// structures for ODEs reconstruction
		short4* ode_poliDEV;
		short4* ode_moniDEV;
		short2* offset_poliDEV;
		short2* offset_moniDEV;

		//MALLOCC
		CUDA_CHECK_RETURN(cudaMalloc((void**) &ode_poliDEV, rowStructPoli * sizeof(short4)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &ode_moniDEV, rowStructMoni * sizeof(short4)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &offset_poliDEV, Nb_species * sizeof(short2)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &offset_moniDEV, rowOffsetMoni * sizeof(short2)));
		//MEMCPY
		CUDA_CHECK_RETURN(cudaMemcpy(ode_poliDEV, ode_poli, rowStructPoli * sizeof(short4), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(ode_moniDEV, ode_moni, rowStructMoni * sizeof(short4), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(offset_poliDEV, offset_poli, Nb_species * sizeof(short2), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(offset_moniDEV, offset_moni, rowOffsetMoni * sizeof(short2), cudaMemcpyHostToDevice));

		// other pointers for the simulation
		int* cs_vectorDEV;
		T* c_vectorDEV;
		T* t_vectorDEV;
		T* atol_vectorDEV;
		T* M_0DEV;
		T* M_feedDEV;

		// Runge-Kutta-Fehlberg pointers
		T *dev_k1, *dev_k2, *dev_k3, *dev_k4, *dev_k5, *dev_k6, *dev_w1, *dev_w2, *dev_u, *dev_R, *dev_delta;
		int *dev_value;

		//pointer solutions
		T *dinamic;
		T *dinamicDEV;
		dinamic = (T*) malloc((Nb_speciesSaving * Nb_times) * sizeof(T));

		//MALLOCC
		CUDA_CHECK_RETURN(cudaMalloc((void**) &cs_vectorDEV, Nb_speciesSaving * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &c_vectorDEV, Nb_reaction * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &t_vectorDEV, Nb_times * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &atol_vectorDEV, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &M_0DEV, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &M_feedDEV, Nb_species * sizeof(T)));

		cudaDeviceSynchronize();
		//MEMCPY
		CUDA_CHECK_RETURN(cudaMemcpy(cs_vectorDEV, cs_vector_point, Nb_speciesSaving * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(c_vectorDEV, c_vector_point, Nb_reaction * sizeof(T), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(t_vectorDEV, t_vector_point, Nb_times * sizeof(T), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(atol_vectorDEV, atol_vector_point, Nb_species * sizeof(T), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(M_0DEV, M_0_point, Nb_species * sizeof(T), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(M_feedDEV, M_feed_point, Nb_species * sizeof(T), cudaMemcpyHostToDevice));


		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_k1, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_k2, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_k3, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_k4, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_k5, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_k6, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_w1, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_w2, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_u, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_R, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_delta, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dinamicDEV, (Nb_speciesSaving * Nb_times) * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_value, Nb_species * sizeof(int)));


		// Backward euler pointers
		T *dev_y_v, *dev_y_v1, *dev_x, *dev_matJac, *dev_matJacTran, *dev_I, *I, *dev_xJ, *dev_xJ1, *dev_xJ2;
		I = (T*) malloc((lenJac) * sizeof(T));
		T* vett = (T*) malloc((lenJac) * sizeof(T));

		T *dev_y1, *dev_y2, *dev_y3, *dev_y4, *dev_y5, *dev_y6;


		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_y_v, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_y_v1, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_x, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_I, (lenJac) * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_xJ, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_xJ1, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_xJ2, Nb_species * sizeof(T)));

		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_matJac, (lenJac) * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_matJacTran, (lenJac) * sizeof(T)));

		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_y1, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_y2, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_y3, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_y4, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_y5, Nb_species * sizeof(T)));
		CUDA_CHECK_RETURN(cudaMalloc((void**) &dev_y6, Nb_species * sizeof(T)));

		int p = 0;
		for(int i=0; i < Nb_species; i++)
		{
			for(int j=0; j < Nb_species; j++)
			{
				if(i==j)
					I[p] = 1;
				else
					I[p] = 0;
				p = p + 1;
			}
		}

		CUDA_CHECK_RETURN(cudaMemcpy(dev_I, I, (lenJac) * sizeof(T), cudaMemcpyHostToDevice));


		/* ********************************************************************************************************************************************* */

		/* ******************** Selection number of blocks and threads ******************** */

		int *blockThread = (int*)malloc(2 * sizeof (int));
		int sm = returnSM();
		// treads and block ODEs
		returnBlockThread(sm, Nb_species, folder, "ODE", 0, blockThread);
		int n_blocks = blockThread[0];
		int n_threads = blockThread[1];

		// treads and block Jacobian
		returnBlockThread(sm, (lenJac), folder, "JAC", 1, blockThread);
		int n_blocksJac = blockThread[0];
		int n_threadsJac = blockThread[1];

		// treads and block sampling species
		returnBlockThread(sm, Nb_speciesSaving, folder, "DINAM-SUBSPEC", 1, blockThread);
		int n_blocksDin = blockThread[0];
		int n_threadsDin = blockThread[1];

		free(blockThread);

		if(verbose)
		{
			cout << "\n\n***************************************************************";
			cout << "\nTHREAD INFORMATIONS:\n";
			cout << "Number of blocks: " << n_blocks << endl;
			cout << "Number of threads: " << n_threads << endl;
			cout << "Number of Jacobian blocks: " << n_blocksJac << endl;
			cout << "Number of Jaobian threads: " << n_threadsJac << endl;
		}

		/* ********************************************************************************************************************************************* */
		/* ******************** Simulation ******************** */

		if(verbose)
		{
			print(left_matrix, right_matrix, cs_vector, c_vector, t_vector, atol_vector, M_0, M_feed, modelkind, volume);
			cout << "\n\n###############################################################";
			cout << "\nStart ode solver\n\n" ;
		}


		left_matrix.clear();
		right_matrix.clear();
		cs_vector.clear();
		c_vector.clear();
		t_vector.clear();
		atol_vector.clear();
		M_0.clear();
		M_feed.clear();

		int i = 0;
		int count = 0;
		T t1 = t_vector_point[0];
		T t = 0;
		T dt = dt_rkf;

		CUDA_CHECK_RETURN(cudaMemcpy(dev_y1, M_0DEV, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
		inizializza<T><<<n_blocks, n_threads>>>(Nb_species, dev_y2);
		inizializza<T><<<n_blocks, n_threads>>>(Nb_species, dev_y3);
		inizializza<T><<<n_blocks, n_threads>>>(Nb_species, dev_y4);
		inizializza<T><<<n_blocks, n_threads>>>(Nb_species, dev_y5);
		inizializza<T><<<n_blocks, n_threads>>>(Nb_species, dev_y6);
		cudaDeviceSynchronize();

		int numStep = 1;
		int ret = 0;
		T newDtBDF = dt_ei;

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);

		if(verbose)
		{
			cout << "****Start simulation***" << endl;
		}

		while(i < (Nb_times))
		{
			if(fabs(t1 - t) < 1e-5)
			{
				saveDinamic<T><<<n_blocksDin, n_threadsDin>>>(Nb_speciesSaving, i, cs_vectorDEV, dev_y1, dinamicDEV);
				cudaDeviceSynchronize();
				i++;
				if(verbose)
				{
					if (i == Nb_times/4)
						cout << "****1/4 of simulation***" << endl;
					else if (i == Nb_times/2 )
						cout << "****1/2 of simulation***" << endl;
					else if (i == 3*Nb_times/4)
						cout << "****3/4 of simulation***" << endl;
				}
				t1 = t_vector_point[i];
			}

			// k1
			generate_ode<T><<<n_blocks, n_threads>>>(ode_poliDEV, offset_poliDEV, ode_moniDEV, offset_moniDEV, c_vectorDEV, dev_y1, dev_k1, Nb_species);
			cudaDeviceSynchronize();
			rk_Fehlberg1<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y1, dev_k1, dev_u, M_feedDEV);
			cudaDeviceSynchronize();

			// k2
			generate_ode<T><<<n_blocks, n_threads>>>(ode_poliDEV, offset_poliDEV, ode_moniDEV, offset_moniDEV, c_vectorDEV, dev_u, dev_k2, Nb_species);
			cudaDeviceSynchronize();
			rk_Fehlberg2<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y1, dev_k1, dev_k2, dev_u, M_feedDEV);
			cudaDeviceSynchronize();

			// k3
			generate_ode<T><<<n_blocks, n_threads>>>(ode_poliDEV, offset_poliDEV, ode_moniDEV, offset_moniDEV, c_vectorDEV, dev_u, dev_k3, Nb_species);
			cudaDeviceSynchronize();
			rk_Fehlberg3<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y1, dev_k1, dev_k2, dev_k3, dev_u, M_feedDEV);
			cudaDeviceSynchronize();

			// k4
			generate_ode<T><<<n_blocks, n_threads>>>(ode_poliDEV, offset_poliDEV, ode_moniDEV, offset_moniDEV, c_vectorDEV, dev_u, dev_k4, Nb_species);
			cudaDeviceSynchronize();
			rk_Fehlberg4<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y1, dev_k1, dev_k2, dev_k3, dev_k4, dev_u, M_feedDEV);
			cudaDeviceSynchronize();

			// k5
			generate_ode<T><<<n_blocks, n_threads>>>(ode_poliDEV, offset_poliDEV, ode_moniDEV, offset_moniDEV, c_vectorDEV, dev_u, dev_k5, Nb_species);
			cudaDeviceSynchronize();
			rk_Fehlberg5<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y1, dev_k1, dev_k2, dev_k3, dev_k4, dev_k5, dev_u, M_feedDEV);
			cudaDeviceSynchronize();

			// k6
			generate_ode<T><<<n_blocks, n_threads>>>(ode_poliDEV, offset_poliDEV, ode_moniDEV, offset_moniDEV, c_vectorDEV, dev_u, dev_k6, Nb_species);
			cudaDeviceSynchronize();
			rk_Fehlberg6<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_k6, M_feedDEV);
			cudaDeviceSynchronize();

			// u and w
			rk_Fehlberg_agg<T><<<n_blocks, n_threads>>>(Nb_species, dev_y1, dev_k1, dev_k2, dev_k3, dev_k4, dev_k5, dev_k6, dev_w1, M_feedDEV, 0);
			cudaDeviceSynchronize();
			rk_Fehlberg_agg<T><<<n_blocks, n_threads>>>(Nb_species, dev_y1, dev_k1, dev_k2, dev_k3, dev_k4, dev_k5, dev_k6, dev_w2, M_feedDEV, 1);
			cudaDeviceSynchronize();

			calculate_delta<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_w1, dev_w2, dev_R, dev_delta, dev_value, atol_vectorDEV);
			cudaDeviceSynchronize();

			//acceptance u as solution
			count = 0;
			thrust::device_vector<short> vec(dev_value, dev_value + Nb_species);
			count =  thrust::count(vec.begin(), vec.end(), 1);

			if(count == Nb_species)
			{
				CUDA_CHECK_RETURN(cudaMemcpy(dev_y6, dev_y5, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
				CUDA_CHECK_RETURN(cudaMemcpy(dev_y5, dev_y4, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
				CUDA_CHECK_RETURN(cudaMemcpy(dev_y4, dev_y3, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
				CUDA_CHECK_RETURN(cudaMemcpy(dev_y3, dev_y2, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
				CUDA_CHECK_RETURN(cudaMemcpy(dev_y2, dev_y1, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
				CUDA_CHECK_RETURN(cudaMemcpy(dev_y1, dev_w1, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
				cudaDeviceSynchronize();


				thrust::device_vector<double> d_x(dev_delta, dev_delta + Nb_species);
				thrust::device_vector<double>::iterator iter = thrust::min_element(d_x.begin(), d_x.end());
				T d_temp = *iter;
				d_temp = dt * d_temp;

				if((t + d_temp) > t1)
				{
					dt = fabs(t1 - t);
				}
				else
				{
					dt = d_temp;
				}

				t = t + dt;
			}
			else

			// repeat RKF or BE
			{
				thrust::device_vector<double> d_x(dev_delta, dev_delta + Nb_species);
				thrust::device_vector<double>::iterator iter = thrust::min_element(d_x.begin(), d_x.end());
				T d_temp = *iter;
				T dt_old = dt;
				dt = dt * d_temp;
				if(dt < isStiff | isnan(dt) | isinf(dt) | fabs(dt_old - dt) < 1e-3)
				{
					dt = newDtBDF;

					if( (t + dt) > t1 )
					{
						dt = fabs(t1 - t);
					}

					CUDA_CHECK_RETURN(cudaMemcpy(dev_y_v, dev_y1, Nb_species*sizeof(T),cudaMemcpyDeviceToDevice));

					int itN = 0;
					T errN = 1000;
					bool luFail = false;
					// Newton-Raphson method
					while(itN < iterN & errN > tollN)
					{
						generate_ode<T><<<n_blocks, n_threads>>>(ode_poliDEV, offset_poliDEV, ode_moniDEV, offset_moniDEV, c_vectorDEV, dev_y_v, dev_x, Nb_species);
						cudaDeviceSynchronize();

						switch (numStep)
						{
						case 6:
							calculate_nonLinear6th<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y_v, dev_y1, dev_y2, dev_y3, dev_y4, dev_y5, dev_y6, dev_x);
							cudaDeviceSynchronize();
							break;
						case 1:
							calculate_nonLinear1th<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y_v, dev_y1, dev_x);
							cudaDeviceSynchronize();
							break;
						case 2:
							calculate_nonLinear2th<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y_v, dev_y1, dev_y2, dev_x);
							cudaDeviceSynchronize();
							break;
						case 3:
							calculate_nonLinear3th<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y_v, dev_y1, dev_y2, dev_y3, dev_x);
							cudaDeviceSynchronize();
							break;
						case 4:
							calculate_nonLinear4th<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y_v, dev_y1, dev_y2, dev_y3, dev_y4, dev_x);
							cudaDeviceSynchronize();
							break;
						case 5:
							calculate_nonLinear5th<T><<<n_blocks, n_threads>>>(Nb_species, dt, dev_y_v, dev_y1, dev_y2, dev_y3, dev_y4, dev_y5, dev_x);
							cudaDeviceSynchronize();
							break;
						default:
							exit(-2);
						}

						if(itN == 0)
						{
							createJac<T><<<n_blocks, n_threads>>>(ode_poliDEV, offset_poliDEV, ode_moniDEV, offset_moniDEV, c_vectorDEV, dev_y_v, dev_matJac, Nb_species);
							cudaDeviceSynchronize();
							switch (numStep)
							{
							case 6:
								aggiorna_Jacobiano6<T><<<n_blocksJac, n_threadsJac>>>((lenJac), dt, dev_matJac, dev_I);
								cudaDeviceSynchronize();
								break;
							case 1:
								aggiorna_Jacobiano1<T><<<n_blocksJac, n_threadsJac>>>((lenJac), dt, dev_matJac, dev_I);
								cudaDeviceSynchronize();
								break;
							case 2:
								aggiorna_Jacobiano2<T><<<n_blocksJac, n_threadsJac>>>((lenJac), dt, dev_matJac, dev_I);
								cudaDeviceSynchronize();
								break;
							case 3:
								aggiorna_Jacobiano3<T><<<n_blocksJac, n_threadsJac>>>((lenJac), dt, dev_matJac, dev_I);
								cudaDeviceSynchronize();
								break;
							case 4:
								aggiorna_Jacobiano4<T><<<n_blocksJac, n_threadsJac>>>((lenJac), dt, dev_matJac, dev_I);
								cudaDeviceSynchronize();
								break;
							case 5:
								aggiorna_Jacobiano5<T><<<n_blocksJac, n_threadsJac>>>((lenJac), dt, dev_matJac, dev_I);
								cudaDeviceSynchronize();
								break;
							default:
								exit(-2);
							}

							transposteJacobian<T><<<n_blocksJac, n_threadsJac>>>(Nb_species*Nb_species, Nb_species, dev_matJac, dev_matJacTran);
							cudaDeviceSynchronize();
						}

						// LU decomposition --> solver linear system
						ret = linearSolverLU(handle, Nb_species, dev_matJacTran, dev_x, dev_xJ);

						if(ret == -1)
						{
							aggiornaNaN<T><<<n_blocks, n_threads>>>(Nb_species, dev_y_v, dev_y1);
							newDtBDF = dt * 0.1;
							luFail = true;
							break;
						}

						aggiorna_Newton<T><<<n_blocks, n_threads>>>(Nb_species, dev_y_v, dev_y_v1, dev_xJ, M_feedDEV);
						cudaDeviceSynchronize();

						divNorm<T><<<n_blocks, n_threads>>>(Nb_species, dev_xJ, dev_y_v1, dev_y_v1, M_feedDEV);
						cudaDeviceSynchronize();

						thrust::device_vector<T> d_x(dev_y_v1, dev_y_v1 + Nb_species);

						square unary_op;
						thrust::plus<T> binary_op;
						T init = 0.0f;
						T norm = sqrt(thrust::transform_reduce(d_x.begin(), d_x.end(), unary_op, init, binary_op));

						errN = norm;
						itN = itN + 1;
					}

					if(luFail == false)
					{
						count = 0;

						thrust::device_vector<T> d_y(dev_y_v, dev_y_v + Nb_species);
						checkEI unary_op;
						thrust::plus<T> binary_op;
						int init1 = 0;
						count  = thrust::transform_reduce(d_y.begin(), d_y.end(), unary_op, init1, binary_op);

						if(count > 0)
						{
							aggiornaNaN<T><<<n_blocks, n_threads>>>(Nb_species, dev_y_v, dev_y1);
							newDtBDF = dt * 0.1;
						}
					}
					
					CUDA_CHECK_RETURN(cudaMemcpy(dev_y6, dev_y5, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
					CUDA_CHECK_RETURN(cudaMemcpy(dev_y5, dev_y4, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
					CUDA_CHECK_RETURN(cudaMemcpy(dev_y4, dev_y3, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
					CUDA_CHECK_RETURN(cudaMemcpy(dev_y3, dev_y2, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
					CUDA_CHECK_RETURN(cudaMemcpy(dev_y2, dev_y1, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));
					CUDA_CHECK_RETURN(cudaMemcpy(dev_y1, dev_y_v, Nb_species * sizeof(T),cudaMemcpyDeviceToDevice));

					t = t + dt;
				}
			}

			if(numStep < bdf)
				numStep++;
		}
		cudaEventRecord( stop, 0 );
		cudaEventSynchronize( stop );
		float tempo = 0;
		cudaEventElapsedTime( &tempo, start, stop );

		if(verbose)
		{
			cout << "****End simulation***" << endl;
		}

		/* ********************************************************************************************************************************************* */
		/* ******************** Saving dynamics ******************** */

		CUDA_CHECK_RETURN(cudaMemcpy(dinamic, dinamicDEV, (Nb_times * Nb_speciesSaving) * sizeof(T), cudaMemcpyDeviceToHost));

		int pos = 0;
		folder1 = folder + slash + "output" + slash + "Solution";
		fd=fopen(folder1.c_str(), "w");
		if( fd==NULL )
		{
			perror("Error open output file Solution.txt");
			exit(-1);
		}

		if(modelkind)
		{
			T avo = 6.022e23;
			for(int i = 0; i < Nb_speciesSaving * Nb_times; i++)
				dinamic[i] = dinamic[i]*(volume * avo);
		}

		for (int i = 0; i < Nb_times; i++)
		{
			fprintf (fd, "%g\t", t_vector_point[i]);
			pos = i * Nb_speciesSaving;
			for(int j = 0; j < Nb_speciesSaving; j++)
			{
				fprintf(fd, "%.8g\t", dinamic[pos]);
				pos++;
			}
			fprintf(fd, "\n");
		}
		fclose(fd);

		if(verbose)
		{
			cout << "\n\nEnd ode solver\n" ;
			cout << "###############################################################\n\n";

		}
		/* ********************************************************************************************************************************************* */
		return tempo;
	}
	/* ******************** Convertion stochastic into deterministic ******************** */
	void conversionToDeterministic(vector<T> &M_0, vector<T> &M_feed, vector<T> &c_vector, vector<vector<short> > A, T volume, int verbose)
	{
		T avo = 6.022e23;
		T conv =  volume*avo;

		if(verbose)
		{
			cout << "\n\n***************************************************************";
			cout << "\nConvertion initial conditions and feeds with volume: " << volume << endl;
		}

		for(int i=0; i < M_0.size(); i++)
		{
			M_0[i] = M_0[i]/conv;
			M_feed[i] = M_feed[i]/conv;
		}

		for(int i=0; i < c_vector.size(); i++)
		{
			//cout << i << "\n";
			if(i == 0 & verbose)
				cout << "Convertion reaction rate constants with volume: " << volume << endl;

			int value = count_nonzeroRow(A, i);
			if(value >= 3)
			{
				cout << "There are one or more reactions with higher order than the second one" << endl;
				exit(-1);
			}
			else
			{
				int om = 0;
				for(int j = 0; j < A[0].size(); j++)
				{
					if(A[i][j] >= 2)
					{
						c_vector[i] =  (c_vector[i] * conv)/2.0;
						om = 1;
						break;
					}
				}
				if(om == 0)
				{
					if(count_nonzeroRow(A, i) > 1)
						c_vector[i] = c_vector[i]*conv;
				}
			}
		}
	}


	/* ******************** Methods to creaye ODEs structures ******************** */

	int count_nonzeroRow(vector<vector<short> > A, int i)
	{
		int count = 0;
		for(int j = 0; j < A[0].size(); j++)
		{
			if(A[i][j] > 0)
				count++;
		}
		return count;
	}

	int count_nonzeroRows(vector<vector<short> > A)
	{
		int count = 0;
		for(int i = 0; i < A.size(); i++)
		{
			int somma = 0;
			for(int j = 0; j < A[0].size(); j++)
			{
				if(A[i][j] > 0)
					somma++;
			}
			if(somma > 0)
				count++;
		}
		return count;
	}

	int count_nonzero(vector<vector<short> > A)
	{
		int count = 0;
		for(int i = 0; i < A.size(); i++)
		{
			for(int j = 0; j < A[0].size(); j++)
			{
				if(A[i][j] != 0)
					count ++;
			}
		}
		return count;
	}

	vector<vector<short> > createTransposte(vector<vector<short> > A, vector<vector<short> > B)
		  {
		vector<short> vect;
		vector<vector<short> > H;

		for(int j = 0; j < A[0].size(); j++)
		{
			vect.clear();
			for(int i = 0; i < A.size(); i++)
			{
				int temp = B[i][j] - A[i][j];
				vect.push_back(temp);
			}
			H.push_back(vect);
		}
		return H;
		  }

	void create_sparse(short4* sparse, vector<vector<short> > matrix, short2* offset, int flag)
	{
		int row = (int) matrix.size();
		int col = (int) matrix[0].size();
		int count = 0;
		int index_old = 0;
		for(int i=0; i < row; i++)
		{
			short index = 0;
			for(int j=0; j < col; j++)
			{
				if(matrix[i][j] != 0)
				{
					sparse[count].x = i;
					sparse[count].y = j;
					sparse[count].z = matrix[i][j];
					if(flag == 1)
						sparse[count].w = j;
					count++;
					index++;
				}
			}
			offset[i].x = index_old;
			offset[i].y = index_old + index;
			index_old = index_old + index;
		}
	}


	/* ******************** Print methods ******************** */

	void print(vector<vector<short> > left_matrix, vector<vector<short> > right_matrix,
			vector<int> cs_vector1, vector<T> c_vector1, vector<T> t_vector1, vector<T> atol_vector1,
			vector<T> M_0, vector<T> M_feed1, int modelkind1, T volume)
	{

		//cout << "\nleft_side" << "\n";
		//printMatrix(left_matrix);

		//cout << "\nright_side" << "\n";
		//printMatrix(right_matrix);

		//cout << "\n\n***************************************************************";
		//cout << "\nSampling times...\n";
		//printArray(t_vector1);

		//cout << "\n\n***************************************************************";
		//cout << "\nConcerned species...\n";
		//printArrayshort(cs_vector1);

		//cout << "\n\n***************************************************************";
		//cout << "\nReaction rate constants...\n";
		//printArray(c_vector1);

		//cout << "\n\n***************************************************************";
		//cout << "\nInitial values...\n";
		//printArray(M_0);

		//cout << "\nM_feed" << "\n";
		//printArray(M_feed1);

		//cout << "\natol_vector" << "\n";
		//printArray(atol_vector1);

		cout << "\n\n***************************************************************";
		cout << "\nThe model has " << left_matrix[0].size() << " species and " << left_matrix.size() << " reactions"<<"\n";

		if(!modelkind1)
		{
			cout << "\n\n***************************************************************";
			cout << "\nmodelkind = deterministic" << "\n";
		}
		else
		{
			cout << "\n\n***************************************************************";
			cout << "\nmodelkind = stochastic" << "\n";
			cout << "\n\n***************************************************************";
			cout << "\nvolume = " << volume <<"\n";
		}

	}

	void printMatrix(vector<vector<short> > matrix)
	{
		for(int i = 0; i < matrix.size(); i++)
		{
			for(int j = 0; j < matrix[0].size(); j++)
			{
				cout << matrix[i][j] << "\t";
			}
			cout << "\n";
		}
	}

	void printArrayshort(vector<short> vector)
	{
		for(int i = 0; i < vector.size(); i++)
		{
			cout << vector[i] << "\t";

		}
		cout << "\n";
	}

	void printArray(vector<T> vector)
	{
		for(int i = 0; i < vector.size(); i++)
		{
			cout << vector[i] << "\t";

		}
		cout << "\n";
	}

	int returnSM()
	{
		cudaDeviceProp  prop;
		CUDA_CHECK_RETURN(cudaGetDeviceProperties(&prop, 0));
		return prop.multiProcessorCount;
	}

	void returnBlockThread(int sm, int thread, string folder, string type, int save, int* info)
	{
		info[0] = 0;
		info[1] = 0;
		string folder1 = folder + slash + "output" + slash + "infoBlockThread";
		FILE* fd;

		if(save == 0)
			fd=fopen(folder1.c_str(), "w");
		else
			fd=fopen(folder1.c_str(), "a");

		if(save == 0)
			fprintf(fd, "%s \n", type.c_str());
		else
			fprintf(fd, "\n\n%s \n", type.c_str());

		fprintf(fd, "Allocating %d threads over %d available streaming multiprocessors\n", thread, sm);
		float threadF = thread;
		int putative_threads_per_block =  minThread(256.0, ceil(threadF/sm));
		fprintf(fd, "Putative TPB: %d \n", putative_threads_per_block);

		float used_blocks;
		float final_tpb;
		if(putative_threads_per_block % 32 == 0)
		{
			final_tpb = putative_threads_per_block * 1.;
			used_blocks = ceil(thread/final_tpb);
		}
		else
		{
			fprintf(fd, "Adapting to warp size\n");
			final_tpb =  (putative_threads_per_block/32+1)*32.;
			used_blocks = ceil(thread/final_tpb);
		}
		info[0] = (int) used_blocks;
		info[1] = (int) final_tpb;
		int spawned = (int) final_tpb*used_blocks;

		fprintf(fd, "Final TPB: %d \n", info[1]);
		fprintf(fd, "Number of used blocks: %d \n", info[0]);
		fprintf(fd, "Total threads spawned: %d \n", spawned);
		fclose(fd);

	}

	int minThread(int n1, int n2)
	{
		if (n1 < n2)
			return n1;
		else
			return n2;
	}


	int linearSolverLU(cusolverDnHandle_t handle, int n, const T* Acopy, const T* b, T* x)
	{
		int bufferSize = 0;
		int *info = NULL;
		T* buffer = NULL;
		T* A = NULL;
		int *ipiv = NULL; // pivoting sequence
		int h_info = 0;

		int ret = 0;


		cusolverDnTTgetrf_bufferSize(handle, n, n, (T*)Acopy, n, &bufferSize);

		cudaMalloc(&info, sizeof(int));
		cudaMalloc(&buffer, sizeof(T)*bufferSize);
		cudaMalloc(&A, sizeof(T)*n*n);
		cudaMalloc(&ipiv, sizeof(int)*n);


		// prepare a copy of A because getrf will overwrite A with L
		cudaMemcpy(A, Acopy, sizeof(T)*n*n, cudaMemcpyDeviceToDevice);
		cudaMemset(info, 0, sizeof(int));

		cusolverDnTTgetrf(handle, n, n, A, n, buffer, ipiv, info);
		cudaDeviceSynchronize();

		cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);

		if ( 0 != h_info )
		{
			ret = -1;
		}
		else
		{
			cudaMemcpy(x, b, sizeof(T)*n, cudaMemcpyDeviceToDevice);
			cusolverDnTTgetrs(handle, CUBLAS_OP_N, n, 1, A, n, ipiv, x, n, info);
			cudaDeviceSynchronize();
		}

		if (info  ) { (cudaFree(info  )); }
		if (buffer) { (cudaFree(buffer)); }
		if (A     ) { (cudaFree(A)); }
		if (ipiv  ) { (cudaFree(ipiv));}

		return ret;
	}

	//cusolverDn<TT>getrf_bufferSize
	inline cusolverStatus_t cusolverDnTTgetrf_bufferSize(cusolverDnHandle_t handle, int n, int n1,  const double* Acopy, int n2, int* bufferSize)
	{
		return cusolverDnDgetrf_bufferSize(handle, n, n, (double*)Acopy, n, bufferSize);
	}

	inline cusolverStatus_t cusolverDnTTgetrf_bufferSize(cusolverDnHandle_t handle, int n, int n1,  const float* Acopy, int n2, int* bufferSize)
	{
		return cusolverDnSgetrf_bufferSize(handle, n, n, (float*)Acopy, n, bufferSize);
	}

	//cusolverDn<TT>getrf
	inline cusolverStatus_t cusolverDnTTgetrf(cusolverDnHandle_t handle, int n, int n1,  double* A, int n2, double* buffer, int* ipiv, int* info)
	{
		return cusolverDnDgetrf(handle, n, n, A, n, buffer, ipiv, info);
	}

	inline cusolverStatus_t cusolverDnTTgetrf(cusolverDnHandle_t handle, int n, int n1,  float* A, int n2, float* buffer, int* ipiv, int* info)
	{
		return cusolverDnSgetrf(handle, n, n, A, n, buffer, ipiv, info);
	}

	//cusolverDn<TT>getrs
	inline cusolverStatus_t cusolverDnTTgetrs(cusolverDnHandle_t handle, cublasOperation_t trans, int n, int n1,  double* A, int n2, int* ipiv, double* x, int n3, int* info)
	{
		return cusolverDnDgetrs(handle, trans, n, n1, A, n, ipiv, x, n3, info);
	}

	inline cusolverStatus_t cusolverDnTTgetrs(cusolverDnHandle_t handle, cublasOperation_t trans, int n, int n1,  float* A, int n2, int* ipiv, float* x, int n3, int* info)
	{
		return cusolverDnSgetrs(handle, trans, n, n1, A, n, ipiv, x, n3, info);
	}
};