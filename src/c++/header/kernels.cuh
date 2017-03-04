#include "cuda.h"
#include "cuda_runtime.h"
#include <iostream>   // std::cout

template<class T>
__global__ void generate_ode(short4* ode_poli, short2* offset_poli, short4* ode_moni, short2* offset_moni, T* c_vector, T* y0, T* y1, const int Nb_species)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < Nb_species)
	{
		y1[gid] = 0;

		int start = offset_poli[gid].x;
		int end = offset_poli[gid].y;

		for(int i = start; i < end; i++)
		{
			int y = ode_poli[i].y;
			int signum = ode_poli[i].z;
			int k_pos = ode_poli[i].w;

			int start_moni = offset_moni[y].x;
			int end_moni = offset_moni[y].y;

			T temp = 0;
			for(int j = start_moni; j < end_moni; j++)
			{
				if(j == start_moni)
				{
					if(ode_moni[j].z == 1)
					{
						temp = temp + signum * c_vector[k_pos] * y0[ode_moni[j].y];
					}
					else
					{
						if(ode_moni[j].z == 2)
						{
							temp = temp + signum * c_vector[k_pos] * y0[ode_moni[j].y] * y0[ode_moni[j].y];
						}
					}
				}
				else
				{
					if(ode_moni[j].z == 1)
					{
						temp = temp * y0[ode_moni[j].y];
					}
					else
					{
						if(ode_moni[j].z == 2)
						{
							temp = temp * y0[ode_moni[j].y] * y0[ode_moni[j].y];
						}
					}
				}
			}
			y1[gid] = y1[gid] + temp;
		}
	}
}

/*********************************************************************************/
//RUNGE KUTTA Fehlberg kernels//

template<class T>
__global__ void rk_Fehlberg1(const int s, T dt, T* u0, T* k1, T* u, T* feed)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		if(feed[gid] != 0.0)
			k1[gid] = feed[gid];
		else
		{
			k1[gid] = k1[gid]*dt;
			u[gid] = u0[gid] + k1[gid] / 4.0;
		}
	}
}

template<class T>
__global__ void rk_Fehlberg2(const int s, T dt, T* u0, const T* __restrict__ k1, T* k2, T* u, T* feed)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		if(feed[gid] != 0.0)
			k2[gid] = feed[gid];
		else
		{
			k2[gid] = k2[gid]*dt;
			u[gid] = u0[gid] + (3.0*k1[gid])/32.0 + (9.0*k2[gid])/32.0;
		}
	}
}

template<class T>
__global__ void rk_Fehlberg3(const int s, T dt, T* u0, const T* __restrict__ k1, const T* __restrict__ k2, T* k3, T* u, T* feed)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		if(feed[gid] != 0.0)
			k3[gid] = feed[gid];
		else
		{
			k3[gid] = k3[gid]*dt;
			u[gid] = u0[gid] + (1932.0*k1[gid])/2197.0 - (7200.0*k2[gid])/2197.0 + (7296.0*k3[gid])/2197.0;
		}
	}
}

template<class T>
__global__ void rk_Fehlberg4(const int s, T dt, T* u0, const T* __restrict__ k1, const T* __restrict__ k2, const T* __restrict__ k3, T* k4, T* u, T* feed)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		if(feed[gid] != 0.0)
			k4[gid] = feed[gid];
		else
		{
			k4[gid] = k4[gid]*dt;
			u[gid] = u0[gid] + (439.0*k1[gid])/216.0 - (8.0*k2[gid]) + (3680.0*k3[gid])/513.0 - (845.0*k4[gid])/4104.0;
		}
	}
}

template<class T>
__global__ void rk_Fehlberg5(const int s, T dt, T* u0, const T* __restrict__ k1, const T* __restrict__ k2, const T* __restrict__ k3, const T* __restrict__ k4, T* k5, T* u, T* feed)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		if(feed[gid] != 0.0)
			k5[gid] = feed[gid];
		else
		{
			k5[gid] = k5[gid]*dt;
			u[gid] = u0[gid] - (8.0*k1[gid])/27.0 + (2.0*k2[gid]) - (3544.0*k3[gid])/2565.0 + (1859.0*k4[gid])/4104.0 - (11.0*k5[gid])/40.0;
		}
	}
}

template<class T>
__global__ void rk_Fehlberg6(const int s, T dt, T* k6, T* feed)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		if(feed[gid] != 0.0)
			k6[gid] = feed[gid];
		else
			k6[gid] = k6[gid]*dt;
	}
}

template<class T>
__global__ void rk_Fehlberg_agg(const int s, T *u0, const T* __restrict__ k1, const T* __restrict__ k2, const T* __restrict__ k3, const T* __restrict__ k4, const T* __restrict__ k5, const T* __restrict__ k6, T *w, T *feed, int caso)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		if(caso == 0)
		{
			if(feed[gid] != 0.0)
				w[gid] = feed[gid];
			else
				w[gid] = u0[gid] + (25.0*k1[gid])/216.0 + (1408.0*k3[gid])/2565.0 + (2197.0*k4[gid])/4104.0- k5[gid]/5.0;
		}
		else
		{
			if(feed[gid] != 0.0)
				w[gid] = feed[gid];
			else
				w[gid] = u0[gid] + (16.0*k1[gid])/135.0 + (6656.0*k3[gid])/12825.0 + (28561.0*k4[gid])/56430.0 - (9.0*k5[gid]/50.0) + (2.0*k6[gid])/55.0;
		}
	}
}

template<class T>
__global__ void calculate_delta(const int s, T dt,  const T* __restrict__ w1, const T* __restrict__ w2,
		T *R, T *delta, int *value, const T* __restrict__ tolerance)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		R[gid] = std::fabs(w1[gid]-w2[gid])/dt;
		T ex = 0.25;
		delta[gid] = 0.84*std::pow(tolerance[gid]/R[gid], ex);
		if(R[gid] <= tolerance[gid] && !isnan(R[gid]) && !isinf(R[gid]) && w1[gid] > 0)
			value[gid] = 1;
		else
			value[gid] = 0;
	}
}


/******************************************************************************/
// Create Jacobian//

template<class T>
__global__ void createJac(short4* ode_poli, short2* offset_poli, short4* ode_moni, short2* offset_moni, T* c_vector, T* y0, T* y1, int Nb_species)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < Nb_species)
	{
		for(int i = 0; i < Nb_species ; i++)
		{
			createJacSub(ode_poli, offset_poli, ode_moni, offset_moni, c_vector, y0, y1, Nb_species, i);
			__syncthreads();
		}
	}
}

template<class T>
__device__ void createJacSub(short4* ode_poli, short2* offset_poli, short4* ode_moni, short2* offset_moni, T* c_vector, T* y0, T* y1, int Nb_species, int index2)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < Nb_species)
	{
		int ii = gid * Nb_species + index2;
		y1[ii] = 0;
		int start = offset_poli[gid].x;
		int end = offset_poli[gid].y;

		for(int i = start; i < end; i++)
		{
			int y = ode_poli[i].y;
			int signum = ode_poli[i].z;
			int k_pos = ode_poli[i].w;

			int start_moni = offset_moni[y].x;
			int end_moni = offset_moni[y].y;

			T temp = 0;
			for(int j = start_moni; j < end_moni; j++)
			{
				if(ode_moni[j].z == 1 & ode_moni[j].y == index2)
				{
					temp = temp + signum * c_vector[k_pos];
					int start_moni1 = offset_moni[y].x;
					int end_moni1 = offset_moni[y].y;
					for(int jj = start_moni1; jj < end_moni1; jj++)
					{
						if(ode_moni[jj].y != index2)
							temp = temp*y0[ode_moni[jj].y];
					}
				}
				else
				{
					if(ode_moni[j].z == 2 & ode_moni[j].y == index2)
					{
						temp = temp + signum * c_vector[k_pos] * 2 * y0[ode_moni[j].y];

						int start_moni1 = offset_moni[y].x;
						int end_moni1 = offset_moni[y].y;
						for(int jj = start_moni1; jj < end_moni1; jj++)
						{
							if(ode_moni[jj].y != index2)
								temp = temp*y0[ode_moni[jj].y];
						}
					}
				}
			}
			if(temp != 0)
				y1[ii] = y1[ii] + temp;
		}
	}
}

template<class T>
__global__ void transposteJacobian(const int s, int rowsA, const T* dev_matJac, T* dev_matJacTran)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		int i = gid/rowsA;
		int j = gid%rowsA;
		dev_matJacTran[gid] = dev_matJac[rowsA*j + i];
	}
}


/******************************************************************************/
// Aggiornamenti sistema non lineare//

template<class T>
__global__ void calculate_nonLinear1th(const int s, T h, T* y_v, T* y1, T* x)
{
	// x è la f(y_n)
	// y_v vettore dell'iterazione di Newton
	// y1 vettore y_(n-1)
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		x[gid] = -(y_v[gid] - y1[gid] -h*x[gid]);
	}
}

template<class T>
__global__ void calculate_nonLinear2th(const int s, T h, T* y_v, T* y1, T* y2, T* x)
{
	// x è la f(y_n)
	// y_v vettore dell'iterazione di Newton
	// y1 => y_(n-1) ...

	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		x[gid] = -(y_v[gid] -(4.0*y1[gid])/3.0 + y2[gid]/3.0 - (2.0*h*x[gid])/3.0);
	}
}

template<class T>
__global__ void calculate_nonLinear3th(const int s, T h, T* y_v, T* y1, T* y2, T* y3, T* x)
{
	// x è la f(y_n)
	// y_v vettore dell'iterazione di Newton
	// y1 => y_(n-1) ...
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
			x[gid] = -(y_v[gid] -
				(18.0*y1[gid])/11.0 +
				(9.0*y2[gid])/11.0 -
				(2.0*y3[gid])/11.0 -
				(6.0*h*x[gid])/11.0);
	}
}

template<class T>
__global__ void calculate_nonLinear4th(const int s, T h, T* y_v, T* y1, T* y2, T* y3, T* y4, T* x)
{
	// x è la f(y_n)
	// y_v vettore dell'iterazione di Newton
	// y1 => y_(n-1) ...
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
			x[gid] = -(y_v[gid] -
				(48.0*y1[gid])/25.0 +
				(36.0*y2[gid])/25.0 -
				(16.0*y3[gid])/25.0 +
				(3.0*y4[gid])/25.0 -
				(12.0*h*x[gid])/25.0);
	}
}

template<class T>
__global__ void calculate_nonLinear5th(const int s, T h, T* y_v, T* y1, T* y2, T* y3, T* y4, T* y5, T* x)
{
	// x è la f(y_n)
	// y_v vettore dell'iterazione di Newton
	// y1 => y_(n-1) ...
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
			x[gid] = -(y_v[gid] -
				(300.0*y1[gid])/147.0 +
				(300.0*y2[gid])/147.0 -
				(200.0*y3[gid])/147.0 +
				(75.0*y4[gid])/147.0 -
				(12.0*y5[gid])/147.0 -
				(60.0*h*x[gid])/147.0);
	}
}

template<class T>
__global__ void calculate_nonLinear6th(const int s, T h, T* y_v, T* y1, T* y2, T* y3, T* y4, T* y5, T* y6, T* x)
{
	// x è la f(y_n)
	// y_v vettore dell'iterazione di Newton
	// y1 => y_(n-1) ...
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		x[gid] = -(y_v[gid] -
				(360.0*y1[gid])/147.0 +
				(450.0*y2[gid])/147.0 -
				(400.0*y3[gid])/147.0 +
				(225.0*y4[gid])/147.0 -
				(72.0*y5[gid])/147.0 +
				(10.0*y6[gid])/147.0 -
				(60.0*h*x[gid])/147.0);
	}
}


/******************************************************************************/
// Aggiornamenti Jacobiano//
template<class T>
__global__ void aggiorna_Jacobiano1(const int s, T h, T* J, T* I)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;

	if(gid < s)
	{
		J[gid] = I[gid] - h*J[gid];
	}
}

template<class T>
__global__ void aggiorna_Jacobiano2(const int s, T h, T* J, T* I)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;

	if(gid < s)
	{
		J[gid] = I[gid] - (2.0*h*J[gid])/3.0;
	}
}

template<class T>
__global__ void aggiorna_Jacobiano3(const int s, T h, T* J, T* I)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;

	if(gid < s)
	{
		J[gid] = I[gid] - (6.0*h*J[gid])/11.0;
	}
}

template<class T>
__global__ void aggiorna_Jacobiano4(const int s, T h, T* J, T* I)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;

	if(gid < s)
	{
		J[gid] = I[gid] - (12.0*h*J[gid])/25.0;
	}
}

template<class T>
__global__ void aggiorna_Jacobiano5(const int s, T h, T* J, T* I)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;

	if(gid < s)
	{
		J[gid] = I[gid] - (60.0*h*J[gid])/137.0;
	}
}

template<class T>
__global__ void aggiorna_Jacobiano6(const int s, T h, T* J, T* I)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;

	if(gid < s)
	{
		J[gid] = I[gid] - (60.0*h*J[gid])/147.0;
	}
}


/******************************************************************************/
// Aggiornamento Newton//
template<class T>
__global__ void aggiorna_Newton(const int s, T* y_v, T* y_v1, T*  delta, T* feed)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		y_v1[gid] = y_v[gid];

		if(feed[gid] == 0)
			y_v[gid] = y_v[gid] + delta[gid];
		else
			y_v[gid] = feed[gid];

		delta[gid] = y_v[gid] - y_v1[gid];
	}
}

template<class T>
__global__ void divNorm(const int s, T* x, T* x_old, T* x_new, T* feed)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{

		if(feed[gid] == 0)
			x_new[gid] = 0;
		else
		{
			if(x_old[gid] != 0)
				x_new[gid] = x[gid] / x_old[gid];
			else
				x_new[gid] = x[gid];
		}
	}

}


/******************************************************************************/
// Altre funzioni//

template<class T>
__global__ void inizializza(const int s, T* x)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		x[gid] = 0.0;
	}
}

template<class T>
__global__ void aggiornaNaN(const int s, T* y_v, T* y1)
{
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	if(gid < s)
	{
		if(isnan(y_v[gid]) || isinf(y_v[gid]) || y_v[gid] < 0)
			y_v[gid] = y1[gid];
	}
}


template<class T>
__global__ void saveDinamic(const int ssp, const int t, const int* __restrict__ subspec, const T* __restrict__ ci, T* dinamic)
{
	int pos;
	int gid = threadIdx.x + blockDim.x * blockIdx.x;
	int p = t * ssp;
	if(gid < ssp)
	{
		pos = subspec[gid];
		dinamic[p + gid] = ci[pos];
	}

}
