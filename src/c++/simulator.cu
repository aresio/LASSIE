#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <cstring>
#include <cmath>
#include "header/Simulation.cuh"

using namespace std;

#ifdef _WIN32
    #include "dirent.h"
#elif _WIN64
	#include "dirent.h"
#elif __linux__
	#include <dirent.h>
#else
    #include <dirent.h>
#endif



__host__ vector<string> openDirectory(string path);

int main(int argc, char* argv[])
{
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	int os = 0;

	#ifdef _WIN32
	os = 0;
	#elif _WIN64
	os = 0;
	#elif __linux__
	os = 1;
	#else
	os = 2;
	#endif


	int sis;

	//matrix reagents
	string A;
	//matrix products
	string B;
	//vector reaction rate
	string c_vector;
	//vector sampling times
	string t_vector;
	//vector initial conditions
	string MX_0;

	//optinal files
	string atol_vector = "NA";
	string alphabet = "NA";
	string be_step = "NA";
	string newton_iter = "NA";
	string newton_tol = "NA";
	string rkf_step = "NA";
	string stiffness_tol = "NA";
	string volume = "NA";
	string modelkind = "NA";
	string cs_vector = "NA";
	string M_feed = "NA";

	// check mandatory parameters
	if(argc == 1)
	{
		cout << "\nError: the first parameter must be the input folder" << "\n";
		cout << "\nError: the second parameter must be the output folder" << "\n";
		cout << "***************************************************************\n\n";
		exit(-3);
	}

	if(argc == 2)
	{
		cout << "\nEror: the second parameter must be the output folder" << "\n";
		cout << "***************************************************************\n\n";
		exit(-2);
	}

	//read model files
	string path = argv[1];

	vector<string> files = openDirectory(path);
	int count[5] = {0};
	for(int ff = 0; ff < files.size(); ff++)
	{
		if(strcmp(files[ff].c_str(), "left_side") == 0)
		{
			A = path + slash + files[ff];
			count[0] = 1;
		}
		else if(strcmp(files[ff].c_str(), "right_side") == 0)
		{
			B = path + slash + files[ff];
			count[1] = 1;
		}
		else if(strcmp(files[ff].c_str(), "c_vector") == 0)
		{
			c_vector = path + slash + files[ff];
			count[2] = 1;
		}
		else if(strcmp(files[ff].c_str(), "t_vector") == 0)
		{
			t_vector = path + slash + files[ff];
			count[3] = 1;
		}
		else if(strcmp(files[ff].c_str(), "MX_0") == 0 || strcmp(files[ff].c_str(), "M_0") == 0)
		{
			MX_0 = path + slash + files[ff];
			count[4] = 1;
		}
		else if(strcmp(files[ff].c_str(), "cs_vector") == 0)
		{
			cs_vector = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "M_feed") == 0)
		{
			M_feed = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "modelkind") == 0)
		{
			modelkind =  path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "atol_vector") == 0)
		{
			atol_vector = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "alphabet") == 0)
		{
			alphabet = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "be_step") == 0)
		{
			be_step = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "newton_iter") == 0)
		{
			newton_iter = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "newton_tol") == 0)
		{
			newton_tol = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "rkf_step") == 0)
		{
			rkf_step = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "stiffness_tol") == 0)
		{
			stiffness_tol = path + slash + files[ff];
		}
		else if(strcmp(files[ff].c_str(), "volume") == 0)
		{
			volume = path + slash + files[ff];
		}
	}

	//check mandatory model files
	bool toExit = false;
	if(count[0] == 0)
	{
		cout << "\n\n***************************************************************\n\n";
		cout << "\nError: there is not left side matrix" << "\n";
		cout << "***************************************************************\n\n";
		toExit = true;
	}
	if(count[1] == 0)
	{
		cout << "\n\n***************************************************************\n\n";
		cout << "\nError: there is not right side matrix" << "\n";
		cout << "***************************************************************\n\n";
		toExit = true;
	}
	if(count[2] == 0)
	{
		cout << "\n\n***************************************************************\n\n";
		cout << "\nError: there is not kinetic constants vector" << "\n";
		cout << "***************************************************************\n\n";
		toExit = true;
	}
	if(count[3] == 0)
	{
		cout << "\n\n***************************************************************\n\n";
		cout << "\nError: there is not vector of time instants" << "\n";
		cout << "***************************************************************\n\n";
		toExit = true;
	}
	if(count[4] == 0)
	{
		cout << "\n\n***************************************************************\n\n";
		cout << "\nError: there is not vector of initial conditions" << "\n";
		cout << "***************************************************************\n\n";
		toExit = true;
	}

	if(toExit)
	{
		exit(-4);
	}

	//creation output folder
	string folder = argv[2];
	string folder1 = "\"" + folder + "\"";

	DIR* dir;
	dir = opendir(folder1.c_str());

	if(dir == NULL)
	{
		if(folder1[0] == '-' || folder1[0] == ' ')
		{
			cout << "\n\n***************************************************************\n\n";
			cout << "\nError: the second parameter must be the output folder" << "\n";
			cout << "***************************************************************\n\n";
			exit(-2);
		}
		if (os == 0)
		{
			folder1 = "MD " + folder1;
			sis = system(folder1.c_str());
		}
		else
		{
			folder1 = "mkdir " + folder1;
			sis = system(folder1.c_str());
		}
	}

  	folder1 = "\"" + folder + slash + "output" + "\"";
  	dir = opendir(folder1.c_str());
	if(dir == NULL)
	{
		if (os == 0)
		{
			folder1 =  "MD " + folder1;
			sis = system(folder1.c_str());
		}
		else
		{
			folder1 =  "mkdir " + folder1;
			sis = system(folder1.c_str());
		}
	}

	string temp;
	string temp1 = "double";
	int verbose = 0;
	string prec = "-double";;
	//check optional parameters
	for(int i = 3; i < argc; i++)
	{
		temp = argv[i];
		if(strcmp(temp.c_str(), "-v") == 0)
		{
			verbose = 1;
		}
		if(strcmp(temp.c_str(), "-p") == 0)
		{
			if (i == argc-1)
			{
				cout << "\nError: you must insert the precison: double or float" << "\n";
				cout << "***************************************************************\n\n";
				exit(-5);
			}

			temp1 = argv[i+1];
			if((strcmp(temp1.c_str(), "double") != 0) && (strcmp(temp1.c_str(), "float") != 0))
			{
				cout << "\nError: the parameter must be double or float" << "\n";
				cout << "***************************************************************\n\n";
				exit(-6);
			}
			else
			{
				prec = "-" + temp1;
			}
		}
	}

	cout << "\n***************************************************************\n";
	cout << "\nLASSIE: A large-scale simulator of mass-action kinetics models\n";

	// stringstream stream;
	// stream << verbose;
	// string s;

	float timeSim = 0;

	//run simulation using single floating point precision
	if(strcmp(prec.c_str(), "-float") == 0)
	{
		Simulation<float>* sim = new Simulation<float>();
		timeSim = sim -> run(A, B, c_vector, cs_vector, t_vector, MX_0, M_feed, modelkind, folder,
				verbose, atol_vector, be_step, newton_iter, newton_tol, rkf_step, stiffness_tol, volume, temp1);
	}
	//run simulation using double floating point precision
	else
	{
		if (strcmp(prec.c_str(), "-double") == 0)
		{
			Simulation<double>* sim = new Simulation<double>();
			timeSim = sim -> run(A, B, c_vector, cs_vector, t_vector, MX_0, M_feed, modelkind, folder,
				verbose, atol_vector, be_step, newton_iter, newton_tol, rkf_step, stiffness_tol, volume, temp1);
		}
	}
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );
	float tempo = 0;
	cudaEventElapsedTime( &tempo, start, stop );
	tempo /= 1000;
	timeSim /= 1000;

	//printf("\nSimulation running time: %f seconds\n", timeSim);
	printf("Total running time: %f seconds\n", tempo);
	cout << "\n***************************************************************\n\n";


	return 0;
}

//method to read all file in a folder
__host__ vector<string> openDirectory(string path)
{
	string str;
	vector<string> files;
	DIR*    dir;
	struct dirent *pdir;

	dir = opendir(path.c_str());
	if(dir == NULL)
	{
		cout << "\n\n***************************************************************\n\n";
		cout << "\nError: the input folder does not exist" << "\n";
		cout << "***************************************************************\n\n";
		exit(-1);
	}

	while ((pdir = readdir(dir)) != NULL )
	{
		str = pdir->d_name;
		if(strcmp(str.c_str(), ".") != 0 & strcmp(str.c_str(), "..") != 0 & strcmp(str.c_str(), ".DS_Store") != 0)
			files.push_back(str);
	}
	return files;
}
