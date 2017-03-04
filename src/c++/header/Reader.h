#ifndef RK4SOLVER_H
#define RK4SOLVER_H
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstring>
using namespace std;

template <class T>
class Reader
{
public:

  Reader()
  {

  }

  int readType(string file)
  {
    f = file;
    string line;
    ifstream infile(f.c_str());
    getline(infile, line);
    istringstream buf(line);
    if(strcmp(line.c_str(), "deterministic") == 0 || strcmp(line.c_str(), "concentration") == 0)
      return 0;
    else if(strcmp(line.c_str(), "stochastic") == 0 || strcmp(line.c_str(), "number") == 0)
      return 1;
    else
      return 2;
  }

  vector<T> readByLine(string file)
  {
    vector<T> v;
    f = file;
    string line;
    ifstream infile(f.c_str());
    while (getline(infile, line))
    {
      istringstream buf(line);
      T x = strtod(line.c_str(), NULL);
      v.push_back(x);
    }
    return v;
  }

  vector<T> readByTab(string file)
  {
    vector<T> v;
    f = file;
    string line;
    ifstream infile(f.c_str());
    while (getline(infile, line))
    {
      istringstream buf(line);
      for(string token; getline(buf, token, '\t'); )
      {
        T x = strtod(token.c_str(), NULL);
        v.push_back(x);
      }
    }
    return v;
  }

  vector<short> readByTabshort(string file)
  {
    vector<short> v;
    f = file;
    string line;
    ifstream infile(f.c_str());
    while (getline(infile, line))
    {
      istringstream buf(line);
      for(string token; getline(buf, token, '\t'); )
      {
        T x = (int) strtod(token.c_str(), NULL);
        v.push_back(x);
      }
    }
    return v;
  }

  vector<int> readByLineshort(string file)
  {
    vector<int> v;
    f = file;
    string line;
    ifstream infile(f.c_str());
    while (getline(infile, line))
    {
      istringstream buf(line);
      T x = (int) strtod(line.c_str(), NULL);
      v.push_back(x);
    }
    return v;
  }

  vector<vector<short> > readMatrix(string file)
  {
    vector<vector<short> > matrix;
    f = file;
    string line;
    ifstream infile(f.c_str());
    vector<short> v1;
    while (getline(infile, line))
    {
      v1.clear();
      istringstream buf(line);
      for(string token; getline(buf, token, '\t'); )
      {
        int x = (int) strtod(token.c_str(), NULL);
        v1.push_back(x);
      }
      matrix.push_back(v1);
    }
    return matrix;

  }

  T readSingleValue(string file)
  {
    f = file;
    string line;
    ifstream infile(f.c_str());
    getline(infile, line);
    istringstream buf(line);
    return strtod(line.c_str(), NULL);
  }



private:
  string f;
};
#endif
