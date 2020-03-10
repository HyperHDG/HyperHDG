#include <fstream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <iostream>

#include <chrono>
#include <experimental/filesystem>

#include <HyperHDG/HyAssert.hxx>

using namespace std;
namespace fs = std::experimental::filesystem;

string hyCythonize( const vector<string>& names )
{
  hy_assert( names.size() >= 2 ,
             "The size of the names vector must be large enough for all needed compile options!" );
  
  cout << "Cythonizing " << names[1] << " ... " << flush;
  
  // Auxiliary string for name of internal Python class
  
  string python_name = names[1];
  replace( python_name.begin(), python_name.end(), '<', '_' );
  replace( python_name.begin(), python_name.end(), '>', '_' );
  replace( python_name.begin(), python_name.end(), ',', '_' );
  
  // Define names of input file and output file, and define system command to cythonize files.
  
  string infileName = "./cython/" + names[0];
  string outfileName = "./build/CythonFiles/" + python_name;
  
  string cythonCommand = "cd ./build/CythonFiles/; cython -3 --cplus " + python_name + ".pyx";
  // TODO: Make automatic adaption to latest python version. Now solved with python3.*. This gives
  // a warning.
  string compileCommand = "g++ \
    -pthread -g  -I\
    /usr/include/python3.* \
    -I. -Iinclude -fwrapv -O2 -Wall -g \
    -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 \
    -fPIC --std=c++17 \
    -c " + outfileName + ".cpp \
    -o " + outfileName + ".o";
  string linkCommand = "g++ \
    -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro \
    -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fstack-protector-strong -Wformat \
    -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 "
    + outfileName + ".o "
    + "-o build/" + python_name + ".so \
    -llapack -lstdc++fs";
  
  // Introduce auxiliary variables needed for copying files and substituting keywords.
  
  string line, word;
  std::istringstream linestream;
  ifstream infile;
  ofstream outfile;
  
  // Check whether or not file needs to be recompiled
  
  bool success = true;
  infile.open(infileName + ".pxd");
  std::getline(infile, line);
  linestream = std::istringstream(line);
  infile.close();
  
  linestream >> word; if (word != "#")     success = false;
  linestream >> word; if (word != "C++:")  success = false;
  linestream >> word;
  
  if (success)
  {
    fs::path so_file  = fs::current_path().string() + "/build/" + python_name + ".so";
    if (exists(so_file))
    {
      fs::path pyx_file = fs::current_path().string() + "/cython/" + names[0] + ".pyx";
      hy_assert( exists(pyx_file) , "File needs to exist!" );
      fs::path pxd_file = fs::current_path().string() + "/cython/" + names[0] + ".pxd";
      hy_assert( exists(pxd_file) , "File needs to exist!" );
      fs::path cxx_file = fs::current_path().string() + "/" + word;
      hy_assert( exists(cxx_file) , "File needs to exist!" );
      
      auto so_time  = fs::last_write_time(so_file);
      auto pyx_time = fs::last_write_time(pyx_file);
      auto pxd_time = fs::last_write_time(pxd_file);
      auto cxx_time = fs::last_write_time(cxx_file);
      
      std::time_t so_t = decltype(so_time)::clock::to_time_t(so_time);
      std::time_t pyx_t = decltype(pyx_time)::clock::to_time_t(pyx_time);
      std::time_t pxd_t = decltype(pxd_time)::clock::to_time_t(pxd_time);
      std::time_t cxx_t = decltype(cxx_time)::clock::to_time_t(cxx_time);
      
      if (so_t > pyx_t && so_t > pxd_t && so_t > cxx_t)
      {
        cout << " DONE without recompilation." << endl;
        return python_name;
      }
    }
  }
  
  // Copy .pxd file to build location.
  
  infile.open(infileName + ".pxd");
  outfile.open(outfileName + ".pxd");
  
  while ( std::getline(infile, line) )
  {
    for (unsigned int i = 0; i < line.size() && line[i] == ' '; ++i)  outfile << " ";
    linestream = std::istringstream(line);
    while ( !linestream.eof() && ! line.empty() )
    {
      linestream >> word;
      if (word == "C++ClassName")     word = "\"" + names[1] + "\"";
      if (word == "CythonClassName")  word = python_name + "_Cython";
      if (word == "PythonClassName")  word = python_name;
      outfile << word << " ";
    }
    outfile << endl;
  }
  
  infile.close();
  outfile.close(); 
  
  // Copy .pyx file to build location.
  
  infile.open(infileName + ".pyx");
  outfile.open(outfileName + ".pyx");
  
  while ( std::getline(infile, line) )
  {
    for (unsigned int i = 0; i < line.size() && line[i] == ' '; ++i)  outfile << " ";
    linestream = std::istringstream(line);
    while ( !linestream.eof() && ! line.empty() )
    {
      linestream >> word;
      if (word == "C++ClassName")     word = "\"" + names[1] + "\"";
      if (word == "CythonClassName")  word = python_name + "_Cython";
      if (word == "PythonClassName")  word = python_name;
      outfile << word << " ";
    }
    outfile << endl;
  }
  
  infile.close();
  outfile.close(); 
  
  // Execute system commands to cythonize, compile, and link class.
  
  system(cythonCommand.c_str());
  system(compileCommand.c_str());
  system(linkCommand.c_str());
  
  cout << " DONE with compilation." << endl;
  
  return python_name;
}


int main()
{
  hyCythonize({ "hyCythonize" , "hyCythonizer" });
  return 0;
}
