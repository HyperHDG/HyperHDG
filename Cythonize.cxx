#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#include <HyperHDG/HyAssert.hxx>

using namespace std;

void hyCython( vector<string> names )
{
  hy_assert( names.size() >= 2 ,
             "The size of the names vector must be large enough for all needed compile options!" );
  
  cout << "Cythonizing " << names[1] << " ... " << flush;
  
  // Define names of input file and output file, and define system command to cythonize files.
  
  string infileName = "./cython/" + names[0];
  string outfileName = "./build/CythonFiles/" + names[2];
  
  string cythonCommand = "cd ./build/CythonFiles/; cython -3 --cplus " + names[2] + ".pyx";
  // TODO: Make automatic adaption to latest python version.
  string compileCommand = "g++ \
    -pthread -g  -I\
    /usr/include/python3.6 \
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
    + "-o build/" + names[2] + ".so \
    -llapack";
  
  // Introduce auxiliary variables needed for copying files and substituting keywords.
  
  string line, word;
  std::istringstream linestream;
  ifstream infile;
  ofstream outfile;
  
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
      if (word == "CythonClassName")  word = names[2] + "_Cython";
      if (word == "PythonClassName")  word = names[2];
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
      if (word == "CythonClassName")  word = names[2] + "_Cython";
      if (word == "PythonClassName")  word = names[2];
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
  
  cout << " DONE." << endl;
}


int main()
{
  vector<string> names;
  names.push_back("Cython");
  names.push_back("xxx");
  names.push_back("importer");
  hyCython(names);
  return 0;
}
