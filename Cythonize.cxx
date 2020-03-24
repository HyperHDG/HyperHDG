#ifndef PYVERMAJ
#define PYVERMAJ -1
#endif

#ifndef PYVERMIN
#define PYVERMIN -1
#endif


#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

// #include <chrono>
#include <experimental/filesystem>

#include <HyperHDG/HyAssert.hxx>

using namespace std;
namespace fs = experimental::filesystem;


/*!*************************************************************************************************
 * \brief   Function serving as C++ counterpart of just-in-time compilation.
 * 
 * In HyperHDG, a Python script can be executed using the libraries functions (which have been
 * implemented in C++ utilizing templates). Thus, we provide a just-in-time compilation framework
 * that allows to use the scripts without prior compilation of the used classes. In contrast, we
 * provide our own import function \c hyImport (cf. hyImport.py), which is the Python counterpart
 * to this function hyCythonize.
 * 
 * If hyImport is called with a vector of strings (file name and class name plus possible template
 * specifications), it uses Cython based interfaces to call hyCythnoize which itself builds
 * Cython based interfaces and compiles the respective files, i.e. a .hxx/.cxx, a .pxd, and a .pyx
 * file to a .so file which is needed by the overall Python script.
 *        
 * \param   names     Vector containing specifying names.
 * \retval  name      Name associated to created .so file.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
string hyCythonize
( 
  const vector<string>& names, const vector<string>& filenames,
  const unsigned int ver_maj, const unsigned int ver_min
)
{
  static_assert( PYVERMAJ != -1 && PYVERMIN != -1 ,
                 "Python verion needs to be set as compile flags!" );

  hy_assert( ver_maj == PYVERMAJ && ver_min == PYVERMIN ,
             "The Python version, this program is executed under needs to be the same as the one "
             << "the program has been compiled with. This seems not to be the case." );

  hy_assert( names.size() >= 2 ,
             "The size of the names vector must be large enough for all needed compile options!" );
  
  cout << "Cythonizing " << names[1] << " ... " << flush;
  
  // Auxiliary string for name of internal Python class
  
  string python_name = names[1];
  python_name.erase(std::remove(python_name.begin(), python_name.end(), ' '), python_name.end());
  replace( python_name.begin(), python_name.end(), '+', '_' );
  replace( python_name.begin(), python_name.end(), '-', '_' );
  replace( python_name.begin(), python_name.end(), '*', '_' );
  replace( python_name.begin(), python_name.end(), ':', '_' );
  replace( python_name.begin(), python_name.end(), '<', '_' );
  replace( python_name.begin(), python_name.end(), '>', '_' );
  replace( python_name.begin(), python_name.end(), ',', '_' );
  
  // Define names of input file and output file.
  
  string infileName = "./cython/" + names[0];
  string outfileName = "./build/CythonFiles/" + python_name;
  
  // String defining current Python version.

  string pyVersion = to_string(PYVERMAJ) + "." + to_string(PYVERMIN);
  
  // Define commands to compile, cythonize, and link files.
  
  string cythonCommand = "cd ./build/CythonFiles/; cython -3 --cplus " + python_name + ".pyx";
  string compileCommand = "g++ \
    -pthread -g  -I/usr/include/python" + pyVersion + " -I. -Iinclude -fwrapv -O2 -Wall -g \
    -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 \
    -fPIC --std=c++17 -DPYVERMAJ=" + to_string(PYVERMAJ) + " -DPYVERMIN=" + to_string(PYVERMIN) +
    " -c " + outfileName + ".cpp -o " + outfileName + ".o";
  string linkCommand = "g++ \
    -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro \
    -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fstack-protector-strong -Wformat \
    -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 "
    + outfileName + ".o " + "-o build/SharedObjects/" + python_name + ".so \
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
    fs::path so_file  = fs::current_path().string() + "/build/SharedObjects/" + python_name + ".so";
    if (exists(so_file))
    {
      fs::path pyx_file = fs::current_path().string() + "/cython/" + names[0] + ".pyx";
      hy_assert( exists(pyx_file) , "File needs to exist!" );
      fs::path pxd_file = fs::current_path().string() + "/cython/" + names[0] + ".pxd";
      hy_assert( exists(pxd_file) , "File needs to exist!" );
      fs::path cxx_file = fs::current_path().string() + "/" + word;
      if ( !exists(cxx_file) )  cxx_file = fs::current_path().string() + "/include/" + word;
      hy_assert( exists(cxx_file) , "File needs to exist!" );
      
      auto so_time  = fs::last_write_time(so_file);
      auto pyx_time = fs::last_write_time(pyx_file);
      auto pxd_time = fs::last_write_time(pxd_file);
      auto cxx_time = fs::last_write_time(cxx_file);
/*      
      std::time_t so_t = decltype(so_time)::clock::to_time_t(so_time);
      std::time_t pyx_t = decltype(pyx_time)::clock::to_time_t(pyx_time);
      std::time_t pxd_t = decltype(pxd_time)::clock::to_time_t(pxd_time);
      std::time_t cxx_t = decltype(cxx_time)::clock::to_time_t(cxx_time);
*/      
      if (so_time > pyx_time && so_time > pxd_time && so_time > cxx_time)
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
      if      (word == "C++ClassName")     word = "\"" + names[1] + "\"";
      else if (word == "CythonClassName")  word = python_name + "_Cython";
      else if (word == "PythonClassName")  word = python_name;
      if (word == "IncludeFiles")
        for (unsigned int i = 0; i < filenames.size(); ++i)
          outfile << "cdef extern from \"<" << filenames[i] << ">\" : pass" << endl;
      else  outfile << word << " ";
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
  
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-result" 
  
  system(cythonCommand.c_str());
  system(compileCommand.c_str());  
  system(linkCommand.c_str());
  
  #pragma GCC diagnostic pop
  
  // Finish program.
  
  cout << " DONE with compilation." << endl;
  
  return python_name;
}
/*!*************************************************************************************************
 * \brief   Main function builds hyCythonize function using itself.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
int main()
{
  hyCythonize({ "hyCythonize" , "hyCythonizer" } , { } , PYVERMAJ , PYVERMIN );
  return 0;
}
