#define HY_PYTHON_VERSION 3

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
 * \brief   Function to capture output of a shell command.
 * 
 * This function is inspired (almost copied, retrieved March 11, 2020) from 
 * https://www.jeremymorgan.com/
 * tutorials/c-programming/how-to-capture-the-output-of-a-linux-command-in-c/
 *         
 * \param   cmd       Shell command.
 * \retval  return    Return string of the shell command.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
string GetStdOutFromCommand(string cmd)
{
  string data;
  FILE *stream;
  const int max_buffer = 256;
  char buffer[max_buffer];
  cmd.append(" 2>&1");

  stream = popen(cmd.c_str(), "r");
  if (stream) {
    while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL)  data.append(buffer);
    pclose(stream);
  }
  return data;
}
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
  
  // Define names of input file and output file.
  
  string infileName = "./cython/" + names[0];
  string outfileName = "./build/CythonFiles/" + python_name;
  
  // Find used python version.
  
  int ver_start, ver_end;
  string pyVersion = GetStdOutFromCommand("python" + to_string(HY_PYTHON_VERSION) + " --version");
  ver_start = pyVersion.find(to_string(HY_PYTHON_VERSION));
  ver_end = pyVersion.find('.', ver_start);
  ver_end = pyVersion.find('.', ver_end+1);
  hy_assert( ver_start != ver_end && ver_end != -1 , "Found python version is invalid." );
  pyVersion = pyVersion.substr(ver_start, ver_end - ver_start);
  
  // Define commands to compile, cythonize, and link files.
  
  string cythonCommand = "cd ./build/CythonFiles/; cython -3 --cplus " + python_name + ".pyx";
  string compileCommand = "g++ \
    -pthread -g  -I/usr/include/python" + pyVersion + " -I. -Iinclude -fwrapv -O2 -Wall -g \
    -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 \
    -fPIC --std=c++17 -c " + outfileName + ".cpp -o " + outfileName + ".o";
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
  
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-result" 
  
  system(cythonCommand.c_str());
  system(compileCommand.c_str());  
  system(linkCommand.c_str());
  
  #pragma GCC diagnostic pop 
  
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
  hyCythonize({ "hyCythonize" , "hyCythonizer" });
  return 0;
}
