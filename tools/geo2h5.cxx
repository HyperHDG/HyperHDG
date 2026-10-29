#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/read_domain.hxx>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Convert a .geo domain file into a .geo.h5 file (cf. include/HyperHDG/read_domain_hdf5.hxx):
//
//     geo2h5 domains/cross.geo [domains/cross.geo.h5] [--no-vtkhdf]
//
// The output name defaults to the input name plus ".h5". The dimensions are template parameters
// of read_domain, so they are read from the header first and dispatched over the supported
// combinations below.

static unsigned int header_value(const std::string& filename, const std::string& keyword)
{
  std::ifstream file(filename);
  hy_check(file.is_open(), "Cannot open " << filename << ".");
  std::string line, word, equal_sign;
  unsigned int value = 0;
  while (std::getline(file, line))
  {
    std::istringstream stream(line);
    if (stream >> word && word == keyword && stream >> equal_sign >> value && equal_sign == "=")
      return value;
  }
  hy_check(false, "The keyword " << keyword << " has not been found in " << filename << ".");
  return 0;
}

template <unsigned int hyEdge_dim, unsigned int space_dim>
static void convert(const std::string& input, const std::string& output, const bool vtkhdf)
{
  const auto domain = read_domain<hyEdge_dim, space_dim>(input);
  write_domain_hdf5(output, domain, vtkhdf);
  std::cout << input << " -> " << output << ": " << domain.n_points << " points, "
            << domain.n_hyNodes << " hypernodes, " << domain.n_hyEdges << " hyperedges of "
            << "dimension " << hyEdge_dim << " in " << space_dim << "d, " << domain.n_properties
            << " properties" << (vtkhdf && hyEdge_dim <= 3 && space_dim <= 3 ? ", VTKHDF view" : "")
            << std::endl;
}

int main(int argc, char** argv)
{
  std::string input, output;
  bool vtkhdf = true;
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if (arg == "--no-vtkhdf")
      vtkhdf = false;
    else if (input.empty())
      input = arg;
    else if (output.empty())
      output = arg;
    else
      input.clear(), i = argc;  // too many arguments
  }
  if (input.empty())
  {
    std::cerr << "usage: geo2h5 <domain.geo> [<domain.geo.h5>] [--no-vtkhdf]" << std::endl;
    return 1;
  }
  if (output.empty())
    output = input + ".h5";

  const unsigned int hyEdge_dim = header_value(input, "HyperEdge_Dim");
  const unsigned int space_dim = header_value(input, "Space_Dim");
  switch (hyEdge_dim * 10 + space_dim)
  {
    case 11:
      convert<1, 1>(input, output, vtkhdf);
      break;
    case 12:
      convert<1, 2>(input, output, vtkhdf);
      break;
    case 13:
      convert<1, 3>(input, output, vtkhdf);
      break;
    case 14:
      convert<1, 4>(input, output, vtkhdf);
      break;
    case 15:
      convert<1, 5>(input, output, vtkhdf);
      break;
    case 22:
      convert<2, 2>(input, output, vtkhdf);
      break;
    case 23:
      convert<2, 3>(input, output, vtkhdf);
      break;
    case 33:
      convert<3, 3>(input, output, vtkhdf);
      break;
    default:
      std::cerr
        << "geo2h5: hyperedge dimension " << hyEdge_dim << " in " << space_dim
        << " space dimensions is not compiled in (add it to the switch in tools/geo2h5.cxx)."
        << std::endl;
      return 1;
  }
  return 0;
}
