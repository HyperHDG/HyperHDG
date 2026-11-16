#pragma once  // Ensure that file is included only once in a single compilation.

#include <nanobench.h>

#include <filesystem>
#include <fstream>
#include <string>

namespace tests::performance
{
/*!*************************************************************************************************
 * \brief   Write the results of a benchmark as JSON to output/performance/<name>.json.
 **************************************************************************************************/
inline void save(const ankerl::nanobench::Bench& bench, const std::string& name)
{
  std::filesystem::create_directories("output/performance");
  std::ofstream file("output/performance/" + name + ".json");
  ankerl::nanobench::render(ankerl::nanobench::templates::json(), bench, file);
}
}  // namespace tests::performance
