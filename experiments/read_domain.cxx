#include <HyperHDG/read_domain.hxx>
#include <HyperHDG/dense_la.hxx>
#include <print>
#include <format>
#include <vector>
#include <cstdint>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::println(stderr, "ERROR: usage: {} <network>", argv[0]);
    return 1;
  }

  auto info = read_domain_geobin<1, 3,  std::vector, Point<3, double>, uint32_t, uint32_t, uint32_t>(argv[1]);

  std::println("{} {} {}", info.points.size(), info.hyNodes_hyEdge.size(), info.hyEdge_properties.size());
}
