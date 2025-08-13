#ifndef STATS_H
#define STATS_H

#include <cstddef>

struct SimpleStats {
  double min, max, sum, avg, stddev;
};

void compute_stats(double* values, size_t n, SimpleStats* stats);

#endif // STATS_H
