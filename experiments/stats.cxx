#include "stats.hxx"
#include <cmath>

void compute_stats(double* values, size_t n, SimpleStats* stats) {
  if (n == 0)
    return;
  stats->sum = stats->avg = stats->stddev = 0;
  stats->max = stats->min = values[0];
  for (size_t p = 0; p < n; p++) {
    stats->sum += values[p];
    stats->max = std::fmax(stats->max,values[p]);
    stats->min = std::fmin(stats->min,values[p]);
  }
  stats->avg = stats->sum / n;
  for (size_t p = 0; p < n; p++) {
    double d = values[p] - stats->avg;
    stats->stddev += d*d;
  }
  // sample variance
  stats->stddev = sqrt(1./(n-1) * stats->stddev);
}
