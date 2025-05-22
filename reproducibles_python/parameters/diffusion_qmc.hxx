#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <cmath>
#include <array>
#include <algorithm>

#include <iostream>

template <unsigned int space_dimT, typename param_float_t = double>
struct TestParametersDiffusionAffine
{
 typedef std::vector<param_float_t> param_time_t;
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t arr)
  {
    return 1. / inverse_diffusion_coeff(point, arr);
  }
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                               const param_time_t arr)
  {
    param_float_t r = 5.;
    unsigned int kl[] = {1, 1};
    std::for_each(arr.begin(), arr.end(), [&r, &kl, point](const param_float_t& n) {
      param_float_t w = pow(kl[0] * kl[0] + kl[1] * kl[1], -1.3);
      w *= sin(kl[0] * M_PI * point[0]) * sin(kl[1] * M_PI * point[1]);
      r += n * w;
      //std::cout << kl[0] << " " << kl[1] << " " << kl[0] * kl[0] + kl[1] * kl[1] << " " << n << std::endl;
      akt(kl);
    });
    return r; 
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t = 0)
  {
    return point[0];
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t mean_indicator(const Point<space_dimT, param_float_t>& p,
                                       const param_time_t = 0)
  {
    if(p[0] >= 0.2 && p[0] <= 0.8) {
      if(p[1] >= 0.2 && p[1] <= 0.8)
        return 1.;
    }
    return 0.;
  }
private:
  static void akt(unsigned int kl[2]) {
    if(kl[0] == 1) {
      kl[0] = kl[1] + 1;
      kl[1] = 1;
    } else {
      kl[0] -= 1;
      kl[1] += 1;
    }
  }
};

template <unsigned int space_dimT, typename param_float_t = double>
struct TestParametersDiffusionLognormal
{
 typedef std::vector<param_float_t> param_time_t;
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t arr)
  {
    return 1. / inverse_diffusion_coeff(point, arr);
  }
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                               const param_time_t arr)
  {
    param_float_t r = 0.;
    unsigned int kl[] = {1, 1};
    std::for_each(arr.begin(), arr.end(), [&r, &kl, point](const param_float_t& n) {
      param_float_t w = pow(kl[0] * kl[0] + kl[1] * kl[1], -1.3);
      w *= sin(kl[0] * M_PI * point[0]) * sin(kl[1] * M_PI * point[1]);
      r += n * w;
      //std::cout << kl[0] << " " << kl[1] << " " << kl[0] * kl[0] + kl[1] * kl[1] << " " << n << std::endl;
      akt(kl);
    });
    return exp(r); 
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t = 0)
  {
    return point[0];
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t mean_indicator(const Point<space_dimT, param_float_t>& p,
                                       const param_time_t = 0)
  {
    if(p[0] >= 0.2 && p[0] <= 0.8) {
      if(p[1] >= 0.2 && p[1] <= 0.8)
        return 1.;
    }
    return 0.;
  }
private:
  static void akt(unsigned int kl[2]) {
    if(kl[0] == 1) {
      kl[0] = kl[1] + 1;
      kl[1] = 1;
    } else {
      kl[0] -= 1;
      kl[1] += 1;
    }
  }
};

template <unsigned int space_dimT, typename param_float_t = double>
struct TestParametersDiffusionGevreyAffine
{
 typedef std::vector<param_float_t> param_time_t;
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t arr)
  {
    return 1. / inverse_diffusion_coeff(point, arr);
  }
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                               const param_time_t arr)
  {
    param_float_t r = 5.;
    unsigned int kl[] = {1, 1};
    std::for_each(arr.begin(), arr.end(), [&r, &kl, point](const param_float_t& n) {
      param_float_t w = pow(kl[0] * kl[0] + kl[1] * kl[1], -1.3);
      w *= sin(kl[0] * M_PI * point[0]) * sin(kl[1] * M_PI * point[1]);
      r += zeta(n) * w;
      //std::cout << kl[0] << " " << kl[1] << " " << kl[0] * kl[0] + kl[1] * kl[1] << " " << n << std::endl;
      akt(kl);
    });
    return r; 
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t = 0)
  {
    return point[0];
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t mean_indicator(const Point<space_dimT, param_float_t>& p,
                                       const param_time_t = 0)
  {
    if(p[0] >= 0.2 && p[0] <= 0.8) {
      if(p[1] >= 0.2 && p[1] <= 0.8)
        return 1.;
    }
    return 0.;
  }
private:
  static constexpr param_float_t omega = 1. / (1.25 - 1.);
  static param_float_t zeta(const param_float_t t) {
    return exp(- pow(t + 0.5, omega));
  }
  static void akt(unsigned int kl[2]) {
    if(kl[0] == 1) {
      kl[0] = kl[1] + 1;
      kl[1] = 1;
    } else {
      kl[0] -= 1;
      kl[1] += 1;
    }
  }
};

template <unsigned int space_dimT, typename param_float_t = double>
struct TestParametersDiffusionGevreyLognormal
{
 typedef std::vector<param_float_t> param_time_t;
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t arr)
  {
    return 1. / inverse_diffusion_coeff(point, arr);
  }
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                               const param_time_t arr)
  {
    param_float_t r = 0.;
    unsigned int kl[] = {1, 1};
    std::for_each(arr.begin(), arr.end(), [&r, &kl, point](const param_float_t& n) {
      param_float_t w = pow(kl[0] * kl[0] + kl[1] * kl[1], -1.3);
      w *= sin(kl[0] * M_PI * point[0]) * sin(kl[1] * M_PI * point[1]);
      r += n * w;
      //std::cout << kl[0] << " " << kl[1] << " " << kl[0] * kl[0] + kl[1] * kl[1] << " " << n << std::endl;
      akt(kl);
    });
    return exp(r); 
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t = 0)
  {
    return point[0];
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_time_t = 0)
  {
    return 0.;
  }
  static param_float_t mean_indicator(const Point<space_dimT, param_float_t>& p,
                                       const param_time_t = 0)
  {
    if(p[0] >= 0.2 && p[0] <= 0.8) {
      if(p[1] >= 0.2 && p[1] <= 0.8)
        return 1.;
    }
    return 0.;
  }
private:
  static constexpr param_float_t omega = 1. / (1.25 - 1.);
  static int sign(const param_float_t t) {
    if(t == 0)
      return 0;
    else if(t > 0)
      return 1;
    else
      return -1;
  }
  static param_float_t zeta(const param_float_t t) {
    param_float_t r = exp( - pow(t, -omega));
    return r * sign(t);
  }
  static void akt(unsigned int kl[2]) {
    if(kl[0] == 1) {
      kl[0] = kl[1] + 1;
      kl[1] = 1;
    } else {
      kl[0] -= 1;
      kl[1] += 1;
    }
  }
};
