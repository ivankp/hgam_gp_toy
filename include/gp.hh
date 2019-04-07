#ifndef GP_HH
#define GP_HH

#include <array>
#include <vector>
#include <iterator>
#include <type_traits>
#include <cmath>

#include "linalg.hh"

// Gaussian Process =================================================

template <typename Xs, typename Ys, typename Us, typename Ts, typename Kernel>
std::vector<std::array<double,2>> GP(
  const Xs& xs, // training points coordinates
  const Ys& ys, // training points values
  const Us& us, // noise variances (add to diagonal)
  const Ts& ts, // test points
  Kernel&& kernel // kernel function
) {
  using std::distance;
  using std::next;
  using namespace linalg;

  const auto x_begin = begin(xs);
  const auto x_end = end(xs);
  const auto nx = distance(x_begin, x_end);
  using nx_t = std::remove_const_t<decltype(nx)>;
  const auto N = utn(nx);

  double* const L = new double[N+nx+nx];
  double* const y = L + N;
  double* const k = y + nx;

  { double* k = L;
    auto u = begin(us);
    // compute covariance matrix K
    for (auto a=x_begin; a!=x_end; ++a) {
      for (auto b=x_begin; ; ++b, ++k) {
        *k = kernel(*a,*b);
        if (b==a) { *k += *u; ++u; ++k; break; }
      }
    }
  }

  cholesky(L,N); // K = L L

  // mean = k* (LL)^-1 y

  { auto y_ = begin(ys);
    for (nx_t i=0; i<nx; ++i) {
      y[i] = *y_;
      ++y_;
    }
  }

  solve_triang(L,y,nx); // Solve L^-1 y

  const auto t_begin = begin(ts);
  const auto t_end = end(ts);
  const auto nt = distance(t_begin, t_end);
  using nt_t = std::remove_const_t<decltype(nt)>;

  std::vector<std::array<double,2>> out(nt);

  { auto t = t_begin;
    for (nt_t i=0; i<nt; ++i) {
      for (nx_t j=0; j<nx; ++j)
        *(k+j) = kernel( *t, *next(x_begin,j) );

      solve_triang(L,k,nx); // Solve k* L^-1

      out[i] = {
        // mean = (k* L^-1) (L^-1 y)
        dot(k, y, nx),
        // var = k** - (k* L^-1) (L^-1 k*)
        std::sqrt(kernel(*t,*t) - dot(k, k, nx))
      };
      ++t;
    }
  }

  delete[] L;

  return out;
}

// log margimal likelihood
// this function is intended for optimization of kernel parameters
// the returned value is -likelihood and ignores constant terms
template <typename Xs, typename Ys, typename Us, typename Kernel,
          typename... H>
double gp_logml_opt(
  const Xs& xs, // training points coordinates
  const Ys& ys, // training points values
  const Us& us, // noise variances (add to diagonal)
  Kernel&& kernel, // kernel function
  const H&... hs // kernel (hyper)parameters
) {
  using std::distance;
  using std::next;
  using namespace linalg;

  const auto x_begin = begin(xs);
  const auto x_end = end(xs);
  const auto nx = distance(x_begin, x_end);
  using nx_t = std::remove_const_t<decltype(nx)>;
  const auto N = utn(nx);

  double* const L = new double[N+nx];
  double* const y = L + N;

  { double* k = L;
    auto u = begin(us);
    // compute covariance matrix K
    for (auto a=x_begin; a!=x_end; ++a) {
      for (auto b=x_begin; ; ++b, ++k) {
        *k = kernel(*a,*b,hs...);
        if (b==a) { *k += *u; ++u; ++k; break; }
      }
    }
  }

  cholesky(L,N); // K = L L

  double log_det = 0;
  for (unsigned i=0, j=0; i<nx; j += ++i+1) {
    log_det += std::log(L[j]);
  }

  { auto y_ = begin(ys);
    for (nx_t i=0; i<nx; ++i) {
      y[i] = *y_;
      ++y_;
    }
  }

  solve_triang(L,y,nx); // Solve L^-1 y

  return 0.5*dot(y,y,nx) + log_det;
}

#endif
