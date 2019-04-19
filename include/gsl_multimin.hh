#include <gsl/gsl_multimin.h>
#include <cstdio>
#include <array>

namespace ivanp {

namespace detail { namespace gsl_multimin {

template <class F, size_t... I>
decltype(auto) apply_impl(
  F&& f, const gsl_vector* args, std::index_sequence<I...>
) {
  return f( gsl_vector_get(args,I)... );
}
template <size_t N, typename F>
decltype(auto) apply(F&& f, const gsl_vector* args) {
  return apply_impl(
    std::forward<F>(f), args,
    std::make_index_sequence<N>{} );
}

template <typename T, size_t... I>
void set_start_step_impl(
  gsl_vector* start, gsl_vector* step,
  const T& start_step,
  std::index_sequence<I...>
) {
  using expander = int[];
  (void)expander{0, ((void)(
    gsl_vector_set( start, I, std::get<0>(std::get<I>(start_step)) ),
    gsl_vector_set( step , I, std::get<1>(std::get<I>(start_step)) )
  ), 0)...};
}
template <size_t N>
void set_start_step(
  gsl_vector* start, gsl_vector* step,
  const std::array<std::array<double,2>,N>& start_step
) {
  set_start_step_impl(start,step,start_step,std::make_index_sequence<N>{});
}

}}

struct gsl_multimin_opts {
  bool verbose = false;
  double tolerance = 1e-3;
  size_t max_iter = 1000;
  const gsl_multimin_fminimizer_type* min_type
    = gsl_multimin_fminimizer_nmsimplex2;
};

template <size_t N, typename F>
auto gsl_multimin (
  const std::array<std::array<double,2>,N>& start_step,
  F&& f,
  const gsl_multimin_opts& opts = { }
) -> std::array<double,N> {
  using namespace detail::gsl_multimin;

  /* Initialize method and iterate */
  gsl_multimin_function minex_func;
  minex_func.n = N;
  minex_func.params = &f;
  minex_func.f = [](const gsl_vector *v, void *p){
    return apply<N>( (*reinterpret_cast<F*>(p)), v );
  };

  gsl_vector *x = gsl_vector_alloc(N);
  gsl_vector *step = gsl_vector_alloc(N);
  set_start_step(x,step,start_step);

  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(opts.min_type,N);
  gsl_multimin_fminimizer_set(s, &minex_func, x, step);

  size_t iter = 0;
  int status;
  double size;

  do {
    ++iter;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, opts.tolerance);

    if (opts.verbose) {
      printf("%5lu f = %10.3e size = %9.3e", iter, s->fval, size);
      for (size_t i=0; i<N; ++i)
        printf(" %10.3e",gsl_vector_get(s->x,i));
      printf("\n");
    }
  } while (status == GSL_CONTINUE && iter < opts.max_iter);

  std::array<double,N> ret;
  for (size_t i=N; i; ) { --i; ret[i] = gsl_vector_get(s->x,i); }

  gsl_vector_free(x);
  gsl_vector_free(step);
  gsl_multimin_fminimizer_free(s);

  return ret;
}

}
