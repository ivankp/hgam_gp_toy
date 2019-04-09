#ifndef INTEGRATION_HH
#define INTEGRATION_HH

template <typename F>
double integrate(unsigned n, double a, double b, F&& f) {
  const double d = (b-a)/n;
  double sum = f(a)/2;
  for (unsigned i=1; i<n; ++i)
    sum += f(a + i*d);
  sum += f(b)/2;
  return sum*d;
}

#endif
