#include <iostream>
#include <vector>

#include <nlohmann/json.hpp>

#include <TMinuit.h>
#include <TFile.h>
#include <TH1.h>

#include "ivanp/enumerate.hh"
#include "ivanp/time_seed.hh"
#include "ivanp/binner.hh"
#include "ivanp/string.hh"
#include "ivanp/root/minuit.hh"

#include "wls.hh"
#include "gp.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::string;
using std::vector;
using nlohmann::json;
using namespace ivanp;

template <typename T>
const T& cref(const json& j) { return j.get_ref<const T&>(); }

using hist = binner<double,
  std::tuple<axis_spec<uniform_axis<double>,0,0>> >;

TH1D* to_root(const char* name, const hist& h) {
  const auto& axis = h.axis();
  const unsigned n = axis.nbins();
  TH1D* _h = new TH1D(name, "", n, axis.min(), axis.max());
  auto* val = dynamic_cast<TArrayD*>(_h)->GetArray() + 1;
  for (unsigned i=0; i<n; ++i)
    val[i] = h[{i}];
  return _h;
}

int main(int argc, char* argv[]) {
  json in;
  std::cin >> in;

  long seed;
  try { seed = in["seed"]; }
  catch (...) { seed = time_seed(); }
  TEST(seed)
  std::mt19937 gen(seed);

  hist h({55,105,160});

  const auto& bkg = in["bkg"];
  const auto& sig = in["sig"];

  auto bkg_f = [
    c = bkg["poly"].get<vector<double>>(),
    e = (bool) bkg["exp"]
  ](double x){
    double xn = 1, p = c[0];
    for (unsigned i=1, n=c.size(); i<n; ++i)
      p += c[i] * (xn *= x);
    return e ? std::exp(p) : p;
  };

  const double bkg_a = bkg["a"];
  const double bkg_b = bkg["b"];
  const size_t bkg_n = bkg["n"];
  // const double sig_a = sig["a"];
  // const double sig_b = sig["b"];
  // const size_t sig_n = sig["n"];

  // Monte Carlo ====================================================
  std::uniform_real_distribution<double>
    dist_x(bkg_a,bkg_b),
    dist_y(0,bkg_f(bkg_a));

  vector<double> mc(bkg_n);
  TEST(bkg_n)
  for (size_t i=0; i<bkg_n; ) {
    const double x = dist_x(gen);
    if (dist_y(gen) < bkg_f(x)) { h(mc[i] = x); ++i; }
  }

  // Weighted Least Squares =========================================
  std::array<double(*)(double),3> fs {
    [](double x){ return 1.;  },
    [](double x){ return x;   },
    [](double x){ return x*x; }
  };

  const auto& axis = h.axis();
  const auto nbins = axis.nbins();
  std::vector<double> A;
  A.reserve(fs.size()*nbins);
  for (auto& f : fs) {
    double a = axis.edge(0), b;
    for (unsigned i=0; i<nbins; ) {
      ++i;
      b = axis.edge(i);
      A.push_back(f(a+(b-a)/2));
      a = b;
    }
  }
  std::vector<double> us;
  us.reserve(nbins);
  for (double y : h) us.push_back(y==0 ? 1 : y);

  std::array<double,fs.size()> wls_cs;
  wls(
    A.data(),
    h.bins().data(), // observed values
    us.data(), // variances
    nbins, // number of values
    fs.size(), // number of parameters
    wls_cs.data() // fit coefficients
  );
  cout << "WLS:\n";
  for (auto c : wls_cs) cout << ' ' << c;
  cout << "\nNormalized coefficients:\n";
  const double norm = 1./(
      55.*wls_cs[0]
    + (14575./2.)*wls_cs[1]
    + (2938375./3.)*wls_cs[2]
  );
  for (auto& c : wls_cs) cout << ' ' << (c*=norm);
  cout << endl;

  // Likelihood fit =================================================
  const auto& fit = in["fit"];
  auto mLogL = make_minuit(2,
    [&](const double* c) -> double {
      long double logl = 0.;
      const unsigned n = mc.size();
      #pragma omp parallel for reduction(+:logl)
      for (unsigned i=0; i<n; ++i) {
        const double x = mc[i];
        const double c0 = (6./330.)-(43725./330.)*c[0]-(5876750./330.)*c[1];
        logl += std::log(c0 + c[0]*x + c[1]*x*x);
      }
      return -2.*logl;
    }
  );
  mLogL.SetPrintLevel(fit["verbose"]);
  for (unsigned i=1; i<wls_cs.size(); ++i) {
    mLogL.DefineParameter(
      i-1,
      cat('c',i).c_str(), // name
      wls_cs[i], // start
      0.1, 0, 0 // step, min, max
    );
  }
  mLogL.Migrad();

  cout << "\nLikelihood fit:\n";
  for (unsigned i=1; i<wls_cs.size(); ++i) {
    double val, err;
    mLogL.GetParameter(i-1,val,err);
    cout << 'c' << i << ": " << val << " +- " << err << '\n';
  }
  cout << endl;

  TFile fout("toy_test.root","recreate");
  to_root("first_toy",h);
  fout.Write();
}
