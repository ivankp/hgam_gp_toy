#include <iostream>
#include <vector>
#include <cmath>

#include <nlohmann/json.hpp>

#include <TMinuit.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>

#include "ivanp/time_seed.hh"
#include "ivanp/binner.hh"
#include "ivanp/string.hh"
#include "ivanp/root/minuit.hh"

#include "wls.hh"
#include "gp.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

constexpr double sqrt_pi_2 = std::sqrt(M_PI_2);
constexpr double sqrt_1_2 = M_SQRT1_2;

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
  const auto& bins = h.bins();
  for (unsigned i=0; i<n; ++i)
    val[i] = bins[i];
  return _h;
}

template <typename F, size_t Np=0>
TF1* mkfcn(const char* name, F&& f, const std::array<double,Np>& p = {}) {
  TF1 *_f = new TF1( name, std::forward<F>(f), 105, 160, p.size() );
  for (size_t i=0; i<p.size(); ++i)
    _f->SetParameter(i,p[i]);
  _f->SetNpx(1000);
  _f->Write();
  return _f;
}

int main(int argc, char* argv[]) {
  json in;
  std::cin >> in;

  long seed;
  try { seed = in["seed"]; }
  catch (...) { seed = time_seed(); }
  TEST(seed)
  std::mt19937 gen(seed);

  const std::array<double,2> range = in["range"];

  const auto& bkg = in["bkg"];
  const size_t bkg_n = bkg["n"];

  const auto& sig = in["sig"];
  const size_t sig_n = sig["n"];

  const auto& fit = in["fit"];

  // Monte Carlo ====================================================
  vector<double> mc(bkg_n+sig_n);
  hist h({fit["nbins"],range[0],range[1]});

  // background -----------------------------------------------------
  auto gen_bkg = [ // polynomial or exp(poly)
    c = bkg["poly"].get<vector<double>>(),
    e = (bool) bkg["exp"]
  ](double x){
    double xn = 1, p = c[0];
    for (unsigned i=1, n=c.size(); i<n; ++i)
      p += c[i] * (xn *= x);
    return e ? std::exp(p) : p;
  };

  TEST(bkg_n)
  { std::uniform_real_distribution<double>
      dist_x(range[0],range[1]),
      dist_y(0,gen_bkg(range[0]));
    for (size_t i=0; i<bkg_n; ) {
      const double x = dist_x(gen);
      if (dist_y(gen) < gen_bkg(x)) { h(mc[i] = x); ++i; }
    }
  }

  // signal ---------------------------------------------------------
  struct DSCB_params_t {
    double mu, s, aL, nL, rL, eL, aH, nH, rH, eH, norm;
    DSCB_params_t(const json& sig, const std::array<double,2>& range)
    : mu(sig.at("muCB")), s(sig.at("sCB")),
      aL(sig.at("aLow" )), nL(sig.at("nLow" )), rL(nL/aL), eL(std::exp(-0.5*aL*aL)),
      aH(sig.at("aHigh")), nH(sig.at("nHigh")), rH(nH/aH), eH(std::exp(-0.5*aH*aH))
    {
      // norm = 1. / ( // infinite tails
      //     eL * rL / (nL-1)
      //   + eH * rH / (nH-1)
      //   + sqrt_pi_2 * (std::erf(aH*sqrt_1_2) - std::erf(aL*sqrt_1_2))
      // );

      const double gL = (mu-range[0])/s + rL - aL;
      const double gH = (range[1]-mu)/s + rL - aL;
      norm = 1. / (
          eL * (rL - std::pow(rL/gL,nL)*gL) / (nL-1)
        + eH * (rH - std::pow(rH/gH,nH)*gH) / (nH-1)
        + sqrt_pi_2 * (std::erf(aH*sqrt_1_2) - std::erf(aL*sqrt_1_2))
      );
    }
  } DSCB_params(sig,range);

  auto gen_sig = [ // Double-sided Crystal Ball
    &cb = DSCB_params
  ](double x) {
    const double t = (x-cb.mu)/cb.s;
    if (t < -cb.aL)
      return cb.eL * std::pow( (cb.rL - cb.aL - t)/cb.rL, -cb.nL );
    if (t > cb.aH)
      return cb.eH * std::pow( (cb.rH - cb.aH + t)/cb.rH, -cb.nH );
    return std::exp(-0.5*t*t);
  };

  TEST(sig_n)
  { std::uniform_real_distribution<double>
      dist_x(range[0],range[1]),
      dist_y(0,1);
    for (size_t i=bkg_n, n=bkg_n+sig_n; i<n; ) {
      const double x = dist_x(gen);
      if (dist_y(gen) < gen_sig(x)) { h(mc[i] = x); ++i; }
    }
  }

  // Weighted Least Squares =========================================
  const auto& axis = h.axis();
  const auto nbins = axis.nbins();
  constexpr size_t nfs = 3;
  std::array<double,nfs> wls_cs;
  {
    std::vector<double> A;
    A.reserve(nfs*nbins);
    for (auto& f : { // unary + decays lambdas to function pointers
      +[](double x){ return 1.;  },
      +[](double x){ return x;   },
      +[](double x){ return x*x; }
    }) {
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

    wls(
      A.data(),
      h.bins().data(), // observed values
      us.data(), // variances
      nbins, // number of values
      nfs, // number of parameters
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
  }

  // Likelihood fit =================================================
  { auto mLogL = make_minuit(2,
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

    cout << "\nLikelihood fit (background only):\n";
    for (unsigned i=0, n=mLogL.GetNumPars(); i<n; ++i) {
      double val, err;
      mLogL.GetParameter(i,val,err);
      cout << val << " +- " << err << '\n';
    }
    cout << endl;
  }

  auto fit_sig = [ // Normalized Double-sided Crystal Ball
    &cb = DSCB_params
  ](double x) {
    const double t = (x-cb.mu)/cb.s;
    double f;
    if (t < -cb.aL)
      f = cb.eL * std::pow( (cb.rL - cb.aL - t)/cb.rL, -cb.nL );
    else if (t > cb.aH)
      f = cb.eH * std::pow( (cb.rH - cb.aH + t)/cb.rH, -cb.nH );
    else
      f = std::exp(-0.5*t*t);
    return f * cb.norm;
  };
  auto fit_fcn = [&](double x, const double* c) {
    const double c0 = (6./330.)-(43725./330.)*c[1]-(5876750./330.)*c[2];
    return (1.-c[0]) * (c0 + c[1]*x + c[2]*x*x) + c[0] * fit_sig(x);
  };
  std::array<double,3> fit_c;
  { auto mLogL = make_minuit(3,
      [&](const double* c) -> double {
        long double logl = 0.;
        const unsigned n = mc.size();
        #pragma omp parallel for reduction(+:logl)
        for (unsigned i=0; i<n; ++i) {
          const double x = mc[i];
          logl += std::log( fit_fcn(x,c) );
        }
        return -2.*logl;
      }
    );
    mLogL.SetPrintLevel(fit["verbose"]);
    mLogL.DefineParameter(
      0, "sig", 0, // num, name, start
      0.01, 0, 0 // step, min, max
    );
    for (unsigned i=1; i<wls_cs.size(); ++i) {
      mLogL.DefineParameter(
        i,
        cat('c',i).c_str(), // name
        wls_cs[i], // start
        0.1, 0, 0 // step, min, max
      );
    }
    // mLogL.FixParameter(0);
    mLogL.Migrad();

    cout << "\nLikelihood fit (signal + background):\n";
    for (unsigned i=0, n=mLogL.GetNumPars(); i<n; ++i) {
      double val, err;
      mLogL.GetParameter(i,val,err);
      fit_c[i] = val;
      cout << val << " +- " << err << '\n';
    }
    cout << endl;
  }

  // Save output ====================================================
  TFile fout("toy_test.root","recreate");
  to_root("first_toy",h);

  mkfcn("gen_bkg",[&](double* x, double*){ return gen_bkg(*x); });
  mkfcn("gen_sig",[&](double* x, double*){ return gen_sig(*x); });
  mkfcn("fit",[&](double* x, double* p){ return fit_fcn(*x,p); },fit_c);

  fout.Write();
}
