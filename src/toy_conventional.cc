#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <boost/optional.hpp>

#include <nlohmann/json.hpp>
#include "json/boost_optional.hh"

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>

#include "ivanp/time_seed.hh"
#include "ivanp/binner.hh"
#include "ivanp/string.hh"
#include "ivanp/timed_counter.hh"
#include "ivanp/root/minuit.hh"
#include "ivanp/math/running_stats.hh"

#include "wls.hh"
#include "gp.hh"
#include "integration.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

constexpr double sqrt_pi_2 = std::sqrt(M_PI_2);
constexpr double sqrt_1_2 = M_SQRT1_2;

using std::cout;
using std::endl;
using std::get;
using std::string;
using std::vector;
using nlohmann::json;
using namespace ivanp;
using linalg::sq;

template <typename T>
const T& cref(const json& j) { return j.get_ref<const T&>(); }

const char* cstr(const json& j) {
  return j.get_ptr<const json::string_t*>()->c_str();
}

using hist = binner<double,
  std::tuple<axis_spec<uniform_axis<double>,0,0>> >;

namespace nlohmann {
template <> struct adl_serializer<hist::axis_type<0>> {
  using axis = hist::axis_type<0>;
  static void from_json(const json& j, axis& a) {
    a = { j["nbins"], j["min"], j["max"] };
  }
};
}

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

  TFile fout(cstr(in["output"]),"recreate");

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

  const size_t tot_n = bkg_n + sig_n;
  const double sig_frac = double(sig_n)/double(tot_n);
  TEST(tot_n)
  TEST(sig_frac)

  // Monte Carlo ====================================================
  const unsigned nbins = in["nbins"];
  hist h({nbins,range[0],range[1]});
  hist sig_yield(in["yield_hist"]);
  running_stats yield_stats;
  vector<double> mc(tot_n), us(nbins);

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
  mkfcn("gen_bkg",[&](double* x, double*){ return gen_bkg(*x); });
  TEST(bkg_n)

  // signal ---------------------------------------------------------
  struct DSCB_t { // Double-sided Crystal Ball
    double mu, s, aL, nL, rL, eL, aH, nH, rH, eH, norm;
    DSCB_t(const json& sig, const std::array<double,2>& range)
    : mu(sig.at("muCB")), s(sig.at("sCB")),
      aL(sig.at("aLow" )), nL(sig.at("nLow" )), rL(nL/aL), eL(std::exp(-0.5*aL*aL)),
      aH(sig.at("aHigh")), nH(sig.at("nHigh")), rH(nH/aH), eH(std::exp(-0.5*aH*aH)),
      norm(1./integrate(200,range[0],range[1],*this))
    { }
    double operator()(double x) const {
      const double t = (x-mu)/s;
      if (t < -aL)
        return eL * std::pow( (rL - aL - t)/rL, -nL );
      if (t > aH)
        return eH * std::pow( (rH - aH + t)/rH, -nH );
      return std::exp(-0.5*t*t);
    }
  } const cb(sig,range);
  mkfcn("gen_sig",[&](double* x, double*){ return cb(*x); });
  TEST(sig_n)

  // Fit functions --------------------------------------------------
  auto fit_bkg = [](double x, const double* c) {
    const double c0 = (6./330.)-(43725./330.)*c[0]-(5876750./330.)*c[1];
    return c0 + c[0]*x + c[1]*x*x;
  };
  auto fit_bkg_sig = [&](double x, const double* c) {
    return (1.-c[0]) * fit_bkg(x,c+1) + c[0] * cb(x) * cb.norm;
  };

  // TOYS ***********************************************************
  for (timed_counter<unsigned> toy(in["ntoys"]); !!toy; ++toy) {
    { std::uniform_real_distribution<double>
        dist_x(range[0],range[1]),
        dist_y(0,gen_bkg(range[0]));
      for (size_t i=0; i<bkg_n; ) {
        const double x = dist_x(gen);
        if (dist_y(gen) < gen_bkg(x)) h(mc[i++] = x);
      }
    }

    { std::uniform_real_distribution<double>
        dist_x(range[0],range[1]),
        dist_y(0,1);
      for (size_t i=bkg_n, n=bkg_n+sig_n; i<n; ) {
        const double x = dist_x(gen);
        if (dist_y(gen) < cb(x)) h(mc[i++] = x);
      }
    }

    vector<double> bin_centers(nbins);
    { double a = h.axis().edge(0), b;
      for (unsigned i=0; i<nbins; ) {
        ++i;
        b = h.axis().edge(i);
        bin_centers[i-1] = a + (b-a)/2;
        a = b;
      }
    }

    // Weighted Least Squares =========================================
    constexpr size_t nfs = 3; // number of function
    std::array<double,nfs> wls_cs;

    std::vector<double> A;
    A.reserve(nfs*nbins);
    for (auto& f : { // unary + decays lambdas to function pointers
      +[](double x){ return 1.;  },
      +[](double x){ return x;   },
      +[](double x){ return x*x; }
    }) {
      for (unsigned i=0; i<nbins; ++i)
        A.push_back(f(bin_centers[i]));
    }

    for (unsigned i=0; i<nbins; ++i)
      us[i] = h[{i}] ?: 1;

    wls(
      A.data(),
      h.bins().data(), // observed values
      us.data(), // variances
      nbins, // number of values
      nfs, // number of parameters
      wls_cs.data() // fit coefficients
    );
    const double norm = 1./( // Normalize coefficients
        55.*wls_cs[0]
      + (14575./2.)*wls_cs[1]
      + (2938375./3.)*wls_cs[2]
    );
    for (auto& c : wls_cs) c*=norm;

    // Likelihood fit ===============================================
    auto mLogL = make_minuit(3,
      [&](const double* c) -> double {
        long double logl = 0.;
        const unsigned n = mc.size();
        #pragma omp parallel for reduction(+:logl)
        for (unsigned i=0; i<n; ++i) {
          logl += std::log( fit_bkg_sig(mc[i],c) );
        }
        return -2.*logl;
      }
    );
    mLogL.SetPrintLevel(-1);
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
    mLogL.Migrad();

    std::array<double,3> fit_c;
    for (const auto& p : mLogL.pars()) fit_c[p.i] = p.val;

    if (toy==0u) {
      to_root("first_toy",h);
      mkfcn("first_fit",
        [&](double* x, double* p){ return fit_bkg_sig(*x,p); }, fit_c);
    }

    const double rat = fit_c[0]/sig_frac;
    sig_yield(rat);
    yield_stats(rat);

    // Reset ========================================================
    for (auto& b : h) b = { };

  } // End TOY loop *************************************************

  { TH1* h = to_root("sig_yield",sig_yield);
    h->SetTitle(cat(std::setprecision(3),
      "fit/given signal fraction : "
      "#mu = ",yield_stats.mean(),", #sigma = ",yield_stats.stdev()
    ).c_str());
  }

  TNamed("args", cat(std::setw(2),json{
    {"seed", seed},
    {"bkg", in["bkg"]},
    {"sig", in["sig"]},
    {"yield", {
      {"ntoys", in["ntoys"]},
      {"mean", yield_stats.mean()},
      {"stdev", yield_stats.stdev()}
    }}
  })).Write();

  fout.Write();
}
