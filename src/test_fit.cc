#include <iostream>
#include <vector>
#include <cmath>

#include <boost/optional.hpp>

#include <nlohmann/json.hpp>
#include "json/boost_optional.hh"

#include <TMinuit.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>

#include "ivanp/time_seed.hh"
#include "ivanp/binner.hh"
#include "ivanp/string.hh"
#include "ivanp/timed_counter.hh"
#include "ivanp/root/minuit.hh"

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
  TNamed("seed",cat(seed).c_str()).Write();
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
  mkfcn("gen_bkg",[&](double* x, double*){ return gen_bkg(*x); });

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
  struct DSCB_t { // Double-sided Crystal Ball
    double mu, s, aL, nL, rL, eL, aH, nH, rH, eH, norm;
    DSCB_t(const json& sig, const std::array<double,2>& range)
    : mu(sig.at("muCB")), s(sig.at("sCB")),
      aL(sig.at("aLow" )), nL(sig.at("nLow" )), rL(nL/aL), eL(std::exp(-0.5*aL*aL)),
      aH(sig.at("aHigh")), nH(sig.at("nHigh")), rH(nH/aH), eH(std::exp(-0.5*aH*aH))
    {
      const double gL = (mu-range[0])/s + rL - aL;
      const double gH = (range[1]-mu)/s + rL - aL;
      norm = 1. / ( s * ( // σ from the Jacobian
          eL * (rL - std::pow(rL/gL,nL)*gL) / (nL-1)
        + eH * (rH - std::pow(rH/gH,nH)*gH) / (nH-1)
        + sqrt_pi_2 * (std::erf(aL*sqrt_1_2) + std::erf(aH*sqrt_1_2))
      ));
      TEST(norm)
      norm = 1./integrate(200,105,160,*this);
      TEST(norm)
    }
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

  // auto gaus = [](double x){
  //   return (1./std::sqrt(M_PI*2*4))*std::exp(-0.5*std::pow((x-125.)/2.,2));
  // };
  // TEST(integrate(100,105,160,gaus))

  TEST(sig_n)
  { std::uniform_real_distribution<double>
      dist_x(range[0],range[1]),
      dist_y(0,1);
    for (size_t i=bkg_n, n=bkg_n+sig_n; i<n; ) {
      const double x = dist_x(gen);
      if (dist_y(gen) < cb(x)) { h(mc[i] = x); ++i; }
    }
  }

  to_root("mc",h);

  const auto& axis = h.axis();
  const auto nbins = axis.nbins();
  vector<double> bin_centers(nbins);
  { double a = axis.edge(0), b;
    for (unsigned i=0; i<nbins; ) {
      ++i;
      b = axis.edge(i);
      bin_centers[i-1] = a + (b-a)/2;
      a = b;
    }
  }

  // Weighted Least Squares =========================================
  constexpr size_t nfs = 3;
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

  std::vector<double> us; // bin variances
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

  // Fit without signal =============================================
  auto fit_bkg = [&](double x, const double* c) {
    const double c0 = (6./330.)-(43725./330.)*c[0]-(5876750./330.)*c[1];
    return c0 + c[0]*x + c[1]*x*x;
  };
  { auto mLogL = make_minuit(2,
      [&](const double* c) -> double {
        long double logl = 0.;
        const unsigned n = mc.size();
        #pragma omp parallel for reduction(+:logl)
        for (unsigned i=0; i<n; ++i) {
          logl += std::log(fit_bkg(mc[i],c));
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
    std::array<double,2> fit_c;
    for (const auto& p : mLogL.pars()) {
      fit_c[p.i] = p.val;
      cout << p;
    }
    cout << endl;

    mkfcn("fit_bkg",
      [&](double* x, double* p){ return fit_bkg(*x,p); }, fit_c);
  }

  // Fit with signal ================================================
  auto fit_bkg_sig = [&](double x, const double* c) {
    return (1.-c[0]) * fit_bkg(x,c+1) + c[0] * cb(x) * cb.norm;
  };
  { auto mLogL = make_minuit(3,
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
    std::array<double,3> fit_c;
    for (const auto& p : mLogL.pars()) {
      fit_c[p.i] = p.val;
      cout << p;
    }
    cout << endl;
    TEST(sig_frac)
    cout << endl;

    mkfcn("fit_bkg_sig",
      [&](double* x, double* p){ return fit_bkg_sig(*x,p); }, fit_c);

    // Likelihood scan ----------------------------------------------
    boost::optional<std::tuple<unsigned,double,double>> scan = fit["scan"];
    if (scan && get<0>(*scan)!=0) {
      cout << "Yield likelihood scan "
        << get<0>(*scan) <<": "<< get<1>(*scan) <<' '<< get<2>(*scan) << endl;
      TH1D* h_logl = new TH1D(
        "logL scan","Yield likelihood scan;Signal fraction;-2logL",
        get<0>(*scan), get<1>(*scan), get<2>(*scan));
      const double step = (get<2>(*scan)-get<1>(*scan))/get<0>(*scan);
      const double start = get<1>(*scan) + 0.5*step;
      double ignore;
      std::array<double,3> fit_c;
      for (timed_counter<int> i(get<0>(*scan)); !!i; ++i) {
        mLogL.SetPrintLevel(-1);
        mLogL.DefineParameter(
          0, "sig", step**i+start, // num, name, start
          0, 0, 0 // step, min, max
        );
        mLogL.FixParameter(0);
        mLogL.Migrad();
        for (unsigned i=0; i<fit_c.size(); ++i)
          mLogL.GetParameter(i,fit_c[i],ignore);
        h_logl->SetBinContent(i+1,mLogL.f(fit_c.data()));
      }
    }
    cout << endl;

    { hist h_fit(axis), h_diff(axis);
      auto& b_fit  = h_fit .bins();
      auto& b_diff = h_diff.bins();
      double a = axis.edge(0), b;
      for (unsigned i=0; i<nbins; ) {
        ++i;
        b = axis.edge(i);
        b_diff[i-1] = h[{i-1}] - (
          b_fit[i-1] = tot_n*integrate(10,a,b,
            [&,c=fit_c.data()](double x){ return fit_bkg_sig(x,c); })
        );
        a = b;
      }
      to_root("binned_fit",h_fit);
      to_root("binned_diff",h_diff);
    }

  }

  // Fit with GP ====================================================
  if (fit["gp"].get<bool>()) {
    vector<double> gp_ys(nbins);
    const auto& bins = h.bins();
    // size_t call = 0;
    auto mLogL = make_minuit(5,
      [&](const double* c) -> double {
        long double logl = 0.;
        const unsigned n = mc.size();
        #pragma omp parallel for reduction(+:logl)
        for (unsigned i=0; i<n; ++i) {
          logl += std::log( fit_bkg_sig(mc[i],c) );
        }

        double a = axis.edge(0), b;
        for (unsigned i=0; i<nbins; ) {
          ++i;
          b = axis.edge(i);
          gp_ys[i-1] = bins[i-1] - tot_n*integrate(10,a,b,
            [&](double x){ return fit_bkg_sig(x,c); });
          a = b;
        }

        // TEST((++call))
        return -2.*(logl - gp_logml_opt(bin_centers, gp_ys, us,
          [](auto a, auto b, double s, double l){
            return s * std::exp(-0.5*sq((a-b)/l));
          }, c[3], c[4]));
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
    mLogL.DefineParameter(
      3, "gp_s", 1, // num, name, start
      0.1, 0.01, 10000 // step, min, max
    );
    mLogL.DefineParameter(
      4, "gp_l", 1, // num, name, start
      0.1, 0.1, 100 // step, min, max
    );
    // mLogL.FixParameter(0);
    mLogL.Migrad();

    cout << "\nLikelihood fit (signal + background):\n";
    std::array<double,3> fit_c;
    for (const auto& p : mLogL.pars()) {
      if (p.i<fit_c.size()) fit_c[p.i] = p.val;
      cout << p;
    }
    cout << endl;
    TEST(sig_frac)
    cout << endl;

    mkfcn("fit_bkg_sig_gp",
      [&](double* x, double* p){ return fit_bkg_sig(*x,p); }, fit_c);

    // const auto gp = GP(xs,ys,us,
    //   generator(0,nt,[a=xs.front(),s=(xs.back()-xs.front())/(nt-1)](auto i){
    //     return a + s*i;
    //   }),
    //   [](auto a, auto b){
    //     return std::exp((-0.5/sq(2))*sq(a-b));
    //   }
    // );
  }

  fout.Write();
}
