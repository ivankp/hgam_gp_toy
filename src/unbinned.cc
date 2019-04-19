#include <iostream>
#include <vector>
#include <cmath>

#include <nlohmann/json.hpp>

#include <TFile.h>
#include <TH1.h>

#include "ivanp/time_seed.hh"
#include "ivanp/string.hh"
#include "ivanp/timed_counter.hh"
#include "ivanp/binner.hh"

#include "gp.hh"
#include "generator.hh"
#include "gsl_multimin.hh"
#include "json/struct.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::get;
using std::string;
using std::vector;
using namespace ivanp;
using linalg::sq;
using nlohmann::json;

JSON_SERIALIZER(gsl_multimin_opts,(verbose)(tolerance)(max_iter))

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

int main(int argc, char* argv[]) {
  json in;
  std::cin >> in;

  TFile fout(in["output"].get<string>().c_str(),"recreate");

  long seed;
  try { seed = in["seed"]; }
  catch (...) { seed = time_seed(); }
  TEST(seed)
  TNamed("seed",cat(seed).c_str()).Write();
  std::mt19937 gen(seed);

  const std::array<double,2> range = in["range"];

  vector<double> mc(in["mc"]["nevents"].get<size_t>());
  hist h({in["mc"]["nbins"],range[0],range[1]});

  auto gen_bkg = [ // polynomial or exp(poly)
    c = in["mc"]["poly"].get<vector<double>>(),
    e = in["mc"]["exp"].get<bool>()
  ](double x){
    double xn = 1, p = c[0];
    for (unsigned i=1, n=c.size(); i<n; ++i)
      p += c[i] * (xn *= x);
    return e ? std::exp(p) : p;
  };
  // mkfcn("gen_bkg",[&](double* x, double*){ return gen_bkg(*x); });

  cout << "Generating " << mc.size() << " events" << endl;
  { std::uniform_real_distribution<double>
      dist_x(range[0],range[1]),
      dist_y(0,gen_bkg(range[0]));
    for (size_t i=0, n=mc.size(); i<n; ) {
      const double x = dist_x(gen);
      if (dist_y(gen) < gen_bkg(x)) { h(mc[i] = x); ++i; }
    }
  }

  to_root("mc",h);

  cout << "GP optimization" << endl;
  const auto hs = gsl_multimin<2>({{
      { 1., 0.1 },
      { 1., 0.1 }
    }}, [&](double s, double l){
      return gp_logml_opt(mc,
        generator(0,mc.size(),[](auto){ return 1; }), // ys
        generator(0,mc.size(),[](auto){ return 1; }), // us
        [](auto a, auto b, double s, double l){ // kernel
          return s * std::exp(-0.5*sq((a-b)/l));
        }, s, l);
    }, in["gp"]["opt"]
  );

  const unsigned nt = in["gp"]["ntest"];
  const auto gp = GP(mc,
    generator(0,mc.size(),[](auto){ return 1; }), // ys
    generator(0,mc.size(),[](auto){ return 1; }), // us
    generator(0,nt,[a=range[0],s=(range[1]-range[0])/(nt-1)](auto i){ // ts
      return a + s*i;
    }),
    [s=hs[0],l=hs[1]](auto a, auto b){ // kernel
      return s * std::exp(-0.5*sq((a-b)/l));
    }
  );

  fout.Write();
}
