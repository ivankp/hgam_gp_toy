#include <iostream>
#include <vector>
#include <cstdlib>

#include "linalg.hh"

using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
  std::vector<double> c(argc-3);
  for (int i=3; i<argc; ++i)
    c[i-3] = atof(argv[i]);

  cout << "old coeffs:\n";
  for (double c : c) cout << c << '\n';
  cout.flush();

  linalg::change_poly_coords(c.data(),c.size(),atof(argv[1]),atof(argv[2]));

  cout << "\nnew coeffs:\n";
  for (double c : c) cout << c << '\n';
  cout.flush();
}
