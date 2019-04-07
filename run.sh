#!/bin/bash

./bin/toy_test << CARD
{
  "seed": null,
  "bkg": {
    "exp": true,
    "poly": [ 1.394e+01, -7.083e-02, 1.349e-04 ],
    "_poly": [ 301204, -3706.58, 11.6843 ],
    "min": 105, "max": 160, "n": 67484
  },
  "sig": {
    "min": 121, "max": 129, "n": 286
  },
  "fit": {
    "verbose": 1,
    "exp": false,
    "nbins": 55
  }
}
CARD
