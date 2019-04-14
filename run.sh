#!/bin/bash

./bin/toy_test << CARD
{
  "seed": null,
  "range": [ 105, 160 ],
  "bkg": {
    "exp": true,
    "poly": [ 1.394e+01, -7.083e-02, 1.349e-04 ],
    "_poly": [ 301204, -3706.58, 11.6843 ],
    "n": 67484
  },
  "sig": {
    "muCB"  : 125.14,
    "sCB"   : 1.82,
    "aLow"  : 1.55,
    "nLow"  : 9.93,
    "aHigh" : 1.82,
    "nHigh" : 8.42,
    "n": 286
  },
  "fit": {
    "verbose": 1,
    "exp": false,
    "nbins": 55,
    "scan": [ 0, -0.01, 0.02 ],
    "gp": true
  }
}
CARD
