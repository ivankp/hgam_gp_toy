#!/bin/bash

./bin/toy_test << CARD
{
  "seed": null,
  "bkg": {
    "exp": true,
    "poly": [ 1.394e+01, -7.083e-02, 1.349e-04 ],
    "_poly": [ 301204, -3706.58, 11.6843 ],
    "a": 105, "b": 160,
    "n": 1e3
  },
  "sig": {
    "a": 121, "b": 129,
    "n": 100
  },
  "fit": {
    "verbose": 1,
    "exp": false
  }
}
CARD
