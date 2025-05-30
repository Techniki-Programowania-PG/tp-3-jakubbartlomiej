from __future__ import annotations

import scikit_build_example as m

def test_sin():
    m.sinPattern(0.16, 0, 6.28, 1000)

def test_cos():
    m.cosPattern(0.16, -3.1415, 3.1415, 1000)

def test_square():
    m.squarePattern(1, -3.1415, 3.1415, 1000)

def test_sawtooth():
    m.sawtoothPattern(1, -3.1415, 3.1415, 1000)