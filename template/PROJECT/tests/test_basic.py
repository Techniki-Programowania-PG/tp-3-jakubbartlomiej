from __future__ import annotations

import scikit_build_example as m

def test_sin():
    m.showPattern(m.sinPattern(0.16, 0, 6.28, 1000))

def test_cos():
    m.showPattern(m.cosPattern(0.16, -3.1415, 3.1415, 1000))

def test_square():
    m.showPattern(m.squarePattern(1, -3.1415, 3.1415, 1000))

def test_sawtooth():
    m.showPattern(m.sawtoothPattern(1, -3.1415, 3.1415, 1000))

def test_dft():
    m.showDFTPattern(m.discreteFourierTransform(m.sinPattern(10, 0, 3.14, 1000)))

def test_inverse_dft():
    m.inverseDiscreteFourierTransform(m.discreteFourierTransform(m.sinPattern(10, 0, 3.14, 1000)))