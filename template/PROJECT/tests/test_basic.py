from __future__ import annotations

import scikit_build_example as m

def test_add():
    assert m.add(1, 2) == 3


def test_sub():
    assert m.subtract(1, 2) == -1

m.sinPattern(1, 0, 3, 1000)