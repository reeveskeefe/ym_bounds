from ym_bounds.contraction import contraction_sequence, check_summability_and_product
from ym_bounds.su3 import dim_su3, casimir_su3, tail_bound_arxiv_style
from ym_bounds.tube_cost import spectral_gap_from_tube_cost
from ym_bounds.string_tension import sigma_phys_from_lattice
from ym_bounds.clustering import exp_clustering_bound

def test_su3_basic():
    assert dim_su3(1,0) == 3  # fundamental
    assert round(casimir_su3(1,0), 6) == round(4/3, 6)

def test_contraction_and_product():
    etas = contraction_sequence(eta0=0.05, A=3.0, steps=10)
    s, prod, s_ok, p_ok = check_summability_and_product(etas, C=0.2)
    assert s_ok and p_ok
    assert prod > 0.0

def test_tail_bound():
    partial, tail = tail_bound_arxiv_style(beta=6.0, p0=5, q0=0)
    assert partial >= 0 and tail >= 0

def test_gap_and_sigma():
    lam, m0 = spectral_gap_from_tube_cost(0.4)
    assert 0 < lam < 1 and m0 == 0.4
    assert sigma_phys_from_lattice(0.045, 0.08) > 0

def test_clustering():
    assert exp_clustering_bound(1.0, 0.3, 1.0) < 1.0