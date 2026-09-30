from adt.core import (canon, marks, observed_prime_power_threshold, odd_part,
                      periodic_cycle_points, periodic_cycles,
                      predicted_prime_power_threshold, solvable_fingerprint, zset_mul)


def test_burnside_product():
    assert canon(zset_mul({2: 1}, {3: 1})) == ((6, 1),)
    assert canon(zset_mul({2: 1}, {4: 1})) == ((4, 2),)


def test_mobius_equivalent_marks():
    assert marks({1: 2, 3: 1}, 4) == (2, 2, 5, 2)


def test_collision_65_119_direct_and_formula():
    assert solvable_fingerprint((5, 13)) == solvable_fingerprint((7, 17))
    assert all(canon(periodic_cycles(65, c)) == canon(periodic_cycles(119, c)) for c in (0, -2))


def test_formula_matches_direct_for_small_primes():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        f0, fm2 = solvable_fingerprint((p,))
        assert f0 == canon(periodic_cycles(p, 0))
        assert fm2 == canon(periodic_cycles(p, -2))


def test_total_periodic_counts_recover_local_odd_parts():
    # N_0=1+m_- and N_-2=(m_-+m_+)/2, hence the two totals recover
    # m_-=odd(p-1) and m_+=odd(p+1), which proves prime separation.
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 97, 101):
        f0, fm2 = solvable_fingerprint((p,))
        n0 = sum(length * count for length, count in f0)
        nm2 = sum(length * count for length, count in fm2)
        assert n0 == 1 + odd_part(p - 1)
        assert nm2 == (odd_part(p - 1) + odd_part(p + 1)) // 2
        assert 2 * nm2 - n0 + 1 == odd_part(p + 1)


def test_multiplier_threshold_all_small_cases():
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        for c in range(-3, 4):
            assert predicted_prime_power_threshold(p, c) == observed_prime_power_threshold(p, c)


def test_multiplier_zero_one_and_other_cases():
    assert periodic_cycle_points(3, -1) == [(0, 2)]
    assert predicted_prime_power_threshold(3, -1) is None
    assert observed_prime_power_threshold(3, -1) is None
    assert periodic_cycle_points(3, 1) == [(2,)]
    assert predicted_prime_power_threshold(3, 1) == 1
    assert observed_prime_power_threshold(3, 1) == 1
    assert marks(periodic_cycles(3, 1), 1) == (1,)
    assert marks(periodic_cycles(9, 1), 1) == (0,)
    assert predicted_prime_power_threshold(5, -3) == 8
    assert observed_prime_power_threshold(5, -3) == 8


def test_all_critical_examples_are_equal_at_every_depth():
    for c in (1, 2):
        assert predicted_prime_power_threshold(41, c) is None
        assert canon(periodic_cycles(41, c)) == canon(periodic_cycles(41 * 41, c))
