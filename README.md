# Certificates for Theorem 1

Computer-assisted verification scripts for the proof that, for the (3,6)-regular LDPC 4-cycle
elementary trapping set (ETS), the **Split** repair scheme yields a strictly lower local error rate
than the **Extend** scheme on the SNR interval `[u0, u_crit]`, for every iteration `l >= 0`.

All scripts are written for **GNU Octave** and use **outward-rounded interval arithmetic** (except
N1, which is a floating-point verification; see below).

## Requirements

- GNU Octave (any recent version)
- The `interval` package for outward-rounded interval arithmetic

One-time installation of the package:

```octave
pkg install -forge interval
```

## How to run

Put all `.m` files in a single directory. Each script is **self-contained**: it carries its own
copies of the `phi`, `phi'`, `phi^{-1}` and `Phi` interval kernels, so no external helper files are
needed.

### Run everything at once

```octave
cd('/path/to/certificates')
pkg load interval
clear functions
run_all
```

This runs all seven certificates in order and prints a summary table with each certificate's
minimum margin and PASS/FAIL result.

### Run individual certificates

```octave
cd('/path/to/certificates')
pkg load interval
clear functions

out1 = certify_A0prime;          % A0"
out2 = verify_left_region;       % C0 and C1 (one script, two margins)
out3 = verify_right_region;      % C2
out4 = certify_u_monotone;       % C2'
out5 = certify_lambda_monotone;  % Certificate 1
out6 = verify_l_monotone;        % N1 (floating-point verification)
```

Every script prints `RESULT: PASS` together with its minimum margin, and aborts with a non-zero
status if the corresponding certificate fails.

## Certificates

| ID | File | Call | Statement | Min. margin | Type |
|---|---|---|---|---|---|
| A0" | `certify_A0prime.m` | `out = certify_A0prime;` | `Delta_0'(u) > 0` on `[u0_lo, u0_hi + 1e-4]` | 0.159190 | interval |
| C0 | `verify_left_region.m` | `out = verify_left_region;` | `Delta_0(u) > 0` on `[u0_lo + 1e-4, 2.0]` | 6.27e-7 | interval |
| C1 | `verify_left_region.m` | (same call) | `Delta_tail(u) > 0` on `[u0_lo, 2.0]` | 9.80e-4 | interval |
| C2 | `verify_right_region.m` | `out = verify_right_region;` | `W_l(u) > R_l(u)` at `u = 2.0`, every `l >= 0` | 0.602231 | interval |
| C2' | `certify_u_monotone.m` | `out = certify_u_monotone;` | `(W-R)'_0(u) > 0` on `[2.0, u_crit]` | 0.022896 | interval |
| Cert. 1 | `certify_lambda_monotone.m` | `out = certify_lambda_monotone;` | `lambda'(h) < 1` on `[u0_lo, 2.4]` | 0.275111 | interval |
| N1 | `verify_l_monotone.m` | `out = verify_l_monotone;` | `(W-R)'_l >= (W-R)'_0` for `l <= 60`, 49-point grid | 0.718357 | float |

## Interval certificates vs. floating-point verification

All certificates except **N1** are **interval certificates**. The algorithm bisects the target
interval, computes an outward-rounded interval enclosure of the target function on each leaf, and
requires the infimum of every enclosure to be positive. Because the enclosures contain the true
values and the leaves cover the whole domain, the conclusion is rigorous and independent of
floating-point rounding.

**N1** is a **floating-point verification** on a finite grid: it checks `(W-R)'_l(u) >= (W-R)'_0(u)`
for `l = 0, ..., 60` and 49 values of `u` in `[2.0, u_crit]`. It is the only numerical
(non-rigorous) step in the proof. Its rounding error is about `2^-52 ~ 2.2e-16`, far below the
observed minimum gap `0.718357`.

## Thresholds and enclosures

- The crossing point `u0` (at `l = 0`) is enclosed by
  `[u0_lo, u0_hi] = [1.2385458589018437, 1.2385858589018435]`.
- The critical point `u_crit` (where `lambda'(h*) = 1`) is certified by the upper bound
  `2.6309118378181706`.
- Certificates A0", C0, C1 and Cert. 1 use the **whole enclosure** `[u0_lo, u0_hi]`, so their
  conclusions hold for the **true** `u0` (not just its midpoint).

## Reproducibility

All scripts are self-contained: each has its own copies of the `phi`, `phi'`, `phi^{-1}` and `Phi`
interval kernels. Only the `interval` package is required. Because interval arithmetic is
outward-rounded, the minimum margins should reproduce up to the last few digits; small differences
may arise from different Octave or `interval` package versions.

If you change the certified enclosure of `u0` or the certified upper bound of `u_crit`, then A0",
C0, C1 and Cert. 1 (which use the `u0` enclosure) and C2, C2' (which use `u_crit`) must be re-run.

## Files

| File | Purpose |
|---|---|
| `README.md` | This document |
| `run_all.m` | Runs all seven certificates and prints a summary table |
| `certify_A0prime.m` | Certificate A0" |
| `verify_left_region.m` | Certificates C0 and C1 |
| `verify_right_region.m` | Certificate C2 |
| `certify_u_monotone.m` | Certificate C2' |
| `certify_lambda_monotone.m` | Certificate 1 |
| `verify_l_monotone.m` | Numerical verification N1 |
| `CITATION.cff` | Citation metadata (fill in author and repository fields before release) |

## Citation

If you use these certificates, please cite the accompanying paper (see `CITATION.cff`).
Note: `CITATION.cff` contains `TODO` placeholders for the author names and the repository URL;
replace them before publishing. Linking this repository to Zenodo and minting a DOI is recommended
for long-term archiving.
