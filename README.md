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
(non-rigorous) step in the proof.

Its two limitations should be stated explicitly:

1. **Grid, not interval.** The claim is checked at 49 discrete values of `u`. The variation of
   `(W-R)'_l` between grid points is not controlled, so this is not a statement about the whole
   interval `[2.0, u_crit]`.
2. **No accumulated-error bound.** The value `2^-52 ~ 2.2e-16` is the rounding scale of a *single*
   operation; it does not bound the accumulation of rounding through the recursion.

Consequently, the part of Lemma 4 concerning `l >= 1` (and its extension to `l > 60`) is
**supported numerically, not proved**. The extension to `l > 60` cannot be obtained from the
contraction argument of Lemma 3, because `h'_l` does not converge: since
`h'_{l+1} = 1 + lambda'(h_l) h'_l` and `lambda'(h_l) -> 1` as `u -> u_crit`, the derivative may grow
linearly instead of decaying.

The `l = 0` part of Lemma 4, certified by C2', **is** rigorous.

### Status of the theorem

| Part | Status |
|---|---|
| Left region `[u0, u_hat]` (A0", C0, C1) | rigorous |
| Right region, C2 at `u = 2.0` (prefix and quantitative tail) | rigorous |
| Right region, C2' (`l = 0`) | rigorous |
| Lemma 4, `l >= 1` (N1) | numerical evidence only |
| Certificate 1 (bisection invariant) | rigorous |

So the result should be described as a **mixed interval-certificate / floating-point verification**,
not as a fully rigorous proof for all `l >= 0`.

## A candidate route to a rigorous `l`-direction tail lemma

Numerical exploration suggests the following statement, which is *not* proved here but is a clean
target:

> **(Open)** For `u` in `[2.0, u_crit]`, the ratio `(W-R)'_l / h'_l` is nondecreasing in `l`.

If this holds, then `h'_l >= h'_0 = 1` gives

```
(W-R)'_l >= ((W-R)'_0 / h'_0) * h'_l >= (W-R)'_0 > 0,
```

which closes Lemma 4 for all `l`. Numerically the ratio rises monotonically from about `0.15` at
`l = 0` to an asymptotic value between `1` and `1.3` (`1.30` at `u = 2.5`, `1.17` at `u = 2.6`,
approaching `1` at `u = u_crit`).

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
| `CITATION.cff` | Citation metadata (repository-level only; author/paper fields are commented out) |

## Citation

If you use these certificates, please cite this repository (see `CITATION.cff`).

Note: `CITATION.cff` currently contains only repository-level metadata. The author names, the
repository URL and the reference to the accompanying paper are left commented out at the end of the
file, to be filled in once the paper is finalized. Linking this repository to Zenodo and minting a
DOI is recommended for long-term archiving.
