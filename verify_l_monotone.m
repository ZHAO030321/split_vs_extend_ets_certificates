function out = verify_l_monotone()
% VERIFY_L_MONOTONE  Numerical verification N1 (Octave, double precision).
%
% Verifies:  (W-R)'_l(u) >= (W-R)'_0(u)  for l = 0..60 and u on a fine grid in [2.0, u_crit],
% hence (with C2' giving (W-R)'_0 > 0) the monotonicity in l of Lemma 4.
%
% This is a FLOATING-POINT verification, NOT an interval certificate: it checks the claim on a
% finite grid and at finite l.  It is the only numerical step in the proof.
%
% Run:
%   out = verify_l_monotone();
%
% Expected: OVERALL PASS, min over (u,l) of (W-R)'_l - (W-R)'_0 ~ 0.718 (> 0).

  uR = 2.6309118378181706;   % certified upper bound of u_crit
  M  = 60;

  [phi, phiinv, dphi, lam] = build_kernels();

  % ---- smoke test ----
  fprintf('smoke: phi(1.0) = %.10f  (expect ~0.6498866)\n', phi(1.0));
  fprintf('smoke: Phi(1.0) = %.10f  (expect ~0.8413447)\n', Phi(1.0));

  % ---- per-u table ----
  us = [2.0, 2.3, 2.5, 2.6, 2.63];
  fprintf('\nper-u check (l = 0..%d):\n', M);
  fprintf('  %-7s %10s %10s %10s %10s   %-9s %s\n', 'u','(W-R)''_0','(W-R)''_10','(W-R)''_30','(W-R)''_60','min incr','monotone');
  allok = true; gmin = +Inf;
  for u = us
    s = wrp(u, M, phi, phiinv, dphi, lam);
    inc = min(diff(s));
    ok = (inc >= -1e-9) && (min(s) > 0);
    allok = allok && ok;
    gmin = min(gmin, min(s));
    fprintf('  %-7.3f %10.4f %10.4f %10.4f %10.4f   %+.2e  %s\n', ...
            u, s(1), s(11), s(31), s(61), inc, yn(ok));
  end
  fprintf('  global min (W-R)'' = %.4f (should be at l = 0)\n', gmin);

  % ---- fine grid: the actual claim ----
  us2 = linspace(2.0, uR, 49);
  worst = +Inf; wu = 0; wk = 0; fails = 0;
  for u = us2
    s = wrp(u, M, phi, phiinv, dphi, lam);
    for k = 2:(M+1)
      d = s(k) - s(1);
      if d < worst, worst = d; wu = u; wk = k-1; end
      if d < -1e-9, fails = fails + 1; end
    end
  end
  fprintf('\nfine grid: %d u values in [2.0, u_crit], l = 0..%d\n', numel(us2), M);
  fprintf('  min over (u,l) of (W-R)''_l - (W-R)''_0 = %.6g  at u=%.6f, l=%d\n', worst, wu, wk);
  fprintf('  #violations = %d\n', fails);
  ok2 = (fails == 0);
  fprintf('  RESULT: %s\n', pf(ok2));
  pass = allok && ok2;
  fprintf('\nOVERALL: %s\n', pf(pass));
  out = struct('min_diff', worst, 'fails', fails, 'cells', numel(us2), 'pass', pass);
end

% ==================== kernels (float) ====================
function [phi, phiinv, dphi, lam] = build_kernels()
  T  = 40.0; NT = 8001; tt = linspace(-T, T, NT);
  XLO = 0.02; XHI = 20.0; NX = 4001;
  xg  = linspace(XLO, XHI, NX);
  phv = zeros(1, NX);
  sech2 = 1 ./ cosh(tt/2);
  for ii = 1:NX
    x = xg(ii);
    integ = exp(-(tt.^2)/(4*x)) .* sech2;
    Ii = trapz(tt, integ);
    phv(ii) = exp(-x/4) / sqrt(4*pi*x) * Ii;
  end
  if any(~isfinite(phv)) || any(diff(phv) >= 0)
    error('phi table invalid (must be finite and strictly decreasing).');
  end
  phi    = @(x) interp1(xg, phv, x, 'pchip');
  phiinv = @(p) interp1(fliplr(phv), fliplr(xg), p, 'pchip');
  dx     = xg(2) - xg(1);
  dphi   = @(x) (phi(x+dx) - phi(x-dx)) / (2*dx);
  lam    = @(x) 2 * phiinv(1 - (1 - phi(x)).^5);
end

function y = Phi(z)
  y = 0.5 * erfc(-z / sqrt(2));
end
function y = ndens(z)
  y = exp(-(z.^2)/2) / sqrt(2*pi);
end

% ==================== derivative recursion ====================
function s = wrp(u, M, phi, phiinv, dphi, lam)
  N  = M + 1;
  h = zeros(1,N); dh = zeros(1,N);
  g = zeros(1,N); r = zeros(1,N); v = zeros(1,N); e = zeros(1,N);
  bta = zeros(1,N); dlt = zeros(1,N);
  dg = zeros(1,N); dr = zeros(1,N); dv = zeros(1,N); de = zeros(1,N);
  dbta = zeros(1,N); ddlt = zeros(1,N);

  h(1) = u; dh(1) = 1;
  phu = phi(u); dpu = dphi(u);
  for k = 1:N
    ph = phi(h(k)); q = 1 - ph;
    g(k) = q^4; r(k) = q^2;
    v(k) = phiinv(1 - q^5);
    e(k) = phiinv(1 - (1 - phu) * q^4);
    bta(k) = u + v(k); dlt(k) = bta(k) + e(k);
    pp = dphi(h(k));
    dlam = 10 * q^4 * pp / dphi(lam(h(k))/2);
    dg(k) = -4 * q^3 * pp * dh(k);
    dr(k) = -2 * q * pp * dh(k);
    dv(k) = (dlam/2) * dh(k);
    de(k) = (dpu * q^4 + (1 - phu) * 4 * q^3 * pp * dh(k)) / dphi(e(k));
    dbta(k) = 1 + dv(k); ddlt(k) = dbta(k) + de(k);
    if k < N
      dh(k+1) = 1 + dlam * dh(k);
      h(k+1)  = u + 2 * v(k);
    end
  end

  aE = zeros(1,N); bE = zeros(1,N); aS = zeros(1,N); bS = zeros(1,N); zz = zeros(1,N);
  daE = zeros(1,N); dbE = zeros(1,N); daS = zeros(1,N); dbS = zeros(1,N); dz = zeros(1,N);
  aE(1) = dlt(1); bE(1) = bta(1); aS(1) = bta(1); bS(1) = bta(1); zz(1) = u;
  daE(1) = ddlt(1); dbE(1) = dbta(1); daS(1) = dbta(1); dbS(1) = dbta(1); dz(1) = 1;
  for k = 1:M
    aE(k+1) = g(k)*bE(k) + dlt(k);  bE(k+1) = g(k)*aE(k) + bta(k);
    aS(k+1) = r(k)*zz(k) + bta(k);  bS(k+1) = g(k)*aS(k) + bta(k);
    zz(k+1) = r(k)*bS(k) + u;
    daE(k+1) = dg(k)*bE(k) + g(k)*dbE(k) + ddlt(k);
    dbE(k+1) = dg(k)*aE(k) + g(k)*daE(k) + dbta(k);
    daS(k+1) = dr(k)*zz(k) + r(k)*dz(k) + dbta(k);
    dbS(k+1) = dg(k)*aS(k) + g(k)*daS(k) + dbta(k);
    dz(k+1)  = dr(k)*bS(k) + r(k)*dbS(k) + 1;
  end

  B = zeros(1,N); C = zeros(1,N); dB = zeros(1,N); dC = zeros(1,N);
  dA1 = zeros(1,N); dA2 = zeros(1,N);
  for k = 1:N
    B(k)  = g(k)*aS(k) + r(k)*zz(k) + bta(k);
    C(k)  = 2*r(k)*bS(k) + u;
    dB(k) = dg(k)*aS(k) + g(k)*daS(k) + dr(k)*zz(k) + r(k)*dz(k) + dbta(k);
    dC(k) = 2*dr(k)*bS(k) + 2*r(k)*dbS(k) + 1;
    dA1(k)= 2*dg(k)*bE(k) + 2*g(k)*dbE(k) + ddlt(k);
    dA2(k)= 2*dg(k)*aE(k) + 2*g(k)*daE(k) + dbta(k);
  end

  s = zeros(1,N);
  for k = 1:N
    Tk = C(k)^2 - B(k)^2;
    R  = exp(-Tk/2) / (C(k) * Phi(C(k)));
    hC = ndens(C(k)) / Phi(C(k));
    dR = R * (B(k)*dB(k) - dC(k)*(C(k) + hC + 1/C(k)));
    s(k) = 2*dB(k) - dA1(k) - dA2(k) - dR;   % (W-R)'_k
  end
end

function s = yn(ok)
  if ok, s = 'YES'; else, s = 'NO'; end
end

function s = pf(ok)
  if ok, s = 'PASS'; else, s = 'FAIL'; end
end
