function out = verify_left_region()
% VERIFY_LEFT_REGION  Self-contained interval certificate for the LEFT region.
%
% Proves (rigorously, with outward-rounded interval arithmetic):
%   Certificate 2 :  Delta_0(u)    > 0  for u in (u0, uhat]
%   Certificate 3 :  Delta_tail(u) > 0  for u in [u0, uhat]
%
% These two, together with the right-region certificate (W_l > R_l on [uhat, u_crit],
% run from your existing certify_critical_ets_final.m) and Lemmas 1-3 of the paper,
% imply Split < Extend for all u in [u0, u_crit] and all iterations l >= 0.
%
% REQUIREMENTS:
%   GNU Octave + the 'interval' package.  (INTLAB/MATLAB would need `infsup` instead of
%   `infsupdec` and a few API changes; this file is written for Octave.)
%
% HOW TO RUN:
%   1) Save this file as  verify_left_region.m
%   2) In Octave:   pkg load interval
%   3)              out = verify_left_region();
%
% A quick smoke test of the kernels runs first; then the bisection certificate.
% Expected result: RESULT: PASS, with
%     min lower bound Delta_tail ~ 2.65e-3  (> 0)
%     min lower bound Delta_0    ~ 6.5e-3   (> 0)
%
% UNTESTED in a live Octave: authored offline; if anything errors, the message pinpoints
% the stage (smoke / seed / bisection) so it can be fixed quickly.

  if exist('infsupdec','file') == 0 && exist('infsupdec','class') == 0
    error('Octave interval package not loaded. Run: pkg load interval');
  end

  % ---------------- configuration ----------------
  u0   = 1.2385458589018437;   % certified LOWER bound of u0 (conservative; true u0 in [u0_lo, u0_hi])
  eps0 = 1.0e-4;               % Delta_0 endpoint gap handled by A0" (certified on [u0_lo, u0_hi+eps0])
  uhat = 2.0;
  hstar = 4.6;               % safe upper bound for hbar(u), u <= uhat
  MAX_U_CELLS = 20000;
  T_TRUNC = 30;
  PI = infsupdec('pi');
  phiCache = containers.Map('KeyType','char','ValueType','any');
  hbarCache = containers.Map('KeyType','char','ValueType','any');

  fprintf('============================================================\n');
  fprintf('Left-region interval certificate (Delta_0, Delta_tail on [u0,uhat])\n');
  fprintf('============================================================\n');

  % ---------------- smoke test of the kernels ----------------
  fprintf('[0/3] smoke test of the kernels ...\n');
  p1 = phi_point_ext(1.0, PI, phiCache, T_TRUNC);
  fprintf('  phi(1.0)      = [%.9g, %.9g]   (expect ~0.6498866)\n', inf(p1), sup(p1));
  p005 = phi_point_ext(0.05, PI, phiCache, T_TRUNC);
  fprintf('  phi(0.05)     = [%.9g, %.9g]   (expect ~0.9756005)  <-- u-branch\n', inf(p005), sup(p005));
  P1 = Phi_point(1.0, PI);
  fprintf('  Phi(1.0)      = [%.9g, %.9g]   (expect ~0.8413447)\n', inf(P1), sup(P1));
  if ~(inf(p1)<0.6499 && sup(p1)>0.6498)
    error('smoke: phi(1.0) enclosure does not overlap (0.6498, 0.6499)');
  end
  if ~(inf(p005)<0.9757 && sup(p005)>0.9755)
    error('smoke: phi(0.05) enclosure does not overlap (0.9755, 0.9757)');
  end
  if ~(inf(P1)<0.8414 && sup(P1)>0.8413)
    error('smoke: Phi(1.0) enclosure does not overlap (0.8413, 0.8414)');
  end
  fprintf('  smoke test PASS\n');

  % ---------------- seed table for phi^{-1} over [0.02, 35] ----------------
  fprintf('[1/3] building phi^{-1} seed table (600 pts over [0.02,35]) ...\n');
  [xTab, phiTab, pTab, xTabRev] = build_seed_table();
  fprintf('  seed table built (phi(0.02)=%.6g, phi(35)=%.3g)\n', phiTab(1), phiTab(end));

  % ---------------- bisect [u0, uhat] and certify ----------------
  fprintf('[2/3] bisecting [u0,uhat] and certifying ...\n');
  pending = [u0, uhat];
  nCells = 0;
  worstDT = +Inf;
  worstD0 = +Inf;

  while ~isempty(pending)
    ab = pending(end,:); pending(end,:) = [];
    a = ab(1); b = ab(2);

    % fixed-point enclosure hbar([a,b]) = [hbar(a), hbar(b)] (hbar increasing in u)
    ha = hbar_single(a, hstar, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab, hbarCache);
    hb = hbar_single(b, hstar, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab, hbarCache);
    hLo = inf(ha); hHi = sup(hb);

    % limiting states (Extend) over u in [a,b], hbar in [hLo,hHi]
    uI = infsupdec(a,b);
    hI = infsupdec(hLo,hHi);
    [gS,rS,vS,eS,bS,dS] = external_stars(uI, hI, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
    den = 1 - gS.^2;
    if inf(den) <= 0, error('1-g_*^2 <= 0 in a cell; abort.'); end
    aE = (dS + gS.*bS)./den;
    bE = (bS + gS.*dS)./den;
    A1inf = 2.*gS.*bE + dS;
    A2inf = 2.*gS.*aE + bS;

    % one-step values B1, C1 over u in [a,b]
    [B1, C1] = B1C1(uI, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);

    % Delta_tail (Certificate 3)
    DT = 2.*lnPhi_iv(B1, PI) + lnPhi_iv(C1, PI) ...
         - lnPhi_iv(A1inf, PI) - lnPhi_iv(A2inf, PI);

    % Delta_0 (Certificate 2), on the part of [a,b] with u >= u0+eps0
    D0lower = +Inf;
    if b > u0 + eps0
      uu0 = infsupdec(max(a, u0+eps0), b);
      [A10,A20,B0,C0] = ell0(uu0, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
      D0 = 2.*lnPhi_iv(B0, PI) + lnPhi_iv(C0, PI) ...
           - lnPhi_iv(A10, PI) - lnPhi_iv(A20, PI);
      D0lower = inf(D0);
    end

    d0Pass = (b <= u0 + eps0) || (D0lower > 0);
    if (inf(DT) > 0) && d0Pass
      % leaf passes: record its (positive) margins
      worstDT = min(worstDT, inf(DT));
      worstD0 = min(worstD0, D0lower);
    else
      mid = (a+b)/2;
      if mid <= a || mid >= b
        error('u-bisection reached floating-point resolution; NOT CERTIFIED.');
      end
      pending(end+1,:) = [a,mid]; %#ok<AGROW>
      pending(end+1,:) = [mid,b]; %#ok<AGROW>
    end
    nCells = nCells + 1;
    if mod(nCells, 20) == 0
      fprintf('  cells processed=%d, pending=%d, worstDT=%.3g, worstD0=%.3g\n', ...
              nCells, size(pending,1), worstDT, worstD0);
    end
    if nCells > MAX_U_CELLS
      error('Exceeded MAX_U_CELLS; NOT CERTIFIED.');
    end
  end

  % ---------------- report ----------------
  pass = (worstDT > 0) && (worstD0 > 0);
  fprintf('[3/3] RESULT\n');
  fprintf('  min lower bound Delta_tail = %.17g  (>0 required)\n', worstDT);
  fprintf('  min lower bound Delta_0    = %.17g  (>0 required)\n', worstD0);
  fprintf('  cells processed = %d\n', nCells);
  fprintf('  RESULT: %s\n', pass_text(pass));
  if ~pass
    error('Left-region certificate FAILED; NOT CERTIFIED.');
  end
  out = struct('Delta_tail_min_lower', worstDT, 'Delta_0_min_lower', worstD0, ...
               'cells', nCells, 'pass', pass);
  fprintf('Left region CERTIFIED: Delta_0>0 on (u0,uhat], Delta_tail>0 on [u0,uhat].\n');
end

% ======================================================================
%  phi kernels  (hybrid: u-representation for x<0.1, t-representation else)
% ======================================================================

function P = phi_point_ext(x, PI, phiCache, T_TRUNC)
  if x < 0.02 || x > 35
    error('phi_point_ext requires 0.02 <= x <= 35; got %.17g', x);
  end
  if x < 0.1
    P = phi_u(x, PI);
    return;
  end
  key = sprintf('%.17g', x);
  if isKey(phiCache, key)
    bb = phiCache(key); P = infsupdec(bb(1), bb(2)); return;
  end
  if x >= 2.5
    NN=2048; rr=2.0; xmin=2.5;
  elseif x >= 1.0
    NN=2048; rr=2.0; xmin=1.0;
  else
    NN=4096; rr=0.875; xmin=0.1;
  end
  TT = T_TRUNC;
  hh = (2*TT)/NN;
  t  = -TT + (0:NN)*hh;
  tI = infsupdec(t); XI = infsupdec(x); hI = infsupdec(hh);
  f = exp(-(tI.^2)./(4.*XI)).*sech(tI./2);
  S = (hI./3).*(f(1)+f(end)+infsupdec(4).*sum(f(2:2:NN))+infsupdec(2).*sum(f(3:2:NN-1)));
  rI = infsupdec(rr); xminI = infsupdec(xmin);
  M4I = infsupdec(24).*exp((rI.^2)./(4.*xminI)) ./ (cos(rI./2).*(rI.^4));
  EI  = (infsupdec(2*TT)./infsupdec(180)).*(hI.^4).*M4I;
  E = sup(EI);
  tailI = infsupdec(8).*exp(infsupdec(-TT/2));
  tailU = sup(tailI);
  J = S + infsupdec(-E,E) + infsupdec(0,tailU);
  prefactor = exp(-XI./4).*(4.*PI.*XI).^(-0.5);
  P = prefactor.*J;
  phiCache(key) = [inf(P), sup(P)];
end

function P = phi_u(x, PI)
  T = 6.0; N = 2000; r = 1.0;
  h  = (2*T)/N;
  u  = -T + (0:N)*h;
  uI = infsupdec(u); xI = infsupdec(x); hI = infsupdec(h); rI = infsupdec(r);
  f = exp(-(uI.^2)) .* sech(sqrt(xI).*uI);
  S = (hI./3).*(f(1)+f(end)+infsupdec(4).*sum(f(2:2:N))+infsupdec(2).*sum(f(3:2:N-1)));
  cosR = cos(sqrt(xI).*rI);
  if inf(cosR) <= 0, error('Cauchy radius invalid (cos(sqrt(x) r) <= 0).'); end
  M4 = infsupdec(24).*exp(rI.^2) ./ (cosR .* (rI.^4));
  E  = sup((infsupdec(2*T)./180).*(hI.^4).*M4);
  tailU = sup(exp(infsupdec(-T*T))./infsupdec(T));
  J = S + infsupdec(-E, E) + infsupdec(0, tailU);
  P = exp(-xI./4) ./ sqrt(PI) .* J;
end

function [xTab, phiTab, pTab, xTabRev] = build_seed_table()
  xTab = linspace(0.02, 35.0, 400);
  phiTab = zeros(size(xTab));
  for j = 1:numel(xTab)
    phiTab(j) = phi_fast(xTab(j));
  end
  if any(~isfinite(phiTab)) || any(diff(phiTab) >= 0)
    error('Seed table is not finite and strictly decreasing.');
  end
  pTab    = fliplr(phiTab);
  xTabRev = fliplr(xTab);
end

function y = phi_fast(x)
  if x < 0.1
    integrand = @(u) exp(-u.^2).*sech(sqrt(x)*u);
    I = integral(integrand, -10, 10, 'RelTol', 2e-11, 'AbsTol', 2e-14);
    y = exp(-x/4)/sqrt(pi)*I;
  else
    integrand = @(t) exp(-(t.^2)./(4*x)).*sech(t./2);
    I = integral(integrand, -Inf, Inf, 'RelTol', 2e-11, 'AbsTol', 2e-14);
    y = exp(-x/4)/sqrt(4*pi*x)*I;
  end
end

function P = phi_inverse_enclosure_ext(Ptarget, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab)
  pLo = max(0, inf(Ptarget));
  pHi = min(1, sup(Ptarget));
  if ~(isfinite(pLo) && isfinite(pHi) && pLo <= pHi)
    error('Invalid target probability interval for phi inverse.');
  end
  pMidRaw = pLo + (pHi-pLo)/2;
  pMid = min(max(pMidRaw, pTab(1)), pTab(end));
  seed = interp1(pTab, xTabRev, pMid, 'pchip');
  phiSeed = interp1(xTab, phiTab, seed, 'pchip');
  jj = find(xTab <= seed, 1, 'last');
  if isempty(jj), jj = 1; end
  if jj >= numel(xTab), jj = numel(xTab)-1; end
  slope = (phiTab(jj+1)-phiTab(jj))/(xTab(jj+1)-xTab(jj));
  if ~(isfinite(seed) && isfinite(phiSeed) && isfinite(slope) && slope < 0)
    error('Inverse seed/slope invalid.');
  end
  INV_PAD = 5.0e-9;
  xLo = seed + (pHi - phiSeed)/slope - INV_PAD;
  xHi = seed + (pLo - phiSeed)/slope + INV_PAD;
  xLo = min(35, max(0.02, xLo));
  xHi = min(35, max(0.02, xHi));
  if xLo >= xHi
    xLo = max(0.02, seed-INV_PAD); xHi = min(35, seed+INV_PAD);
  end
  if xLo >= xHi
    error('Initial inverse bracket collapsed at a domain endpoint.');
  end
  passLo = false; passHi = false; stepLo = INV_PAD; stepHi = INV_PAD;
  for it = 1:35
    if ~passLo
      Plo = phi_point_ext(xLo, PI, phiCache, T_TRUNC);
      passLo = inf(Plo) >= pHi;
      if ~passLo
        if xLo <= 0.02, error('Could not certify inverse lower bracket at x=0.02.'); end
        xLo = max(0.02, xLo - stepLo); stepLo = 2*stepLo;
      end
    end
    if ~passHi
      Phi = phi_point_ext(xHi, PI, phiCache, T_TRUNC);
      passHi = sup(Phi) <= pLo;
      if ~passHi
        if xHi >= 35, error('Could not certify inverse upper bracket at x=35.'); end
        xHi = min(35, xHi + stepHi); stepHi = 2*stepHi;
      end
    end
    if passLo && passHi
      P = infsupdec(xLo, xHi); return;
    end
  end
  error('Could not certify inverse bracket. p=[%.17g, %.17g], seed=%.17g', pLo, pHi, seed);
end

function P = phi_range_iv(X, PI, phiCache, T_TRUNC)
  xl = inf(X); xu = sup(X);
  if xl == xu
    P = phi_point_ext(xl, PI, phiCache, T_TRUNC);
  else
    Pl = phi_point_ext(xl, PI, phiCache, T_TRUNC);
    Pu = phi_point_ext(xu, PI, phiCache, T_TRUNC);
    P = infsupdec(inf(Pu), sup(Pl));   % phi decreasing
  end
end

% ======================================================================
%  normal CDF Phi and ln Phi (rigorous)
% ======================================================================

function P = Phi_point(x, PI)
  if x <= 0
    P = infsupdec('0.5'); return;
  end
  if x >= 6.0
    xI  = infsupdec(x);
    phiX = exp(-(xI.^2)./2) ./ sqrt(infsupdec(2).*PI);
    Qlo = phiX .* xI ./ (xI.^2 + 1);
    Qhi = phiX ./ xI;
    P = infsupdec(1 - sup(Qhi), 1 - inf(Qlo));
    return;
  end
  N = 2000;
  h = x / N;
  hI = infsupdec(h);
  t  = (0:N) * h;
  tI = infsupdec(t);
  g  = exp(-(tI.^2)./2) ./ sqrt(infsupdec(2).*PI);
  S  = (hI./3) .* (g(1) + g(end) + infsupdec(4).*sum(g(2:2:N)) + infsupdec(2).*sum(g(3:2:N-1)));
  xI  = infsupdec(x);
  M4  = max(infsupdec(6), abs(xI.^4 - 6.*xI.^2 + 3)) ./ sqrt(infsupdec(2).*PI);
  E   = sup((xI./180) .* (hI.^4) .* M4);
  P = infsupdec('0.5') + S + infsupdec(-E, E);
end

function L = lnPhi_iv(X, PI)
  xl = inf(X); xu = sup(X);
  Pl = Phi_point(xl, PI);
  Pu = Phi_point(xu, PI);
  P = infsupdec(inf(Pl), sup(Pu));   % Phi increasing
  L = log(P);                        % ln increasing
end

% ======================================================================
%  state recursion / certificate helpers
% ======================================================================

function hb = hbar_single(u, hstar, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab, hbarCache)
  key = sprintf('%.15g', u);
  if isKey(hbarCache, key)
    bb = hbarCache(key); hb = infsupdec(bb(1), bb(2)); return;
  end
  lo = u; hi = hstar;
  for it = 1:25
    mid = (lo+hi)/2;
    Fmid = infsupdec(mid) - lambda_point(mid, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
    if inf(Fmid) >= u
      hi = mid;
    elseif sup(Fmid) <= u
      lo = mid;
    else
      break;
    end
  end
  hb = infsupdec(lo, hi);
  hbarCache(key) = [inf(hb), sup(hb)];
end

function L = lambda_point(x, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab)
  ph = phi_point_ext(x, PI, phiCache, T_TRUNC);
  pV = 1 - (1-ph).^5;
  L = 2.*phi_inverse_enclosure_ext(pV, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
end

function [gS,rS,vS,eS,bS,dS] = external_stars(uI, hI, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab)
  phI = phi_range_iv(hI, PI, phiCache, T_TRUNC);
  gS = (1 - phI).^4;
  rS = sqrt(gS);
  vS = phi_inverse_enclosure_ext(1 - (1 - phI).^5, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
  phU = phi_range_iv(uI, PI, phiCache, T_TRUNC);
  eS = phi_inverse_enclosure_ext(1 - (1 - phU).*(1 - phI).^4, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
  bS = uI + vS;
  dS = bS + eS;
end

function [B1, C1] = B1C1(uI, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab)
  ph0 = phi_range_iv(uI, PI, phiCache, T_TRUNC);
  g0 = (1-ph0).^4; r0 = sqrt(g0);
  v0 = phi_inverse_enclosure_ext(1-(1-ph0).^5, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
  beta0 = uI + v0;
  h1 = uI + 2.*v0;                       % lambda(u) = 2*VN(u)
  ph1 = phi_range_iv(h1, PI, phiCache, T_TRUNC);
  g1 = (1-ph1).^4; r1 = sqrt(g1);
  v1 = phi_inverse_enclosure_ext(1-(1-ph1).^5, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
  beta1 = uI + v1;
  aS1 = r0.*uI + beta0;
  bS1 = g0.*beta0 + beta0;
  z1  = r0.*beta0 + uI;
  B1 = g1.*aS1 + r1.*z1 + beta1;
  C1 = 2.*r1.*bS1 + uI;
end

function [A1,A2,B,C] = ell0(uI, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab)
  ph0 = phi_range_iv(uI, PI, phiCache, T_TRUNC);
  g0 = (1-ph0).^4; r0 = sqrt(g0);
  v0 = phi_inverse_enclosure_ext(1-(1-ph0).^5, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
  e0 = v0;                               % Ex(u,0) = VN(u)
  beta0 = uI + v0; delta0 = beta0 + e0;
  A1 = 2.*g0.*beta0 + delta0;
  A2 = 2.*g0.*delta0 + beta0;
  B  = g0.*beta0 + r0.*uI + beta0;
  C  = 2.*r0.*beta0 + uI;
end

function s = pass_text(p)
  if p, s = 'PASS'; else, s = 'FAIL'; end
end
