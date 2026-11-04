function out = verify_right_region()
% VERIFY_RIGHT_REGION  Interval certificate for the RIGHT region (single point uhat).
%
% Proves (outward-rounded interval arithmetic): at u = uhat = 2.0,
%   W_l > R_l   for every iteration l >= 0,
% where  W = 2B - A1 - A2,  T = C^2 - B^2,  R = exp(-T/2)/(C*Phi(C)).
%
% By Lemma 1 (comparison), W_l > R_l implies Delta_l > 0, i.e. Split < Extend at uhat.
% Combined with the u-monotonicity of W_l - R_l (Lemma 4), this covers the whole right
% region [uhat, u_crit].
%
% Requires GNU Octave + 'interval' package. Run:
%   pkg load interval
%   out = verify_right_region();
%
% Self-contained: includes the phi kernels (hybrid u/t representation) and the normal-CDF
% kernel. The finite prefix is l = 0..M-1 (orbit converges fast at u=2.0); the tail
% l >= M uses the monotonicity of the states.
%
% UNTESTED in a live Octave: authored offline; review before use.

  if exist('infsupdec','file') == 0 && exist('infsupdec','class') == 0
    error('Octave interval package not loaded. Run: pkg load interval');
  end

  uhat  = 2.0;         % left end of the right region (hardest point: W-R is smallest here)
  M     = 500;         % finite prefix l = 0..M-1  (orbit converges fast at u=2.0)
  hstar = 4.6;         % safe orbit upper bound (>= hbar(u) for all u <= u_crit)
  T_TRUNC = 30;
  PI = infsupdec('pi');
  phiCache = containers.Map('KeyType','char','ValueType','any');

  fprintf('============================================================\n');
  fprintf('Right-region interval certificate: W_l > R_l at u = %.10g\n', uhat);
  fprintf('============================================================\n');

  % ---- smoke test (same kernels as the left region) ----
  fprintf('[0/4] smoke test of the kernels ...\n');
  p1 = phi_point_ext(1.0, PI, phiCache, T_TRUNC);
  P1 = Phi_point(1.0, PI);
  fprintf('  phi(1.0) = [%.9g, %.9g]  (expect ~0.6499)\n', inf(p1), sup(p1));
  fprintf('  Phi(1.0) = [%.9g, %.9g]  (expect ~0.8413)\n', inf(P1), sup(P1));
  if ~(inf(p1)<0.6499 && sup(p1)>0.6498), error('smoke phi(1.0) failed'); end
  if ~(inf(P1)<0.8414 && sup(P1)>0.8413), error('smoke Phi(1.0) failed'); end
  fprintf('  smoke test PASS\n');

  % ---- seed table ----
  fprintf('[1/4] building phi^{-1} seed table ...\n');
  [xTab, phiTab, pTab, xTabRev] = build_seed_table();
  fprintf('  done\n');

  % ---- orbit propagation (single point u = uhat) ----
  fprintf('[2/4] propagating orbit and states, l = 0..%d ...\n', M);
  uI = infsupdec(uhat);
  N = M+1;
  h = infsupdec(zeros(N,1));  g = h; r = h; v = h; e = h; beta = h; delta = h;
  h(1) = uI;
  for k = 1:N
    phI = phi_range_iv(h(k), PI, phiCache, T_TRUNC);
    g(k) = (1-phI).^4;
    r(k) = sqrt(g(k));
    v(k) = phi_inverse_enclosure_ext(1-(1-phI).^5, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
    phU = phi_range_iv(uI, PI, phiCache, T_TRUNC);
    e(k) = phi_inverse_enclosure_ext(1-(1-phU).*(1-phI).^4, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
    beta(k)  = uI + v(k);
    delta(k) = beta(k) + e(k);
    if k < N
      hNext = uI + 2.*v(k);               % h_{k+1} = u + lambda(h_k) = u + 2*VN(h_k)
      hLo = max(inf(h(k)), inf(hNext));
      hHi = min(sup(hNext), hstar);
      if hLo > hHi, error('orbit enclosure inconsistent at ell=%d', k); end
      h(k+1) = infsupdec(hLo, hHi);
    end
  end

  % ---- internal states and W, R ----
  aE = infsupdec(zeros(N,1)); bE = aE; aS = aE; bS = aE; z = aE;
  aE(1) = delta(1); bE(1) = beta(1);
  aS(1) = beta(1); bS(1) = beta(1); z(1) = uI;
  for k = 1:M
    aE(k+1) = g(k).*bE(k) + delta(k);
    bE(k+1) = g(k).*aE(k) + beta(k);
    aS(k+1) = r(k).*z(k) + beta(k);
    bS(k+1) = g(k).*aS(k) + beta(k);
    z(k+1)  = r(k).*bS(k) + uI;
  end
  A1 = 2.*g.*bE + delta;
  A2 = 2.*g.*aE + beta;
  B  = g.*aS + r.*z + beta;
  C  = 2.*r.*bS + uI;
  T  = C.^2 - B.^2;
  R  = infsupdec(zeros(N,1));
  for k = 1:N
    R(k) = exp(-T(k)./2) ./ (C(k) .* Phi_iv(C(k), PI));   % exact R
  end

  % ---- finite prefix: W_l > R_l for l = 0..M-1 ----
  fprintf('[3/4] checking finite prefix W_l > R_l ...\n');
  worstWR = +Inf;  idxW = 0;
  for k = 1:M
    Wk = 2.*B(k) - A1(k) - A2(k);
    margin = inf(Wk) - sup(R(k));
    if margin < worstWR
      worstWR = margin; idxW = k-1;
    end
  end
  prefixPass = (worstWR > 0);
  fprintf('  min over l=0..%d of lower(W_l) - upper(R_l) = %.17g  (>0 required)\n', M-1, worstWR);
  fprintf('  attained at l = %d\n', idxW);
  fprintf('  finite-prefix W>R: %s\n', pass_text(prefixPass));

  % ---- tail: l >= M ----
  fprintf('[4/4] checking tail l >= %d ...\n', M);
  hbar = h(N);                       % orbit limit ~ hbar(uhat) = 2.336
  phI = phi_range_iv(hbar, PI, phiCache, T_TRUNC);
  gS = (1-phI).^4; rS = sqrt(gS);
  vS = phi_inverse_enclosure_ext(1-(1-phI).^5, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
  phU = phi_range_iv(uI, PI, phiCache, T_TRUNC);
  eS = phi_inverse_enclosure_ext(1-(1-phU).*(1-phI).^4, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
  bSstar = uI + vS; dSstar = bSstar + eS;
  den = 1 - gS.^2;
  if inf(den) <= 0, error('1-g_*^2 <= 0; tail NOT CERTIFIED.'); end
  aEstar = (dSstar + gS.*bSstar)./den;
  bEstar = (bSstar + gS.*dSstar)./den;
  A1star = 2.*gS.*bEstar + dSstar;
  A2star = 2.*gS.*aEstar + bSstar;
  BM = B(N); CM = C(N);
  tailW = 2.*BM - A1star - A2star;       % W_l >= 2 B_M - A1_* - A2_*  for l >= M
  if inf(CM) <= 0, error('C_M <= 0; tail NOT CERTIFIED.'); end
  tailR = 2 ./ inf(CM);                  % R_l < 2/C_l <= 2/C_M
  tailWLower = inf(tailW);
  tailPass = tailWLower > tailR;
  fprintf('  tailW = 2*B_M - A1_* - A2_*  lower = %.17g\n', tailWLower);
  fprintf('  tailR = 2/C_M                upper = %.17g\n', tailR);
  fprintf('  tail W>R: %s\n', pass_text(tailPass));

  pass = prefixPass && tailPass;
  fprintf('\n============================================================\n');
  fprintf('RIGHT-REGION CERTIFICATE (u = %.10g)\n', uhat);
  fprintf('  finite prefix W_l > R_l, l=0..%d : %s\n', M-1, pass_text(prefixPass));
  fprintf('  tail W_l > R_l, l>=%d           : %s\n', M, pass_text(tailPass));
  fprintf('  RESULT: %s\n', pass_text(pass));
  fprintf('============================================================\n');
  if ~pass, error('Right-region certificate FAILED; NOT CERTIFIED.'); end
  out = struct('uhat', uhat, 'M', M, 'prefixMargin', worstWR, 'prefixIndex', idxW, ...
               'tailWLower', tailWLower, 'tailRUpper', tailR, 'pass', pass);
  fprintf('CERTIFIED: W_l > R_l at u = %.10g for every l >= 0.\n', uhat);
end

% ======================================================================
%  phi kernels (hybrid u/t representation; domain [0.02, 35])
% ======================================================================
function P = phi_point_ext(x, PI, phiCache, T_TRUNC)
  if x < 0.02 || x > 35, error('phi_point_ext: 0.02 <= x <= 35 required, got %.17g', x); end
  if x < 0.1, P = phi_u(x, PI); return; end
  key = sprintf('%.17g', x);
  if isKey(phiCache, key), bb = phiCache(key); P = infsupdec(bb(1), bb(2)); return; end
  if x >= 2.5, NN=2048; rr=2.0; xmin=2.5;
  elseif x >= 1.0, NN=2048; rr=2.0; xmin=1.0;
  else, NN=4096; rr=0.875; xmin=0.1; end
  TT = T_TRUNC; hh = (2*TT)/NN; t = -TT + (0:NN)*hh;
  tI = infsupdec(t); XI = infsupdec(x); hI = infsupdec(hh);
  f = exp(-(tI.^2)./(4.*XI)).*sech(tI./2);
  S = (hI./3).*(f(1)+f(end)+infsupdec(4).*sum(f(2:2:NN))+infsupdec(2).*sum(f(3:2:NN-1)));
  rI = infsupdec(rr); xminI = infsupdec(xmin);
  M4I = infsupdec(24).*exp((rI.^2)./(4.*xminI)) ./ (cos(rI./2).*(rI.^4));
  EI = (infsupdec(2*TT)./infsupdec(180)).*(hI.^4).*M4I; E = sup(EI);
  tailU = sup(infsupdec(8).*exp(infsupdec(-TT/2)));
  J = S + infsupdec(-E,E) + infsupdec(0,tailU);
  P = exp(-XI./4).*(4.*PI.*XI).^(-0.5) .* J;
  phiCache(key) = [inf(P), sup(P)];
end

function P = phi_u(x, PI)
  T = 6.0; N = 2000; r = 1.0;
  h = (2*T)/N; u = -T + (0:N)*h;
  uI = infsupdec(u); xI = infsupdec(x); hI = infsupdec(h); rI = infsupdec(r);
  f = exp(-(uI.^2)) .* sech(sqrt(xI).*uI);
  S = (hI./3).*(f(1)+f(end)+infsupdec(4).*sum(f(2:2:N))+infsupdec(2).*sum(f(3:2:N-1)));
  cosR = cos(sqrt(xI).*rI);
  if inf(cosR) <= 0, error('Cauchy radius invalid'); end
  M4 = infsupdec(24).*exp(rI.^2) ./ (cosR .* (rI.^4));
  E = sup((infsupdec(2*T)./180).*(hI.^4).*M4);
  tailU = sup(exp(infsupdec(-T*T))./infsupdec(T));
  J = S + infsupdec(-E,E) + infsupdec(0,tailU);
  P = exp(-xI./4) ./ sqrt(PI) .* J;
end

function [xTab, phiTab, pTab, xTabRev] = build_seed_table()
  xTab = linspace(0.02, 35.0, 400);
  phiTab = zeros(size(xTab));
  for j = 1:numel(xTab), phiTab(j) = phi_fast(xTab(j)); end
  if any(~isfinite(phiTab)) || any(diff(phiTab) >= 0), error('Seed table invalid'); end
  pTab = fliplr(phiTab); xTabRev = fliplr(xTab);
end

function y = phi_fast(x)
  if x < 0.1
    I = integral(@(u) exp(-u.^2).*sech(sqrt(x)*u), -10, 10, 'RelTol', 2e-11, 'AbsTol', 2e-14);
    y = exp(-x/4)/sqrt(pi)*I;
  else
    I = integral(@(t) exp(-(t.^2)./(4*x)).*sech(t./2), -Inf, Inf, 'RelTol', 2e-11, 'AbsTol', 2e-14);
    y = exp(-x/4)/sqrt(4*pi*x)*I;
  end
end

function P = phi_inverse_enclosure_ext(Ptarget, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab)
  pLo = max(0, inf(Ptarget)); pHi = min(1, sup(Ptarget));
  if ~(isfinite(pLo) && isfinite(pHi) && pLo <= pHi), error('Invalid p interval'); end
  pMid = min(max(pLo + (pHi-pLo)/2, pTab(1)), pTab(end));
  seed = interp1(pTab, xTabRev, pMid, 'pchip');
  phiSeed = interp1(xTab, phiTab, seed, 'pchip');
  jj = find(xTab <= seed, 1, 'last');
  if isempty(jj), jj = 1; end
  if jj >= numel(xTab), jj = numel(xTab)-1; end
  slope = (phiTab(jj+1)-phiTab(jj))/(xTab(jj+1)-xTab(jj));
  if ~(isfinite(seed) && isfinite(phiSeed) && isfinite(slope) && slope < 0), error('seed/slope invalid'); end
  INV_PAD = 5.0e-9;
  xLo = seed + (pHi - phiSeed)/slope - INV_PAD;
  xHi = seed + (pLo - phiSeed)/slope + INV_PAD;
  xLo = min(35, max(0.02, xLo)); xHi = min(35, max(0.02, xHi));
  if xLo >= xHi, xLo = max(0.02, seed-INV_PAD); xHi = min(35, seed+INV_PAD); end
  if xLo >= xHi, error('bracket collapsed'); end
  passLo = false; passHi = false; stepLo = INV_PAD; stepHi = INV_PAD;
  for it = 1:35
    if ~passLo
      Plo = phi_point_ext(xLo, PI, phiCache, T_TRUNC);
      passLo = inf(Plo) >= pHi;
      if ~passLo
        if xLo <= 0.02, error('inverse lower bracket fail at x=0.02'); end
        xLo = max(0.02, xLo - stepLo); stepLo = 2*stepLo;
      end
    end
    if ~passHi
      Phi = phi_point_ext(xHi, PI, phiCache, T_TRUNC);
      passHi = sup(Phi) <= pLo;
      if ~passHi
        if xHi >= 35, error('inverse upper bracket fail at x=35'); end
        xHi = min(35, xHi + stepHi); stepHi = 2*stepHi;
      end
    end
    if passLo && passHi, P = infsupdec(xLo, xHi); return; end
  end
  error('inverse bracket not certified, p=[%.17g, %.17g]', pLo, pHi);
end

function P = phi_range_iv(X, PI, phiCache, T_TRUNC)
  xl = inf(X); xu = sup(X);
  if xl == xu, P = phi_point_ext(xl, PI, phiCache, T_TRUNC);
  else
    Pl = phi_point_ext(xl, PI, phiCache, T_TRUNC);
    Pu = phi_point_ext(xu, PI, phiCache, T_TRUNC);
    P = infsupdec(inf(Pu), sup(Pl));   % phi decreasing
  end
end

% ======================================================================
%  normal CDF Phi (rigorous)
% ======================================================================
function P = Phi_point(x, PI)
  if x <= 0, P = infsupdec('0.5'); return; end
  if x >= 6.0
    xI = infsupdec(x);
    phiX = exp(-(xI.^2)./2) ./ sqrt(infsupdec(2).*PI);
    Qlo = phiX .* xI ./ (xI.^2 + 1);
    Qhi = phiX ./ xI;
    P = infsupdec(1 - sup(Qhi), 1 - inf(Qlo)); return;
  end
  N = 2000; h = x/N; hI = infsupdec(h); t = (0:N)*h; tI = infsupdec(t);
  g = exp(-(tI.^2)./2) ./ sqrt(infsupdec(2).*PI);
  S = (hI./3).*(g(1)+g(end)+infsupdec(4).*sum(g(2:2:N))+infsupdec(2).*sum(g(3:2:N-1)));
  xI = infsupdec(x);
  M4 = max(infsupdec(6), abs(xI.^4 - 6.*xI.^2 + 3)) ./ sqrt(infsupdec(2).*PI);
  E = sup((xI./180).*(hI.^4).*M4);
  P = infsupdec('0.5') + S + infsupdec(-E, E);
end

function P = Phi_iv(X, PI)
  xl = inf(X); xu = sup(X);
  Pl = Phi_point(xl, PI); Pu = Phi_point(xu, PI);
  P = infsupdec(inf(Pl), sup(Pu));   % Phi increasing
end

function s = pass_text(p)
  if p, s = 'PASS'; else, s = 'FAIL'; end
end
