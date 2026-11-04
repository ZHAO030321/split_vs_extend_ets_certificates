function out = certify_lambda_monotone()
% CERTIFY_LAMBDA_MONOTONE  Certificate 1: lambda'(h) < 1 on [u0, 2.4].
%
% This certifies that F(h) = h - lambda(h) is strictly increasing on [u0, 2.4], which is
% all that the fixed-point bisection needs: for u in [u0, 2.0] the orbit limit satisfies
% hbar(u) <= hbar(2.0) = 2.336, so the bisection invariant "hi >= hbar" only requires
% F increasing on [u0, 2.336] (lambda'(2.4) = 0.433, margin 0.567 below 1).  No behaviour
% near the critical point h_* is needed.
%
% Method: bisection in h. For each cell [a,b] in [u0, 2.4]:
%   lambda'(h) = 10*(1-phi(h))^4 * phi'(h) / phi'(lambda(h)/2),
% evaluated with interval arithmetic, and checked to be < 1.
%
% phi' is computed by the t-representation for x >= 0.1 and by the u-representation
% (well-conditioned) for x in [0.02, 0.1]; the latter is needed because
% lambda(u0)/2 = 0.0226 < 0.1.
%
% Requires GNU Octave + 'interval' package. Run:
%   pkg load interval
%   out = certify_lambda_monotone();
%
% UNTESTED in a live Octave: authored offline; review before use.

  if exist('infsupdec','file') == 0 && exist('infsupdec','class') == 0
    error('Octave interval package not loaded. Run: pkg load interval');
  end

  u0    = 1.2385458589018437;   % certified LOWER bound of u0 (conservative; true u0 in [u0_lo, u0_hi])
  hstar = 2.4;                  % safe upper bound on hbar(2.0)=2.336; the bisection only needs F increasing on [u0, 2.336]
  epsEnd = 0.0;                 % no critical-point gap needed
  MAX_H_CELLS = 40000;
  T_TRUNC = 30;
  T_DER = 48;
  PI = infsupdec('pi');
  phiCache = containers.Map('KeyType','char','ValueType','any');

  fprintf('============================================================\n');
  fprintf('Certificate 1: lambda''(h) < 1 on [%.10g, %.10g]\n', u0, hstar-epsEnd);
  fprintf('============================================================\n');

  % ---- seed table ----
  fprintf('[1/3] building phi^{-1} seed table ...\n');
  [xTab, phiTab, pTab, xTabRev] = build_seed_table();
  fprintf('  done\n');

  % ---- bisect [u0, hstar-epsEnd], certify lambda' < 1 ----
  fprintf('[2/3] bisecting and certifying lambda'' < 1 ...\n');
  pending = [u0, hstar-epsEnd];
  nCells = 0;
  worstGap = +Inf;   % min over leaves of (1 - sup(lambda'))  (worst margin below 1)

  while ~isempty(pending)
    ab = pending(end,:); pending(end,:) = [];
    a = ab(1); b = ab(2);
    hI = infsupdec(a,b);

    % phi(h) over [a,b]
    ph = phi_range_iv(hI, PI, phiCache, T_TRUNC);
    % phi'(h) at the endpoints (phi' increasing; hull of endpoint enclosures)
    Pp_a = phi_prime_iv(a, PI, phiCache, T_DER);
    Pp_b = phi_prime_iv(b, PI, phiCache, T_DER);
    Pp = infsupdec(min(inf(Pp_a),inf(Pp_b)), max(sup(Pp_a),sup(Pp_b)));
    % lambda(h)/2 at the endpoints (lambda increasing)
    la = lambda_point(a, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
    lb = lambda_point(b, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
    argLo = inf(la)/2; argHi = sup(lb)/2;
    % phi'(lambda(h)/2) over [argLo, argHi]
    Pd_a = phi_prime_iv(argLo, PI, phiCache, T_DER);
    Pd_b = phi_prime_iv(argHi, PI, phiCache, T_DER);
    Pd = infsupdec(min(inf(Pd_a),inf(Pd_b)), max(sup(Pd_a),sup(Pd_b)));

    Lp = infsupdec(10).*(1-ph).^4 .* Pp ./ Pd;

    if sup(Lp) < 1
      worstGap = min(worstGap, 1 - sup(Lp));
    else
      mid = (a+b)/2;
      if mid <= a || mid >= b
        error('h-bisection reached resolution; NOT CERTIFIED.');
      end
      pending(end+1,:) = [a,mid]; %#ok<AGROW>
      pending(end+1,:) = [mid,b]; %#ok<AGROW>
    end
    nCells = nCells + 1;
    if mod(nCells,50)==0
      fprintf('  cells=%d pending=%d worstGap=%.3g\n', nCells, size(pending,1), worstGap);
    end
    if nCells > MAX_H_CELLS, error('Exceeded MAX_H_CELLS'); end
  end

  pass = (worstGap > 0);
  fprintf('[3/3] RESULT\n');
  fprintf('  min over leaves of (1 - sup(lambda'')) = %.17g  (>0 required)\n', worstGap);
  fprintf('  cells processed = %d\n', nCells);
  fprintf('  RESULT: %s\n', pass_text(pass));
  if ~pass, error('Certificate 1 FAILED'); end
  out = struct('u0', u0, 'hstar', hstar, 'epsEnd', epsEnd, 'worstGap', worstGap, ...
               'cells', nCells, 'pass', pass);
  fprintf('CERTIFIED: lambda''(h) < 1 on [%.10g, %.10g].\n', u0, hstar-epsEnd);
end

% ==================== phi kernels ====================
function P = phi_point_ext(x, PI, phiCache, T_TRUNC)
  if x < 0.02 || x > 35, error('phi_point_ext: 0.02<=x<=35 required, got %.17g', x); end
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
    P = infsupdec(inf(Pu), sup(Pl));
  end
end

function L = lambda_point(x, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab)
  ph = phi_point_ext(x, PI, phiCache, T_TRUNC);
  pV = 1 - (1-ph).^5;
  L = 2.*phi_inverse_enclosure_ext(pV, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
end

% ==================== phi' kernels ====================
function Pp = phi_prime_iv(x, PI, phiCache, T_DER)
% Enclosure of phi'(x) at a point x; dispatch t-representation (x>=0.1) / u-rep (x<0.1).
  if x < 0.1
    Pp = phi_prime_u(x, PI);
  else
    Pp = phi_prime_t(x, PI, phiCache, T_DER);
  end
end

function Pp = phi_prime_u(x, PI)
% phi'(x) via the u-representation, for x in [0.02, 0.1]:
%   phi(x) = exp(-x/4)/sqrt(pi) * J(x),  J(x) = int e^{-u^2} sech(sqrt(x) u) du,
%   phi'(x) = exp(-x/4)/sqrt(pi) * [ (-1/4) J(x) - (1/(2 sqrt(x))) K(x) ],
%   K(x) = int u e^{-u^2} sech(sqrt(x) u) tanh(sqrt(x) u) du.
  T = 6.0; N = 2000; r = 1.0;
  h = (2*T)/N; u = -T + (0:N)*h;
  uI = infsupdec(u); xI = infsupdec(x); hI = infsupdec(h); rI = infsupdec(r);
  % J
  fJ = exp(-(uI.^2)) .* sech(sqrt(xI).*uI);
  SJ = (hI./3).*(fJ(1)+fJ(end)+infsupdec(4).*sum(fJ(2:2:N))+infsupdec(2).*sum(fJ(3:2:N-1)));
  cosR = cos(sqrt(xI).*rI);
  if inf(cosR) <= 0, error('Cauchy radius invalid (J)'); end
  M4J = infsupdec(24).*exp(rI.^2) ./ (cosR .* (rI.^4));
  EJ = sup((infsupdec(2*T)./180).*(hI.^4).*M4J);
  tailJ = sup(exp(infsupdec(-T*T))./infsupdec(T));
  J = SJ + infsupdec(-EJ, EJ) + infsupdec(0, tailJ);
  % K (|tanh| <= 1 on the disk since sqrt(x) r <= 0.32 < pi/4)
  fK = uI .* exp(-(uI.^2)) .* sech(sqrt(xI).*uI) .* tanh(sqrt(xI).*uI);
  SK = (hI./3).*(fK(1)+fK(end)+infsupdec(4).*sum(fK(2:2:N))+infsupdec(2).*sum(fK(3:2:N-1)));
  M4K = infsupdec(24).*infsupdec(T+r).*exp(rI.^2) ./ (cosR .* (rI.^4));
  EK = sup((infsupdec(2*T)./180).*(hI.^4).*M4K);
  tailK = sup(exp(infsupdec(-T*T)));
  K = SK + infsupdec(-EK, EK) + infsupdec(-tailK, tailK);
  % phi'(x)
  Pp = exp(-xI./4) ./ sqrt(PI) .* ( (-infsupdec(1)./4).*J - (infsupdec(1)./(2.*sqrt(xI))).*K );
end

function Pp = phi_prime_t(x, PI, phiCache, T_DER)
% phi'(x) via the t-representation, for x >= 0.1 (adapted from the u_crit code).
  key = sprintf('%.17g_d', x);
  if isKey(phiCache, key), bb = phiCache(key); Pp = infsupdec(bb(1), bb(2)); return; end
  if x >= 2.5, NN=2048; rr=2.0; xmin=2.5;
  elseif x >= 1.0, NN=2048; rr=2.0; xmin=1.0;
  else, NN=4096; rr=0.875; xmin=0.1; end
  TT = T_DER; hh = (2*TT)/NN; t = -TT + (0:NN)*hh;
  tI = infsupdec(t); XI = infsupdec(x); hI = infsupdec(hh);
  f0 = exp(-(tI.^2)./(4.*XI)).*sech(tI./2);
  f1 = (tI.^2)./(4.*XI.^2).*f0;
  S0 = (hI./3).*(f0(1)+f0(end)+infsupdec(4).*sum(f0(2:2:NN))+infsupdec(2).*sum(f0(3:2:NN-1)));
  S1 = (hI./3).*(f1(1)+f1(end)+infsupdec(4).*sum(f1(2:2:NN))+infsupdec(2).*sum(f1(3:2:NN-1)));
  rI = infsupdec(rr); xminI = infsupdec(xmin);
  cosR = cos(rI./2);
  if inf(cosR) <= 0, error('Cauchy radius invalid (t-derivative)'); end
  M40 = infsupdec(24).*exp((rI.^2)./(4.*xminI)) ./ (cosR.*(rI.^4));
  Q1 = (exp(infsupdec(-1))./xminI + (rI.^2)./(4.*xminI.^2)) .* exp((rI.^2)./(4.*xminI));
  M41 = infsupdec(24).*Q1 ./ (cosR.*(rI.^4));
  EI0 = (infsupdec(2*TT)./infsupdec(180)).*(hI.^4).*M40;
  EI1 = (infsupdec(2*TT)./infsupdec(180)).*(hI.^4).*M41;
  E0 = sup(EI0); E1 = sup(EI1);
  tail0U = sup(infsupdec(8).*exp(infsupdec(-TT/2)));
  tail1I = (infsupdec(2*TT^2+8*TT+16)./(xminI.^2)) .* exp(infsupdec(-TT/2));
  tail1U = sup(tail1I);
  I0 = S0 + infsupdec(-E0,E0) + infsupdec(0,tail0U);
  I1 = S1 + infsupdec(-E1,E1) + infsupdec(0,tail1U);
  prefactor = exp(-XI./4).*(4.*PI.*XI).^(-0.5);
  Pp = prefactor .* (I1 - (infsupdec(1)./4 + infsupdec(1)./(2.*XI)) .* I0);
  phiCache(key) = [inf(Pp), sup(Pp)];
end

function s = pass_text(p)
  if p, s = 'PASS'; else, s = 'FAIL'; end
end
