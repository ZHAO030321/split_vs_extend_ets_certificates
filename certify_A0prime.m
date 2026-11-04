function out = certify_A0prime()
% CERTIFY_A0PRIME
% Interval certificate for (L1) of the LEFT region:
%
%     Delta_0'(u) > 0   for all u in [u0_lo, u0_hi + eps0],   eps0 = 1e-4,
%
% where [u0_lo, u0_hi] is the certified enclosure of the crossing point u0
% (u0_lo = 1.2385458589018437, u0_hi = 1.2385858589018435).  Using the whole
% enclosure (not its midpoint) makes the certificate valid for the TRUE u0.
%
% (L1) is needed because verify_left_region.m handles the endpoint gap
% [u0_lo, u0_lo+eps0] by *passing it without checking Delta_0*: it relies on
% A0".  Certifying Delta_0' > 0 on the WHOLE segment, together with
% Delta_0(u0) = 0, gives Delta_0 > 0 on (u0, u0_hi + eps0], which overlaps
% verify_left_region's Delta_0 certificate on [u0_lo + eps0, uhat].
% All computations use outward-rounded interval arithmetic (GNU Octave +
% 'interval').
%
% Delta_0'(u) = 2 h(B0) B0' + h(C0) C0' - h(A10) A10' - h(A20) A20',
%   h(x) = n(x)/Phi(x)  (n = standard normal density),
%   and B0,C0,A10,A20 are the ell=0 decision values with their u-derivatives.
%
% Run:
%   pkg load interval
%   out = certify_A0prime();
%
% Expected: RESULT: PASS, min lower bound ~ 0.1594 (> 0).
% UNTESTED in a live Octave: authored offline; review before use.

  if exist('infsupdec','file') == 0 && exist('infsupdec','class') == 0
    error('Octave interval package not loaded. Run: pkg load interval');
  end

  u0_lo = 1.2385458589018437;  % certified LOWER bound of u0
  u0_hi = 1.2385858589018435;  % certified UPPER bound of u0
  eps0  = 1.0e-4;              % endpoint gap width (same as verify_left_region.m)
  MAX_CELLS = 20000;
  T_TRUNC = 30;
  T_DER   = 48;
  PI = infsupdec('pi');
  phiCache = containers.Map('KeyType','char','ValueType','any');

  fprintf('============================================================\n');
  fprintf('A0" certificate:  Delta_0''(u) > 0 on [u0_lo, u0_hi + eps0]\n');
  fprintf('u0 in [%.16g, %.16g],  eps0 = %.3g\n', u0_lo, u0_hi, eps0);
  fprintf('============================================================\n');

  [xTab, phiTab, pTab, xTabRev] = build_seed_table();
  fprintf('seed table built\n');

  pending = [u0_lo, u0_hi + eps0];
  nCells = 0;
  worst  = +Inf;

  while ~isempty(pending)
    ab = pending(end,:); pending(end,:) = [];
    a = ab(1); b = ab(2);
    uI = infsupdec(a,b);

    Dp = D0prime_iv(uI, PI, phiCache, T_TRUNC, T_DER, xTab, pTab, xTabRev, phiTab);
    g  = inf(Dp);

    if g > 0
      worst = min(worst, g);
    else
      mid = (a+b)/2;
      if mid <= a || mid >= b
        error('bisection reached floating-point resolution; NOT CERTIFIED.');
      end
      pending(end+1,:) = [a,mid]; %#ok<AGROW>
      pending(end+1,:) = [mid,b]; %#ok<AGROW>
    end
    nCells = nCells + 1;
    if nCells > MAX_CELLS, error('Exceeded MAX_CELLS; NOT CERTIFIED.'); end
  end

  pass = (worst > 0);
  fprintf('RESULT\n');
  fprintf('  min lower bound of Delta_0'' = %.17g   (> 0 required)\n', worst);
  fprintf('  cells processed = %d\n', nCells);
  fprintf('  RESULT: %s\n', pass_text(pass));
  if ~pass, error('A0" certificate FAILED; NOT CERTIFIED.'); end
  fprintf('A0" CERTIFIED: Delta_0'' > 0 on [u0_lo, u0_hi + eps0].\n');
  out = struct('u0_lo', u0_lo, 'u0_hi', u0_hi, 'eps0', eps0, 'min_lower', worst, 'cells', nCells, 'pass', pass);
end

% ======================================================================
%  Delta_0'(u)  (interval enclosure)
% ======================================================================
function Dp = D0prime_iv(uI, PI, phiCache, T_TRUNC, T_DER, xTab, pTab, xTabRev, phiTab)
  % ell = 0: h0 = u, so v0 = phi^{-1}(1-(1-phi(u))^5) = lambda(u)/2
  phu  = phi_range_iv(uI, PI, phiCache, T_TRUNC);
  dphu = phi_prime_iv_range(uI, PI, phiCache, T_DER);
  v0   = phi_inverse_enclosure_ext(1 - (1-phu).^5, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
  dphv = phi_prime_iv_range(v0, PI, phiCache, T_DER);          % phi'(lambda(u)/2)
  dlam = infsupdec(10) .* (1-phu).^4 .* dphu ./ dphv;          % lambda'(u)
  v0p  = dlam ./ 2;                                            % v0' = lambda'(u)/2

  g0  = (1-phu).^4;               g0p = infsupdec(-4) .* (1-phu).^3 .* dphu;
  r0  = (1-phu).^2;               r0p = infsupdec(-2) .* (1-phu) .* dphu;

  b0  = uI + v0;                  b0p = infsupdec(1) + v0p;         % beta_0 = u + v0
  d0  = uI + 2.*v0;               d0p = infsupdec(1) + 2.*v0p;      % delta_0 = u + 2 v0

  B0  = b0.*(1+g0) + r0.*uI;
  B0p = b0p.*(1+g0) + b0.*g0p + r0p.*uI + r0;
  C0  = 2.*r0.*b0 + uI;
  C0p = 2.*r0p.*b0 + 2.*r0.*b0p + infsupdec(1);
  A10 = d0 + 2.*g0.*b0;
  A10p= d0p + 2.*g0p.*b0 + 2.*g0.*b0p;
  A20 = b0 + 2.*g0.*d0;
  A20p= b0p + 2.*g0p.*d0 + 2.*g0.*d0p;

  hB  = Phi_dens_iv(B0,  PI) ./ Phi_iv(B0,  PI);
  hC  = Phi_dens_iv(C0,  PI) ./ Phi_iv(C0,  PI);
  hA1 = Phi_dens_iv(A10, PI) ./ Phi_iv(A10, PI);
  hA2 = Phi_dens_iv(A20, PI) ./ Phi_iv(A20, PI);

  Dp  = 2.*hB.*B0p + hC.*C0p - hA1.*A10p - hA2.*A20p;
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
    P = infsupdec(inf(Pu), sup(Pl));   % phi decreasing
  end
end

% ==================== phi' kernels ====================
function Pp = phi_prime_iv(x, PI, phiCache, T_DER)
  if x < 0.1, Pp = phi_prime_u(x, PI);
  else, Pp = phi_prime_t(x, PI, phiCache, T_DER); end
end

function Pp = phi_prime_iv_range(X, PI, phiCache, T_DER)
  xl = inf(X); xu = sup(X);
  P1 = phi_prime_iv(xl, PI, phiCache, T_DER);
  P2 = phi_prime_iv(xu, PI, phiCache, T_DER);
  Pp = infsupdec(min(inf(P1),inf(P2)), max(sup(P1),sup(P2)));
end

function Pp = phi_prime_u(x, PI)
  T = 6.0; N = 2000; r = 1.0;
  h = (2*T)/N; u = -T + (0:N)*h;
  uI = infsupdec(u); xI = infsupdec(x); hI = infsupdec(h); rI = infsupdec(r);
  fJ = exp(-(uI.^2)) .* sech(sqrt(xI).*uI);
  SJ = (hI./3).*(fJ(1)+fJ(end)+infsupdec(4).*sum(fJ(2:2:N))+infsupdec(2).*sum(fJ(3:2:N-1)));
  cosR = cos(sqrt(xI).*rI);
  if inf(cosR) <= 0, error('Cauchy radius invalid (J)'); end
  M4J = infsupdec(24).*exp(rI.^2) ./ (cosR .* (rI.^4));
  EJ = sup((infsupdec(2*T)./180).*(hI.^4).*M4J);
  tailJ = sup(exp(infsupdec(-T*T))./infsupdec(T));
  J = SJ + infsupdec(-EJ, EJ) + infsupdec(0, tailJ);
  fK = uI .* exp(-(uI.^2)) .* sech(sqrt(xI).*uI) .* tanh(sqrt(xI).*uI);
  SK = (hI./3).*(fK(1)+fK(end)+infsupdec(4).*sum(fK(2:2:N))+infsupdec(2).*sum(fK(3:2:N-1)));
  M4K = infsupdec(24).*infsupdec(T+r).*exp(rI.^2) ./ (cosR .* (rI.^4));
  EK = sup((infsupdec(2*T)./180).*(hI.^4).*M4K);
  tailK = sup(exp(infsupdec(-T*T)));
  K = SK + infsupdec(-EK, EK) + infsupdec(-tailK, tailK);
  Pp = exp(-xI./4) ./ sqrt(PI) .* ( (-infsupdec(1)./4).*J - (infsupdec(1)./(2.*sqrt(xI))).*K );
end

function Pp = phi_prime_t(x, PI, phiCache, T_DER)
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

% ==================== Phi and normal density ====================
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
  P = infsupdec(inf(Pl), sup(Pu));
end

function P = Phi_dens_iv(X, PI)
  % standard normal density over an interval: general sign-safe enclosure
  xl = inf(X); xu = sup(X);
  pl = exp(-(infsupdec(xl).^2)./2) ./ sqrt(infsupdec(2).*PI);
  pu = exp(-(infsupdec(xu).^2)./2) ./ sqrt(infsupdec(2).*PI);
  P = infsupdec(min(inf(pl),inf(pu)), max(sup(pl),sup(pu)));
end

function s = pass_text(p)
  if p, s = 'PASS'; else, s = 'FAIL'; end
end
