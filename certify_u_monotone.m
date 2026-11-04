function out = certify_u_monotone()
% CERTIFY_U_MONOTONE  Lemma 4 (l = 0 part): W_0(u) - R_0(u) is strictly increasing in u.
%
% Certifies, by interval arithmetic, that the l = 0 u-derivative
%     (W-R)'_0 = W_0' - R_0'  >  0   for every u in [uL, u_crit],
% where W_0 = 2 B_0 - A1_0 - A2_0,  R_0 = exp(-T_0/2)/(C_0 Phi(C_0)), T_0 = C_0^2 - B_0^2,
% and the l = 0 states and their u-derivatives are explicit functions of u (no recursion).
%
% Combined with Lemma 6 ((W-R)'_l is nondecreasing in l), this gives Lemma 4 for all l:
% the minimum of (W-R)' over l is at l = 0, so certifying l = 0 suffices.
%
% Requires GNU Octave + 'interval' package. Run:
%   pkg load interval
%   out = certify_u_monotone();
%
% UNTESTED in a live Octave: authored offline; review before use.

  if exist('infsupdec','file') == 0 && exist('infsupdec','class') == 0
    error('Octave interval package not loaded. Run: pkg load interval');
  end

  uL = 2.0;
  uR = 2.6309118378181706;   % certified upper bound of u_crit
  MAX_U_CELLS = 40000;
  T_TRUNC = 30;
  T_DER = 48;
  PI = infsupdec('pi');
  phiCache = containers.Map('KeyType','char','ValueType','any');

  fprintf('============================================================\n');
  fprintf('Lemma 4 (l=0): (W-R)''_0 > 0 on [%.10g, %.10g]\n', uL, uR);
  fprintf('============================================================\n');

  [xTab, phiTab, pTab, xTabRev] = build_seed_table();
  fprintf('seed table built\n');

  pending = [uL, uR];
  nCells = 0;
  worstGap = +Inf;

  while ~isempty(pending)
    ab = pending(end,:); pending(end,:) = [];
    a = ab(1); b = ab(2);
    uI = infsupdec(a,b);

    % ---- l = 0 states and their u-derivatives (explicit) ----
    phU  = phi_range_iv(uI, PI, phiCache, T_TRUNC);          % phi(u)
    ppU  = phi_prime_iv_range(uI, PI, phiCache, T_DER);      % phi'(u)
    q    = 1 - phU;
    g0   = q.^4;
    dg0  = infsupdec(-4).*q.^3 .* ppU;                       % g0' = -4(1-phi)^3 phi'
    r0   = q.^2;
    dr0  = infsupdec(-2).*q .* ppU;                          % r0' = -2(1-phi) phi'
    v0   = phi_inverse_enclosure_ext(1 - q.^5, PI, phiCache, T_TRUNC, xTab, pTab, xTabRev, phiTab);
    ppv  = phi_prime_iv_range(v0, PI, phiCache, T_DER);      % phi'(v0)
    dv0  = infsupdec(5).*q.^4 .* ppU ./ ppv;                 % v0' = 5(1-phi)^4 phi'/phi'(v0)
    beta0  = uI + v0;
    dbeta0 = infsupdec(1) + dv0;
    delta0  = beta0 + v0;                                    % e0 = v0, so delta0 = u + 2 v0
    ddelta0 = infsupdec(1) + 2.*dv0;

    A1_0 = 2.*g0.*beta0 + delta0;
    A2_0 = 2.*g0.*delta0 + beta0;
    B_0  = (1+g0).*beta0 + r0.*uI;
    C_0  = 2.*r0.*beta0 + uI;

    dA1_0 = 2.*dg0.*beta0 + 2.*g0.*dbeta0 + ddelta0;
    dA2_0 = 2.*dg0.*delta0 + 2.*g0.*ddelta0 + dbeta0;
    dB_0  = (1+g0).*dbeta0 + dg0.*beta0 + dr0.*uI + r0;
    dC_0  = 2.*dr0.*beta0 + 2.*r0.*dbeta0 + infsupdec(1);

    dW0 = 2.*dB_0 - dA1_0 - dA2_0;
    T0  = C_0.^2 - B_0.^2;
    R0  = exp(-T0./2) ./ (C_0 .* Phi_iv(C_0, PI));
    hC  = Phi_dens_iv(C_0, PI) ./ Phi_iv(C_0, PI);           % phi(C)/Phi(C)
    dR0 = R0 .* ( B_0.*dB_0 - dC_0.*( C_0 + hC + 1./C_0 ) );
    dWR0 = dW0 - dR0;

    gap = inf(dWR0);
    if gap > 0
      worstGap = min(worstGap, gap);
    else
      mid = (a+b)/2;
      if mid <= a || mid >= b, error('u-bisection reached resolution; NOT CERTIFIED.'); end
      pending(end+1,:) = [a,mid]; %#ok<AGROW>
      pending(end+1,:) = [mid,b]; %#ok<AGROW>
    end
    nCells = nCells + 1;
    if mod(nCells,20)==0
      fprintf('  cells=%d pending=%d worstGap=%.3g\n', nCells, size(pending,1), worstGap);
    end
    if nCells > MAX_U_CELLS, error('Exceeded MAX_U_CELLS'); end
  end

  pass = (worstGap > 0);
  fprintf('RESULT\n');
  fprintf('  min over cells of lower((W-R)''_0) = %.17g  (>0 required)\n', worstGap);
  fprintf('  cells = %d\n', nCells);
  fprintf('  RESULT: %s\n', pass_text(pass));
  if ~pass, error('Lemma 4 (l=0) certificate FAILED'); end
  out = struct('uL', uL, 'uR', uR, 'worstGap', worstGap, 'cells', nCells, 'pass', pass);
  fprintf('CERTIFIED: (W-R)''_0 > 0 for u in [%.10g, %.10g].\n', uL, uR);
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

% ==================== normal CDF ====================
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
  xl = inf(X); xu = sup(X);
  Pl = exp(-(infsupdec(xu).^2)./2) ./ sqrt(infsupdec(2).*PI);
  Pu = exp(-(infsupdec(xl).^2)./2) ./ sqrt(infsupdec(2).*PI);
  P = infsupdec(inf(Pl), sup(Pu));
end

function s = pass_text(p)
  if p, s = 'PASS'; else, s = 'FAIL'; end
end
