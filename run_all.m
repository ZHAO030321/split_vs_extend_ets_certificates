function run_all()
% RUN_ALL  Run all seven certificates and print a summary table.
%
% Usage (from the directory containing the certificate scripts):
%   pkg load interval
%   clear functions
%   run_all
%
% A0", C0, C1, C2, C2' and Certificate 1 are rigorous interval certificates.
% N1 is a floating-point verification (the only non-rigorous step).

  here = fileparts(mfilename('fullpath'));
  if ~isempty(here), addpath(here); end

  if exist('infsupdec','file') == 0 && exist('infsupdec','class') == 0
    error('Run "pkg load interval" first.');
  end

  names = {'A0"', 'C0', 'C1', 'C2', 'C2''', 'Cert. 1', 'N1'};
  files = {'certify_A0prime.m', 'verify_left_region.m', 'verify_left_region.m', ...
           'verify_right_region.m', 'certify_u_monotone.m', ...
           'certify_lambda_monotone.m', 'verify_l_monotone.m'};
  types = {'interval','interval','interval','interval','interval','interval','float'};

  n = numel(names);
  pass = false(1, n);
  marg = nan(1, n);

  fprintf('============================================================\n');
  fprintf(' Running all seven certificates\n');
  fprintf('============================================================\n');

  % ---------- 1) A0" ----------
  hdr('A0"', 1);
  try
    r = certify_A0prime;
    pass(1) = logical(r.pass); marg(1) = r.min_lower;
  catch err
    fprintf('\n  ERROR: %s\n', err.message);
  end

  % ---------- 2,3) C0 and C1 ----------
  hdr('C0 and C1 (one script)', 2);
  try
    r = verify_left_region;
    pass(2) = logical(r.pass); marg(2) = r.Delta_0_min_lower;
    pass(3) = logical(r.pass); marg(3) = r.Delta_tail_min_lower;
  catch err
    fprintf('\n  ERROR: %s\n', err.message);
  end

  % ---------- 4) C2 ----------
  hdr('C2', 4);
  try
    r = verify_right_region;
    pass(4) = logical(r.pass); marg(4) = r.prefixMargin;
  catch err
    fprintf('\n  ERROR: %s\n', err.message);
  end

  % ---------- 5) C2' ----------
  hdr('C2''', 5);
  try
    r = certify_u_monotone;
    pass(5) = logical(r.pass); marg(5) = r.worstGap;
  catch err
    fprintf('\n  ERROR: %s\n', err.message);
  end

  % ---------- 6) Certificate 1 ----------
  hdr('Certificate 1', 6);
  try
    r = certify_lambda_monotone;
    pass(6) = logical(r.pass); marg(6) = r.worstGap;
  catch err
    fprintf('\n  ERROR: %s\n', err.message);
  end

  % ---------- 7) N1 ----------
  hdr('N1 (floating-point)', 7);
  try
    r = verify_l_monotone;
    pass(7) = logical(r.pass); marg(7) = r.min_diff;
  catch err
    fprintf('\n  ERROR: %s\n', err.message);
  end

  % ---------- summary ----------
  fprintf('\n============================================================\n');
  fprintf(' SUMMARY\n');
  fprintf('============================================================\n');
  fprintf('  %-8s %-27s %-9s %-13s %s\n', 'ID', 'File', 'Type', 'MinMargin', 'Result');
  fprintf('  %s\n', repmat('-', 1, 70));
  for k = 1:n
    if isnan(marg(k))
      ms = '-';
    else
      ms = sprintf('%.6g', marg(k));
    end
    fprintf('  %-8s %-27s %-9s %-13s %s\n', names{k}, files{k}, types{k}, ms, pf(pass(k)));
  end

  allpass = all(pass);
  fprintf('\n  OVERALL: %s\n', pf(allpass));
  if ~allpass
    fprintf('\n  Some certificate did not report pass=true. Check that script''s own output.\n');
  end
end

function hdr(name, idx)
  fprintf('\n------------------------------------------------------------\n');
  fprintf('[%d/7] %s\n', idx, name);
  fprintf('------------------------------------------------------------\n');
end

function s = pf(ok)
  if ok, s = 'PASS'; else, s = 'FAIL'; end
end
