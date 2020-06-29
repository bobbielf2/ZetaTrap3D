function [A, An] = Lap3dSLPmat(t,s,if_jac)
% LAP3DSLPMAT.  dense Laplace SLP Nystrom eval matrix from sources to targets
%
% [A An] = Lap3dSLPmat(t,s)
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*N) nodes, and, if An requested, nx (3*N) normals
%  if_jac - whether or not to include the jacobian (speed) s.sp
% Outputs:
%  A - (M*N) matrix getting potentials from density values
%  An - (M*N) matrix getting t.nx directional derivs from density values

% Modified from the same func in Alex Barnett's BIE3D package

if nargin<3 || isempty(if_jac), if_jac = 1; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
ir = 1./sqrt(rr);
if if_jac
    A = bsxfun(@times, ir, s.w*(1/4/pi));           % including src quadr wei
else % don't include jacobian
    A = bsxfun(@times, ir, (s.w./s.sp)*(1/4/pi));   % including src quadr wei, excl jacobian
end
if size(s.x,2) == size(t.x,2) && norm(s.x - t.x) < 1e-14
    A(diagind(A)) = 0;
end
if nargout>1                  % targ deriv wanted
  ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
  if if_jac
      An = bsxfun(@times, ddottn./rr.*ir, s.w*(-1/4/pi));       % monopole deriv, incl src quad wei
  else % don't include jacobian
      An = bsxfun(@times, ddottn./rr.*ir, s.w./s.sp*(-1/4/pi)); % monopole deriv, incl src quad wei, excl jacobian
  end
  if size(s.x,2) == size(t.x,2) && norm(s.x - t.x) < 1e-14
      An(diagind(An)) = 0;
  end
end

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
%
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);
