function [A, An] = Helm3dSLPmat(t,s,k,if_jac)
% HELM3DSLPMAT.  dense Helmholtz SLP Nystrom eval matrix from sources to targets
% 
% [A, An] = Helm3dSLPmat(t,s,k)
% 
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*N) nodes, and, if An requested, nx (3*N) normals
%  k - wavenumber, satisfies imag(k)>0
%  if_jac - whether or not to include the jacobian (speed) s.sp
% Outputs:
%  A - (M*N) matrix getting potentials from density values
%  An - (M*N) matrix getting t.nx directional derivs from density values
%
% Bowei Wu 3/13/20

if nargin<4 || isempty(if_jac), if_jac = 1; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
if if_jac
    A = bsxfun(@times, exp(1i*k*r)./r, s.w*(1/4/pi));          % including src quadr wei
else
    A = bsxfun(@times, exp(1i*k*r)./r, (s.w./s.sp)*(1/4/pi));  % including src quadr wei, excl jacobian
end
if size(s.x,2) == size(t.x,2) && norm(s.x - t.x) < 1e-14
    A(diagind(A)) = 0;
end
if nargout>1                  % targ deriv wanted
  ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
  if if_jac
      An = bsxfun(@times, (1-1i*k*r).*exp(1i*k*r).*ddottn./(sqrt(rr).*rr), s.w*(-1/4/pi));          % monopole deriv, incl src quad wei
  else
      An = bsxfun(@times, (1-1i*k*r).*exp(1i*k*r).*ddottn./(sqrt(rr).*rr), (s.w./s.sp)*(-1/4/pi));  % monopole deriv, incl src quad wei, excl jacobian
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
