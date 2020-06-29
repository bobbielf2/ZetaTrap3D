function [A, An] = Helm3dDLPmat(t,s,k)
% HELM3DDLPMAT.  dense Helmholtz DLP Nystrom eval matrix from sources to targets
% 
% [A, An] = Helm3dDLPmat(t,s,k)
% 
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*N) nodes, and, if An requested, nx (3*N) normals
%  k - wavenumber, satisfies imag(k)>0
% Outputs:
%  A - (M*N) matrix getting potentials from density values
%  An - (M*N) matrix getting t.nx directional derivs from density values
%
% Bowei Wu 3/13/20


d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
ddotsn = bsxfun(@times,d1,s.nx(1,:))+bsxfun(@times,d2,s.nx(2,:))+bsxfun(@times,d3,s.nx(3,:));
A = bsxfun(@times, (1-1i*k*r).*exp(1i*k*r).*ddotsn./(r.*rr), s.w*(1/4/pi));  % including src quadr wei
if size(s.x,2) == size(t.x,2) && norm(s.x - t.x) < 1e-14
    A(diagind(A)) = 0;
end
if nargout>1                  % targ deriv wanted
  ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
  tndotsn = bsxfun(@times,s.nx(1,:),t.nx(1,:)')+bsxfun(@times,s.nx(2,:),t.nx(2,:)')+bsxfun(@times,s.nx(3,:),t.nx(3,:)');
  An = bsxfun(@times, exp(1i*k*r)./(r.*rr).*((k^2-3*(1-1i*k*r)./rr).*ddotsn.*ddottn+(1-1i*k*r).*tndotsn), s.w*(1/4/pi));  % dipole deriv, incl src quad wei
  if size(s.x,2) == size(t.x,2) && norm(s.x - t.x) < 1e-14
      An(diagind(An)) = 0;
  end
end

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
%
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);
