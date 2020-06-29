function [A, An] = Lap3dDLPmat(t,s)
% LAP3DDLPMAT.  dense Laplace DLP Nystrom eval matrix from sources to targets
%
% [A An] = Lap3dDLPmat(t,s)
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, nx (3*N) unit normals, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*N) nodes, and, if An requested, nx (3*N) normals
% Outputs:
%  A - (M*N) matrix getting potentials from density values
%  An - (M*N) matrix getting t.nx directional derivs from density values

% Modified from the same func in Alex Barnett's BIE3D package

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
ny = (1/4/pi) * s.nx;     % apply prefactor here, cheaper (don't use s.nx now!)
ddotsn = bsxfun(@times,d1,ny(1,:))+bsxfun(@times,d2,ny(2,:))+bsxfun(@times,d3,ny(3,:));  % M*N
A =  bsxfun(@times, ddotsn ./ (sqrt(rr).*rr), s.w);  % including src quadr wei
if size(s.x,2) == size(t.x,2) && norm(s.x - t.x) < 1e-14
    A(diagind(A)) = 0;
end
if nargout>1                  % targ deriv wanted ... not the fastest
  ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
  tndotsn = bsxfun(@times,t.nx(1,:)',ny(1,:)) + bsxfun(@times,t.nx(2,:)',ny(2,:)) + bsxfun(@times,t.nx(3,:)',ny(3,:));
  An = bsxfun(@times, (1./(sqrt(rr).*rr)).*(-3./rr.*ddottn.*ddotsn + tndotsn), s.w);  % dipole deriv, incl src quad wei
  if size(s.x,2) == size(t.x,2) && norm(s.x - t.x) < 1e-14
    An(diagind(An)) = 0;
  end
end

function i = diagind(A)
% DIAGIND  Return indices of diagonal of square matrix
%
% Example usage:   A = randn(3,3); A(diagind(A)) = 0;
N = size(A,1); i = sub2ind([N,N], 1:N, 1:N); i = i(:);
