% DEMO
% zeta trapezoidal quadrature for Laplace potentials on a surface patch 
% in 3D, using ONLY the first fundamental form (instead of requiring 
% higher-order geometric derivatives as in the first version)

% Bowei Wu 11/29/2020

% generate random quartic patch by sampling a point on a random surface
m = randi(4); % petal number of generating curve
n = randi(5); % twist number along toroidal direction
a = max(0.5*(1+0.3*randn()),0);
s = wobblytorus(m,n,a); % toroidal surface
u = (rand()*2-1)*pi;
v = (rand()*2-1)*pi;
% derivatives at the sampled point (for surface patch construction only)
[~,ru,rv,ruu,ruv,rvv,ruuu,ruuv,ruvv,rvvv,...
    ruuuu,ruuuv,ruuvv,ruvvv,rvvvv] = s.get_x_xp_xpp_xppp(u,v);

% compute first fundamental form
E=dot(ru,ru); F=dot(ru,rv); G=dot(rv,rv);

% define density function (with compact support)
rpat = 0.5; % radius of the patch
a = randn(); b = randn();
f = @(u,v) (a*cos(a+u)+b*sin(b+v)).* exp(-40*((u.^2+v.^2)/rpat^2).^4); 

% convergence test
ka = 10*(randn()+1i*randn()); % random wavenumber
hh = 0.5.^([3:8,10]);
val_d = []; val_s = []; val_sn = [];
for h = hh
    % parametric grid & density func
    n = ceil(rpat/h);
    [u,v] = ndgrid(h*(-n:n)); u = u(:).'; v = v(:).';
    sig = f(u,v).';
    
    % define geometry (quartic surface patch centered at 0)
    [rvec,r_u,r_v] = quarticSurf(u,v,ru,rv,ruu,ruv,rvv,ruuu,ruuv,...
                                ruvv,rvvv,ruuuu,ruuuv,ruuvv,ruvvv,rvvvv);	
	% plot surface
   	if h == 2^-5, figure(1); subplot(2,1,1), PlotPatch(n,rvec,sig,ru,rv); end 
    
    % define kernels via punctured trapezoidal rule
    r2 = dot(rvec,rvec);  	% r^2
    r = sqrt(r2);          	% r
    ker_s = exp(1i*ka*r)./r*h^2;	% SLP
    ker_d0 = exp(1i*ka*r).*(1-1i*ka*r)./(r2.*r)*h^2; % DLP & SLPn
    ind = find(r2==0);     	% location index of singular pt
    ker_s(ind) = 1i*ka*h^2;         % constant term correction
    ker_d0(ind) = 0;                % punctured trapezoidal rule
    
    % zeta correction
    ord = 5;                    % desired order of convergence
    Q = E*u.^2+2*F*u.*v+G*v.^2; % 1st fundamental form
    r2mQ = r2-Q;                % smooth factors for product integration
    if 0    % use specific r, 1/r, 1/r^3 corrections
        % SLP: 1/r+1i*ka-ka^2*r/2+...
        ker_s = zeta_correct_r1(ker_s,ind,h,r2mQ,1,ord,E,F,G);      % 1/r
        ker_s = zeta_correct_r(ker_s,ind,h,r2mQ,-ka^2/2,ord,E,F,G);	% -k^2/2*r
        % DLP & SLPn: 1/r^3+ka^2/2r+...
        offset = 1; % smoothness offset, kernel has an extra O(h^2) smoothness
        ker_d0 = zeta_correct_r3(ker_d0,ind,h,r2mQ,1,ord,E,F,G,offset);    	% 1/r^3
        ker_d0 = zeta_correct_r1(ker_d0,ind,h,r2mQ,ka^2/2,ord,E,F,G,offset);% k^2/2/r
    else    % use the generic r^p correction
        % SLP: 1/r+1i*ka-ka^2*r/2+...
        p = -1; c = 1;
        ker_s = zeta_correct_rp(ker_s,p,ind,h,r2mQ,c,ord,E,F,G);  	% 1/r
        p = 1; c = -ka^2/2;
        ker_s = zeta_correct_rp(ker_s,p,ind,h,r2mQ,c,ord,E,F,G);    % -k^2/2*r
        % DLP & SLPn: 1/r^3+ka^2/2r+...
        offset = 1; % smoothness offset, kernel has an extra O(h^2) smoothness
        p = -3; c = 1;
        ker_d0 = zeta_correct_rp(ker_d0,p,ind,h,r2mQ,c,ord,E,F,G,offset);	% 1/r^3
        p = -1; c = ka^2/2;
        ker_d0 = zeta_correct_rp(ker_d0,p,ind,h,r2mQ,c,ord,E,F,G,offset);  % k^2/2/r
    end

    % multiply back smooth components for each kernel
    nJ = cross(r_u,r_v);% normal vec r_u x r_v (not normalized)
    J = vecnorm(nJ);    % jacobian
    nJ0 = cross(ru,rv); % normal vec at 0
    n0 = nJ0/norm(nJ0); % unit normal at 0
    ker_s = J.*ker_s;                   % SLP: jacobian
    ker_d = -dot(rvec,nJ).*ker_d0;      % DLP: -r.(r_u x r_v)
    ker_sn=sum(rvec.*n0,1).*J.*ker_d0;  % SLPn: (r.n0)*jacobian
    
    
    val_s = [val_s; ker_s*sig];     % SLP vals
    val_d = [val_d; ker_d*sig];     % DLP vals
    val_sn = [val_sn; ker_sn*sig];  % SLP normal gradient 
end

% plot err
err_s = abs(val_s(1:end-1)-val_s(end));
err_d = abs(val_d(1:end-1)-val_d(end));
err_sn= abs(val_sn(1:end-1)-val_sn(end));
hh = hh(1:end-1).';
subplot(2,1,2)
loglog(hh,hh.^ord/hh(end)^ord*err_d(end)*2,'k--',hh,err_s,'o-',hh,err_d,'v-',hh,err_sn,'*-');
title(sprintf('Helmholtz $\\kappa=%.2f%+.2fi$',real(ka),imag(ka)),'interpreter','latex')
xlabel('$h$','interpreter','latex')
ylabel('error','interpreter','latex','rotation',90)
legend({['$O(h^',num2str(ord),')$'],'SLP','DLP','SLP grad'},'interpreter','latex','location','nw')

% display convergence order
fprintf(['Numerical orders of convergence:\n at h =',repmat(' %.1e,',1,numel(hh)),'\n'],hh)
fprintf([' DLP:\t', repmat('\t%3.1f,  ',1,numel(err_d)-1),'\n'],log2(err_d(1:end-1)./err_d(2:end)))
fprintf([' SLP:\t', repmat('\t%3.1f,  ',1,numel(err_s)-1),'\n'],log2(err_s(1:end-1)./err_s(2:end)))
fprintf([' SLP grad: ', repmat('\t%3.1f,  ',1,numel(err_sn)-1),'\n'],log2(err_sn(1:end-1)./err_sn(2:end)))

function [rvec,r_u,r_v] = quarticSurf(u,v,...
                    ru,rv,ruu,ruv,rvv,ruuu,ruuv,ruvv,rvvv,...
                  	ruuuu,ruuuv,ruuvv,ruvvv,rvvvv)
% construct a quartic surface centered at (0,0) from the given derivatives
% Define quartic surface, r
d1 = ru.*u+rv.*v;
d2 = (ruu.*u.^2+2*ruv.*u.*v+rvv.*v.^2)/2;
d3 = (ruuu.*u.^3+3*ruuv.*(u.^2.*v)+3*ruvv.*(u.*v.^2)+rvvv.*v.^3)/6;
d4 = (ruuuu.*u.^4+4*ruuuv.*(u.^3.*v)+6*ruuvv.*(u.^2.*v.^2)+4*ruvvv.*(u.*v.^3)+rvvvv.*v.^4)/24;
rvec = d1 + d2 + d3 + d4;
% Define derivatives of the surface, dr/du, dr/dv
dr2u = ruu.*u+ruv.*v;
dr2v = ruv.*u+rvv.*v;
dr3u = (ruuu.*u.^2+2*ruuv.*(u.*v)+ruvv.*v.^2)/2;
dr3v = (ruuv.*u.^2+2*ruvv.*(u.*v)+rvvv.*v.^2)/2;
dr4u = (ruuuu.*u.^3+3*ruuuv.*(u.^2.*v)+3*ruuvv.*(u.*v.^2)+ruvvv.*v.^3)/6;
dr4v = (ruuuv.*u.^3+3*ruuvv.*(u.^2.*v)+3*ruvvv.*(u.*v.^2)+rvvvv.*v.^3)/6;
r_u = ru+dr2u+dr3u+dr4u;
r_v = rv+dr2v+dr3v+dr4v;
end

function ker = zeta_correct_r(ker,ind,h,r2mQ,c,ord,E,F,G,offset)
% zeta correction for r kernel on surface
% Input:
%   ker = kernel matrix via punctured trapezoidal rule
%   ind = index location of the singular point
%   h = mesh spacing
%   r2mQ = r^2-Q
%   c = canstant prefactor for the kernel correction
%   ord = desired order of convergence
%   E,F,G = first fundamental form coeffs
%   offset = stencil offset parameter, assume an O(h^(2*offset)) extra
%            smoothness for the non-singular part of the kernel
if nargin<10, offset=0; end
M = ceil((ord-5)/2)-offset;
if M<0, return; end     % no correction needed
siz = sqrt(numel(ker))*[1,1];   % mesh size
[sub1,sub2] = ind2sub(siz,ind); % subscripts of singular point
for m = 0:2*M
    l1 = ceil(3*m/2)+offset;% stencil inner layer = l1
    l2 = M+m+offset;     	% stencil outer layer = l2+1
    Qpow = m-0.5;           % power of the quadratic form
    fac = c*binom(1/2,m)*h^(3-2*m)*r2mQ.^m;	% compute smooth factors
    tau = zeta_weights(l1,l2,Qpow,E,F,G);   % compute correction weights
    ker = zeta_correct(ker,sub1,sub2,l1,l2,fac,tau); % apply correction
end
end

function ker = zeta_correct_r1(ker,ind,h,r2mQ,c,ord,E,F,G,offset)
% zeta correction for 1/r kernel on surface
% Input:
%   ker = kernel matrix via punctured trapezoidal rule
%   ind = index location of the singular point
%   h = mesh spacing
%   r2mQ = r^2-Q
%   c = canstant prefactor for the kernel correction
%   ord = desired order of convergence
%   E,F,G = first fundamental form coeffs
%   offset = stencil offset parameter, assume an O(h^(2*offset)) extra
%            smoothness for the non-singular part of the kernel
if nargin<10, offset=0; end
M = ceil((ord-3)/2)-offset;
if M<0, return; end
siz = sqrt(numel(ker))*[1,1];   % mesh size
[sub1,sub2] = ind2sub(siz,ind); % subscripts of singular point
for m = 0:2*M
    l1=ceil(3*m/2)+offset;  % stencil inner layer = l1
    l2=M+m+offset;       	% stencil outer layer = l2+1
    Qpow = m+0.5;           % power of the quadratic form
    fac = c*binom(-1/2,m)*h^(1-2*m)*r2mQ.^m;% compute smooth factors
    tau = zeta_weights(l1,l2,Qpow,E,F,G);   % compute correction weights
    ker = zeta_correct(ker,sub1,sub2,l1,l2,fac,tau); % apply correction
end
end

function ker = zeta_correct_r3(ker,ind,h,r2mQ,c,ord,E,F,G,offset)
% zeta correction for 1/r^3 kernel on surface
% Input:
%   ker = kernel matrix via punctured trapezoidal rule
%   ind = index location of the singular point
%   h = mesh spacing
%   r2mQ = r^2-Q
%   c = canstant prefactor for the kernel correction
%   ord = desired order of convergence
%   E,F,G = first fundamental form coeffs
%   offset = stencil offset parameter, assume an O(h^(2*offset)) extra
%            smoothness for the non-singular part of the kernel
if nargin<10, offset=0; end
M = ceil((ord-1)/2)-offset;
if M<0, return; end
siz = sqrt(numel(ker))*[1,1];   % mesh size
[sub1,sub2] = ind2sub(siz,ind); % subscripts of singular point
for m = 0:2*M
    l1=ceil(3*m/2)+offset;	% stencil inner layer = l1
    l2=M+m+offset;        	% stencil outer layer = l2+1
    Qpow = m+1.5;           % power of the quadratic form
    fac = c*binom(-3/2,m)*h^(-1-2*m)*r2mQ.^m;   % compute smooth factors
    tau = zeta_weights(l1,l2,Qpow,E,F,G);       % compute correction weights
    ker = zeta_correct(ker,sub1,sub2,l1,l2,fac,tau); % apply correction
end
end

function ker = zeta_correct_rp(ker,p,ind,h,r2mQ,c,ord,E,F,G,offset)
% zeta correction for r^p kernel on surface
% Input:
%   ker = kernel matrix via punctured trapezoidal rule
%   p   = power in r^p
%   ind = index location of the singular point
%   h = mesh spacing
%   r2mQ = r^2-Q
%   c = canstant prefactor for the kernel correction
%   ord = desired order of convergence
%   E,F,G = first fundamental form coeffs
%   offset = stencil offset parameter, assume an O(h^(2*offset)) extra
%            smoothness for the non-singular part of the kernel
if mod(p,2)==0 && p >= 0, warning('p=%d,r^p is smooth'); return; end
if nargin<11, offset=0; end
% calculate M, such that 2M+1 terms in the kernel expansion needs correction
M = ceil((ord-p)/2)-2-offset;
if M < 0, return; end
siz = sqrt(numel(ker))*[1,1];   % mesh size
[sub1,sub2] = ind2sub(siz,ind); % subscripts of singular point
for m = 0:2*M               % correct 2M+1 terms in the kernel expansion
    l1=ceil(3*m/2)+offset;	% stencil inner layer = l1
    l2=M+m+offset;        	% stencil outer layer = l2+1
    Qpow = m-p/2;           % power of the quadratic form in the m-th term
    fac = c*binom(p/2,m)*h^(2+p-2*m)*r2mQ.^m; 	% compute smooth factors
    tau = zeta_weights(l1,l2,Qpow,E,F,G);       % compute correction weights
    ker = zeta_correct(ker,sub1,sub2,l1,l2,fac,tau); % apply correction
end
end

function ker = zeta_correct(ker,sub1,sub2,l1,l2,fac,tau)
% Generic zeta correction for a kernel matrix

nt = numel(tau);              % correction weights DoF
siz = sqrt(numel(ker))*[1,1]; % mesh size
% 1. inner layer (symmetric about axes)
u=l1:-1:0; v=0:l1;
taus=tau(1:l1+1);
if l1>1, u=[u,-l1+1:-1]; v=[v,1:l1-1]; taus=[taus,tau(2:l1)]; end % off-axis copies
if l1>0, u=[u,-u]; v=[v,-v]; taus = kron([1,1],taus); end % central symmetry copies
inds = sub2ind(siz,sub1+u,sub2+v);      % kernel indices
ker(inds) = ker(inds)+taus.*fac(inds); 
% 2. intermediate layers (symmetric aboutt center)
if l2 > l1
    u=[]; v=[]; % assemble stencil layers l, l1<l<=l2
    for l=l1+1:l2, u=[u,(l:-1:-l+1)]; v=[v,0:l-1,l:-1:1]; end
    inds=sub2ind(siz,sub1+[u,-u],sub2+[v,-v]);	% kernel indices
    taus = kron([1,1],tau(l1+2:nt-l2));         % tau copies
    ker(inds) = ker(inds) + taus.*fac(inds);
end
% 3. outer layer (anti-symmetric about axes)
if l2 > 0
    u=[l2:-1:1,-l2:-1]'; v=[1:l2,1:l2]';
    inds=sub2ind(siz,sub1+[u;-u],sub2+[v;-v]);	% kernel indices
    taus=kron([1,-1,1,-1],tau(nt-l2+1:nt));     % tau copies
    ker(inds) = ker(inds) + taus.*fac(inds);
end
end

function tau = zeta_weights(l1,l2,Qpow,E,F,G)
% Compute zeta correction weights

% form monomial powers (a,b) for each row u^a * v^b
% such that a+b=2*l, for all l1<=l<=l2
a = []; b = [];
for p=2*(l1:l2), a=[a;(p:-1:0)']; b=[b;(0:p)']; end
% Construct fitting matrix A on the (u,v) reduced stencil 
n = numel(a); A = zeros(n);
% 1. inner layer (symmetric about axes)
u=l1:-1:0; v=0:l1;
A(:,1:l1+1)=u.^a.*v.^b;
if l1>0, A(:,1:l1+1)=2*A(:,1:l1+1); end % on-axis pts
if l1>1, A(:,2:l1)=A(:,2:l1)+2*(-u(2:l1)).^a.*v(2:l1).^b; end % off-axis pts
% 2. intermediate layers (symmetric aboutt center)
if l2 > l1
    u=[]; v=[]; % assemble stencil layers l, l1<l<=l2
    for l=l1+1:l2, u=[u,l:-1:-l+1]; v=[v,0:l-1,l:-1:1]; end 
    A(:,l1+2:n-l2)=2*u.^a.*v.^b;
end
% 3. outer layer (anti-symmetric about axes)
if l2>0, u=l2:-1:1; v=1:l2; A(:,n-l2+1:n)=2*(u.^a.*v.^b-(-u).^a.*v.^b); end
% Compute Wigner limits
W = [];
for l = l1:l2
    s = Qpow-l;
    fallpow = 1/prod(-s-l+1:-s); % Knuth falling factorial
    Z = epstein_zeta_deriv(2*s,l,E,F,G); % l-th derivatives of zeta fun
    W = [W;-fallpow*Z];
end
tau = double(sym(A)\W).'; % solve for zeta weights (row vector)
end

function c = binom(n,k)
% generalized binomial coeff, n = any real number, k = positive integer
c = prod(n-k+1:n)/gamma(k+1); 
end

function sq = epstein_zeta_deriv(s,d,E,F,G)
% Compute all possible d-th parametric partial derivative of 
% the Epstein zeta function, \partial^d Z(s;E,F,G)
if d == 0
    sq = epstein_zeta(s,E,F,G,0,0,0);
elseif d == 1
    [~, d1_E] = epstein_zeta(s,E,F,G,1,0,0);    % TODO: combine these
    [~, d1_G] = epstein_zeta(s,E,F,G,0,0,1);
    [~, d1_F] = epstein_zeta(s,E,F,G,0,0.5,0);
    sq = [d1_E; d1_F; d1_G];
elseif d == 2
    [~, ~, d2_E] = epstein_zeta(s,E,F,G,1,0,0); % TODO: combine these
    [~, ~, d2_G] = epstein_zeta(s,E,F,G,0,0,1);
    [~, ~, d2_F] = epstein_zeta(s,E,F,G,0,0.5,0);
    [~,~,d2_EpF] = epstein_zeta(s,E,F,G,1,.5,0);
    [~,~,d2_GpF] = epstein_zeta(s,E,F,G,0,.5,1);
    sq = [d2_E; (d2_EpF-d2_E-d2_F)/2; d2_F; (d2_GpF-d2_F-d2_G)/2; d2_G];
elseif d == 3
    [~,~,~, d3_E] = epstein_zeta(s,E,F,G,1,0,0);
    [~,~,~, d3_G] = epstein_zeta(s,E,F,G,0,0,1);
    [~,~,~, d3_F] = epstein_zeta(s,E,F,G,0,0.5,0);
    [~,~,~,d3_EpF] = epstein_zeta(s,E,F,G,1,.5,0);
    [~,~,~,d3_EmF] = epstein_zeta(s,E,F,G,1,-.5,0);
    [~,~,~,d3_GpF] = epstein_zeta(s,E,F,G,0,.5,1);
    [~,~,~,d3_GmF] = epstein_zeta(s,E,F,G,0,-.5,1);
    sq = [d3_E; (d3_EpF-d3_EmF-2*d3_F)/6; (d3_EpF+d3_EmF-2*d3_E)/6;...
          d3_F; (d3_GpF+d3_GmF-2*d3_G)/6; (d3_GpF-d3_GmF-2*d3_F)/6; d3_G];
elseif d == 4
    [~,~,~,~,   d4_E] = epstein_zeta(s,E,F,G,1,0,0);
    [~,~,~,~,   d4_G] = epstein_zeta(s,E,F,G,0,0,1);
    [~,~,~,~,   d4_F] = epstein_zeta(s,E,F,G,0,0.5,0);
    [~,~,~,~, d4_EpF] = epstein_zeta(s,E,F,G,1,.5,0);
    [~,~,~,~,d4_Ep2F] = epstein_zeta(s,E,F,G,1,1,0);
    [~,~,~,~, d4_EmF] = epstein_zeta(s,E,F,G,1,-.5,0);
    [~,~,~,~, d4_GpF] = epstein_zeta(s,E,F,G,0,.5,1);
    [~,~,~,~,d4_Gp2F] = epstein_zeta(s,E,F,G,0,1,1);
    [~,~,~,~, d4_GmF] = epstein_zeta(s,E,F,G,0,-.5,1);
    sq = [d4_E;...
        (6*d4_EpF-2*d4_EmF-d4_Ep2F-3*d4_E+12*d4_F)/24;...
        (d4_EpF+d4_EmF-2*d4_E-2*d4_F)/12;...
        (d4_Ep2F-3*d4_EpF-d4_EmF+3*d4_E-12*d4_F)/24;...
        d4_F;...
        (d4_Gp2F-3*d4_GpF-d4_GmF+3*d4_G-12*d4_F)/24;...
        (d4_GpF+d4_GmF-2*d4_G-2*d4_F)/12;...
        (6*d4_GpF-2*d4_GmF-d4_Gp2F-3*d4_G+12*d4_F)/24;
        d4_G];
else
    error('%d-th derivatives not yet implemented!',d)
end
end

function PlotPatch(n,rvec,sig,ru,rv)
scatter3(0,0,0,100,'k','filled');
hold on;
quiver3(0,0,0,ru(1),ru(2),ru(3),0.5,'linewidth',2)  % r_u vector
quiver3(0,0,0,rv(1),rv(2),rv(3),0.5,'linewidth',2)  % r_v vector
surf(reshape(rvec(1,:),2*n+1,2*n+1),reshape(rvec(2,:),2*n+1,2*n+1),reshape(rvec(3,:),2*n+1,2*n+1),reshape(sig,2*n+1,2*n+1))
hold off
axis equal
colorbar
title('random quartic patch','interpreter','latex')
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$z$','interpreter','latex')
legend({'target','$\mathbf{r}_u$','$\mathbf{r}_v$'},'interpreter','latex')
end
