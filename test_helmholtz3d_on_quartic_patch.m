% O(h^5) corrected trapezoidal rule for Helmholtz double-layer potential,
% single-layer potential, and its normal gradient, on a surface patch in 3D
% See the manuscript [1], Example 1, for more details.
%
% [1] Wu, B., & Martinsson, P.G. (2020, arXiv:2007.02512). Corrected
%     Trapezoidal Rules for Boundary Integral Equations in Three
%     Dimensions.

% generate random quartic patch
m = randi(4); % petal number of generating curve
n = randi(5); % twist number along toroidal direction
a = max(0.5*(1+0.3*randn()),0);
% sample a random point on a surface
s = wobblytorus(m,n,a);
u = (rand()*2-1)*pi;
v = (rand()*2-1)*pi;
[~,ru,rv,ruu,ruv,rvv,ruuu,ruuv,ruvv,rvvv,...
    ruuuu,ruuuv,ruuvv,ruvvv,rvvvv] = s.get_x_xp_xpp_xppp(u,v);

% pre-compute correction weights (for all 3 types of potentials)
[Cd0,Cd01,Cd02,Cd1,Cd2,Cd3,Cd4,Cd5,...
    Cs0,Cs01,Cs02,Cs1,Cs2,Cs3,Cs4,Cs5,...
    Csn01,Csn1,Csn2] = wigner_type_lims(...
                                ru,rv,ruu,ruv,rvv,ruuu,ruuv,ruvv,rvvv,...
                                ruuuu,ruuuv,ruuvv,ruvvv,rvvvv);

% define density function (with compact support)
rpat = 0.5; % radius of the patch
a = randn(); b = randn();
f = @(u,v) (a*cos(a+u)+b*sin(b+v)).* exp(-40*((u.^2+v.^2)/rpat^2).^4); 

% convergence test
ka = 10*(randn()+1i*randn()); % random wavenumber
hh = 0.5.^([2:8,10]);
val_d = []; val_s = []; val_sn = [];
for h = hh
    % parametric grid & density func
    n = ceil(rpat/h);
    [u,v] = ndgrid(h*(-n:n)); u = u(:).'; v = v(:).';
    sig = f(u,v).';
    
    % define geometry (quartic surface)
    d1 = ru.*u+rv.*v;
    d2 = (ruu.*u.^2+2*ruv.*u.*v+rvv.*v.^2)/2;
    d3 = (ruuu.*u.^3+3*ruuv.*(u.^2.*v)+3*ruvv.*(u.*v.^2)+rvvv.*v.^3)/6;
    d4 = (ruuuu.*u.^4+4*ruuuv.*(u.^3.*v)+6*ruuvv.*(u.^2.*v.^2)+4*ruvvv.*(u.*v.^3)+rvvvv.*v.^4)/24;
    rvec = d1 + d2 + d3 + d4;
    r2 = dot(rvec,rvec);
    
    % define r.(ru x rv), etc.
    d2u = ruu.*u+ruv.*v; d2v = ruv.*u+rvv.*v;
    d3u = (ruuu.*u.^2+2*ruuv.*(u.*v)+ruvv.*v.^2)/2;
    d3v = (ruuv.*u.^2+2*ruvv.*(u.*v)+rvvv.*v.^2)/2;
    d4u = (ruuuu.*u.^3+3*ruuuv.*(u.^2.*v)+3*ruuvv.*(u.*v.^2)+ruvvv.*v.^3)/6;
    d4v = (ruuuv.*u.^3+3*ruuvv.*(u.^2.*v)+3*ruvvv.*(u.*v.^2)+rvvvv.*v.^3)/6;
    dru = ru+d2u+d3u+d4u;
    drv = rv+d2v+d3v+d4v;
    nj  = cross(dru,drv);
    rnj = dot(rvec,nj);
    j   = vecnorm(nj);
    nj0 = cross(ru,rv);
    j0  = norm(nj0);
    rn0= sum(rvec .* (nj0/j0));
    
    if h == 2^-5, PlotPatch(n,rvec,sig,ru,rv); end % plot surface
    
    % define kernels
    r = sqrt(r2);
    ker_d = -rnj./(r2.*r).*exp(1i*ka*r).*(1-1i*ka*r)*h^2; % Helm DLP
    ker_s = exp(1i*ka*r)./r*h^2; % Helm SLP
    ker_sn= rn0./(r2.*r).*exp(1i*ka*r).*(1-1i*ka*r)*h^2; % Helm SLP normal gradient
    
    % 9-point correction of 2-periodic trapezoidal rule, O(h^5) accuracy
    
    % locate indices of 9-point stencil
    ind = find(r2==0);                              %( 0, 0)
    siz = [2*n+1,2*n+1];
    [ind1,ind2] = ind2sub(siz,ind);
    indup = sub2ind(siz,ind1+1,ind2);               %(+1, 0)
    indum = sub2ind(siz,ind1-1,ind2);               %(-1, 0)
    indvp = sub2ind(siz,ind1,ind2+1);               %( 0,+1)
    indvm = sub2ind(siz,ind1,ind2-1);               %( 0,-1)
    induvp = sub2ind(siz,ind1+[1,-1],ind2+[1,-1]);  %(+1,+1) and (-1,-1)
    induvm = sub2ind(siz,ind1+[1,-1],ind2+[-1,1]);  %(+1,-1) and (-1,+1)
    
    % correct DLP
    ker_d(ind) = ((Cd0-Cd3-Cd4)+(Cd01+Cd02*ka^2/2)*h^2)*h;  %( 0, 0)
    ker_d(indup) = ker_d(indup)+(Cd3+Cd1*h)*h/2;            %(+1, 0)
    ker_d(indum) = ker_d(indum)+(Cd3-Cd1*h)*h/2;            %(-1, 0)
    ker_d(indvp) = ker_d(indvp)+(Cd4+Cd2*h)*h/2;            %( 0,+1)
    ker_d(indvm) = ker_d(indvm)+(Cd4-Cd2*h)*h/2;            %( 0,-1)
    ker_d(induvp) = ker_d(induvp)+Cd5*h/4;                  %(+1,+1) and (-1,-1)
    ker_d(induvm) = ker_d(induvm)-Cd5*h/4;                  %(+1,-1) and (-1,+1)
    
    % do the same for SLP
    ker_s(ind) = ((Cs0-Cs3-Cs4)+1i*ka*h+(Cs01+Cs02*ka^2/2)*h^2)*h;	%( 0, 0)
    ker_s(indup) = ker_s(indup)+(Cs3+Cs1*h)*h/2;                    %(+1, 0)
    ker_s(indum) = ker_s(indum)+(Cs3-Cs1*h)*h/2;                    %(-1, 0)
    ker_s(indvp) = ker_s(indvp)+(Cs4+Cs2*h)*h/2;                    %( 0,+1)
    ker_s(indvm) = ker_s(indvm)+(Cs4-Cs2*h)*h/2;                    %( 0,-1)
    ker_s(induvp) = ker_s(induvp)+Cs5*h/4;                          %(+1,+1) and (-1,-1)
    ker_s(induvm) = ker_s(induvm)-Cs5*h/4;                          %(+1,-1) and (-1,+1)
    
    % do the same for SLP-grad
    ker_sn(ind) = ((Cd0-Cd3-Cd4)+(Csn01+Cd02*ka^2/2)*h^2)*h/j0;	%( 0, 0)
    ker_sn(indup) = ker_sn(indup)+(Cd3+Csn1*h)*h/2/j0;         	%(+1, 0)
    ker_sn(indum) = ker_sn(indum)+(Cd3-Csn1*h)*h/2/j0;          %(-1, 0)
    ker_sn(indvp) = ker_sn(indvp)+(Cd4+Csn2*h)*h/2/j0;          %( 0,+1)
    ker_sn(indvm) = ker_sn(indvm)+(Cd4-Csn2*h)*h/2/j0;          %( 0,-1)
    ker_sn(induvp) = ker_sn(induvp)+Cd5*h/4/j0;                 %(+1,+1) and (-1,-1)
    ker_sn(induvm) = ker_sn(induvm)-Cd5*h/4/j0;                 %(+1,-1) and (-1,+1)
    
    % include Jacobian for SLP and SLP-grad
    ker_s = j.*ker_s;
    ker_sn= j.*ker_sn;
    
    val_d = [val_d; ker_d*sig];     % DLP vals
    val_s = [val_s; ker_s*sig];     % SLP vals
    val_sn= [val_sn; ker_sn*sig];   % SLP normal gradient 
end

% plot err
err_d = abs(val_d(1:end-1)-val_d(end));
err_s = abs(val_s(1:end-1)-val_s(end));
err_sn= abs(val_sn(1:end-1)-val_sn(end));
hh = hh(1:end-1).';
subplot(1,2,2)
loglog(hh,hh.^5/hh(end)^5*err_d(end)*2,'k--',hh,err_d,'o-',hh,err_s,'v-',hh,err_sn,'*-');
title(sprintf('Helmholtz $\\kappa=%.2f%+.2fi$',real(ka),imag(ka)),'interpreter','latex')
xlabel('$h$','interpreter','latex')
ylabel('error','interpreter','latex','rotation',90)
legend({['$O(h^',num2str(5),')$'],'DLP','SLP','SLP grad'},'interpreter','latex','location','nw')

% display convergence order
fprintf(['convergence order of DLP:\n', repmat('  %.1f  ',1,numel(err_d)-1),'\n'],log2(err_d(1:end-1)./err_d(2:end)))
fprintf(['convergence order of SLP:\n', repmat('  %.1f  ',1,numel(err_s)-1),'\n'],log2(err_s(1:end-1)./err_s(2:end)))
fprintf(['convergence order of SLP normal graident:\n', repmat('  %.1f  ',1,numel(err_sn)-1),'\n'],log2(err_sn(1:end-1)./err_sn(2:end)))

function [Cd0,Cd01,Cd02,Cd1,Cd2,Cd3,Cd4,Cd5,...
    Cs0,Cs01,Cs02,Cs1,Cs2,Cs3,Cs4,Cs5,...
    Csn01,Csn1,Csn2] = wigner_type_lims(...
                                ru,rv,ruu,ruv,rvv,ruuu,ruuv,ruvv,rvvv,...
                                ruuuu,ruuuv,ruuvv,ruvvv,rvvvv)
% Compute Wigner-type limits for the O(h^5)  correction of Helmholtz DLP, 
% SLP, and normal gradient of SLP.
%
% Cdx are Wigner-type limits for the DLP
% * Cd0 is coefficient of the O(h) component of D0
% * Cd01 and Cd02 are coefficients of the O(h^3) components of D0
% * Cd1--Cd5 are coefficients of the O(h^3) factors D1--D5
%
% Csx are Wigner-type limits for the SLP
% * Cs0 is coefficient of the O(h) component of D0
% * Cs01 and Cs02 are coefficients of the O(h^3) components of D0
% * Cs1--Cs5 are coefficients of the O(h^3) factors D1--D5
%
% Csnx are Wigner-type limits for the normal gradient of SLP
% * Csnx == Cdx, except Csn01, Csn1, Csn2 are different

% compute coefficients of expansion of the geometry
E = dot(ru,ru); F = dot(ru,rv); G = dot(rv,rv);     % 1st fund form
[Lb2,Lb3,Lb4,La3,La4,Lc3,Lc4] = wigner_type_coeff(...
                                ru,rv,ruu,ruv,rvv,ruuu,ruuv,ruvv,rvvv,...
                                ruuuu,ruuuv,ruuvv,ruvvv,rvvvv);

% construct the \square^{(1)}, \square^{(2)} & \square^{(3)} operators

% first compute all second & third derivatives of Epstein zeta Z(s;E,F,G)
s = -1;
[~, d1_E,   d2_E,   d3_E,   d4_E] = epstein_zeta(s,E,F,G,1,0,0);
[~, d1_G,   d2_G,   d3_G,   d4_G] = epstein_zeta(s,E,F,G,0,0,1);
[~, d1_F,   d2_F,   d3_F,   d4_F] = epstein_zeta(s,E,F,G,0,0.5,0);
[~,~, d2_EpF, d3_EpF, d4_EpF] = epstein_zeta(s,E,F,G,1,.5,0);
[~,~,      ~,      ~,d4_Ep2F] = epstein_zeta(s,E,F,G,1,1,0);
[~,~,      ~, d3_EmF, d4_EmF] = epstein_zeta(s,E,F,G,1,-.5,0);
[~,~, d2_GpF, d3_GpF, d4_GpF] = epstein_zeta(s,E,F,G,0,.5,1);
[~,~,      ~,      ~,d4_Gp2F] = epstein_zeta(s,E,F,G,0,1,1);
[~,~,      ~, d3_GmF, d4_GmF] = epstein_zeta(s,E,F,G,0,-.5,1);
sq2 = [d2_E, (d2_EpF-d2_E-d2_F)/2, d2_F, (d2_GpF-d2_F-d2_G)/2, d2_G];
sq3 = [d3_E, (d3_EpF-d3_EmF-2*d3_F)/6, (d3_EpF+d3_EmF-2*d3_E)/6,...
    d3_F, (d3_GpF+d3_GmF-2*d3_G)/6, (d3_GpF-d3_GmF-2*d3_F)/6, d3_G];
sq4 = [d4_E, (6*d4_EpF-2*d4_EmF-d4_Ep2F-3*d4_E+12*d4_F)/24,...
    (d4_EpF+d4_EmF-2*d4_E-2*d4_F)/12,...
    (d4_Ep2F-3*d4_EpF-d4_EmF+3*d4_E-12*d4_F)/24,d4_F,...
    (d4_Gp2F-3*d4_GpF-d4_GmF+3*d4_G-12*d4_F)/24,...
    (d4_GpF+d4_GmF-2*d4_G-2*d4_F)/12,...
    (6*d4_GpF-2*d4_GmF-d4_Gp2F-3*d4_G+12*d4_F)/24,d4_G];

% D0
[Cs02,Cd02,~,~,~,Cs0,Cd0] = epstein_zeta(s,E,F,G,Lb2(1),Lb2(2)/2,Lb2(3));
Cs0 = -Cs0;         % O(h) err
Cd02 = -Cd02;      	% O(h^3) i*kappa error
Lb6 = conv(La4,Lb2)+conv(La3,Lb3);
Lc6 = conv(La4,Lb2)+conv(La3,Lc3);
La6 = conv(La3,La3);
Lb8 = conv(La6,Lb2);
Cd01 = (dot(Lb4,sq2)+2*dot(Lb6,sq3)+dot(Lb8,sq4)); % O(h^3) err, DLP
Cs01 = -(2*dot(La4,sq2)+dot(La6,sq3));             % O(h^3) err, SLP
Csn01= (dot(Lc4,sq2)+2*dot(Lc6,sq3)+dot(Lb8,sq4)); % SLP-grad O(h^3) err

% D1
Lb31 = conv([1,0],Lb3);
Lb61 = conv(conv([1,0],La3),Lb2);
La31 = conv([1,0],La3);
Lc31 = conv([1,0],Lc3);
Cd1 = 2*(dot(Lb31,sq2)+dot(Lb61,sq3)); %DLP
Cs1 = -2*dot(La31,sq2);                %SLP
Csn1= 2*(dot(Lc31,sq2)+dot(Lb61,sq3)); %SLP-grad

% D2
Lb32 = conv([0,1],Lb3);
Lb62 = conv(conv([0,1],La3),Lb2);
La32 = conv([0,1],La3);
Lc32 = conv([0,1],Lc3);
Cd2 = 2*(dot(Lb32,sq2)+dot(Lb62,sq3)); %DLP
Cs2  = -2*dot(La32,sq2);               %SLP
Csn2 =2*(dot(Lc32,sq2)+dot(Lb62,sq3)); %SLP-grad

% D3
Cd3 = 2*dot([Lb2,0,0],sq2); %DLP
Cs3 = -2*d1_E;              %SLP

% D4
Cd4 = 2*dot([0,0,Lb2],sq2); %DLP
Cs4 = -2*d1_G;              %SLP

% D5
Cd5 = 2*dot([0,Lb2,0],sq2); %DLP
Cs5 = -2*d1_F;              %SLP
end

function [Lb2,Lb3,Lb4,La3,La4,Lc3,Lc4] = wigner_type_coeff(...
    ru,rv,ruu,ruv,rvv,ruuu,ruuv,ruvv,rvvv,...
    ruuuu,ruuuv,ruuvv,ruvvv,rvvvv)
% Calculate coeffs in the expansion of geometry functions, 
% such as r.(ru x rv), in terms of u & v.
% Calculation done via vectorized convolutions.

% pre-define vectors coefficients of geometry expansion
nj = cross(ru,rv);
D1 = [ru,rv]; D2 = [ruu,2*ruv,rvv]/2; D3 = [ruuu,3*ruuv,3*ruvv,rvvv]/6;
D4 = [ruuuu,4*ruuuv,6*ruuvv,4*ruvvv,rvvvv]/24;
D1u = ru; D2u = [ruu,ruv]; D3u = [ruuu,2*ruuv,ruvv]/2; D4u = [ruuuu,3*ruuuv,3*ruuvv,ruvvv]/6;
D1v = rv; D2v = [ruv,rvv]; D3v = [ruuv,2*ruvv,rvvv]/2; D4v = [ruuuv,3*ruuvv,3*ruvvv,rvvvv]/6;

% O(u,v)^2 term coeffs ( = 2nd fundamental form * jacobian)
Lb2 = [dot(ruu,nj), 2*dot(ruv,nj), dot(rvv,nj)];

% O(u,v)^3 term coeffs
Lb3 =(detconv(D1,D2u,D2v)+detconv(D2,D1u,D2v)+detconv(D2,D2u,D1v)+...
 detconv(D1,D1u,D3v)+detconv(D1,D3u,D1v)+detconv(D3,D1u,D1v)) * -2;
La3 = dotconv([ru,rv],[ruu,2*ruv,rvv]);
Lc3 = dotconv(nj,D3*2);

% O(u,v)^4 term coeffs
Lb4 = detconv(D4,D1u,D1v)+detconv(D1,D4u,D1v)+detconv(D1,D1u,D4v)+...   %(1,1,4)
 detconv(D1,D2u,D3v)+detconv(D1,D3u,D2v)+...                         %(1,2,3)
 detconv(D2,D1u,D3v)+detconv(D2,D3u,D1v)+...
 detconv(D3,D1u,D2v)+detconv(D3,D2u,D1v); Lb4=Lb4 * -4;
La4 = dotconv([ru,rv],[ruuu,3*ruuv,3*ruvv,rvvv])/3+dotconv([ruu,2*ruv,rvv],[ruu,2*ruv,rvv])/4;
Lc4 = dotconv(nj,D4*4);

% assisting functions
    function c = dotconv(a,b)
        % convolution of two row vectors a,b whose elements are 3-vectors,
        % the usual multiplications in convolution is replace with dot
        % products
        c = conv(a(1,:),b(1,:))+conv(a(2,:),b(2,:))+conv(a(3,:),b(3,:));
    end
    function c = crossconv(a,b)
        % convolution of two row vectors a,b whose elements are 3-vectors,
        % the usual multiplications in convolution is replace with cross
        % products
        c = zeros(3,size(a,2)+size(b,2)-1);
        c(1,:) = conv(a(2,:),b(3,:))-conv(a(3,:),b(2,:));
        c(2,:) = conv(a(3,:),b(1,:))-conv(a(1,:),b(3,:));
        c(3,:) = conv(a(1,:),b(2,:))-conv(a(2,:),b(1,:));
    end
    function d = detconv(a,b,c)
        % convolution of three row vectors a,b,c whose elements are
        % 3-vectors, the usual multiplications in convolution is replace
        % with determinant of the 3 vector-elements being multiplied
        d = dotconv(a,crossconv(b,c));
    end
end

function PlotPatch(n,rvec,sig,ru,rv)
figure; subplot(1,2,1)
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