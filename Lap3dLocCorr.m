function [As, Ad, Asn] = Lap3dLocCorr(s,ord)
% LAP3DLOCCORR. Laplace surface potentials Nystrom local-corrections 
% for the 2-periodic trapezoidal discretization of a toroidal surface.
%
% [As,Ad,Asn] = Lap3dLocCorr(s,ord)

if ord <= 1
    [As, Asn] = Lap3dSLPmat(s,s);
    Ad = Lap3dDLPmat(s,s);
    return
elseif ord <= 3
    [As,Ad,Asn] = CorrTrap1pt(s);
elseif ord <= 5
    % 9-point stencil location for each node
    [Nu,Nv] = ndgrid(1:s.Nu,1:s.Nv);
    sub1=sub2ind([s.Nu,s.Nv],Nu,Nv);                    % self
    sub2=sub2ind([s.Nu,s.Nv],Nu([2:end,1],:),Nv);       % south
    sub3=sub2ind([s.Nu,s.Nv],Nu([end,1:end-1],:),Nv);   % north
    sub4=sub2ind([s.Nu,s.Nv],Nu,Nv(:,[2:end,1]));       % east
    sub5=sub2ind([s.Nu,s.Nv],Nu,Nv(:,[end,1:end-1]));   % west
    sub6=sub2ind([s.Nu,s.Nv],Nu([2:end,1],:),Nv(:,[2:end,1]));          % south-east
    sub7=sub2ind([s.Nu,s.Nv],Nu([end,1:end-1],:),Nv(:,[end,1:end-1]));  % north-west
    sub8=sub2ind([s.Nu,s.Nv],Nu([2:end,1],:),Nv(:,[end,1:end-1]));      % south-west
    sub9=sub2ind([s.Nu,s.Nv],Nu([end,1:end-1],:),Nv(:,[2:end,1]));      % north-east
    
    sub = [sub1(:),sub2(:),sub3(:),sub4(:),sub5(:),sub6(:),sub7(:),sub8(:),sub9(:)];
    
    % Corrected quadr
    [As,Ad,Asn] = CorrTrap9pt(s,sub);
else
    error('ord: order of correction cannot bigger than 5')
end
end

function [As,Ad,Asn] = CorrTrap1pt(s)
% 2-periodic trapezoidal rule 1-point correction

% calculate coefficient (fundamental forms)
E = dot(s.xu,s.xu); F = dot(s.xu,s.xv); G = dot(s.xv,s.xv);
e = dot(s.xuu,s.nxj); f = dot(s.xuv,s.nxj); g = dot(s.xvv,s.nxj);
% O(h) errs
[Cs0,Cd0] = epstein_zeta(1,E,F,G,e,f,g); Cs0 = -Cs0;

% Laplace BIOs
[As, Asn] = Lap3dSLPmat(s,s);
Ad = Lap3dDLPmat(s,s);

% correction
h = 2*pi/s.Nu;
c = 1/4/pi;
N = size(As,1);
ind = sub2ind([N,N], 1:N, 1:N); ind = ind(:);
As(ind)  = Cs0.*s.sp*(h*c);
Ad(ind)  = Cd0*(h*c);
Asn(ind) = Cd0*(h*c);
end

function [As,Ad,Asn] = CorrTrap9pt(s,sub)
% 2-periodic trapezoidal rule 9-point local correction weights

% calculate coefficients (discrepancy factors)
[Cd0,Cd01,Cd1,Cd2,Cd3,Cd4,Cd5,...
    Cs0,Cs01,Cs1,Cs2,Cs3,Cs4,Cs5,...
    Csn01,Csn1,Csn2] = wigner_type_lims(s);

% Laplace BIOs
if_jac = 0; % don't include jacobian in SLP or its normal gradient
[As, Asn] = Lap3dSLPmat(s,s,if_jac);
Ad = Lap3dDLPmat(s,s);

% correction
N = size(s.x,2);
h = 2*pi/s.Nu;
hc = h/4/pi;
jac = s.sp(:);

% self  ( 0, 0)
ind = sub2ind([N,N],(1:N)',sub(:,1));
As(ind)  = ((Cs0-Cs3-Cs4)+Cs01*h^2)*hc;
Ad(ind)  = ((Cd0-Cd3-Cd4)+Cd01*h^2)*hc;
Asn(ind) = ((Cd0-Cd3-Cd4)+Csn01*h^2)./jac*hc;

% south (+1, 0)
ind = sub2ind([N,N],(1:N)',sub(:,2));
As(ind)  = As(ind)  + (Cs3+Cs1*h)*hc/2;
Ad(ind)  = Ad(ind)  + (Cd3+Cd1*h)*hc/2;
Asn(ind) = Asn(ind) + (Cd3+Csn1*h)./jac*hc/2;

% north (-1, 0)
ind = sub2ind([N,N],(1:N)',sub(:,3));
As(ind)  = As(ind)  + (Cs3-Cs1*h)*hc/2;
Ad(ind)  = Ad(ind)  + (Cd3-Cd1*h)*hc/2;
Asn(ind) = Asn(ind) + (Cd3-Csn1*h)./jac*hc/2;

% east  ( 0,+1)
ind = sub2ind([N,N],(1:N)',sub(:,4));
As(ind)  = As(ind)  + (Cs4+Cs2*h)*hc/2;
Ad(ind)  = Ad(ind)  + (Cd4+Cd2*h)*hc/2;
Asn(ind) = Asn(ind) + (Cd4+Csn2*h)./jac*hc/2;

% west  ( 0,-1)
ind = sub2ind([N,N],(1:N)',sub(:,5));
As(ind)  = As(ind)  + (Cs4-Cs2*h)*hc/2;
Ad(ind)  = Ad(ind)  + (Cd4-Cd2*h)*hc/2;
Asn(ind) = Asn(ind) + (Cd4-Csn2*h)./jac*hc/2;

% s-e   (+1,+1) and n-w (-1,-1)
Csn5 = Cd5./jac;
ind = sub2ind([N,N],(1:N)',sub(:,6));
As(ind)  = As(ind)  + Cs5*hc/4;
Ad(ind)  = Ad(ind)  + Cd5*hc/4;
Asn(ind) = Asn(ind) + Csn5*hc/4;
ind = sub2ind([N,N],(1:N)',sub(:,7));
As(ind)  = As(ind)  + Cs5*hc/4;
Ad(ind)  = Ad(ind)  + Cd5*hc/4;
Asn(ind) = Asn(ind) + Csn5*hc/4;

% s-w   (+1,-1) and n-e (-1,+1)
ind = sub2ind([N,N],(1:N)',sub(:,8));
As(ind)  = As(ind)  - Cs5*hc/4;
Ad(ind)  = Ad(ind)  - Cd5*hc/4;
Asn(ind) = Asn(ind) - Csn5*hc/4;
ind = sub2ind([N,N],(1:N)',sub(:,9));
As(ind)  = As(ind)  - Cs5*hc/4;
Ad(ind)  = Ad(ind)  - Cd5*hc/4;
Asn(ind) = Asn(ind) - Csn5*hc/4;

% include Jacobian
As = As .* s.sp; 
Asn= Asn.* s.sp; 
end

function [Cd0,Cd01,Cd1,Cd2,Cd3,Cd4,Cd5,...
    Cs0,Cs01,Cs1,Cs2,Cs3,Cs4,Cs5,...
    Csn01,Csn1,Csn2] = wigner_type_lims(s)
% Compute Wigner-type limits for the O(h^5)  correction of Laplace DLP, 
% SLP, and normal gradient of SLP.
%
% Cdx are Wigner-type limits for the DLP
% * Cd0 and Cd01 are coefficients of the O(h) and O(h^3) components of D0
% * Cd1--Cd5 are coefficients of the O(h^3) factors D1--D5
%
% Csx are Wigner-type limits for the SLP
% * Cs0 and Cs01 are coefficients of the O(h) and O(h^3) components of D0
% * Cs1--Cs5 are coefficients of the O(h^3) factors D1--D5
% 
% Csnx are Wigner-type limits for the normal gradient of SLP
% * Csnx == Cdx, except Csn01, Csn1, Csn2 are different

% compute coefficients of expansion of the geometry
E = dot(s.xu,s.xu); F = dot(s.xu,s.xv); G = dot(s.xv,s.xv); % 1st fund form
[Lb2,Lb3,Lb4,Lb5,Lb6,Lb8,La3,La4,La6,Lc3,Lc4,Lc6] = wigner_type_coeff(s);

% compute \square^{(1)}, \square^{(2)} & \square^{(3)}
% all second & third derivatives of Epstein zeta Z(s;E,F,G)
s = -1;
[~, d1_E, d2_E, d3_E, d4_E] = epstein_zeta(s,E,F,G,1,0,0);
[~, d1_G, d2_G, d3_G, d4_G] = epstein_zeta(s,E,F,G,0,0,1);
[~, d1_F, d2_F, d3_F, d4_F] = epstein_zeta(s,E,F,G,0,0.5,0);
[~,~, d2_EpF, d3_EpF, d4_EpF] = epstein_zeta(s,E,F,G,1,.5,0);
[~,~,      ~,      ~,d4_Ep2F] = epstein_zeta(s,E,F,G,1,1,0);
[~,~,      ~, d3_EmF, d4_EmF] = epstein_zeta(s,E,F,G,1,-.5,0);
[~,~, d2_GpF, d3_GpF, d4_GpF] = epstein_zeta(s,E,F,G,0,.5,1);
[~,~,      ~,      ~,d4_Gp2F] = epstein_zeta(s,E,F,G,0,1,1);
[~,~,      ~, d3_GmF, d4_GmF] = epstein_zeta(s,E,F,G,0,-.5,1);
sq2 = [d2_E; (d2_EpF-d2_E-d2_F)/2; d2_F; (d2_GpF-d2_F-d2_G)/2; d2_G];
sq3 = [d3_E; (d3_EpF-d3_EmF-2*d3_F)/6; (d3_EpF+d3_EmF-2*d3_E)/6;...
    d3_F; (d3_GpF+d3_GmF-2*d3_G)/6; (d3_GpF-d3_GmF-2*d3_F)/6; d3_G];
sq4 = [d4_E; (6*d4_EpF-2*d4_EmF-d4_Ep2F-3*d4_E+12*d4_F)/24;...
    (d4_EpF+d4_EmF-2*d4_E-2*d4_F)/12;...
    (d4_Ep2F-3*d4_EpF-d4_EmF+3*d4_E-12*d4_F)/24;d4_F;...
    (d4_Gp2F-3*d4_GpF-d4_GmF+3*d4_G-12*d4_F)/24;...
    (d4_GpF+d4_GmF-2*d4_G-2*d4_F)/12;...
    (6*d4_GpF-2*d4_GmF-d4_Gp2F-3*d4_G+12*d4_F)/24; d4_G];

% D0
[Cs0,Cd0] = epstein_zeta(1,E,F,G,Lb2(1,:),Lb2(2,:)/2,Lb2(3,:));
Cs0 = -Cs0(:); Cd0 = Cd0(:);                            % O(h) errs
Cd01 = (dot(Lb4,sq2)+2*dot(Lb6,sq3)+dot(Lb8,sq4)).';    % DLP O(h^3) err
Cs01 =-(2*dot(La4,sq2)+dot(La6,sq3)).';                 % SLP O(h^3) err
Csn01= (dot(Lc4,sq2)+2*dot(Lc6,sq3)+dot(Lb8,sq4)).';    % conj DLP O(h^3) err
% D1
Cs1 = -2*dot(La3,sq2(1:4,:)).';  	%SLP: dot((1,0)*La3,sq2)
Cd1 = 2*(dot(Lb3,sq2(1:4,:))...     %DLP: dot((1,0)*Lb3,sq2)
         +dot(Lb5,sq3(1:6,:))).';	%    +dot((1,0)*La3*Lb2,sq2)
Csn1= 2*(dot(Lc3,sq2(1:4,:))...     %conj DLP: dot((1,0)*Lc3,sq2)
         +dot(Lb5,sq3(1:6,:))).';	%         +dot((1,0)*La3*Lc2,sq2)
% D2
Cs2 = -2*dot(La3,sq2(2:5,:)).';  	%SLP: dot((0,1)*La3,sq2)
Cd2 = 2*(dot(Lb3,sq2(2:5,:))...     %DLP: dot((0,1)*Lb3,sq2)
         +dot(Lb5,sq3(2:7,:))).';	%    +dot((0,1)*La3*Lb2,sq2)
Csn2= 2*(dot(Lc3,sq2(2:5,:))...     %conj DLP: dot((0,1)*Lc3,sq2)
         +dot(Lb5,sq3(2:7,:))).';	%         +dot((0,1)*La3*Lc2,sq2)
% D3
Cd3 = 2*dot(Lb2,sq2(1:3,:)).';      %DLP: dot((1,0,0)*Lb2,sq2)
Cs3 =-2*d1_E(:);                    %SLP
% D4
Cd4 = 2*dot(Lb2,sq2(3:5,:)).';      %DLP: dot((0,0,1)*Lb2,sq2)
Cs4 =-2*d1_G(:);                    %SLP
% D5
Cd5 = 2*dot(Lb2,sq2(2:4,:)).';      %DLP: dot((0,1,0)*Lb2,sq2)
Cs5 =-2*d1_F(:);                    %SLP
end

function [Lb2,Lb3,Lb4,Lb5,Lb6,Lb8,La3,La4,La6,Lc3,Lc4,Lc6] = wigner_type_coeff(s)
% Expand r.(s.xu x s.xv) in terms of u & v

npt = size(s.x,2); % num of pts

% O(u,v)^2 term coeffs ( = 2nd fundamental form * jacobian)
Lb2 = [dot(s.nxj,s.xuu); 2*dot(s.nxj,s.xuv); dot(s.nxj,s.xvv)];

% O(u,v)^3 term coeffs
L30 = 2*dot(s.nxj,s.xuuu)/3-dot(s.xu,cross(s.xuu,s.xuv));
L21 = 2*dot(s.nxj,s.xuuv)-dot(s.xu,cross(s.xuu,s.xvv))-dot(s.xv,cross(s.xuu,s.xuv));
L12 = 2*dot(s.nxj,s.xuvv)-dot(s.xu,cross(s.xuv,s.xvv))-dot(s.xv,cross(s.xuu,s.xvv));
L03 = 2*dot(s.nxj,s.xvvv)/3-dot(s.xv,cross(s.xuv,s.xvv));
Lb3 = [L30;L21;L12;L03];

L30 = dot(s.xu,s.xuu);
L21 = dot(s.xu,s.xuv)*2+dot(s.xv,s.xuu);
L12 = dot(s.xv,s.xuv)*2+dot(s.xu,s.xvv);
L03 = dot(s.xv,s.xvv);
La3 = [L30;L21;L12;L03];

L30 = dot(s.nxj,s.xuuu)/3;
L21 = dot(s.nxj,s.xuuv);
L12 = dot(s.nxj,s.xuvv);
L03 = dot(s.nxj,s.xvvv)/3;
Lc3 = [L30;L21;L12;L03];

% O(u,v)^4 term coeffs
L40 = dot(s.nxj,s.xuuuu)/2-dot(s.xu,cross(s.xuu,s.xuuv))+...
    4*dot(s.xu,cross(s.xuv,s.xuuu))/3-dot(s.xv,cross(s.xuu,s.xuuu))/3;
L31 = 2*dot(s.nxj,s.xuuuv)-2*dot(s.xu,cross(s.xuu,s.xuvv))+2*dot(s.xu,cross(s.xuv,s.xuuv))+...
    4*dot(s.xu,cross(s.xvv,s.xuuu))/3-2*dot(s.xv,cross(s.xuu,s.xuuv))+2*dot(s.xv,cross(s.xuv,s.xuuu))/3;
L22 = 3*dot(s.nxj,s.xuuvv)-dot(s.xu,cross(s.xuu,s.xvvv))+...
    3*dot(s.xu,cross(s.xvv,s.xuuv))-3*dot(s.xv,cross(s.xuu,s.xuvv))+dot(s.xv,cross(s.xvv,s.xuuu));
L13 = 2*dot(s.nxj,s.xuvvv)+2*dot(s.xv,cross(s.xvv,s.xuuv))-2*dot(s.xv,cross(s.xuv,s.xuvv))-...
    4*dot(s.xv,cross(s.xuu,s.xvvv))/3+2*dot(s.xu,cross(s.xvv,s.xuvv))-2*dot(s.xu,cross(s.xuv,s.xvvv))/3;
L04 = dot(s.nxj,s.xvvvv)/2+dot(s.xv,cross(s.xvv,s.xuvv))-...
    4*dot(s.xv,cross(s.xuv,s.xvvv))/3+dot(s.xu,cross(s.xvv,s.xvvv))/3;
Lb4 = [L40;L31;L22;L13;L04];

L40 = dot(s.xuuu,s.xu)/3+dot(s.xuu,s.xuu)/4;
L31 = dot(s.xuuu,s.xv)/3+dot(s.xuu,s.xuv)+dot(s.xu,s.xuuv);
L22 = dot(s.xuuv,s.xv)+dot(s.xuu,s.xvv)/2+dot(s.xu,s.xuvv)+dot(s.xuv,s.xuv);
L13 = dot(s.xuvv,s.xv)+dot(s.xuv,s.xvv)+dot(s.xu,s.xvvv)/3;
L04 = dot(s.xvvv,s.xv)/3+dot(s.xvv,s.xvv)/4;
La4 = [L40;L31;L22;L13;L04];

L40 = dot(s.nxj,s.xuuuu)/6;
L31 = 2*dot(s.nxj,s.xuuuv)/3;
L22 = dot(s.nxj,s.xuuvv);
L13 = 2*dot(s.nxj,s.xuvvv)/3;
L04 = dot(s.nxj,s.xvvvv)/6;
Lc4 = [L40;L31;L22;L13;L04];

% O(u,v)^5 term coeffs
% Lb5 = conv(La3,Lb2)
Lb5 = zeros(6,npt);
Lb5(1:4,:) = Lb5(1:4,:) + La3.*Lb2(1,:);
Lb5(2:5,:) = Lb5(2:5,:) + La3.*Lb2(2,:);
Lb5(3:6,:) = Lb5(3:6,:) + La3.*Lb2(3,:);

% O(u,v)^6 term coeffs
% La6 = conv(La3,La3)
La6 = zeros(7,npt);
La6(1:4,:) = La6(1:4,:) + La3.*La3(1,:);
La6(2:5,:) = La6(2:5,:) + La3.*La3(2,:);
La6(3:6,:) = La6(3:6,:) + La3.*La3(3,:);
La6(4:7,:) = La6(4:7,:) + La3.*La3(4,:);
% Lb6 = conv(La4,Lb2)+conv(La3,Lb3)
Lb6 = zeros(7,npt);
Lb6(1:5,:) = Lb6(1:5,:) + La4.*Lb2(1,:);
Lb6(2:6,:) = Lb6(2:6,:) + La4.*Lb2(2,:);
Lb6(3:7,:) = Lb6(3:7,:) + La4.*Lb2(3,:);
Lc6 = Lb6;
Lb6(1:4,:) = Lb6(1:4,:) + La3.*Lb3(1,:);
Lb6(2:5,:) = Lb6(2:5,:) + La3.*Lb3(2,:);
Lb6(3:6,:) = Lb6(3:6,:) + La3.*Lb3(3,:);
Lb6(4:7,:) = Lb6(4:7,:) + La3.*Lb3(4,:);
% Lc6 = conv(La4,Lb2)+conv(La3,Lc3)
Lc6(1:4,:) = Lc6(1:4,:) + La3.*Lc3(1,:);
Lc6(2:5,:) = Lc6(2:5,:) + La3.*Lc3(2,:);
Lc6(3:6,:) = Lc6(3:6,:) + La3.*Lc3(3,:);
Lc6(4:7,:) = Lc6(4:7,:) + La3.*Lc3(4,:);

% O(u,v)^8 term coeffs
% Lb8 = conv(La6,Lb2);
Lb8 = zeros(9,npt);
Lb8(1:7,:) = Lb8(1:7,:) + La6.*Lb2(1,:);
Lb8(2:8,:) = Lb8(2:8,:) + La6.*Lb2(2,:);
Lb8(3:9,:) = Lb8(3:9,:) + La6.*Lb2(3,:);
end