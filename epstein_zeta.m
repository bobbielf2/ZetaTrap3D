function [S,Sd1,Sd2,Sd3,Sd4,Ssp2,Sd1sp2] = epstein_zeta(s,E,F,G,L,M,N)
% Evaluate the Epstein zeta function Z(s) defined by
%       Z(s) = \sum_{(i,j)\neq(0,0)} (E*i^2+2*F*i*j+G*j^2)^(-s/2)
% for Re(s)>2, and analytically continued to the whole complex plane except
% a simple pole at s=2. See the manuscript [1], Appendix E, for more details.
%
% Output:
%   S      = Z(s)
%   Sd1    = (L*{d/dE} + M*{d/dF} + N*{d/dG}) Z(s)
%   Sd2    = (L*{d/dE} + M*{d/dF} + N*{d/dG})^2 Z(s)
%   Sd3    = (L*{d/dE} + M*{d/dF} + N*{d/dG})^3 Z(s)
%   Sd4    = (L*{d/dE} + M*{d/dF} + N*{d/dG})^4 Z(s)
%   Ssp2   = Z(s+2)
%   Sd1sp2 = (L*{d/dE} + M*{d/dF} + N*{d/dG}) Z(s+2)
%
% [1] Wu, B., & Martinsson, P.G. (2020, arXiv:2007.02512). Corrected
%     Trapezoidal Rules for Boundary Integral Equations in Three
%     Dimensions.
%
% Bowei Wu, 2020/3/30
%           4/7 add 1st derivative
%           4/16 add higher derivatives

if nargin == 0, test_epstein_zeta; return; end % no input for unit testing

if numel(s) > 1
    assert(numel(E)==1,'either s or {E,F,G} must be scalar')
    assert(numel(F)==1,'either s or {E,F,G} must be scalar')
    assert(numel(G)==1,'either s or {E,F,G} must be scalar')
else
    assert(isequal(size(E),size(F)),'E and F must have the same size')
    assert(isequal(size(E),size(G)),'E and G must have the same size')
    if s==0 && nargout==1, S=-ones(size(E)); return; end % handle exception: Z(0)=-1.
end

% Preprocessing: make quadratic form has unit determinant
J = sqrt(E.*G-F.^2);
E = E./J; G = G./J; F = F./J;
if nargout > 1, L = L./J; M = M./J; N = N./J; end

% constants for computing derivatives
if nargout > 1
H   = (G.*L+E.*N-2*F.*M)/2;
K   = L.*N-M.^2;
Hsq = K - 2*H.^2;
Hsq2= -2*H.*K-4*H.*Hsq;
Hsq3= (4*H.^2-2*Hsq).*K-4*Hsq.^2-4*H.*Hsq2;
end

% Define the quadratic form
QA  = @(i,j) E*i^2+2*F*i*j+G*j^2;
if nargout > 1, QB = @(i,j) L*i^2+2*M*i*j+N*j^2; end

% Determine summation cutoffs based on decay of incomplete gamma func
lambda = min(min((E+G)/2 - sqrt((E-G).^2+4*F.^2)/2)); % min eigenvalue of Q
n = floor(sqrt(33/pi./lambda))+3;

% summation, exclude origin
S=0; 
if nargout > 1, Sd1=0; Sd2=0; Sd3=0; Sd4=0; end
if nargout > 5, Ssp2=0; Sd1sp2 = 0; end
s1 = s/2; s2 = 1-s1;
for i = 0:n-1     % right-half ij-plane & j-axis
    for j = 1:ceil(sqrt(n^2-i^2))
        % first quadrant + pos j-axis
        if nargout > 1
            x = pi*QA(i,j); y = pi*QB(i,j); ex = exp(-x);
            xsq = y-H.*x; xsq2 = -H.*y-Hsq.*x-H.*xsq;
            xsq3 = (-Hsq+H.^2).*y-(Hsq2.*x+2*Hsq.*xsq+H.*xsq2);
            xsq4 = (-Hsq2+3*H.*Hsq-H.^3).*y-(Hsq3.*x+3*Hsq2.*xsq+3*Hsq.*xsq2+H.*xsq3);
            [g1,g1p1,g1p2,g1p3,g1p4] = incgamma(s1,x);
            [g2,g2p1,g2p2,g2p3,g2p4] = incgamma(s2,x);
            gp1=g1p1+g2p1; gp2=g1p2+g2p2; gp3=g1p3+g2p3; gp4=g1p4+g2p4;
            S   = S+g1+g2;
            Sd1 = Sd1-gp1.*xsq;
            Sd2 = Sd2+gp2.*xsq.^2-gp1.*xsq2;
            Sd3 = Sd3-gp3.*xsq.^3+3*gp2.*xsq.*xsq2-gp1.*xsq3;
            Sd4 = Sd4+gp4.*xsq.^4-6*gp3.*xsq.^2.*xsq2+...
                gp2.*(4*xsq.*xsq3+3*xsq2.^2)-gp1.*xsq4;
        else
            x = pi*QA(i,j);
            S = S+incgamma(s1,x)+incgamma(s2,x);
        end
        if nargout > 5
            Ssp2 = Ssp2+g1p1-(g2.*x-ex)./s1;
            Sd1sp2 = Sd1sp2-(g1p2+g2).*xsq;
        end
        
        % fourth quadrant + pos i-axis
        if nargout > 1
            x = pi*QA(j,-i); y = pi*QB(j,-i); ex = exp(-x);
            xsq = y-H.*x; xsq2 = -H.*y-Hsq.*x-H.*xsq;
            xsq3 = (-Hsq+H.^2).*y-Hsq2.*x-2*Hsq.*xsq-H.*xsq2;
            xsq4 = (-Hsq2+3*H.*Hsq-H.^3).*y-(Hsq3.*x+3*Hsq2.*xsq+3*Hsq.*xsq2+H.*xsq3);
            [g1,g1p1,g1p2,g1p3,g1p4] = incgamma(s1,x);
            [g2,g2p1,g2p2,g2p3,g2p4] = incgamma(s2,x);
            gp1=g1p1+g2p1; gp2=g1p2+g2p2; gp3=g1p3+g2p3; gp4=g1p4+g2p4;
            S   = S+g1+g2;
            Sd1 = Sd1-gp1.*xsq;
            Sd2 = Sd2+gp2.*xsq.^2-gp1.*xsq2;
            Sd3 = Sd3-gp3.*xsq.^3+3*gp2.*xsq.*xsq2-gp1.*xsq3;
            Sd4 = Sd4+gp4.*xsq.^4-6*gp3.*xsq.^2.*xsq2+...
                gp2.*(4*xsq.*xsq3+3*xsq2.^2)-gp1.*xsq4;
        else
            x = pi*QA(j,-i);
            S = S+incgamma(s1,x)+incgamma(s2,x);
        end
        if nargout > 5
            Ssp2 = Ssp2+g1p1-(g2.*x-ex)./s1;
            Sd1sp2 = Sd1sp2-(g1p2+g2).*xsq;
        end
    end
end
% Postprocessing: symmetry add-back & determinant scale-back
if isreal(s1)
    Cs1 = (pi./J).^s1 ./ gamma(s1); % scaling constant
else
    Cs1 = (pi./J).^s1 ./ igamma(s1,0);
end
S  = (2*S - 1 ./s1 - 1 ./s2).*Cs1;
if nargout > 1
    Sd1=2*Sd1.*Cs1 - s1.*H.*S;
    Sd2=2*Sd2.*Cs1 - (s1.*Hsq+(s1.*H).^2).*S - 2*s1.*H.*Sd1;
    Sd3=2*Sd3.*Cs1 - (s1.*Hsq2+3*s1.^2.*H.*Hsq+(s1.*H).^3).*S - 3*(s1.*Hsq+(s1.*H).^2).*Sd1 - 3*s1.*H.*Sd2;
    Sd4=2*Sd4.*Cs1 - (s1.*Hsq3+3*(s1.*Hsq).^2+4*s1.^2.*H.*Hsq2+6*s1.^3.*H.^2.*Hsq+(s1.*H).^4).*S-...
        (4*s1.*Hsq2+12*s1.^2.*H.*Hsq+4*(s1.*H).^3).*Sd1-6*(s1.*Hsq+(s1.*H).^2).*Sd2-4*s1.*H.*Sd3;
end
% handle exception: Z(0)=-1, dZ(0)=0, Z(2)=Inf, dZ(2)=Inf.
if numel(s)==1 && s==0
    S(:)=-1;
    if nargout > 1, Sd1(:)=0; Sd2(:)=0; Sd3(:)=0; Sd4(:)=0; end
elseif numel(s)>1
    ind = s==0;
    S(ind)=-1;
    if nargout > 1, Sd1(ind)=0; Sd2(ind)=0; Sd3(ind)=0; Sd4(ind)=0; end
end

if nargout > 5
    Cs1p1 = Cs1.*(pi./J./s1);
    Ssp2 = (2*Ssp2 - 1 ./(s1+1) - 1 ./(-s1)).*Cs1p1;
    Sd1sp2 = 2*Sd1sp2.*Cs1p1 - (s1+1).*H.*Ssp2;
    % handle exception: Z(0)=-1, dZ(0)=0, Z(2)=Inf, dZ(2)=Inf.
    if numel(s)==1
        if s==-2
            Ssp2(:)=-1; Sd1sp2(:)=0;
        elseif s==0
            Ssp2(:)=Inf; Sd1sp2(:)=Inf;
        end
    elseif numel(s)>1
        Ssp2(s==-2)=-1; Sd1sp2(s==-2)=0;
        Ssp2(s==0)=Inf; Sd1sp2(s==0)=Inf;
    end
end
end

function test_epstein_zeta
% unit testing

% gradient checking
ru = randn(3,1); rv = randn(3,1);
E = dot(ru,ru); F = dot(ru,rv); G = dot(rv,rv);
s = -1.5;
fprintf('s=%.1f, (E,F,G) = (%f,%f,%f), angle=pi*%f\n',s,E,F,G,acos(F/sqrt(E*G))/pi)

% 1. compute your func output
[~,     ~, Sd2_EG, Sd3_EpG, S4d_EpG] = epstein_zeta(s,E,F,G,1,0,1);
[~,     ~,      ~, Sd3_EmG, S4d_EmG] = epstein_zeta(s,E,F,G,1,0,-1);
[~, Sd1_E,  Sd2_E,   Sd3_E,   Sd4_E] = epstein_zeta(s,E,F,G,1,0,0);
[S, Sd1_F,  Sd2_F,       ~,   Sd4_F] = epstein_zeta(s,E,F,G,0,1,0);
[~, Sd1_G,  Sd2_G,   Sd3_G,   Sd4_G] = epstein_zeta(s,E,F,G,0,0,1);
[~,     ~,      ~, Sd3_EpF] = epstein_zeta(s,E,F,G,1,1,0);
[~,     ~,      ~, Sd3_EmF] = epstein_zeta(s,E,F,G,1,-1,0);
Sd2_EG = (Sd2_EG - Sd2_E - Sd2_G)/2;
Sd3_EEG = (Sd3_EpG - Sd3_EmG -2*Sd3_G)/6;
Sd3_EFF = (Sd3_EpF + Sd3_EmF -2*Sd3_E)/6;
Sd4_EEGG = (S4d_EpG+S4d_EmG-2*Sd4_E-2*Sd4_G)/12;

% 2. finite diff approx
epsi = 1e-4;
% (d/dF) and (d^2/dF^2)
S_Fp = epstein_zeta(s,E,F+epsi,G,1,0,0);
S_Fm = epstein_zeta(s,E,F-epsi,G,1,0,0);
Sd1_F_approx = (S_Fp - S_Fm)/(2*epsi);
Sd2_FF_approx = (S_Fp - 2*S + S_Fm)/epsi^2;
% (d^2/dEdG), (d^3/dE^2dG) and (d^3/dEdF^2)
[S_EpGp,Sd1_E_EpGp] = epstein_zeta(s,E+epsi,F,G+epsi,1,0,0);
[S_EpGm,Sd1_E_EpGm] = epstein_zeta(s,E+epsi,F,G-epsi,1,0,0);
[S_EmGp,Sd1_E_EmGp] = epstein_zeta(s,E-epsi,F,G+epsi,1,0,0);
[S_EmGm,Sd1_E_EmGm] = epstein_zeta(s,E-epsi,F,G-epsi,1,0,0);
[~,~,Sd2_E_Gp] = epstein_zeta(s,E,F,G+epsi,1,0,0);
[~,~,Sd2_E_Gm] = epstein_zeta(s,E,F,G-epsi,1,0,0);
[~,Sd1_G_Ep] = epstein_zeta(s,E+epsi,F,G,0,0,1);
[~,Sd1_G_Em] = epstein_zeta(s,E-epsi,F,G,0,0,1);
[~,Sd1_E_Fp] = epstein_zeta(s,E,F+epsi,G,1,0,0);
[~,Sd1_E_Fm] = epstein_zeta(s,E,F-epsi,G,1,0,0);
Sd2_EG_approx = ((S_EpGp - S_EpGm)-(S_EmGp - S_EmGm))/(4*epsi^2);
Sd3_EEG_approx = ((Sd1_E_EpGp - Sd1_E_EpGm)-(Sd1_E_EmGp - Sd1_E_EmGm))/(4*epsi^2);
Sd3_EEG_approx2 = (Sd2_E_Gp - Sd2_E_Gm)/(2*epsi);
Sd3_EEG_approx3 = (Sd1_G_Ep - 2*Sd1_G + Sd1_G_Em)/epsi^2;
Sd3_EFF_approx = (Sd1_E_Fp - 2*Sd1_E + Sd1_E_Fm)/epsi^2;


% 3. verify func output against FDM approx
fprintf('\nverify your derivatives...\n')
fprintf('(d/dF)Z(s):\t\t[FDM approx, your output]=[%f, %f], the same?\n',Sd1_F_approx,Sd1_F)
fprintf('(d^2/dF^2)Z(s):\t\t[FDM approx, your output]=[%f, %f], the same?\n',Sd2_FF_approx,Sd2_F)
fprintf('(d^3/dE^2dG)Z(s):\t3 different FDM approxs: [%f, %f, %f], your output: %f\n',Sd3_EEG_approx,Sd3_EEG_approx2,Sd3_EEG_approx3,Sd3_EEG)

% 4. verify derivative equivalence: (d^2/dF^2)Z(s)/4 = (d^2/dEdG)Z(s) and (d^3/dE^2dG)Z(s) = (d^3/dEdF^2)Z(s)/4
% this test also implicitly verify correctness of the computed Z(s)
fprintf('\nderivatives equivalence check... (implicitly verify correctness of Z(s))\n')
fprintf('FDM approx:\t[(d^2/dF^2)Z(s)/4, (d^2/dEdG)Z(s)] = [%f, %f], the same?\n',Sd2_FF_approx/4,Sd2_EG_approx)
fprintf('your output:\t[(d^2/dF^2)Z(s)/4, (d^2/dEdG)Z(s)] = [%f, %f] the same? (should be the same as prev line, too)\n',Sd2_F/4,Sd2_EG)
fprintf('FDM approx:\t[(d^3/dE^2dG)Z(s), (d^3/dEdF^2)Z(s)/4] = [%f, %f] the same?\n',Sd3_EEG_approx,Sd3_EFF_approx/4)
fprintf('your output:\t[(d^3/dE^2dG)Z(s), (d^3/dEdF^2)Z(s)/4] = [%f, %f] the same? (should be the same as prev line, too)\n',Sd3_EEG,Sd3_EFF/4)
fprintf('your output:\t[(d^4/dE^2dG^2)Z(s), (d^4/dF^4)Z(s)/16] = [%f, %f] the same?\n',Sd4_EEGG,Sd4_F/16)
end
