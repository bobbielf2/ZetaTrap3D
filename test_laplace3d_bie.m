% Solve Laplace BVPs exterior to a toroidal surface 
% using 9-point zeta-correction of the 2-periodic trapezoidal rule
% See the manuscript [1], Example 2, for more details.
%
% [1] Wu, B., & Martinsson, P.G. (2020, arXiv:2007.02512). Corrected
%     Trapezoidal Rules for Boundary Integral Equations in Three
%     Dimensions.

% define toroidal surface
m = 1; % petal number of generating curve
n = 5; % twist number along toroidal direction
a = 0.2;
s = wobblytorus(m,n,a);

% define exact solution
zs = [0.3,-0.9,0; 0.5,1.0,-0; -1,-0.1,0].'; % source inside torus
zt = [2.4;-1.1;0.7]*0.8;                    % target test pt outside torus
strength = randn(3,1)+1;                    % source strength
figure(1); subplot(1,2,1)
plot3(zs(1,:),zs(2,:),zs(3,:),'k.','markersize',20); hold on;
plot3(zt(1),zt(2),zt(3),'r.','markersize',20); showsurf(s); hold off;
legend({'source','target'})
uexac = Lap3dSLPmat(struct('x',zt),struct('x',zs,'w',1))*strength;

% convergence test (WARNING: form full matrix, slow for large N!)
ord = 5;    % pick convergence order, ord = no bigger than 5
bvp = 'nd'; % n = Neumann, d = Dirichlet, nd = both
NN = round(2.^(3:0.5:6));
errN = []; errD = [];
for Nv = NN     % num of nodes in the v-direction
    Nu = Nv;    % assume Nu = Nv (could be different)
    s = quadr_doubleptr(s, [Nu,Nv]);
    
    % Diri & Neu data on grid
    [f, g] = Lap3dSLPmat(s,struct('x',zs,'w',1));
    f = f*strength; g = g*strength;
    
    % precompute locally corrected potential matrix
    [As, Ad, Asn] = Lap3dLocCorr(s,ord);
    
    % BIE solve
    if strcmp(bvp,'n') || strcmp(bvp,'nd')      % solve ext Neu BVP. Ansatz: u = S*tau
        AN = Asn-0.5*eye(size(Asn));            % define Neumann BIO: A=-0.5+D'
        tauN = AN\g;
        uN = Lap3dSLPmat(struct('x',zt),s)*tauN;
        err = abs(uN - uexac)/abs(uexac);
        fprintf('Ns=[%d,%d], N=%d:  \tNeu err = %.3g\n',s.Nu,s.Nv,s.N,err);
        errN = [errN;err];
    end
    if strcmp(bvp,'d') || strcmp(bvp,'nd')      % solve ext Diri BVP. Ansatz: u = D*tau
        AD = Ad+0.5*eye(size(Ad));              % define Dirichlet BIO: A=0.5+D (w\ 1D nullspace)
        AD = AD+(1./vecnorm(s.x-zs(:,1)))'.*s.w;% cancel nullspace with an interior source
        tauD = AD\f;
        uD = Lap3dDLPmat(struct('x',zt),s)*tauD + s.w*tauD/norm(zt-zs(:,1));
        err = abs(uD - uexac)/abs(uexac);
        fprintf('Ns=[%d,%d], N=%d:  \tDiri err = %.3g\n',s.Nu,s.Nv,s.N,err);
        errD = [errD;err];
    end
end

subplot(1,2,2)
Ntot = NN.^2; % total num of points (assume Nu = Nv)
errref = NN.^-ord/NN(end)^-ord*errD(end)*3; % reference error
loglog(Ntot,errD,'o-',Ntot,errN,'o-',Ntot,errref,'k--')
legend({'Dirichlet','Neumann',['O(h^{',num2str(ord),'})']})
title('Laplace'); xlabel('N (total number of points)')
