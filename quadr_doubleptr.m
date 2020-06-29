function s = quadr_doubleptr(s, Ns)
% construct doubly periodic trapezoidal rule on a torus-like surface
% This function is from Alex Barnett's BIE3D package

if nargin == 0, test_quadr_doubleptr; return; end

Nu = Ns(1); Nv = Ns(2);
s.Nu=Nu; s.Nv=Nv; s.N = Nu*Nv;
s.u = (0:Nu-1)'/Nu*2*pi;
s.v = (0:Nv-1)'/Nv*2*pi;
[u,v] = ndgrid(s.u,s.v); % use 'ndgrid' instead of 'meshgrid' to avoid annoying dims swapping
u = u(:)'; v = v(:)'; % turn into row vectors

% 3-by-n arrays of points and partials and normals
[s.x, s.xu, s.xv, s.xuu, s.xuv, s.xvv,...
    s.xuuu,s.xuuv,s.xuvv,s.xvvv,...
    s.xuuuu,s.xuuuv,s.xuuvv,s.xuvvv,s.xvvvv] = s.get_x_xp_xpp_xppp(u,v);
s.nxj = cross(s.xu,s.xv); % unit normal times jacobian
s.sp = vecnorm(s.nxj); %length of normal (jacobian), or "speeds"
s.nx = s.nxj ./ s.sp; % unit normal
s.w  = (2*pi/Nu)*(2*pi/Nv)*s.sp; % quadr weights
% hu   = vecnorm(s.xu)*(2*pi/Nu);  % local node spacing in u-direc (h_1)
% hv   = vecnorm(s.xv)*(2*pi/Nv);  % local node spacing in v-direc (h_2)
% s.hmax = max([hu;hv]);
% s.hmin = min([hu;hv]);

function test_quadr_doubleptr

% Define & show torus-like surface
disp('torus-like surface')
m = 3; % petal number of generating curve
n = 5; % twist number along toroidal direction
a = 0.3;
s = wobblytorus(m,n,a);
Ns = [2,1]*50;
s = quadr_doubleptr(s, Ns);
showsurf(s)

disp('Gauss'' law flux convergence...')
ns = 10;
s_th = rand(1,ns)*2*pi; s_z = a/2 * (rand(1,ns)*2-1); % random src pts inside torus about unit circle (cylindrical coord)
s_str = rand(1,ns)*2 - 1; % random src strengths
err = [];
for nn = 10:5:80
    Nu = 2*nn; Nv = nn;
    s = quadr_doubleptr(s, [Nu,Nv]);
    flux = 0;
    for j = 1:ns
        z = [cos(s_th(j)); sin(s_th(j)); s_z(j)]; % src location
        hold on; plot3(z(1),z(2),z(3),'k.','markersize',20); 
        d = bsxfun(@minus,s.x,z); r = sqrt(sum(d.^2,1));
        ddotn = sum(d.*s.nx,1);
        flux = flux + s_str(j)*sum(s.w.*ddotn./r.^3)/(4*pi);      % surf flux of monopole at z
    end
    fprintf('N=[%d,%d], N=%d:  \terr = %.3g\n',Nu,Nv,s.N,flux - sum(s_str));
    err = [err; abs(flux - sum(s_str))];
end
hold off;
keyboard




