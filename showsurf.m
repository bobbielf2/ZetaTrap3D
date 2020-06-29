function showsurf(s)
% show torus-like surface
% Modified from the same function in Alex Barnett's BIE3D package

if nargin == 0, test_showsurf; return; end
if ~isfield(s,'x'), Ns = [1,1]*60; s = quadr_doubleptr(s, Ns); end

Nr = [s.Nu,s.Nv];   % reshaping sizes
                    % note: assume grid generated by 'ndgrid'. If used 'meshgrid' then need to swap dim(x) and dim(y)!
x = reshape_wrapAround(s.x(1,:),Nr);
y = reshape_wrapAround(s.x(2,:),Nr);
z = reshape_wrapAround(s.x(3,:),Nr);

mesh(x,y,z,'FaceAlpha',0); axis equal; hold on


function x = reshape_wrapAround(x0,Nr)
% reshape and wrap around to complete the torus-like surface for plotting
x = reshape(x0,Nr);
x = [x; x(1,:)]; x = [x, x(:,1)];

function test_showsurf

s = wobblytorus;
Ns = [2,1]*50;
s = quadr_doubleptr(s, Ns);
showsurf(s)