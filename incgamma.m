function [gs,gsp1,gsp2,gsp3,gsp4] = incgamma(s,x)
% My customed fast implementation of the (scaled) upper incomplet gamma 
% function for all real s, including s<=0. (Still slow for complex s.)
%
% Input:
%   (s,x) are the two arguments of the incomplete gamma function defined as
%       incgamma(s,x) = \int_1^\infty t^{s-1}*exp(-x*t) dt
%                     = x^{-s} \int_x^\infty t^{s-1}*exp(-t) dt
%   x>0 and at least one of s and x must be a scalar.
% Output:
%   gs   == igamma(s,x)*x^(-s)
%   gsp1 == igamma(s+1,x)*x^(-(s+1))
%   gsp2 == igamma(s+2,x)*x^(-(s+2))
%   gsp3 == igamma(s+3,x)*x^(-(s+3))
%   gsp4 == igamma(s+3,x)*x^(-(s+4))
%
% Note: this implementation is Fast in the following sense:
% When isreal(s), we use gammainc(x,s,'upper') or expint(x) (FAST) 
% instead of igamma(s,x) (SLOW); in particular:
%   - When s<0, use recursion, since gammainc(x,s) only take s>0
%   - When s=0, use expint(x) == igamma(0,x)
% When s is complex, we fall back to igamma(s,x) (Symbolic Toolbox
% required), so still slow.

% Bowei Wu 5/2020. Updated 11/2020.

if nargin == 0, test_incgamma; return; end % unit test

if isreal(s)
    if numel(s) == 1 && s <= 0
        k = -floor(s);
        ss = s + k;
        gg = zeros([size(x),k+1]);
        if ss == 0  % s is integer
            % expint(x) is equivalent to (but MUCH faster than) igamma(0,x) from the symbolic toolbox
            gg(:,:,1) = expint(x).*x.^(-ss);
        else
            % gammainc(x,s,'upper').*gamma(s) is equivalent to (but MUCH faster than) igamma(s,x) from the symbolic toolbox
            gg(:,:,1) = gammainc(x,ss,'upper').*gamma(ss).*x.^(-ss);
        end
        ex = exp(-x);
        for i = 1:k
            ss = ss - 1;
            gg(:,:,i+1) = (gg(:,:,i).*x-ex)./ss;
        end
        gs = gg(:,:,end);
        if nargout > 1
            if k > 3
                gsp1 = gg(:,:,k);
                gsp2 = gg(:,:,k-1);
                gsp3 = gg(:,:,k-2);
                gsp4 = gg(:,:,k-3);
            elseif k == 3
                gsp1 = gg(:,:,k);
                gsp2 = gg(:,:,k-1);
                gsp3 = gg(:,:,k-2);
                gsp4 = ((s+3).*gsp3+ex)./x;
            elseif k == 2
                gsp1 = gg(:,:,k);
                gsp2 = gg(:,:,k-1);
                gsp3 = ((s+2).*gsp2+ex)./x;
                gsp4 = ((s+3).*gsp3+ex)./x;
            elseif k == 1
                gsp1 = gg(:,:,k);
                gsp2 = ((s+1).*gsp1+ex)./x;
                gsp3 = ((s+2).*gsp2+ex)./x;
                gsp4 = ((s+3).*gsp3+ex)./x;
            elseif k == 0
                gsp1 = ex./x;
                gsp2 = ((s+1).*gsp1+ex)./x;
                gsp3 = ((s+2).*gsp2+ex)./x;
                gsp4 = ((s+3).*gsp3+ex)./x;
            end
        end
    else    % handle s>0 case or numel(s)>1 and numel(x)==1 case
        % TODO: this currently only works if all(s>=0)
        gs = gammainc(x,s,'upper').*gamma(s).*x.^(-s);
        if numel(s)>1, gs(s==0) = expint(x); end % correct where s=0
        if nargout > 1
            ex = exp(-x);
            gsp1 = (s.*gs+ex)./x;
            gsp2 = ((s+1).*gsp1+ex)./x;
            gsp3 = ((s+2).*gsp2+ex)./x;
            gsp4 = ((s+3).*gsp3+ex)./x;
        end
    end
else
    gs = igamma(s,x).*x.^(-s); % slow
    if nargout > 1
        ex = exp(-x);
        gsp1 = (s.*gs+ex)./x;
        gsp2 = ((s+1).*gsp1+ex)./x;
        gsp3 = ((s+2).*gsp2+ex)./x;
        gsp4 = ((s+3).*gsp3+ex)./x;
    end
end
end

function test_incgamma
% unittest
x = 2.^(-5:5);
for s = -5:0.5:5
    g = cell(1,5);
    [g{1},g{2},g{3},g{4},g{5}] = incgamma(s,x);
    
    % verify results
    err_abs = zeros(1,5);
    err_rel = zeros(1,5);
    g_exact = cell(1,5);
    for i = 1:5
        ss = s+i-1;
        g_exact{i} = igamma(ss,x).*x.^(-ss);
        err_abs(i) = norm(g_exact{i}-g{i});
        err_rel(i) = norm((g_exact{i}-g{i})./g_exact{i});
    end
    fprintf('s = %.1f, x = [',s)
    fprintf('%.2f ',x)
    fprintf([']:\n\t',...
              'g(s,x):\t\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+1,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+2,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+3,x):\tabs.err = %.1e \trel.err = %.1e\n\t',...
              'g(s+4,x):\tabs.err = %.1e \trel.err = %.1e\n'],...
              [err_abs;err_rel]);
end
end