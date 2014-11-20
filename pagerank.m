function [urls, adjG, ptG, sp0, sp1, n] = pagerank(root, N, d)
%   root = root node for random surfer
%   N    = index size in pages
%   d    = damping factor
%
%   ex. [u,ag,pg,r0,r1,n] = pagerank('http://www.harvard.edu',50,0.85);
%
%   compute pagerank using 
%   power iteration method AND eig

[urls, adjG] = surfer(root, N);
adjG = full(adjG) * 1.0;

% Construct probability transition matrix
ptG = [];
for i=1:N,
    outlinks = adjG(:,i);
    if nnz(adjG(:,i)) > 0
        L = (nnz(outlinks));
        P_outlinks = (1/L) * outlinks;
        P_outlinks = arrayfun(@(p) tran_prob(p,d,N), P_outlinks);
    else
        P_outlinks = (1/N)*ones(N,1);
    end
    ptG = [ptG, P_outlinks];
end

% (1) MATLAB eig method
tol = 0.00001;
[eigvec,eigval] = eig(ptG);
[deig, col] = max(max(abs(eigval)));
if (eigval(1) <= 1 + tol) && (eigval(1) >= 1 - tol)
    sp0 = eigvec(:,1);
    sp0 = abs(sp0) / sum(abs(sp0));
else
    sp0 = eigvec(col,:);
    sp0 = transpose(sp0);
end


% (2) Power iteration method
sp1 = zeros(N,1);
sp1(1) = 1;
tol = 0.0000001;
n = [0];
while(1)
    iter = ptG * sp1;
    if (norm(iter) <= norm(sp1) + tol) && (norm(iter) >= norm(sp1) - tol)
        sp1 = iter;
        n = n + 1;
        break;
    end
    sp1 = iter;
    n(1) = n(1) + 1;
end
end

function pr = tran_prob(p, d, N)
%   p = initial probability (1/num_outlinks(page))
%   d = damping factor
%   N = index size in pages
%
%   compute transition probability from page i to page j
%   using random surfer model

pr = d*p + (1 - d)*(1/N);
end


function [t_rank, t_index] = report_ten_urls(r, u)
% report ranks / indices of top 10 ranking urls
% from pagerank
[t_rank, t_index] = sort(r, 'descend')
for i=1:10
    u(t_index(i))
    r(t_index(i))
end
end