function [urls, adjG, ptG, sp0, sp1, n, E] = ppagerank(root, N, d)
%   root = root node for random surfer
%   N    = index size in pages
%   d    = damping factor
%
%   ex. [u,ag,pg,r0,r1,n,E] = ppagerank('http://www.harvard.edu',50,0.85);
%
%   compute personalized pagerank by generating
%   a normal distribution of user preferences

[urls, adjG] = surfer(root, N);
adjG = full(adjG) * 1.0;

% Create user preference distribution
pd = makedist('Normal');
E = abs(random(pd, N, 1));
E = E / sum(E);

% Construct probability transition matrix
ptG = [];
for j=1:N,
    outlinks = adjG(:,j);
    if nnz(adjG(:,j)) > 0
        L = (nnz(outlinks));
        P_outlinks = (1/L) * outlinks;
        for i=1:N,
            P_outlinks(i) = d*P_outlinks(i) + (1 - d)*E(i);
            if i == j,
                P_outlinks(i) = 0;
            end
        end
    else
        for i=1:N,
            P_outlinks(i) = E(j);
        end
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
    sp0 = abs(sp0) / sum(abs(sp0));
end


% (2) Power iteration method
sp1 = zeros(N,1);
sp1(1) = 1;
tol = 0.00000001;
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


function [t_rank, t_index] = report_ten_urls(r, u, E)
% report ranks / indices of top 10 ranking urls
% from pagerank
[t_rank, t_index] = sort(r, 'descend')
for i=1:10
    u(t_index(i))
    r(t_index(i))
    E(t_index(i))
end
end