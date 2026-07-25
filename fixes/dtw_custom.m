function [dist, wi, wj] = dtw_custom(A, B, metric)
% Simple DTW between two multivariate time series
% A: features x timeA
% B: features x timeB
% metric: 'absolute' (L1) or 'euclidean' (L2)

if nargin < 3
    metric = 'euclidean';
end

nA = size(A, 2);
nB = size(B, 2);

% Build local cost matrix
C = zeros(nA, nB);
for i = 1:nA
    for j = 1:nB
        diff = A(:,i) - B(:,j);
        if strcmp(metric, 'absolute')
            C(i,j) = sum(abs(diff));
        else
            C(i,j) = sqrt(sum(diff.^2));
        end
    end
end

% Accumulated cost matrix
D = inf(nA, nB);
D(1,1) = C(1,1);

% Fill first row and column
for i = 2:nA
    D(i,1) = D(i-1,1) + C(i,1);
end
for j = 2:nB
    D(1,j) = D(1,j-1) + C(1,j);
end

% Fill rest
for i = 2:nA
    for j = 2:nB
        D(i,j) = C(i,j) + min([D(i-1,j), D(i,j-1), D(i-1,j-1)]);
    end
end

% Traceback to get warp path
wi = zeros(1, nA+nB);
wj = zeros(1, nA+nB);
i = nA; j = nB;
k = 1;
wi(k) = i; wj(k) = j;

while i > 1 || j > 1
    if i == 1
        j = j - 1;
    elseif j == 1
        i = i - 1;
    else
        [~, idx] = min([D(i-1,j-1), D(i-1,j), D(i,j-1)]);
        if idx == 1
            i = i-1; j = j-1;
        elseif idx == 2
            i = i-1;
        else
            j = j-1;
        end
    end
    k = k+1;
    wi(k) = i; wj(k) = j;
end

wi = wi(1:k);
wj = wj(1:k);
dist = D(nA, nB);
end