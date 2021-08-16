function [freq, time] = Sparse_MUSIC(Y, W,K)

% [freq, pow] = Sparse_MUSIC(Y, G)
%
% FDCA for sparse array
%
% Input
% Y: observation
% G: EADF matrix
%
% Output:
% freq: frequency estimate
% pow: power estimate
%
% Reference: 
%
% Written by Pan Jie, panjie@yzu.edu.cn
N = size(W,1);
degrad = pi/180;
Z = W*Y;
raw = fliplr(Z(1:((N-1)/2+1)));
col = Z(((N-1)/2+1):end);
R = toeplitz(col);
tic
[U S V] = svd(R);
En = U(:,K+1:end);
M = (N-1);
GG=En*En';
MM=M+1;
a = zeros(2*MM-1,1);
for j=-(MM-1):(MM+1)
    a(j+MM) = sum( diag(GG,j) );
end

ra=roots([a]);
rb=ra(abs(ra)<1);
[~,I]=sort(abs(abs(rb)-1));
w=angle(rb(I(1:K)));
freq=sort((w/degrad-1));
pow = 0;
time = toc;
end
