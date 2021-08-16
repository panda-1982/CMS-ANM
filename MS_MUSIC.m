function [freq, pow, time] = MS_MUSIC(Y, G,K)

% [freq, pow] = MS_MUSIC(Y, G)
%
% Root-MUSIC with Manifold separation
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
N = size(Y,2);
degrad = pi/180;
Rhat = Y * Y' / N;
tic
[U S V] = svd(Rhat);
En = U(:,K+1:end);
M = size(G,2);
GG=G'*En*En'*G;
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
