function [freq_csr,time] = DualAST_sparse_positive(z,beam,beam_ka,K,lambda,n1)


n_ka = length(z);
n = size(beam,1);
m = size(beam,2);
G_ka_t = beam_ka';
coef = zeros(1,2*m-1);
for j = -m+1:m-1
    coef(j+m) = sum(diag(beam'*beam,j));
end
coef = coef/trace(beam'*beam);
G = eye(m);
G = G(((m-1)/2-n1/2):((m-1)/2+n1/2),:);
cvx_solver sdpt3
cvx_tic
cvx_begin sdp quiet

variable H(n1+1,n1+1) hermitian; % H(w) is of degree m so that |H(w)|^2 is of degree n-1
variable q(n_ka,1) complex;              % Q(w) is a dual certificate of degree n-1
variable q_temp(2*m-1,1) complex;
dual variable Q;
q_temp == G_ka_t*q;
H>=0:Q;
trace(G'*H*G) + real(q_temp(m)) == coef(m);
for j = 1:(m-1)
     sum(diag(G'*H*G,j)) + q_temp(m-j) ==coef(m+j); % 1 - real(Q(w)) = |H(w)|^2
end
for j = 1:(m-1)
     sum(diag(G'*H*G,-j)) + q_temp(m+j) == coef(m-j); % 1 - real(Q(w)) = |H(w)|^2
end
maximize( real(z'*q) - lambda/2*(q'*q)); %LASSO penalty
cvx_end
time = cvx_toc;
time = time(5);
qReal = q_temp; %Coefficients of the real part of Q(w)

G_ka_r = beam_ka;
degrad = pi/180;
[U SS V] = svd(Q);
En = U(:,K+1:end);
% [U SS V]=svd(H);
% En = U(:,1:n1+1-K);
GG=G'*En*En'*G;
M2 = m-1;
MM=M2+1;
r = zeros(2*MM-1,1);
for j=-(MM-1):(MM+1)
    r(j+MM) = sum( diag(GG,j) );
end

ra = roots(r);
rb=ra(abs(ra)<1);
[~,I]=sort(abs(abs(rb)-1));
r_roots_detected=(rb(I(1:K)));
T_est = angle(r_roots_detected)/degrad;
freq_csr  = sort(T_est+179);
end