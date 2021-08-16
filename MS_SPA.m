function [freq, pow,time] = MS_SPA(Y, beam,K, Omega)

% [freq, pow] = SPA(Y, mode_noise, Omega)
% 
% SPA implements the sparse and parametric approach (SPA) for linear array 
% signal processing.
% 
% Input
% Y: observation
% mode_noise: 1 if equal noise variances;
%             2 if different noise variances,
%             default value: 2
% Omega: index set of SLA, valid only in the SLA case, sorted ascendingly
% 
% Output:
% freq: frequency estimate
% pow: power estimate
% sigma: noise variance estimate
% 
% Reference: Z. Yang, L. Xie, and C. Zhang, "A discretization-free sparse 
%   and parametric approach for linear array signal processing", 
%   IEEE Trans. Signal Processing, 2014
% 
% Written by Zai Yang, April 2013

if nargin < 4 || max(Omega) == length(Omega)  % ULA
    Omega = [];
else
    if length(Omega) ~= size(Y,1)
        error('Dimension of Omega not match!');
    end
end
if nargin < 2
    mode_noise = 2;
end

[M, N] = size(Y);

if N > 1
    Rhat = Y * Y' / N;
end
degrad = pi/180;
mode_noise = 1;
cvx_quiet true
cvx_precision default

% solve the SDP (or estimate R)
if N >= M && cond(Rhat) < 1e8       % nonsingular
    [u, sig,time] = SPA_nonsing(Rhat, beam,mode_noise, Omega);
    
elseif N == 1  % single snapshot
    [u, sig] = SPA_SMV(Y, mode_noise, Omega);
    
else           % multisnapshots and singular
    [u, sig] = SPA_singu(Rhat, mode_noise, Omega);
    
end
    R_cov = toeplitz(u);
    [U SS V] = svd(R_cov);
    En = U(:,K+1:end);
    GG=En*En';
    M = size(beam,2);
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
% [freq, pow] = PostProc(u,K);
end
% postprocessing and parameter estimation





function [u, sig,time] = SPA_nonsing(Rhat, beam,mode_noise, Omega)
% SPA when Rhat is nonsingular

if isempty(Omega)         % ULA
    M = size(Rhat, 1);
    N = size(beam, 2);
    Rhatinv = inv(Rhat);
%     TconjR = zeros(M,1);
%     TconjR(1) = trace(Rhatinv);
%     for j = 2:M
%         TconjR(j) = 2*sum(diag(Rhatinv,j-1));
%     end
    
    Rhathalf = sqrtm(Rhat);
    
    switch mode_noise
            
        case 1,  % equal noise variances
            cvx_tic;
            cvx_solver sdpt3
            cvx_begin sdp            
              variable x(M,M) hermitian,
              variable u(N) complex,
              
              [x Rhathalf; Rhathalf beam*toeplitz(u)*beam'] >= 0,
              toeplitz(u) >= 0,
              minimize trace(x) + trace(real(Rhatinv*beam*toeplitz(u)*beam'));
            cvx_end
            sig = 0;
            time = cvx_toc;
            time = time(5);
        case 2,  % different noise variances
            cvx_solver sdpt3
            cvx_begin sdp
              variable x(M,M) hermitian,
              variable u(M) complex,
              variable sig(M), 
              
              sig >= 0,
              toeplitz(u) >= 0,
              [x Rhathalf; Rhathalf toeplitz(u)+diag(sig)] >= 0,
              
              minimize trace(x) + real(TconjR'*u) + real(diag(Rhatinv)')*sig;
            cvx_end
            
        otherwise,
            error('error!');
    end
    
    return
end


%% SLA

M = max(Omega);
Mbar = length(Omega);

Rhatinv = inv(Rhat);
Rhatinv_expand = zeros(M);
Rhatinv_expand(Omega, Omega) = Rhatinv;
TconjR = zeros(M,1);
TconjR(1) = trace(Rhatinv_expand);
for j = 2:M
    TconjR(j) = 2*sum(diag(Rhatinv_expand,j-1));
end

S = zeros(Mbar, M);
S(:, Omega) = eye(Mbar);

Rhathalf = sqrtm(Rhat);

switch mode_noise
        
    case 1,  % equal noise variances
        cvx_solver sdpt3
        cvx_begin sdp
          variable x(Mbar,Mbar) hermitian,
          variable u(M) complex,
          
          toeplitz(u) >= 0,
          [x Rhathalf; Rhathalf S*toeplitz(u)*S'] >= 0,
            
          minimize trace(x) + real(TconjR'*u);
        cvx_end
        
        sig = 0;
        
    case 2,  % different noise variances
        cvx_solver sdpt3
        cvx_begin sdp
          variable x(Mbar,Mbar) hermitian,
          variable u(M) complex,
          variable sig(Mbar), 
          
          sig >= 0,
          toeplitz(u) >= 0,
          [x Rhathalf; Rhathalf S*toeplitz(u)*S'+diag(sig)] >= 0,
            
          minimize trace(x) + real(TconjR'*u) + real(diag(Rhatinv)')*sig;
        cvx_end
        
    otherwise,
        error('error!');
end

end
















function [freq, amp] = PostProc(u,K)
degrad = pi/180;
num_K = K;
% [freq, amp] = PostProc(u)
% PostProc implements the postprocessing procedure given the solution u.
% The result is sorted in the descending order of amp

prec = 1e-5; % can be tuned appropriately

eigval = sort(real(eig(toeplitz(u))));
if eigval(1) > 0
    sigma_in_u = eigval(1);
    u(1) = u(1) - sigma_in_u;
    eigval = eigval - sigma_in_u;
else
    sigma_in_u = 0;
end

K = sum(eigval > prec*eigval(end));

M = length(u);

% solve h
h = toeplitz([conj(u(2)); u(1:M-1)], conj(u(2:K+1))) \ (-u(1:M));

% solve the zeros of H(z) in the form of exp(-1i * 2 * pi * theta)
r = roots([1; h]);
% r = r ./ abs(r);

% solve the amper vector
mat = flipud(vander([r; zeros(M-K,1)]).');
amp = mat(:, 1:K) \ u;
[amp, idx] = sort(real(amp), 'descend');
% Khat = sum(amp>0);
Khat = num_K;
amp = amp(1:Khat);

% determine the frequency
freq = sort( phase(r(idx(1:Khat))) / degrad -1 );
% freq = mod(freq, 1);
% freq = sort(asin(- phase(r(idx(1:Khat)))/pi)/degrad);
end