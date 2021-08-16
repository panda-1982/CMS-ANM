function [freq, pow,time] = MST_GL_ANM(Y, beam_w, K)

% function [freq, pow] = MST_GL_ANM(Y, beam_w, K)
%
% MST_GL_ANM implements the Atomic Norm Minimization via generalized linear spectrum estimation for
% DOA estimation
%
% Input
% Y: observation
%beam_w; EADF matrix
%
% Output:
% freq: frequency estimate
% pow: power estimate
% sigma: noise variance estimate
%
% Written by Pan Jie, 2019   E-mail:panjie@yzu.edu.cn


[M, N] = size(Y);



    Rhat = Y * Y' / N;


cvx_quiet true
cvx_precision default

    [freq, pow,time] = ANM_GL(Rhat, beam_w, K,N);
   
end



function [freq, pow,time] = ANM_GL(Rhat,beam_w, K,N)
% SPA when Rhat is nonsingular
    M = size(beam_w, 2);
    B = size(beam_w,1);
    Rhat = B*Rhat/sum(real(diag(Rhat)));
    delta = svd(Rhat);
    delta = delta(end);
%     eta =sqrt(delta);
    eta = sqrt(N/7000)*B*delta(end);
    degrad = pi/180;
%     cvx_solver sdpt3
    cvx_tic
    cvx_begin sdp
    variable x(B,B) hermitian,
    variable u(M) complex,
    variable R(M,B) complex
    
    [x R'; R toeplitz(u)] >= 0,
    norm(beam_w*R-Rhat,'fro')<=eta
    minimize trace(x) + real(trace(toeplitz(u)));
    cvx_end
    time = cvx_toc;
    time = time(5);
    sig = 0;
    R_cov = beam_w*toeplitz(u)*beam_w';
    [U SS V] = svd(R_cov);
    En = U(:,K+1:end);
    GG=beam_w'*En*En'*beam_w;
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
end







