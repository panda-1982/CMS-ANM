clear all
clear;
index=1;
mode = 1;
rmse_csr_sparse_fd= 0;
rmse_csr_sparse_cfd = 0;
rmse_csr_sparse_full = 0;
rmse_music_25 = 0;
time_csr_sparse_fd = 0;
time_csr_sparse_cfd = 0;
time_csr_sparse_full = 0;
time_music_20 = 0;
time_music_25 = 0;
time_music = 0;
time_GLANM = 0;
time_GLSPICE = 0;
snr_index = 1;
mse_csr_sparse_fd = zeros(1,5);
mse_csr_sparse_cfd = zeros(1,5);
mse_csr_sparse_full = zeros(1,5);
mse_music_20 = zeros(1,5);
mse_music_25 = zeros(1,5);
rmse_GL = 0;
rmse_SPA = 0;
rmse_music = 0;
mse_SPA = zeros(1,5);
mse_GL= zeros(1,5);
mse_music = zeros(1,5);
alpha = 0;

x=[1.608 -1.793 -0.683 -1.136 -0.566 1.789 -0.816 2.187 -0.591];
y=[1.049 -1.256 -1.916 1.035 1.615 -1.308 1.915 0.355 -1.489];
array=[x;y];
m=size(x,2); % 阵元数
K=6; % 信号数
dt=[40  80 105 140 170  230]';
st=90*ones(K,1);

N=1000;
degrad=pi/180;
c =  3e8;
lamda=0.3;
fh = c/lamda;
k0=2*pi/lamda;
radius=2*lamda;
M1=20;
M2= 40;
T=361;
A1=zeros(m,T);
for t=-180:179
    A1(:,t+181)=sparse_steerMultiTar(array,t-1,fh,fh);
end
g=zeros(m,360);
g_ka=zeros(m^2-m+1,360);
for t=1:m
    g(t,:)=ifftshift(ifft(A1(t,:),360));
end
G=g(:,181-M1:181+M1);
A2=zeros(m,K);
k=(0:m-1)';
for t=1:K
    A2(:,t)=sparse_steerMultiTar(array,dt(t),fh,fh);
end
G_ka_origin = kron(conj(G),G);
G_ka_r = zeros(m^2,2*M2+1);
for index1 = -M1:M1
    for index2 = -M1:M1
        G_ka_r(:,index2-index1+M2+1)= G_ka_r(:,index2-index1+M2+1) + G_ka_origin(:,(index1+M1)*(2*M1+1)+index2+M1+1);
    end
end
MM1 = 2*M1+1;



load W_25;
W_25 = W;



iter_num = 100;
delta_index = 1;
delta =  3;
load sig07;
for snr = -10:1:10
%     snr = 10;
    parfor num = 1:iter_num
        warning off;
        Yt = A2*S(1:K,1:N);
        Y = awgn(Yt,snr,'measured');
        R = Y*Y'/N;
        R =R/sum(real(diag(R)));
        delta = svd(R);
        delta = delta(end);
        eta = 1*m*delta;
        X = R(:);
        [freq_csr_sparse,time_sparse30] = DualAST_sparse_positive(X,G,G_ka_r,K, eta,36);
        rmse_csr_sparse_fd=rmse_csr_sparse_fd+(sum(abs(freq_csr_sparse-dt).^2));
        time_csr_sparse_fd = time_csr_sparse_fd+time_sparse30;
        [freq_csr_sparse,time_sparse25] = DualAST_sparse_compressed_positive(X,G,G_ka_r,K, eta,23);
        rmse_csr_sparse_cfd=rmse_csr_sparse_cfd+(sum(abs(freq_csr_sparse-dt).^2));
        time_csr_sparse_cfd = time_csr_sparse_cfd+time_sparse25;
        [freq_csr_sparse,Q,time_pANM] = DualAST_MST_positive(X,G,G_ka_r,K, eta);
        rmse_csr_sparse_full=rmse_csr_sparse_full+(sum(abs(freq_csr_sparse-dt).^2));
        time_csr_sparse_full=time_csr_sparse_full+time_pANM;
        [freq, coef,time_rootmusic] = MS_MUSIC(Y,G,K);
        freq_music = 180+freq;
        rmse_music=rmse_music+(sum(abs(freq_music-dt).^2));
        time_music = time_music+time_rootmusic;
        [freq, coef,time_gl] = MST_GL_ANM(Y,G,K);
        freq_GL = 180+freq;
        rmse_GL=rmse_GL+(sum(abs(freq_GL-dt).^2));
        time_GLANM = time_GLANM+time_gl;
        [freq, coef,time_spa] = MS_SPA(Y, G, K);
        freq_SPA = 180+freq;
        rmse_SPA=rmse_SPA+(sum(abs(freq_SPA-dt).^2));
        time_GLSPICE = time_GLSPICE+time_spa;
        [freq, time_sMUSIC25] = Sparse_MUSIC(X,W_25,K);
        freq_music = sort(180-freq-2);
        rmse_music_25=rmse_music_25+(sum(abs(freq_music-dt).^2));
        time_music_25=time_music_25+time_sMUSIC25;
    end
    mse_music(delta_index) = sqrt(rmse_music/(iter_num*K));
    mse_SPA(delta_index) = sqrt(rmse_SPA/(iter_num*K));
    mse_GL(delta_index) = sqrt(rmse_GL/(iter_num*K));
    mse_music_25(delta_index) = sqrt(rmse_music_25/(iter_num*K));
    mse_csr_sparse_cfd(delta_index) = sqrt(rmse_csr_sparse_cfd/(iter_num*K));
    mse_csr_sparse_fd(delta_index) = sqrt(rmse_csr_sparse_fd/(iter_num*K));
    mse_csr_sparse_full(delta_index) = sqrt(rmse_csr_sparse_full/(iter_num*K));
    sqrt(rmse_GL/(iter_num*K))
    sqrt(rmse_csr_sparse_cfd/(iter_num*K))
    sqrt(rmse_csr_sparse_fd/(iter_num*K))
    sqrt(rmse_csr_sparse_full/(iter_num*K))

    rmse_GL = 0;
    rmse_SPA = 0;
    rmse_music = 0;
    rmse_music_25 = 0;
    rmse_csr_sparse_fd = 0;
    rmse_csr_sparse_cfd = 0;
    rmse_csr_sparse_full = 0;
    delta_index = delta_index + 1;
end
semilogy(-10:1:10,mse_music, '-x');
hold on;
semilogy(-10:1:10,mse_GL, '--h');
hold on;
semilogy(-10:1:10,mse_SPA, '-.^');
hold on;
semilogy(-10:1:10,mse_music_25, '-*');
hold on;
semilogy(-10:1:10,mse_csr_sparse_cfd,'-s');
hold on; 
semilogy(-10:1:10,mse_csr_sparse_fd,'-p');
hold on;
semilogy(-10:1:10,mse_csr_sparse_full,'-o');
hold on;
load crb_csr_snr;
semilogy(-10:1:10,crb,'-k');
legend('Root-MUSIC','GLE-ANM','GL-SPICE','FDCA','Compressed-FD-CMS-ANM','FD-CMS-ANM','CMS-ANM','CRB')
ylabel('RMSE (\circ)')
xlabel('SNR(dB)')















