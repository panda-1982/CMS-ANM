%多频率点阵列流形形成函数

function s=sparse_steerMultiTar(array,theta,f0,Fre)

% ------------------------------------------------------------------------------------
% create steer vector of narrow-band signal
% equal-spaced linear array
% ------------------------------------------------------------------------------------
%      Nsensor = the number of sensor element
%      d       = the normalized element distance by wavelength
%      theta   = the spatial angle 
%      f0      =
%      subband =frequency bins
% ------------------------------------------------------------------------------------

degrad=pi/180;
u=[cos(theta*degrad);sin(theta*degrad)];
s=exp(j*2*pi*Fre/f0*array'*u);