img = imread('cameraman.tif');
%% define noise
PS = ones(size(img)); %default is white noise
PS = fftshift(PS);
isoddPS_y = (size(PS,1)~=2*(floor(size(PS,1)/2)));
isoddPS_x = (size(PS,2)~=2*(floor(size(PS,2)/2)));
PS = PS(1:end-isoddPS_y, 1:end-isoddPS_x);          % ensures even dimensions for the power spectrum
PS = fftshift(PS);
delta = real(ifft2(sqrt(PS)));      
delta = fftshift(delta);
PS = abs(fft2(delta)).^2;
PS = fftshift(PS);    
% noise, to be used only with translation variant transforms (such as orthogonal wavelet)
delta = real(ifft2(sqrt(PS).*exp(j*angle(fft2(randn(size(PS)))))));  
%% Pyramid Construction, Denoising and Pyramid Reconstruction
[pi,pind] = buildWpyr(img,'auto','db8');
[piN,pindN] = buildWpyr(delta,'auto','db8');
pi_filtered = [];
pind = flip(pind);
pindN = flip(pindN);
n=0;
for i=1:size(pind,1)
    pyramid{i} = reshape(pi( n + [ 1 : pind(i,1)*pind(i,2) ] ),pind(i,1),pind(i,2));
    noise{i} = reshape(piN( n + [ 1 : pindN(i,1)*pindN(i,2) ] ),pindN(i,1),pindN(i,2));
    filtered_pyramid = denoi_BLS_GSM_band(pyramid{i},[3,3],noise{i},0,1,1);
    pi_filtered = [pi_filtered;filtered_pyramid(:)];
end
pindN = flip(pindN);
img_denoised = reconWpyr(pi_filtered,pindN);