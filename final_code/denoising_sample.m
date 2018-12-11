%% Define image
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
filt = 'db8';
imfilter = wfilters(filt);
im_sz = size(img);
max_ht = maxPyrHt(im_sz, size(imfilter,2));
ht = max_ht;

[pi,pind] = wavedec2(img,ht,filt);
[piN,pindN] = wavedec2(delta,ht,filt);
pi = pi.';
piN = piN.';
pi_filtered = [];
n=0;
n1=0;
%A
pyramid{1} = reshape(pi( n + [ 1 : pind(1,1)*pind(1,2) ] ),pind(1,1),pind(1,2));
noise{1} = reshape(piN( n1 + [ 1 : pindN(1,1)*pindN(1,2) ] ),pindN(1,1),pindN(1,2));
n = n + pind(1,1)*pind(1,2);
n1 = n1 + pindN(1,1)*pindN(1,2);
filtered_pyramid = denoi_BLS_GSM_band(pyramid{1},[3,3],noise{1},0,1,1);
pi_filtered = [pi_filtered;filtered_pyramid(:)];

for i=2:size(pind,1)-1
    %H
    pyramid{i} = reshape(pi( n + [ 1 : pind(i,1)*pind(i,2) ] ),pind(i,1),pind(i,2));
    noise{i} = reshape(piN( n1 + [ 1 : pindN(i,1)*pindN(i,2) ] ),pindN(i,1),pindN(i,2));
    n = n + pind(i,1)*pind(i,2);
    n1 = n1 + pindN(i,1)*pindN(i,2);
    filtered_pyramid = denoi_BLS_GSM_band(pyramid{i},[3,3],noise{i},0,1,1);
    pi_filtered = [pi_filtered;filtered_pyramid(:)];
    %V
    pyramid{i} = reshape(pi( n + [ 1 : pind(i,1)*pind(i,2) ] ),pind(i,1),pind(i,2));
    noise{i} = reshape(piN( n1 + [ 1 : pindN(i,1)*pindN(i,2) ] ),pindN(i,1),pindN(i,2));
    n = n + pind(i,1)*pind(i,2);
    n1 = n1 + pindN(i,1)*pindN(i,2);
    filtered_pyramid = denoi_BLS_GSM_band(pyramid{i},[3,3],noise{i},0,1,1);
    pi_filtered = [pi_filtered;filtered_pyramid(:)];
    %D
    pyramid{i} = reshape(pi( n + [ 1 : pind(i,1)*pind(i,2) ] ),pind(i,1),pind(i,2));
    noise{i} = reshape(piN( n1 + [ 1 : pindN(i,1)*pindN(i,2) ] ),pindN(i,1),pindN(i,2));
    n = n + pind(i,1)*pind(i,2);
    n1 = n1 + pindN(i,1)*pindN(i,2);
    filtered_pyramid = denoi_BLS_GSM_band(pyramid{i},[3,3],noise{i},0,1,1);
    pi_filtered = [pi_filtered;filtered_pyramid(:)];
end
img_denoised = waverec2(pi_filtered,pind,filt);
%% Plot results
figure,imshowpair(img,img_denoised,'montage');
disp(snr(im2double(img),delta));
disp(snr(im2double(img_denoised),delta));