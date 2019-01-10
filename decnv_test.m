imRaw = fitsread('omiCet_im.fits');
kernelRaw = fitsread('63Cet_im.fits');

kSize = 32; %Width of kernel
showSize = 64; %Image size to show

cnt = size(imRaw,1) / 2;
hs = kSize/2;
ss = showSize / 2;
kernel = kernelRaw(cnt-hs : cnt + hs - 1, cnt-hs : cnt + hs - 1);


figure(1)
subplot(1, 3, 1)
im = imRaw;
imagesc(im(cnt-ss:cnt+ss-1, cnt-ss:cnt+ss-1));
title('Raw image')
axis('equal')
subplot(1, 3, 2)
im = kernelRaw;
imagesc(im(cnt-ss:cnt+ss-1, cnt-ss:cnt+ss-1));
title('Kernel')
axis('equal')

subplot(1, 3, 3)

% imRaw = edgetaper(imRaw, kernel);
im = deconvlucy(imRaw, kernel, 10);
% im = deconvreg(imRaw, kernel, 10e-3);
% im = deconvwnr(imRaw, kernel, 1);

imagesc(im(cnt-ss:cnt+ss-1, cnt-ss:cnt+ss-1))
axis('equal')
title('Deconvolved')
colormap('gray')