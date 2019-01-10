% Already have the .mat file loaded

%%% NoPol data
s = size(allNoPolImages);
imsz = s(1);
l1 = s(3);
l2 = s(4);

nims = l1*l2;
allNPims = zeros(imsz, imsz, nims);
meanNPfluxes = zeros(nims, 1); % Check none are blank
count = 1;
for k = 1:l1
    for l = 1:l2
        im = allNoPolImages(:, :, k, l);;
        allNPims(:,:,count) = im;
        meanNPfluxes(count) = mean(im(:));
        count = count+1;
    end
end

summedNPFlat = mean(allNPims, 3);
% subsum1 = mean(allims(:, :, 1:nims/2), 3);
% subsum2 = mean(allims(:, :, nims/2+1:end), 3);
% diff = subsum1-subsum2;

figure(1)
imagesc(summedNPFlat)
figure(2)
plot(meanNPfluxes)

disp('Press a key...')
pause


%%% Polariser data
s = size(allImages);
imsz = s(1);
l1 = s(3);
l2 = s(4);
l3 = s(5);

nims = l1*l2*l3;
allPims = zeros(imsz, imsz, nims);
meanPfluxes = zeros(nims, 1); % Check none are blank
count = 1;
for k = 1:l1
    for l = 1:l2
        for m = 1:l3
            im = allImages(:, :, k, l, m);
            allPims(:,:,count) = im;
            meanPfluxes(count) = mean(im(:));
            count = count+1;
        end
    end
end

summedPFlat = mean(allPims, 3);
% subsum1 = mean(allims(:, :, 1:nims/2), 3);
% subsum2 = mean(allims(:, :, nims/2+1:end), 3);
% diff = subsum1-subsum2;

figure(1)
imagesc(summedPFlat)
figure(2)
plot(meanPfluxes)


fitswrite(summedNPFlat, 'summedNPFlat.fits')
fitswrite(summedPFlat, 'summedPFlat.fits')

