close all
clear all
clc

Im1 = imread('origpic0000.jpg');
Im1_gray = double(rgb2gray(Im1));

Im2 = imread('origpic0051.jpg');
Im2_gray = double(rgb2gray(Im2));

%% Frame size
[M,N] = size(Im1_gray);


cnt = 1;
for width = 8:8:N
    MSE(cnt) = mean(mean((Im1_gray(:,N-width+1:N) - Im2_gray(:,1:width)).^2));
    cnt = cnt + 1
end

figure
plot(MSE)

width_inx = find(MSE == min(MSE));

crop_im1 = Im1_gray(:,N-(width_inx-1)*8+1:N);
crop_im2 = Im2_gray(:,1:8*(width_inx-1));
figure
% The first frame
subplot(2,2,1)
imshow(uint8(Im1_gray));
% The last frame
subplot(2,2,2)
imshow(uint8(Im2_gray));
subplot(2,2,3)
imshow(uint8(crop_im1))
subplot(2,2,4)
imshow(uint8(crop_im2))

cnt = 1;
for n = 50:-1:10
    img_name = ['origpic00',num2str(n),'.jpg']
    img_tmp = double(rgb2gray(imread(img_name)));
    for ii = 1:1:(N - width_inx*8)
        crop_tmp = img_tmp(:,ii+1:ii+(width_inx - 1)*8);
        MSE_tmp(ii) = mean(mean((crop_tmp - crop_im2).^2));
    end
    shift_vec(cnt) = find(MSE_tmp == min(MSE_tmp));
    rgb_tmp = double(imread(img_name));
    img_save = rgb_tmp(:, shift_vec(cnt) + 1: shift_vec(cnt) + (width_inx - 1)*8,:);
    save_name = ['croppic00',num2str(n), '.jpg'];
    imwrite(uint8(img_save), save_name, 'JPEG')
    cnt = cnt + 1;
    
end





