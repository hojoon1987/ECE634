close all
clear all
clc

search_step = 8;
%% Clean all the intermediate results
system('rm *crop*');
system('rm *mp4');

%% Load first fame and last frame of the video sequence
img_lst = dir('*.jpg');
img_name_lst = {img_lst.name};
N_img = length(img_name_lst);

img_name_tmp = img_name_lst{1};
Im1 = imread(img_name_tmp);
Im1_gray = double(rgb2gray(Im1));

img_name_tmp = img_name_lst{N_img};
Im2 = imread(img_name_tmp);
Im2_gray = double(rgb2gray(Im2));

%% Frame size
[M,N] = size(Im1_gray);

cnt = 1;
for width = search_step:search_step:N
    MSE(cnt) = mean(mean((Im1_gray(:,N-width+1:N) - Im2_gray(:,1:width)).^2));
    cnt = cnt + 1;
end

figure
plot(MSE)
title('The MSE of first and last frame assuming horizontal translation only', 'FontSize', 15)
xlabel('The shift pixel', 'FontSize', 15)
ylabel('The MSE value', 'FontSize', 15)

width_inx = find(MSE == min(MSE))

crop_im1 = Im1_gray(:,N-(width_inx-1)*search_step+1:N);
crop_im2 = Im2_gray(:,1:search_step*(width_inx-1));
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
for n = N_img:-1:1
    img_name = img_name_lst{n};
    img_tmp = double(rgb2gray(imread(img_name)));
    for ii = 1:1:(N - width_inx*search_step)
        crop_tmp = img_tmp(:,ii+1:ii+(width_inx - 1)*search_step);
        MSE_tmp(ii) = mean(mean((crop_tmp - crop_im2).^2));
    end
    shift_vec(cnt) = find(MSE_tmp == min(MSE_tmp));
    rgb_tmp = double(imread(img_name));
    img_save = rgb_tmp(:, shift_vec(cnt) + 1: shift_vec(cnt) + (width_inx - 1)*search_step,:);
    img_compare = img_tmp(:, shift_vec(cnt) + 1: shift_vec(cnt) + (width_inx - 1)*search_step);
    save_name = ['croppic00',num2str(n-1), '.jpg'];
    if n < 11
        save_name = ['croppic000', num2str(n-1), '.jpg'];
    else
    end
    imwrite(uint8(img_save), save_name, 'JPEG')
    diff_img = (abs(double(img_compare) - double(crop_im2)));
    figure
    imshow(uint8(diff_img));
    cnt = cnt + 1;
end



%% Reconstruct the video using 'ffmpeg'
% Manually input in terminal
% system('ffmpeg -r 10 -start_number 0 -i croppic_%4d.jpg -vcodec libx264 -r 30 -pix_fmt yuv420p stable_frame.mp4')



