close all
clear all
clc

search_step = 2;
%% Clean all the intermediate results
system('rm *crop*');
system('rm *mp4');
system('rm sift*');
system('rm step*');

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

width_inx = find(MSE == min(MSE));

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
    old_img_compare = crop_im2;
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
    %     diff_img = (abs(double(img_compare) - double(crop_im2)));
    %     figure
    %     imshow(uint8(diff_img));
    cnt = cnt + 1;
end

crop_img = dir('crop*');
crop_img_lst = {crop_img.name};

%figure
for jj = 1:1:N_img-1
    % Select the area to crop
    crop_name1 = crop_img_lst{jj};
    crop_name2 = crop_img_lst{jj + 1};
    crop_tmp1 = imread(crop_name1);
    crop_tmp2 = imread(crop_name2);
    % Prepare the image for the SIFT operator
    Ia = rgb2gray(crop_tmp1);
    Ib = rgb2gray(crop_tmp2);
    [fa, da] = vl_sift(single(Ia));
    [fb, db] = vl_sift(single(Ib));
    % [matches_raw, scores] = vl_ubcmatch(da, db);
    [matches_select, scores] = vl_ubcmatch(da, db);
    new_img = [Ia,Ib];
    
    % cnt = 1;
    % for ii = 1:1:size(matches_raw,2)
    %     if scores(ii) > mean(scores)
    %         matches_select(:,cnt) = matches_raw(:,ii);
    %         cnt = cnt + 1;
    %     else
    %     end
    % end
    
    %figure
    
    %% Save the SIFT results
    %figure
    %imshow(uint8(new_img))
    %hold on
    clear dx_tmp
    clear dy_tmp
    for ii  = 1:1:size(matches_select,2)
        indx1 = matches_select(1,ii);
        indx2 = matches_select(2,ii);
        x1 = fa(1,indx1);
        y1 = fa(2,indx1);
        x2 = fb(1,indx2);
        y2 = fb(2,indx2);
        dx_tmp(ii) = x2 - x1;
        dy_tmp(ii) = y2 - y1;
        
        %    line([x1,x2 + size(Ia,2)],[y1,y2]);
    end
    %figure
    %scatter(dx_tmp,dy_tmp)
    dx(jj) = median(dx_tmp);
    dy(jj) = median(dy_tmp);
    %hold on
    %scatter(dx(jj),dy(jj), 'r')
    %     save_name = ['siftpic00',num2str(jj-1), '.jpg'];
    %     if jj < 11
    %         save_name = ['siftpic000', num2str(jj-1), '.jpg'];
    %     else
    %     end
    %     saveas(gcf, save_name, 'jpg');
    %     close(gcf)
end
%close all
%dx
%dy
%% Reconstruct the video using 'ffmpeg'
% Manually input in terminal
% system('ffmpeg -r 10 -start_number 0 -i croppic%4d.jpg -vcodec libx264 -r 30 -pix_fmt yuv420p stable_frame.mp4')

%% Calculate the difference across each cropped frame: Temporal Step = 1
step_temporal = 1;
cnt = 1;
for jj = 1:1:N_img - step_temporal
    % Select the area to crop
    crop_name1 = crop_img_lst{jj};
    crop_name2 = crop_img_lst{jj + step_temporal};
    crop_tmp1 = rgb2gray(imread(crop_name1));
    crop_tmp2 = rgb2gray(imread(crop_name2));
    diff_img_1(:,:,cnt) = abs(crop_tmp2 - crop_tmp1);
    %figure
    %imshow(uint8(diff_img_1(:,:,cnt)));
    cnt = cnt + 1;
    save_name = ['step1pic00',num2str(jj-1), '.jpg'];
    if jj < 11
        save_name = ['step1pic000', num2str(jj-1), '.jpg'];
    else
    end
    imwrite(uint8(abs(crop_tmp2 - crop_tmp1)), save_name, 'JPEG')
end
%close all

%% Calculate the difference across each cropped frame: Temporal Step = 1
step_temporal = 3;
cnt = 1;
for jj = 1:1:N_img - step_temporal
    % Select the area to crop
    crop_name1 = crop_img_lst{jj};
    crop_name2 = crop_img_lst{jj + step_temporal};
    crop_tmp1 = rgb2gray(imread(crop_name1));
    crop_tmp2 = rgb2gray(imread(crop_name2));
    diff_img_3(:,:,cnt) = abs(crop_tmp2 - crop_tmp1);
    %figure
    %imshow(uint8(diff_img_3(:,:,cnt)));
    cnt = cnt + 1;
    save_name = ['step3pic00',num2str(jj-1), '.jpg'];
    if jj < 11
        save_name = ['step3pic000', num2str(jj-1), '.jpg'];
    else
    end
    imwrite(uint8(abs(crop_tmp2 - crop_tmp1)), save_name, 'JPEG')
end
%close all

%% Calculate the difference across each cropped frame: Temporal Step = 1
step_temporal = 5;
cnt = 1;
for jj = 1:1:N_img - step_temporal
    % Select the area to crop
    crop_name1 = crop_img_lst{jj};
    crop_name2 = crop_img_lst{jj + step_temporal};
    crop_tmp1 = rgb2gray(imread(crop_name1));
    crop_tmp2 = rgb2gray(imread(crop_name2));
    diff_img_5(:,:,cnt) = abs(crop_tmp2 - crop_tmp1);
    %figure
    %imshow(uint8(diff_img_5(:,:,cnt)));
    cnt = cnt + 1;
    save_name = ['step5pic00',num2str(jj-1), '.jpg'];
    if jj < 11
        save_name = ['step5pic000', num2str(jj-1), '.jpg'];
    else
    end
    imwrite(uint8(abs(crop_tmp2 - crop_tmp1)), save_name, 'JPEG')
end
%close all

%% Find the minimum length of diff_img_n matrices
min_size = size(diff_img_5, 3);
diff_combined = diff_img_1(:,:,1:min_size) + diff_img_3(:,:,1:min_size) ...
    + diff_img_5(:,:,1:min_size);




