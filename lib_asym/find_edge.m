function [tEPlus, tEMinus, sE] = find_edge(im, bgim)

[height, width] = size(im);
% time domain plus edge
tmp1 = im-bgim;
tmp1 = tmp1(:); tmp1(tmp1<=0) = 0;
tmp1 = reshape(tmp1, height, width);
tmp1_1 = (double(tmp1)-double(min(tmp1(:))))/(double(max(tmp1(:))-min(tmp1(:))));
tEPlus = edge(tmp1_1);

% time domain minus edge
tmp2 = bgim-im;
tmp2 = tmp2(:); tmp2(tmp2<=0) = 0;
tmp2 = reshape(tmp2, height, width);
tmp2_1 = (double(tmp2)-double(min(tmp2(:))))/(double(max(tmp2(:))-min(tmp2(:))));
tEMinus = edge(tmp2_1);

% % spatial domain edge
% tmp3 = (double(im)-double(min(im(:))))/(double(max(im(:))-min(im(:))));
% sE = edge(tmp3);

% edge after background substraction
tmp4 = im-bgim;
tmp4 = (double(tmp4)-double(min(tmp4(:))))/(double(max(tmp4(:))-min(tmp4(:))));
sE = edge(tmp4);


% % verification
% display_image(tmp1, min(tmp1(:)), max(tmp1(:)),3);
% figure(4)
% imshow(tEPlus)
% 
% display_image(tmp2, min(tmp2(:)), max(tmp2(:)),5);
% figure(6)
% imshow(tEMinus)
% 
% figure(7)
% imshow(sE)