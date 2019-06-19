function out = normalize_img(in)
% normalize input image and save it with uint8
tmp = in-min(in(:));
tmp = tmp/max(tmp(:));
out = uint8(floor(tmp*255));