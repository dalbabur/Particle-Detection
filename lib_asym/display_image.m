function display_image(in, minv, maxv, n)
% linear mapping uint16 to double(0~1) for displaying
in = double(in);
in = in-double(minv);
in = in/double(maxv-minv);
figure(n)
imshow(in);