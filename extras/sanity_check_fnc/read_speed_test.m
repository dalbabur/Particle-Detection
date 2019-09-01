
cine_file = '3.89uL-0.01vv-H362C-2.cine';
t = zeros(1,100);

for i = 1:100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileNamePath=[cine_folder '\' cine_file];
f=fopen(fileNamePath,'r');

% read size of annotation portion of image object:
fseek(f, int64(info.pImage(frameNum)),-1);
annotationSize = fread(f,1,'*uint32');
% file pointer is now positioned immediately after the
% annotation size field

% position file pointer at last 4 bytes of annotation block:
fseek(f,annotationSize-8,0);

% read actual size of image data:
ImageSize = double(fread(f,1,'*uint32')); %#ok<*NASGU>

% file pointer is now positioned at start of image data!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
a = fread(f,ImageSize/1.25,'ubit10=>uint16', 0, 'b');
t1 = toc;
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileNamePath=[cine_folder2 '\' cine_file];
f=fopen(fileNamePath,'r');

% read size of annotation portion of image object:
fseek(f, int64(info.pImage(frameNum)),-1);
annotationSize = fread(f,1,'*uint32');
% file pointer is now positioned immediately after the
% annotation size field

% position file pointer at last 4 bytes of annotation block:
fseek(f,annotationSize-8,0);

% read actual size of image data:
ImageSize = double(fread(f,1,'*uint32')); %#ok<*NASGU>

% file pointer is now positioned at start of image data!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear a
tic
a = fread(f,ImageSize/1.25,'ubit10=>uint16', 0, 'b');
t2 = toc;
fclose(f);

t(i) = t2/t1;
end