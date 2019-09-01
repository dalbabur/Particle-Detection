%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Droplet/Cell/Bead Detection
%
% Diego Alba 3/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idata] = cineRead2(cine_folder, cine_file,frameNums,info,LinLUT,...
    window_height,window_length,window_origin)
%
% Modiefied version of the cineRead function by Diego Alba, 3/12/2019
% Reads a small area of the frame specified by frameNum from the Phantom
% camera cine file specified by cine_file. The area is defined by
% window_height,window_lenght, and window_origin, all in pixels.
%
% NOTE: this function skips all of the checks implemented in the original;
%       may not be as robust.


% TESTED PARAMETERS
% cine_file ='3.89uL-0.01vv-H362C.cine';
% cine_folder='C:\Users\Diego\Documents\MATLAB\JHU\HUR\code\data';
% frameNum = 1251;
% info=cineInfo(cine_folder,cine_file);
% importfile('LinLUT.mat'); %a conversion from packed 10bit data to real 10bit data
% pixels = 1:16;
% window_origin = 280+350;
% window_length = 150;

resize = 1;
fileNamePath=[cine_folder '\' cine_file];
f=fopen(fileNamePath,'r');

% frames = frameNums(1):skip_frames:frameNums(2);
idata = zeros(length(window_height)*window_length,length(frameNums));

for p = 1:length(window_height) % read one row at a time
    for q = frameNums
        % ensure file pointer location:
        
        % read size of annotation portion of image object:
        fseek(f, int64(info.pImage(q)),-1);
        annotationSize = fread(f,1,'*uint32');
        % file pointer is now positioned immediately after the
        % annotation size field
        
        % position file pointer at last 4 bytes of annotation block:
        fseek(f,annotationSize-8,0);
        
        % read actual size of image data:
        ImageSize = double(fread(f,1,'*uint32')); %#ok<*NASGU>
        
        % file pointer is now positioned at start of image data!
        boi = ftell(f);
        
        fseek(f,boi+(window_height(p)-1)*info.Width*1.25 + window_origin,'bof');
        test = fread(f,window_length,'ubit10=>uint16', 0, 'b');
        if ~isempty(test)
            idata(1+(p-1)*window_length:p*window_length,q==frameNums) = test;
        else
            ftell(f)
            feof(f)
        end
    end
end

% postporcessing as in cineRead
idata=idata+1;
idata=LinLUT(idata);
idata=reshape(idata,window_length,length(window_height),length(frameNums));
idata=fliplr(rot90(idata,-1));
idata = imresize(idata,[length(window_height),window_length*resize],'bicubic');
fclose(f);
end




