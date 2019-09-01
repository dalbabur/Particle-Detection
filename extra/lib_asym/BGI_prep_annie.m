function BGI=BGI_prep_annie(CineFileLocation,CineFileName,FirstFrame,LastFrame, LinLUT,Resize)

%Background Image Generation
BGICineFileName=CineFileName; %Background Cine file name
info_BGI=cineInfo(CineFileLocation,BGICineFileName);
%BGINumberImages=info_BGI.NumFrames;

Sum_BGI_frame=0;%this is the correct one
%Sum_BGI_frame=cineRead(CineFileLocation,BGICineFileName,FirstFrame,info_BGI,LinLUT);
%Sum_BGI_frame=imadjust(Sum_BGI_frame);
%figure,imshow(imadjust(Sum_BGI_frame));
NumFrame=(LastFrame-FirstFrame+1);
for i=1:NumFrame
    %Sum_BGI_frame=Sum_BGI_frame+imadjust(cineRead(CineFileLocation,BGICineFileName,i,info_BGI,LinLUT));
    Sum_BGI_frame=Sum_BGI_frame+cineRead(CineFileLocation,BGICineFileName,FirstFrame+i-1,info_BGI,LinLUT);
end
AVE_BGI_frame=Sum_BGI_frame/(NumFrame);
%figure,imshow(AVE_BGI_frame);
%break
AVE_BGI_frame=imresize(AVE_BGI_frame,Resize,'bicubic');
BGI=AVE_BGI_frame; %imadjust(AVE_BGI_frame);
BGICineFileName_New=[substr(BGICineFileName,0,-5),'_BGI.tif']; 
%Save the Averaged Background Image
imwrite(AVE_BGI_frame,BGICineFileName_New);

end
