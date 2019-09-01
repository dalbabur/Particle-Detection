function BGI=BGI_prep_auto(CineFileLocation,CineFileName,NumFrame, LinLUT,Resize)

%Background Image Generation
BGICineFileName=CineFileName; %Background Cine file name
info_BGI=cineInfo(CineFileLocation,BGICineFileName);
%BGINumberImages=info_BGI.NumFrames;

Sum_BGI_frame=cineRead(CineFileLocation,BGICineFileName,1,info_BGI,LinLUT)/(NumFrame);
for i=1:NumFrame
    Sum_BGI_frame=Sum_BGI_frame+cineRead(CineFileLocation,BGICineFileName,i,info_BGI,LinLUT)/(NumFrame);
end
AVE_BGI_frame=imresize(Sum_BGI_frame,Resize,'bicubic');
BGI=imadjust(AVE_BGI_frame);
BGICineFileName_New=[substr(BGICineFileName,0,-5),'_BGI.tif']; 
%Save the Averaged Background Image
imwrite(AVE_BGI_frame,BGICineFileName_New);

end
