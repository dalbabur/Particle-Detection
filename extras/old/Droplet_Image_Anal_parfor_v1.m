%%%%%
%Droplet images generated by Phantom V1210%
%Author: Claire Hur%
%7/24/13%
%This script finds the shape properties of droplets using Image toolbox of
%MATLAB and other m-Files obtained from MATLAB FIleExchange. In order to
%run this script, several sub function files are required. The list of
%files are as following:
%cineInfo.m, BGI_prep.m, LinLUT.m, cineRead.m, substr.m, dlmcell.m
%User has to change variables to match individual image conditions.
%The list of variable required to modify are:
%Directory of Cine file, CineFileName, thresh, FilterSize, AreaThresh.
%Prior to run the code, Background Image cine file has to be prepared,
%named as CineFileName_BGI.cine (less than 50 frames) and save in the same directory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code is modified on 7/24/13 to use parallel computing tool box in
%MATLAB R2012a. It requires FileExchange file called "parfor_progress.m"
%to see the progress of parallel computing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all
fclose('all');

%Change the directory to where the cine file is
% cine_folder='C:\Users\Diego\Documents\MATLAB\JHU\HUR\code\data';
cine_folder='Z:\Encapsulation\Cell-Encapsulation-Videos-DA\3um';
CineFileName='15.56uL-0.045vv-3um-H363C-long1280-40x.cine'; %Cine file name

mkdir(['data\processed\',CineFileName(1:end-5)])
cd(['data\processed\',CineFileName(1:end-5)])

% parpool('local', 2) % CHANGED THIS FEOM 11 to 2
warning off

thresh=0.01; %imfill thresh hold value. It might vary with image quality
Type='canny'; %Edge detection type


FilterSize=15;
AreaThresh=300; %Debri removal threshhold [px] 10uL/min:200
DefThresh=0.4;  %Deformation threshhold (remove artificial particles/debris)
SolidityThresh=0.7; %The proportion of the pixels in the convex hull that are also in the reasion. Computed as Area/ConvexArea
IntensityThresh=50; %IntensityCheck function threshold
Alpha= 5;

importfile('LinLUT.mat'); %a conversion between packed 10bit data to real 10bit data
Result=[];
reSize=3;
% centroids=[];
% NewXc= [];
% NewYc= [];
% NewA= [];
% NewB=[];
% NewTheta=[];

%wait=2;

%Parameters that need to be changed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Change the directory to where the cine file is
%cd('Y:\Water droplets in FC40\2% Krytox in FC40 in 50X143\60%\876uL_min')

skip=1000;%Number of frames to skip
CineFileLocation=cine_folder
info=cineInfo(CineFileLocation,CineFileName);
NumberImages=info.NumFrames
FirstFrameToAnal=NumberImages-1;
StartFrame=info.startFrame;
LastFrame=info.endFrame;
ActualFirstIndex=abs(StartFrame-LastFrame+FirstFrameToAnal);
ImagesToProcess=abs(ActualFirstIndex+StartFrame-info.endFrame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Save Image Anlaysis Condition
NewHeader={'EdgeFilter' 'Edge Thresh','Median Filter Size' 'Area Thresh' 'Solidity Thresh' 'Deformability Thresh' 'Alpha' 'Skipped Frames' 'First Frame' 'Intensity Thresh'};
SecondLine={Type, thresh,FilterSize,AreaThresh,SolidityThresh,DefThresh,Alpha, skip,FirstFrameToAnal,IntensityThresh};
Condition=vertcat(NewHeader,SecondLine);
ConditionFileName=['Image_Adjustment_Condition_',substr(CineFileName,0,-5),'.csv'];
dlmcell(ConditionFileName,Condition,',');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Average Background Image Generation
AVE_BGI_frame=BGI_prep_auto(CineFileLocation,CineFileName,100,LinLUT,reSize);

% DIEGO
image_processing

%%
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Batch processing
NN=floor(ImagesToProcess/skip);
zzz=[1:1:NN];
Ind_pre=ActualFirstIndex+skip.*zzz-(skip-1);

%Sample image
SampleX=cineRead(CineFileLocation,CineFileName,Ind_pre(1),info,LinLUT);
SampleX=imresize(SampleX,reSize,'bicubic');
SampleX=imadjust(SampleX);
[mrows,ncols]=size(SampleX);
segmentedImageSequence=zeros(mrows,ncols,NN,class(SampleX));

% Store results in a struct ...
%resAsStruct = repmat( struct( 'answer', [] ), [NN, 12] );
resAsCell={};
%parfor_progress(NN); % Initialize
parfor zz=1:NN
    try
        tmpRow = [];tmpRowCell={};
        NewXc= []; NewYc= [];
        NewA= []; NewB= [];
        NewTheta=[];
        %Image Subtraction
        Index=Ind_pre(zz);
        CurrentFrameIndex=num2str(abs(Index+info.startFrame))
        
        switch length(int2str(abs(Index+info.startFrame)))
            case 1
                CurrentFrameIndex=['0000',int2str(abs(Index+info.startFrame))];
            case 2
                CurrentFrameIndex=['000',int2str(abs(Index+info.startFrame))];
            case 3
                CurrentFrameIndex=['00',int2str(abs(Index+info.startFrame))];
            case 4
                CurrentFrameIndex=['0',int2str(abs(Index+info.startFrame))];
            case 5
                CurrentFrameIndex=[int2str(abs(Index+info.startFrame))];
        end
        
        %         SubtractedTitle=['Subtracted_',substr(CineFileName,0,-5),'_',CurrentFrameIndex,'.tif'];
        %         FolderEx2=exist('Subtracted','dir');
        %         if FolderEx2==7
        %             cd('Subtracted');
        %             SkipThisLoop=exist(SubtractedTitle,'file')
        %             cd ..
        %         else
        %             SkipThisLoop=0;
        %         end
        %
        %         if SkipThisLoop == 0
        X=cineRead(CineFileLocation,CineFileName,Index,info,LinLUT);
        X=imresize(X,reSize,'bicubic');
        X=imadjust(X);
        
        New_X=AVE_BGI_frame-X;
        %figure,imshow(imadjust(New_X));
        
        %Image filter to get rid of small defects
        FilteredImage = medfilt2(New_X,[FilterSize FilterSize]);
        % figure, imshow(FilteredImage);
        % pause(wait);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Sharpen th image
        sharpFilter = fspecial('disk',Alpha);
        FilteredImage = imfilter(FilteredImage, sharpFilter, 'replicate');
        % figure, imshow(FilteredImage);
        % pause(wait);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        BW = edge(FilteredImage,Type,thresh);
        % imshow(BW);
        %  pause(wait);
        BW2 = imfill(BW,4,'holes');
        %  imshow(BW2);
        %  pause(wait);
        
        
        cc = bwconncomp(BW2, 8); %find objects with 8-connected pixels
        PossibleDroplets=cc.NumObjects %Number of droplets in the image
        
        if cc.NumObjects~=0
            graindata = regionprops(cc, 'all');
            centroids=cat(1,graindata.Centroid);
            Perimeters=cat(1,graindata.Perimeter);
            Orientations=cat(1,graindata.Orientation);
            MajorAxis=cat(1,graindata.MajorAxisLength);
            MinorAxis=cat(1,graindata.MinorAxisLength);
            EquivDiameter=cat(1,graindata.EquivDiameter);
            Area=cat(1,graindata.Area);
            FilledArea=cat(1,graindata.FilledArea);
            Solidity=cat(1,graindata.Solidity);
            %pause(wait);
            
            
            %Shape Fators
            Deformation=zeros(size(MajorAxis));
            Circularity=zeros(size(MajorAxis));
            AspectRatio=zeros(size(MajorAxis));
            AreaRatio=zeros(size(MajorAxis));
            for ii=1:cc.NumObjects
                Deformation(ii)= (MajorAxis(ii)- MinorAxis(ii))/(MajorAxis(ii)+ MinorAxis(ii));
                Circularity(ii)=4*pi*Area(ii)/Perimeters(ii)^2;
                AspectRatio(ii)= MinorAxis(ii)/MajorAxis(ii);
                AreaRatio(ii)=FilledArea(ii)/Area(ii);
            end
            Deformation=Deformation';
            Circularity=Circularity';
            AspectRatio=AspectRatio';
            
            
            %Save results in text file
            ou=size(Result,1);
            u=length(MajorAxis);
            if u(1)>0
                
                for yy=1:u(1)
                    %Debri Check
                    %Yes=IntensityCheck(X,centroids(yy,1),centroids(yy,2),MajorAxis(yy),Orientations(yy),IntensityThresh);
                    if Circularity(yy)<1 && Deformation(yy)<DefThresh && Solidity(yy)> SolidityThresh %&& Yes==1 && FilledArea(yy) >AreaThresh
                        tmpRow = abs(Index+info.startFrame);                     %image number
                        tmpRow=[tmpRow,centroids(yy,1)];          %particle center X
                        tmpRow=[tmpRow,centroids(yy,2)];          %particle center Y
                        tmpRow=[tmpRow,Orientations(yy)];         %particle tilting Angle [degree,-90 ~ 90]
                        tmpRow=[tmpRow,MajorAxis(yy)];              %Major Axis a
                        tmpRow=[tmpRow,MinorAxis(yy)];              %Minor Axis b
                        tmpRow=[tmpRow,Deformation(yy)];            %Deformation (a-b)/(a+b)
                        tmpRow=[tmpRow,Circularity(yy)];            %Circularity 4*pi*A/P^2
                        tmpRow=[tmpRow,AspectRatio(yy)];            %Aspect Ratio Minor/Major
                        tmpRow=[tmpRow,EquivDiameter(yy)];          %Equivalent Diameter [px]
                        tmpRow=[tmpRow,FilledArea(yy)];                   %Droplet Area
                        tmpRow=[tmpRow,Solidity(yy)];                   %Solidity
                        
                        %New values for overlay plotting
                        NewXc= [NewXc,centroids(yy,1)];
                        NewYc= [NewYc,centroids(yy,2)];
                        NewA= [NewA,MajorAxis(yy)];
                        NewB= [NewB,MinorAxis(yy)];
                        NewTheta=[NewTheta,Orientations(yy)];
                        tmpRowCell=[tmpRowCell;{tmpRow}];
                        %                 else
                        %                     NewXc= [];
                        %                     NewYc= [];
                        %                     NewA= [];
                        %                     NewB= [];
                        %                     NewTheta=[];
                    end
                end
                
            else u(1)==0
                fprintf('There are No Droplets\n');
            end
            
            NewXc= NewXc';
            NewYc= NewYc';
            NewA= NewA';
            NewB= NewB';
            NewTheta=NewTheta';
            
            ActualDroplets=size(NewXc,1)%Including debris and not-focused droplets
            
            
            OriginalTitle=['Original_',substr(CineFileName,0,-5),'_',CurrentFrameIndex,'.tif'];
            OverlayTitle=['Overlay_',substr(CineFileName,0,-5),'_',CurrentFrameIndex,'.tif'];
            SubtractedTitle=['Subtracted_',substr(CineFileName,0,-5),'_',CurrentFrameIndex,'.tif'];
            
            % SAVE ORIGINAL
            imwrite(X,OriginalTitle);     	% Write Original to file
            FolderEx=exist('Original','dir');
            if FolderEx~=7
                mkdir('Original');
            end
            movefile(OriginalTitle,['.\Original\',OriginalTitle]);
            
            % CALCULATE OVERLAY
            [indX,map]=gray2ind(X);       	% Save as index format
            rgbX = ind2rgb(indX,map);      	% Save as RGB format
            rad = 2;                        % radius of center pixels
            %Create Overlay image for Centroid (all droplets)
            coords={};
            if size(NewXc,1)~=0
                for jj=1:size(NewXc,1)
                    if (NewXc(jj)- rad)>1
                        if (NewYc(jj)- rad)>1
                            % find ellipse coordinates and plot
                            [x,y]=ellipse2(NewA(jj)/2,NewB(jj)/2,NewTheta(jj)/(-180)*pi,NewXc(jj,1),NewYc(jj,1),'b');
                            coords{jj}=[x,y];
                            x=round(x);
                            y=round(y);
                            BMW = roipoly(rgbX,x,y);
                            [y,x] = find(BMW==1);
                            for r=1:length(x)
                                rgbX(y(r),x(r),1)=0;
                                rgbX(y(r),x(r),2)=0;
                                rgbX(y(r),x(r),3)=255;
                            end
                            % plot centers with radius rad
                            for j=-rad:rad
                                for k=-rad:rad
                                    NewXc(jj) = round(NewXc(jj));
                                    NewYc(jj) = round(NewYc(jj));
                                    rgbX(NewYc(jj)+j,NewXc(jj)+k,1)=255;
                                    rgbX(NewYc(jj)+j,NewXc(jj)+k,2)=0;
                                    rgbX(NewYc(jj)+j,NewXc(jj)+k,3)=0;
                                end
                            end
                        end
                    end
                end
                
            end
        end
        
        % if g~=0
        % SAVE OVERLAY
        imwrite(rgbX,OverlayTitle);     	% Write Overlay to file
        FolderEx1=exist('Overlay','dir');
        if FolderEx1~=7
            mkdir('Overlay');
        end
        movefile(OverlayTitle,['.\Overlay\',OverlayTitle]);
        
        % SAVE SUBTRACTED
        imwrite(BW2,SubtractedTitle);    % Write Subtracted to file
        FolderEx2=exist('Subtracted','dir');
        if FolderEx2~=7
            mkdir('Subtracted');
        end
        movefile(SubtractedTitle,['.\Subtracted\',SubtractedTitle]);
        
        %
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        resAsCell=[resAsCell;{tmpRowCell}];
        %parfor_progress; % Count
        
        u=0;
        %         end
    catch ME
        disp(zz);
        disp(CurrentFrameIndex);
        throw(ME)
    end
    
end
%parfor_progress(0); % Clean up
% Save data
%Column Header
Header = {'Image#' 'Xc' 'Yc' 'Orientation (degree)' 'Major' 'Minor' 'Deformation' 'Circularity' 'Aspect Ratio' 'Equivalent Diameter' 'Area' 'Solidity'};
ResultArray=cell2mat(resAsCell{1});
for ss=2:size(resAsCell,1)
    NextImageCell=cell2mat(resAsCell{ss});
    ResultArray=vertcat(ResultArray,NextImageCell);
end
CombinedResult=vertcat(Header,num2cell(ResultArray));
OutputFileName=['Droplet_Property_',substr(CineFileName,0,-5),'.csv'];
dlmcell(OutputFileName,CombinedResult,',');

cd('..')
cd('..')

toc
close all
delete(gcp('nocreate'))

