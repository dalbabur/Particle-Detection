%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shape Factor Determination Code
%Written by Claire Hur on 2/2/13
%Modified to to parfor on 8/21/13
%This need SplitVec.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all
warning off
fclose('all');

tic
matlabpool local 11

M = importdata('Analysis Sequence4.csv', ',', 1);
NF=size(M.data,1); %number of flow rate
ParentDirectList=M.textdata((2:end),1);
for dd=22:NF
    ParentDir=ParentDirectList{dd};
    FlowRateNum=M.data(dd,1);
    FlowRateText=[num2str(FlowRateNum),'uL_min'];
    main=[ParentDir,'\',FlowRateText];
    
    if M.data(dd,2)==1
        NumMovie=1;
    else
        NumMovie=M.data(dd,2);
    end
    
    
    for ww=1:NumMovie
        cd(main);
%         if NumMovie==1
%             FileName=['Droplet_Property_',FlowRateText,'.txt'];
%         else
            FileName=['Droplet_Property_',FlowRateText,'_',num2str(ww),'.txt'];
        %end
        Data=importdata(FileName);
        FlowCond=[substr(FileName,17,-4)];
        FlowRate=sscanf(substr(FlowCond,0,-2),'%f');
        VIS=['off'];
        
        
        
        %Assign valiable name for data set
        ImageNum=Data.data(:,1);
        Xc=Data.data(:,2);
        Yc=Data.data(:,3);
        Major=Data.data(:,5);
        D=Data.data(:,10);
        Area1=Data.data(:,11);
        
        
        switch FlowRate
                        case 10
                            DefThresh=0.015; %Deformation Threshold
                            %DiffThresh=0.75;%difference btw edge and ellipse fitting
                            InSideRatio=0.5;
                            AreaThresh1=500;
            case 50
                DefThresh=0.05;
                %DiffThresh=0.7;
                InSideRatio=0.5;
                AreaThresh1=500;
                
            case 100
                DefThresh=0.05;
                %DiffThresh=0.86;
                InSideRatio=0.7;
                AreaThresh1=500;
                %AreaThresh=800;
                %Edge_over_In_thr=3;
            case 165
                DefThresh=0.1;
                % DiffThresh=0.6;
                InSideRatio=0.5;
                AreaThresh1=500;
                %             case 165
                %                 DefThresh=0.3;
                %                 %DiffThresh=0.85;
                %                 InSideRatio=0.5;
                %                 AreaThresh1=500;
            case 248
                DefThresh=0.35;
                %DiffThresh=0.88;
                InSideRatio=0.5;
                AreaThresh1=500;
            
            case 331
                DefThresh=0.4;
                %DiffThresh=0.6;
                InSideRatio=0.5;
                AreaThresh1=500;
                
            case 879
                DefThresh=0.4;
                %DiffThresh=0.6;
                InSideRatio=0.5;
                AreaThresh1=500;
                
        end
        
        CropFactor=1.4;
        %InSideRatio=0.7;
        
        %Filter=10;
        
        
        TotalParticle=length(ImageNum);
        
        Dir=[main,'\Original'];
        Dir2=[main,'\Overlay'];
        Dir3=[main,'\Analysis'];
        cd(main);
        
        %Save Cropping Anlaysis Condition
        NewHeader={'InSide Diameter for Intensity' 'Def Thresh' 'Croppint Factor' };
        SecondLine={ InSideRatio, DefThresh, CropFactor };
        Condition=vertcat(NewHeader,SecondLine);
        ConditionFileName=['Crop_Anal_Condition_',FlowCond,'.csv'];
        dlmcell(ConditionFileName,Condition,',');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        resAsCell={};
        %Create cropped images
        parfor ii=1:TotalParticle
            try
                tmpRow=[]; tmpRowCell={};CP_sort={};
                InnerEdge={};OuterEdge={};
                %Elliminate the small contiminating particles
                if Area1(ii) < AreaThresh1
                    continue
                end
                close all
                Contour=[];
                ContourFile=[];
                CP={};
                InnerEdge={};
                OueterEdge={};
                MidEdge={};
                X=[];R2=[];
                
                Area1(ii);
                %     Particle_to_Anal=TotalParticle-ii
                %     Particle_found=size(Result,1)
                %
                CurrentImage=ImageNum(ii);
                switch length(int2str(CurrentImage))
                    case 1
                        CurrentFrameIndex=['0000',int2str(CurrentImage)];
                    case 2
                        CurrentFrameIndex=['000',int2str(CurrentImage)];
                    case 3
                        CurrentFrameIndex=['00',int2str(CurrentImage)];
                    case 4
                        CurrentFrameIndex=['0',int2str(CurrentImage)];
                    case 5
                        CurrentFrameIndex=[int2str(CurrentImage)];
                end
                
                
                ImageTitle=['Original_',FlowCond,'_',num2str(CurrentFrameIndex),'.tif'];
                cd(Dir);
                X=imread(ImageTitle);
                if size(X,3)>1
                    X=rgb2gray(X);
                end
                
                
                HalfWidth=(CropFactor*Major(ii))/2;
                X_corner=Xc(ii)-HalfWidth;
                Y_corner=Yc(ii)-HalfWidth;
                %Elliminate the error caused by particles located at the edges of
                %images
                if X_corner < 0
                    if Y_corner < 0
                        New_width=min(abs(X_corner),abs(Y_corner));
                        X_corner=0;
                        Y_corner=0;
                    else
                        New_width=abs(X_corner);
                        X_corner=0;
                    end
                    
                else if X_corner >0
                        if (X_corner+2*HalfWidth)-size(X,2) <0
                            if Y_corner < 0
                                New_width=abs(Y_corner);
                                Y_corner=0;
                            else
                                New_width=0;
                            end
                        else
                            continue
                        end
                    end
                end
                %Define cropping map & crop the image
                map=[X_corner Y_corner 2*HalfWidth-New_width 2*HalfWidth-New_width];
                
                
                X_crop=imcrop(X,map);
                %figure('Position',[200,1500,400,300],'visible',VIS),imshow(X_crop);
                %pause(0.2);
                %    h=fspecial('disk',Filter);
                %
                X_crop1=mat2gray(X_crop);
                X_crop=X_crop1;
                ImageShift=size(X_crop1,1)/2;
                
                %Overlay Crop
                cd(Dir2)
                OverTitle=['Overlay_',FlowCond,'_',num2str(CurrentFrameIndex),'.tif'];
                X_Over=imread(OverTitle);
                X_factor=size(X,1)/size(X_Over,1);
                Y_factor=size(X,2)/size(X_Over,2);
                SizeFac=(X_factor+Y_factor)/2;
                X_Over=imresize(X_Over,SizeFac);
                %pause(0.2);
                
                if map(1) > size(X_Over,2)
                    cd(Dir);
                    continue
                end
                X_Over_crop=imcrop(X_Over,map);
                X_Over_crop=mat2gray(X_Over_crop);
                %figure('Position',[600,1500,400,300],'visible',VIS),imshow(X_Over_crop);
                cd(Dir);
                
                %Crop the background image
                cd(main)
                BGITitle=[FlowCond,'_BGI','.tif'];
                X_BGI=imread(BGITitle);
                X_BGI_crop=imcrop(X_BGI,map);
                X_BGI_crop=mat2gray(X_BGI_crop);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                X_crop1=X_BGI_crop-X_crop1;
                
                %     se2 = strel('disk',5);
                %     X_crop1 = imtophat(X_crop1,se2); %Getting rid of illumination variation
                X_crop1=imadjust(X_crop1);
                %     se3 = strel('disk',10);
                %     X_crop1 = imclose(X_crop1,se3);
                
                %Find standard deviation of region in the image (75% of the major/2)
                pcimg=imgpolarcoord(X_crop1);
                %     figure,imshow(pcimg);
                %      pause(0.2);
                
                InSide=ceil(D(ii)/2*InSideRatio);
                pcimg_crop=pcimg(1:InSide,1:end);
                %     imshow(pcimg_crop);
                %     pause(0.2)
                %[N1 intensity1]=hist(pcimg_crop);
                %Int_sum=sum(sum(pcimg_crop));
                Int_std=std2(pcimg_crop);
                Int_mean=mean2(pcimg_crop);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                for jj=1:size(pcimg,2)
                    x1=[jj jj];
                    y1min=ceil(0.8*InSide);
                    %    y1min=1;
                    y1=[y1min size(pcimg,1)];
                    h=improfile(pcimg, x1,y1); %find the intensity profile
                    Yaxis_px=[y1min:1:size(pcimg,1)];
                    
                    ImageAxis=[1:1:size(pcimg,2)];
                    
                    h_new=(max(h)-h);
                    A2=h_new-1.5*Int_std+Int_mean;
                    %A2_all{jj}=A2;
                    %     plot(A2);
                    
                    CP{jj}=find(A2(1:end-1).*A2(2:end)<0); %points where intensity profile crosses zero.
                    CP_sort{jj}=sort(CP{jj});
                    if y1min ~=1
                        InnerEdge{jj}=min(CP{jj})+y1min-1;
                        OuterEdge{jj}=min(CP_sort{jj}(2:end))+y1min-1;
                    else
                        InnerEdge{jj}=min(CP{jj});
                        OuterEdge{jj}=min(CP_sort{jj}(2:end));
                    end
                    
                    
                    %If there are no edge found, match the average point.
                    if isempty(InnerEdge{jj})
                        InnerEdge{jj}=mean2(cell2mat(InnerEdge(1:jj-1)));
                    end
                    if isempty(OuterEdge{jj})
                        OuterEdge{jj}=mean2(cell2mat(OuterEdge(1:jj-1)));
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    MidEdge{jj} = (InnerEdge{jj}+OuterEdge{jj})/2; %Middle point
                end
                
                
                
                
                %Check if there is any empty matrix
                InnerEdge1=cell2mat(InnerEdge);
                OuterEdge1=cell2mat(OuterEdge);
                MidEdge1=cell2mat(MidEdge);
                
                
                T_inner=isnan(InnerEdge1);
                T_outer=isnan(OuterEdge1);
                
                %     if mean(abs(OuterEdge1-size(pcimg,1))/size(pcimg,1)) < 0.1
                %         continue
                %     end
                
                if sum(T_inner)/360==1||sum(T_outer)/360==1
                    continue
                end
                
                
                InnerEdge1 = smooth(ImageAxis,InnerEdge1,0.15,'rloess');%Inner Edge
                OuterEdge1 = smooth(ImageAxis,OuterEdge1,0.15,'rloess');%Outer Edge
                MidEdge1 = smooth(ImageAxis,MidEdge1,0.15,'rloess');%Middle Edge
                
                gap=abs(OuterEdge1-InnerEdge1);
                gap_mean=mean(gap);
                
                %Middle Edge
                xx=MidEdge1'.*cos(ImageAxis./180*pi);
                yy=MidEdge1'.*sin(ImageAxis./180*pi);
                
                %Inner Edge
                xx_in=InnerEdge1'.*cos(ImageAxis./180*pi);
                yy_in=InnerEdge1'.*sin(ImageAxis./180*pi);
                
                %Outer Edge
                xx_out=OuterEdge1'.*cos(ImageAxis./180*pi);
                yy_out=OuterEdge1'.*sin(ImageAxis./180*pi);
                
                %pause
                
                %Make contour array
                Contour(:,1)=xx'+ImageShift; %Middle Edge Contour
                Contour(:,2)=yy'+ImageShift; %Middle Edge Contour
                Contour(:,3)=xx_in'+ImageShift; %Inner Edge Contour
                Contour(:,4)=yy_in'+ImageShift; %Inner Edge Contour
                Contour(:,5)=xx_out'+ImageShift; %Inner Edge Contour
                Contour(:,6)=yy_out'+ImageShift; %Inner Edge Contour
                
                ContourArray=num2cell(Contour);
                
                %Ellipse fit
                AAA=ellipse_fit(yy+ImageShift,xx+ImageShift);
                if isnan(AAA)
                    continue
                end
                
                [mj mn X0 Y0 phi]=ellipse_fit(yy+ImageShift,xx+ImageShift);
                
                R=mj*mn./sqrt((mn.*sin(ImageAxis./180*pi+phi)).^2+(mj.*cos(ImageAxis./180*pi+phi)).^2);
                XX=R.*cos(ImageAxis./180*pi);
                YY=R.*sin(ImageAxis./180*pi);
                
                %Major axis line
                xxx=[X0+mj*sin(phi),X0-mj*sin(phi)];
                yyy=[Y0+mj*cos(phi), Y0-mj*cos(phi)];
                
                %figure([left bottom width height])
                %         figure('Position',[1000,1500,400,300],'visible',VIS),imshow(X_crop);
                %         hold on
                %         gg=plot(xx+ImageShift,yy+ImageShift,'y');
                %         hold off
                
                r=size(X_crop,1);
                c=size(X_crop,2);
                
                if ~isreal(XX)|| ~isreal(YY)
                    continue
                end
                
                %Darker Edge intensity compared to inner area
                IntSum=0;
                IntSumIn=0;
                countIn=0;
                countEdge=0;
                IntSumIn_ave=0;
                IntSum_ave=0;
                for ff=1:size(pcimg,2)
                    for fj=1:size(pcimg,1)
                        if fj < R(ff)-gap_mean/2 || fj > R(ff)+ gap_mean/2
                            countIn=countIn+1;
                            IntSumIn=IntSumIn+pcimg(fj,ff);
                        else
                            IntSum=IntSum+pcimg(fj,ff);
                            countEdge=countEdge+1;
                        end
                    end
                end
                
                IntSumIn_ave=IntSumIn/countIn;
                IntSum_ave=IntSum/countEdge;
                Edge_over_In=IntSum_ave/IntSumIn_ave;
                
                %Turn ellipse coordinate to BW image
                closeBW=poly2mask(abs(XX+ImageShift),abs(YY+ImageShift),r,c);
                %     figure,imshow(closeBW);
                %     pause;
                %Ellipse properties
                ecc=sqrt(1-mn^2/mj^2); %eccentricity
                [K,E]=ellipke(ecc^2);
                Circumference=4*mj*E;
                Def=(mj-mn)/(mj+mn);%deformability
                Area=pi*mj*mn;
                OrienT=abs(pi/2-phi)*180/pi;%orientation (degree)
                
                
                %Calculation of area differences to elliminate wigles
                Edge1=MidEdge1';
                Diff_sum=0;
                for ut=1:size(Edge1,2)
                    Diff_sum=Diff_sum+abs(Edge1(ut)-R(ut));
                end
                
                R2(ii)=1-Diff_sum/Area;
                
                [s,MSS,Mss_ID] = mkdir('Analysis');
                if s==0
                    cd('Analysis')
                else
                    cd('Analysis')
                end
                
                OriginalTitle=['Cropped_',FlowCond,'_','p',num2str(ii),'_',num2str(CurrentImage),'.tif'];
                BGITitle=['BGI_',FlowCond,'_','p',num2str(ii),'_',num2str(CurrentImage),'.tif'];
                
                
                imwrite(X_crop, OriginalTitle);
                FolderEx=exist('Original','dir');
                if FolderEx~=7
                    mkdir('Original');
                end
                movefile(OriginalTitle,['.\Original\',OriginalTitle]);
                
                
                imwrite(X_BGI_crop,BGITitle);
                FolderEx1=exist('BGI','dir');
                if FolderEx1~=7
                    mkdir('BGI');
                end
                movefile(BGITitle,['.\BGI\',BGITitle]);
                
                
                
                %creat data array
                tmpRow=CurrentImage;                     %image number
                tmpRow=[tmpRow,ii];                     %particle number
                tmpRow=[tmpRow,X0+X_corner];          %particle center X
                tmpRow=[tmpRow,Y0+Y_corner];          %particle center Y
                tmpRow=[tmpRow,OrienT];         %particle tilting Angle [degree,-90 ~ 90]
                tmpRow=[tmpRow,2*mj];              %Major Axis a
                tmpRow=[tmpRow,2*mn];              %Minor Axis b
                tmpRow=[tmpRow,Def];%Deformation (a-b)/(a+b)
                tmpRow=[tmpRow,Circumference];            %Perimeter
                tmpRow=[tmpRow,Area];            %Droplet Area
                tmpRow=[tmpRow,2*sqrt(Area/pi)];          %Equivalent Diameter [px]
                tmpRow=[tmpRow,4*pi*Area/Circumference^2];                   %Circularity 4*pi*A/P^2
                tmpRow=[tmpRow,R2(ii)];                     %Coefficient of determination
                tmpRow=[tmpRow,X_corner];                     %X_corner
                tmpRow=[tmpRow,Y_corner];                     %Y_corner
                
                tmpRowCell=[tmpRowCell;{tmpRow}];
                
                
                
                %Save Contour info
                ContourFile=['Contour_',FlowCond,'_','p',num2str(ii),'_',num2str(CurrentImage),'.txt'];
                dlmcell(ContourFile,ContourArray,',');
                FolderEx2=exist('Contour','dir');
                if FolderEx2~=7
                    mkdir('Contour');
                end
                movefile(ContourFile,['.\Contour\',ContourFile]);
                resAsCell=[resAsCell;{tmpRowCell}];
                
            catch ME
                disp(ii);
                disp(CurrentFrameIndex);
                throw(ME)
            end
        end
        
        cd(Dir3);
        %Save Results in txt file
        %Column Header
        Header = {'Image#' 'P#' 'Xc' 'Yc' 'Orientation (degree)' 'Major' 'Minor' 'Deformation' 'Perimeter' 'Area' 'Equivalent Diameter' 'Circularity' 'R^2' 'X_corner' 'Y_corner'};
        ResultArray=cell2mat(resAsCell{1});
        for ss=2:size(resAsCell,1)
            NextImageCell=cell2mat(resAsCell{ss});
            ResultArray=vertcat(ResultArray,NextImageCell);
        end
        CombinedResult=vertcat(Header,num2cell(ResultArray));
        
        %save Droplet_Property.csv CombinedResult -ascii -tabs;
        OutputFileName=['Polished_Droplet_Property_',FlowCond,'.txt'];
        dlmcell(OutputFileName,CombinedResult,',');
        
        open(OutputFileName)
    end
end
toc
close all
matlabpool close
