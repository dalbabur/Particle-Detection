%% Code to measure cell particles. Developed by Srivathsan Kalyan, Hur Lab @ JHU
cd('D:\Documents\Research\Temp Files From Drobo\T-Cell Electroporation\Gen 3.b\L1210 Inlet');
if exist ('currentpath','var')==0
    clear;close all;clc;
else
    clearvars -except currentpath
    cd(currentpath)
end
[name,currentpath]=uigetfile('.tif');
cd(currentpath)
info = imfinfo(name);
framenumber=input('What frame?');
A=imread(name,framenumber);
imshow(A)
rect=getrect;
SE3 = strel('disk', 5);
B=imtophat(A,SE3);
% imshow(B)
contrastAdjusted = imadjust(B);
% imshow(contrastAdjusted)
% imshow(B)
C=(A-B);
[i,j] = find(C >= 200);
% imshow(C)
[centers,radii] = imfindcircles(C,[1 5],'ObjectPolarity','dark', ...
          'Sensitivity',0.95);
index=[];
if isempty(centers)
    C=imadjust(C);
    imshow(C);
    rect=getrect;
    [centers,radii] = imfindcircles(C,[2 13],'ObjectPolarity','dark', ...
        'Sensitivity',.1);
    index=[];
end
for o=1:size(centers,1)
    if centers (o,1) < rect(1) || centers(o,1) > (rect(3)+rect(1))
        index(o,:)= [0,0];
    elseif centers(o,2) < rect(2) || centers(o,2) > (rect(4)+rect(2))
        index(o,:)= [0,0];
    else 
        index(o,:)=[1,1];
    end
end
newcenters=[];
newradius=[];
for p=1:size(centers,1)
    if index(p,1)==1
        newcenters=[newcenters;centers(p,:)];
        newradius= [newradius;radii(p,:)];
    end
end
% newradius=newradius+1;
imshow(A);
hold on
h=viscircles(newcenters,newradius);
numeration=1:size(newcenters,1);
texty=cell(1,size(newcenters,1));
for i=1:size(newcenters,1)
    texty{i}=num2str(numeration(i));
end
hold on
for i=1:size(newcenters,1)
    text(newcenters(i,1),newcenters(i,2),texty{i},'color','g')
end
text(rect(1)-15,rect(2)-15,[texty{1},' - ',texty{end}]);
savename=[name(1:end-4),' ',num2str(framenumber),'.tif'];
saveas(gcf,savename)
prompt  = {'How many more frames to analyze?'};
def     = {'0'};
dlgtitle   ='Frames';
lines   = 1;
answer1  = inputdlg(prompt,dlgtitle,lines,def);
imgnum=answer1{1};
numberofframesleft=str2double(imgnum);
runninglistcenters=[];
for n=1:numberofframesleft
    prompt  = {'Which frame to analyze next?'};
    def     = {'2'};
    dlgtitle   ='Frames';
    lines   = 1;
    figure(n+1)
    answer  = inputdlg(prompt,dlgtitle,lines,def);
    framenum=str2double(answer{1});
    currentframe=imread(name,framenum);
    imshow(currentframe)
    rect=getrect;
    [centers,radii] = imfindcircles(currentframe,[1 8],'ObjectPolarity','dark', ...
        'Sensitivity',.9);
    if isempty(centers)
        currentframe=imadjust(currentframe);
        imshow(currentframe);
        rect=getrect;
        [centers,radii] = imfindcircles(currentframe,[1 8],'ObjectPolarity','dark', ...
            'Sensitivity',.9);
    end
    for o=1:size(centers,1)
        if centers (o,1) < rect(1) || centers(o,1) > (rect(3)+rect(1))
            index(o,:)= [0,0];
        elseif centers(o,2) < rect(2) || centers(o,2) > (rect(4)+rect(2))
            index(o,:)= [0,0];
        else
            index(o,:)=[1,1];
        end
    end
    runninglistcenters(n)=size(newcenters,1);
    for p=1:size(centers,1)
        if index(p,1)==1
            newcenters=[newcenters;centers(p,:)];
            newradius= [newradius;radii(p,:)];
        end
    end
    imshow(currentframe);
    hold on
    h=viscircles(newcenters(runninglistcenters(n)+1:end,:),newradius(runninglistcenters(n)+1:end));
    numeration=runninglistcenters(n)+1:size(newcenters,1);
    texty=cell(1,size(newcenters,1)-runninglistcenters(n));
    for i=1:size(newcenters,1)-runninglistcenters(n)
        texty{i}=num2str(numeration(i));
    end
    hold on
    for i=runninglistcenters(n)+1:size(newcenters,1)
        text(newcenters(i,1),newcenters(i,2),texty{i-runninglistcenters(n)},'color','g')
    end
    text(rect(1)-15,rect(2)-15,[texty{1},' - ',texty{end}]);
    savename=[name(1:end-4),' ',answer{1},'.tif'];
    if exist(savename)~=0
        savename=[name(1:end-4),' ',answer{1},'p2','.tif'];
    end
    saveas(gcf,savename)
end
addendum1=[1:size(newcenters,1)]';
addendum2=ones(size(newcenters,1),1);
allnewcenters=[addendum1,newcenters,newradius,addendum2];
if isempty(n)
    n=0;
end
runninglistcenters(n+1)=size(newcenters,1);
finalbeadnumber=0;
for j=1:numberofframesleft+1
    currentcircles=newcenters(finalbeadnumber+1:runninglistcenters(j),:);
    currentradii=newradius(finalbeadnumber+1:runninglistcenters(j),:);
    figure(j)
    l=size(currentcircles,1);
    for bead=1:size(currentcircles,1)
        xlim([currentcircles(bead,1)-30,currentcircles(bead,1)+30])
        ylim([currentcircles(bead,2)-30,currentcircles(bead,2)+30])
        if bead==1
            pause
        end
        answer3 = questdlg('Is the bead completely circled?','Program Menu', ...
            'Yes','No','Get Bead Center + Radius','Yes'); %check if user is satisfied with the cell
        switch answer3
            case 'Yes'
                allnewcenters(bead+finalbeadnumber,5)=1;
            case 'No'
                allnewcenters(bead+finalbeadnumber,5)=0;
            case 'Get Bead Center + Radius'
                allnewcenters(bead+finalbeadnumber,5)=3;
                [x,y]=getpts;
                centerx=(x(2)+x(1))/2;
                centery=(y(2)+y(1))/2;
                allnewcenters(bead+finalbeadnumber,2:3)=[centerx,centery];
                distance=sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2)/2;
                allnewcenters(bead+finalbeadnumber,4)=distance;
                if distance >= 10
                    allnewcenters(bead+finalbeadnumber,5)=0;
                end
        end
        
    end
    finalbeadnumber=runninglistcenters(j);
end
excelname=['Measured beads for ',name(1:end-4),'.xlsx'];
xlswrite(excelname,allnewcenters)
close all