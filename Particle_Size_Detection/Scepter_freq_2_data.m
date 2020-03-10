% Counts is the counts vector - the 2nd column in the scepter csv output file 
% Bin is the Bin vector - the result of converting the volume (3rd column to radius, by the formula
% R = (3V/(4 pi))^(1/3) * 10 , 10 being the conversion factor.
if exist ('path','var')==0
    clear;close all;clc;
else
    clearvars -except path,
    cd(path)     
end
[name,path]=uigetfile;
min_cell_size=input('What is the expected minimum cell size');
L=xlsread(name);
Count=L(:,2);
Bin=L(:,3);
RBin= (3.*Bin/(4*3.14)).^(1/3)*20;
data=[];
data(1:Count(1))= RBin(1);
for i=2:length(Count)
    l=length(data);
    data(l+1:l+(Count(i)))=RBin(i);
end
cutoff=find(data < min_cell_size);
data=data(cutoff(end):end);
datah=histcounts(data,[5.5:1:20.5],'Normalization','Probability');