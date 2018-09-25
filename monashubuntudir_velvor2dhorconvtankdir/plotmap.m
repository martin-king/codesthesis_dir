clear all; close all;

%fid=fopen('./temp2.dat','r','ieee-le');
% read into 'data' in one very tall column
%[data,count]=fread(fid,inf,'real*4');
%fclose(fid);
ascii=1;

if ascii==1

load('-ascii','plotmap.dat');
temp2=plotmap;
clear plotmap;
temp2=reshape(temp2,130,130);
%theta=linspace(0,1.,113);
%r=linspace(0.5,1.0,98);

%x=[];
%y=[];
%for angle=1:113
% x=[x; r.*cos(theta(angle))];
% y=[y; r.*sin(theta(angle))];
%end
x=linspace(0.,1,130);
y=linspace(0.,1,130);
[XX,YY]=meshgrid(x,y);

%x=x';
%y=y';
%subplot(1,2,1)
figure(2),clf
temp2=temp2;
pcolor(XX,YY,flipud(temp2));
shading interp;
colormap jet;
%axis off;
%subplot(1,2,2)
%contour(x,y,temp2,10,'LineWidth',2);

end
