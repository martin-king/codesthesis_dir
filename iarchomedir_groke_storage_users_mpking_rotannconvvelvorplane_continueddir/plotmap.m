clear all;

%fid=fopen('./temp2.dat','r','ieee-le');
% read into 'data' in one very tall column
%[data,count]=fread(fid,inf,'real*4');
%fclose(fid);
ascii=1;

if ascii==1

load('-ascii','plotmap.dat');
temp2=plotmap;
clear plotmap;
temp2=reshape(temp2,242,1793);
theta=linspace(0,2*pi,1793);
r=linspace(0.5,1.0,242);

x=[];
y=[];
for angle=1:1793
 x=[x; r.*cos(theta(angle))];
 y=[y; r.*sin(theta(angle))];
end
x=x';
y=y';
%subplot(1,2,1)
%pcolor(x,y,temp2);
%shading interp;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ;
%colormap jet;
contour(x,y,temp2,[0.65:0.05:1.0]','-k');
hold on;
contour(x,y,temp2,[0.:0.05:0.65]','--k');
contour(x,y,temp2,[0.65]','-k','LineWidth',2);
axis off;
%subplot(1,2,2)
%contour(x,y,temp2,10,'LineWidth',2);

end
