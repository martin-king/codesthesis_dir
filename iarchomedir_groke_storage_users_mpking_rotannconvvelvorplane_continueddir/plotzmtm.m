load('-ascii','zmtmtempvarr9.dat');

r=linspace(0.5,1.0,242);

figure(1),clf
plot(r,zmtmtempvarr9(:,1),'k+');

figure(2),clf
plot(r,zmtmtempvarr9(:,3),'k+');

figure(3),clf
plot(r,zmtmtempvarr9(:,2),'k+');

