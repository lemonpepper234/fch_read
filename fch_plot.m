figure;
xmin=-6;
xmax=6;
ymin=-6;
ymax=6;
data = load ('Au13opt6d10f121.dat');
data1=load('sphere_integrating121.txt');
data2=load('sphere_integrating122.txt');
data3=load('sphere_integrating127.txt');
contourf(linspace(xmin,xmax,100),linspace(ymin,ymax,100),data);
axis square;

figure;
surf(linspace(xmin,xmax,100),linspace(ymin,ymax,100),data,EdgeColor="none");
axis square;

figure;
plot(linspace(xmin,xmax,30),data1);
hold on
plot(linspace(xmin,xmax,30),data2);
hold on
plot(linspace(xmin,xmax,30),data3);
axis square;

legend('P-orbit', 'D-orbit', 'S-orbit');