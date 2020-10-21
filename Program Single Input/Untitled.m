clear all;
close all;

t = 0:0.001:4*pi;
x = sin(2*pi*t);
plot(x)
axis([0 2000, -2 2])
text(1200,1.2,'S1')
text(200,1.2,'S1')
text(700,-1.2,'S3')