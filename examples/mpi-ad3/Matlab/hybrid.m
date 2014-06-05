clear;
N = 128;
k=[1:N/2]*2*pi;
dx=1/N;
L2=(-2+2*cos(k*dx))/(dx*dx);
L4=(-30+32*cos(k*dx)-2*cos(2*k*dx))/(12.0*dx*dx);

Nj = 1000;
for j = 1:Nj
    nu=0.1*j;
    stp = 1 + nu*L2./(1-nu*L2);
    rho22(j) = max(abs(stp));
    stp = 1 + nu*L4./(1-nu*L4);
    rho44(j) = max(abs(stp));
    stp = 1 + nu*L4./(1-nu*L2);
    rho42(j) = max(abs(stp));
    plot(stp);
    pause
end

plot(rho22,'-k'); hold on;
plot(rho44,'-b') 
plot(rho42,'-r'); hold off;