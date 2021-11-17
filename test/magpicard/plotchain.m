function fin_sol=plotchain(Nsteps,Np,fname)
%Nsteps=64;
%Np=20;
sol=zeros(2,Nsteps,Np,Np);
F=zeros(2,Nsteps,Np,Np);
%  Read the solution
for k = 1:Nsteps
  file=strcat(fname,'/P_step_',sprintf("%5.5i ",k));
  D=load(file);
  sol(1,k,:,:)=reshape(D(1,:),Np,Np);
  sol(2,k,:,:)=reshape(D(2,:),Np,Np);
end

%  Read the Function values (Facke matrix)
for k = 1:Nsteps
  file=strcat(fname,'/F_step_',sprintf("%5.5i ",k));
  D=load(file);
  F(1,k,:,:)=reshape(D(1,:),Np,Np);
  F(2,k,:,:)=reshape(D(2,:),Np,Np);
end

figure(10); clf
subplot(2,1,1)
for j = 1:Np
    Dg=squeeze(sol(1,:,j,j));
    plot(Dg); hold on;
end
hold off;
subplot(2,1,2);
%figure(2); clf
for j = 1:Np-1
    Dg=squeeze(sol(2,:,j+1,j));
    plot(Dg); hold on;
end
hold off;
% $$$ figure(3); clf
% $$$ for j = 1:Nsteps
% $$$     trce(j)=trace(squeeze(sol(1,j,:,:)));
% $$$ end
% $$$ semilogy(abs(trce-trce(1)))
% $$$ hold off;


fin_sol=diag(squeeze(sol(1,end,:,:)));