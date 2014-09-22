%  Routine to plot error versus iteration by Vcycle
clear
Niter=10
figure(1);clf
figure(2);clf
Nx=128; Nstep=64; Nproc=Nstep;
%Nx=32; Nstep=128; Nproc=Nstep;

N_Vrack=[1,2,3,4,10]

for NN = 1:length(N_Vrack)
   N_V=N_Vrack(NN)
fbase=['pfasst_V_'];
fspec=['Niter',num2str(Niter,'%02d'),'_Nx',num2str(Nx,'%03d'),'_Nstep',num2str(Nstep,'%03d'),'_Nproc',num2str(Nproc,'%03d'),'_',num2str(N_V,'%03d'),'_',num2str(Nstep,'%03d')];
fname=['../Dat/',fbase,fspec,'.m']
q=load(fname);
q_end_ind = find(q(:,1)==3 );
q128=q(q_end_ind,:);

iter = q128(:,4);
pde_err = q128(:,6);
ode_err = q128(:,7);
res = q128(:,8);

switch(N_V)
  case 1,
   clr = '-bx'
  case 2,
   clr = '-rx'
  case 3,
   clr = '-gx'
  case 4,
   clr = '-kx'
  case 5,
   clr = '-mx'
  case 6,
   clr = '-bo'
  case 7,
   clr = '-ro'
  case 8,
   clr = '-go'
  case 9,
   clr = '-ko'
  case 10,
   clr = '-mo'
end
% $$$ xup = 1e-4
% $$$ figure(2)
% $$$ subplot(1,3,1)
% $$$ semilogy(iter,ode_err,clr,'MarkerSize',10); 
% $$$ axis([0 10,1e-7,xup])
% $$$ title('ODE Err','FontSize',14), hold on;
% $$$ xlabel('Iteration','FontSize',14)
% $$$ set(gca,'FontSize',14)
% $$$ subplot(1,3,2)
% $$$ semilogy(iter,pde_err,clr,'MarkerSize',10); 
% $$$ axis([0 10,1e-7,xup])
% $$$ xlabel('Iteration','FontSize',14)
% $$$ title('PDE Err','FontSize',14), hold on;
% $$$ set(gca,'FontSize',14)
% $$$ subplot(1,3,3)
% $$$ semilogy(iter,res,clr,'MarkerSize',10); 
% $$$ axis([0 10,1e-11,xup])
% $$$ xlabel('Iteration','FontSize',14)
% $$$ title('Residual','FontSize',14);hold on;
% $$$ set(gca,'FontSize',14)
end
%legend('NumV=1','NumV=2','NumV=3','NumV=4','NumV=10')