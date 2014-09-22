%  Routine to plot error versus iteration
clear
Niter=10
figure(1);clf
figure(2);clf
for k = 5:7
Nx=2^k
Nstep=Nx
Nproc=Nx
fbase=['pfasst_heat_'];
fspec=['Niter',num2str(Niter,'%02d'),'_Nx',num2str(Nx,'%03d'),'_Nstep',num2str(Nstep,'%03d'),'_Nproc',num2str(Nproc,'%03d'),'_',num2str(Nstep,'%03d')];
fname=['../Dat/',fbase,fspec,'.m']
q=load(fname);
q_end_ind = find(q(:,1)==3 )
q128=q(q_end_ind,:)

iter = q128(:,4);
pde_err = q128(:,6);
ode_err = q128(:,7);
res = q128(:,8);

  switch k
   case 5,
      clr = '-ro'
   case 6,
      clr = '-b*'
   case 7,
      clr = '-kx'
  end

% $$$ figure(2)
% $$$ subplot(1,3,1)
% $$$ semilogy(iter,ode_err,clr,'MarkerSize',10); 
% $$$ title('ODE Err','FontSize',14), hold on;
% $$$ xlabel('Iteration','FontSize',14)
% $$$ subplot(1,3,2)
% $$$ semilogy(iter,pde_err,clr,'MarkerSize',10); 
% $$$ xlabel('Iteration','FontSize',14)
% $$$ title('PDE Err','FontSize',14), hold on;
% $$$ subplot(1,3,3)
% $$$ semilogy(iter,res,clr,'MarkerSize',10); 
% $$$ xlabel('Iteration','FontSize',14)
% $$$ title('Residual','FontSize',14);hold on;
% $$$   set(gca,'FontSize',14)
 end
% $$$ legend('N=32','N=64','N=128')
