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

figure(log2(Nx));clf
subplot(1,3,1)
semilogy(iter,ode_err,'-ro');title('ODEerr'); 
subplot(1,3,2)
semilogy(iter,pde_err,'-bo');title('PDEerr') 
subplot(1,3,3)
semilogy(iter,res,'-ko'); title('Res') 
figure(2)
subplot(1,3,1)
semilogy(iter,ode_err,'-ro'); title('ODEerr'), hold on;
subplot(1,3,2)
semilogy(iter,pde_err,'-bo'); title('PDEerr'), hold on;
subplot(1,3,3)
semilogy(iter,res,'-ko'); title('Res'); hold on;
end