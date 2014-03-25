%  Routine to plot PDE err, ODE err, Res for heat
clear
Niter=10
figure(1);clf

Nx=128
Nproc=1
for M = 1:4
  switch M
   case 1,
      fbase=['serial_heat1_'];
   case 2,
      fbase=['serial_heat2_'];
   case 3,
      fbase=['serial_heat3_'];
   case 4,
      fbase=['serial_heat4_'];
  end
  for k = 1:13
    Nstep=2^(k-1);
    dt(k) = 1/Nstep;
    Nt(k)=Nstep;
    fspec=['Niter',num2str(Niter,'%02d'),'_Nx',num2str(Nx,'%03d'),'_Nstep',num2str(Nstep,'%03d'),'_Nproc',num2str(Nproc,'%03d'),'_',num2str(Nstep,'%03d')];
    fname=['../Dat/',fbase,fspec,'.m']
    q=load(fname);
    q_end_ind = find(q(:,1)==3 );
    q128=q(q_end_ind,:);
    
    pde_err(k) = q128(end,6);
    ode_err(k) = q128(end,7);
    res(k) = q128(end,8);
  end

  subplot(1,4,M)
% $$$   loglog(Nt,ode_err,'-ro'); hold on;
% $$$   loglog(Nt,pde_err,'-b*');
% $$$   loglog(Nt,res,'-kx'); 
  semilogy(log2(Nt),ode_err,'-ro', 'MarkerSize',10); hold on;
  semilogy(log2(Nt),pde_err,'-b*', 'MarkerSize',10);
  semilogy(log2(Nt),res,'-kx', 'MarkerSize',10); 
  axis([0 12,1e-18,1])
  legend('ODE Err','PDE Err','Residual','Location','NorthEast')
  xlabel('Log_2(Nsteps)','FontSize',14)
  switch M
   case 1,
      title('Order 1','FontSize',14);
   case 2,
      title('Order 2','FontSize',14);
   case 3,
      title('Order 4','FontSize',14);
   case 4,
      title('Order 8','FontSize',14);
  end
  set(gca,'FontSize',14)

end