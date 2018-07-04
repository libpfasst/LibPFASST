%  Quick and dirty plotter for results 
Niter=10
Nlevel=2
Nproc=32
Nlev=2
Nstep=32
Nblock=Nstep/Nproc

res=zeros(Niter,Nproc);
for k = 0:Nproc-1
   fname=['residual_',num2str(k,'%3.3i'),'.dat'];
   foob=load(fname);
   resk=reshape(foob(:,4),Niter,Nblock,Nlev);
   res(:,k+1)=squeeze(resk(:,end,end));
end
figure(2);clf
for k = 1:Niter
   semilogy([1:Nproc],res(k,:),'-*'); hold on;
end
xlabel('processor')
ylabel('residual')

