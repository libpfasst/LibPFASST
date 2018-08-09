%  Quick and dirty plotter for results 
load('result-size.dat')
Nstep=result_size(1)
Niter=result_size(2)
Nproc=result_size(3)
Nlevel=result_size(4)


Nblock=Nstep/Nproc


for j = 1:Nblock
    res=zeros(Niter,Nproc);
    for k = 0:Nproc-1
        fname=['dat/residual_',num2str(k,'%3.3i'),'.dat'];
        foob=load(fname);
        resk=reshape(foob(:,4),Niter,Nblock,Nlevel);
        res(:,k+1)=squeeze(resk(:,j,end));
    end
    figure(j);clf
    for k = 1:Niter
        semilogy([1:Nproc],res(k,:)+eps,'-*'); hold on;
    end
    xlabel('processor')
    ylabel('residual')
    
end