function Uz = forward_Uz(pex,G,nx,ny,nz,dz,z,npts,dx,dy)
[k,kx,ky] = wave2d(npts,dx,dy);
pf=zeros(npts,npts,nz);anof=zeros(npts,npts,nz);
anoex=zeros(npts,npts);ano=zeros(ny,nx,nz);Uz=zeros(ny,nx);
ydiff=floor((npts-ny)/2); xdiff=floor((npts-nx)/2); 
for i = 1:nz
    pf(:,:,i)=fftshift(fft2(pex(:,:,i)));
end
for K=1:nz
    h1=z(K)-dz/2; 
    h2=z(K)+dz/2;
    anof(:,:,K)=pf(:,:,K).*((2*pi*G+eps)./(k+eps)).*(exp(-k*h1)-exp(-k*h2)); 
end
for K=1:nz
    anoex(:,:,K)=ifft2(ifftshift(anof(:,:,K)));
    ano(:,:,K)=real(anoex((ydiff+1):(ydiff+ny),(xdiff+1):(xdiff+nx),K))*1e5;
    Uz(:,:)=Uz(:,:)+(ano(:,:,K));
end