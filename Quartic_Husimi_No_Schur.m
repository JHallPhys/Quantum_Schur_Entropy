clear all
close all

% Model Parameters
N =151; % Basis size 
alpha=1;
n_efn=10;
hbar_eff=0.14;
F=0.3
Qmax=2;
Pmax=2;
alpha=1;
gamma=0.01;
alpha=alpha*hbar_eff^2;
F=F/hbar_eff;
Qmax=Qmax/hbar_eff;
Pmax=Pmax/hbar_eff;

NT=100;

D=N;
beta=0;
% Calculate Husimi distribution etc
[a,Q,P]=init_number_basis(N,1,1,1); % Get the operators in the number basis
H0=0.5*(P^2)+0.25*alpha*(Q^4)-1i*gamma*P^2;
U=get_umatrix_number_basis(H0,NT,F,Q,N);
[phin,En] = eig(U); % psi are the eigenfns and En matrix of eigs
Es=diag(En);
[~,ind]=sort(real(log(Es)),'descend');
Es=Es(ind);
phin=phin(:,ind);
figure(1)
plot(imag(log(Es)),real(log(Es)),'k.')

figure
plot(1:1:N,real(log(Es)),'k.')
% return
[q,p,z,dz]=init_husimi_grid(Qmax,Pmax,D); % Initialise husimi grid
Hus=get_husimi(D,phin,n_efn,z); % get the Husimi distribution

% Plot things

figure(1) % Eigenvalues
clf
plot(En,'k.')

Hus_av=zeros(D,D);

for n = 1:n_efn
Hus_av=Hus_av+Hus(:,:,n);

figure % Husimi
clf
imagesc(q,p,Hus(:,:,n))
set(gca,'YDir','normal')
colorbar
colormap(viridis)

end


figure % Average husimi
clf
imagesc(q,p,Hus_av) % again scaling of the q p 
set(gca,'YDir','normal')
colorbar
colormap(viridis)
