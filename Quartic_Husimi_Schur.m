clear all
close all

% Model Parameters
N =101; % Basis size 
n_efn=20;
hbar_eff=0.12;
F=0.3
Qmax=1.5;
Pmax=1.5;
alpha=1;
gamma=0.001;
set_efn='G'; % Invarient Subspace: Gain ('G') Loss ('L')
NT=101; % dt=2*pi/NT in integrator
D=N; % Resolution of the Husimi array
alpha=alpha*hbar_eff^2;
F=F/hbar_eff;
Qmax=Qmax/hbar_eff;
Pmax=Pmax/hbar_eff;
% Calculate Husimi distribution etc
[a,Q,P]=init_number_basis(N,1,1,1); % Get the operators in the number basis
H0=0.5*(P^2)+0.25*alpha*(Q^4)-0.5*1i*gamma*P^2;
U=get_umatrix_number_basis(H0,NT,F,Q,N);
[phin,En] = schur(U); % Get the Schur vectors
[phin,Es]=get_schur_ordered(N,En,phin,set_efn) ; %Reorder the schur vectors
[q,p,z,dz]=init_husimi_grid(Qmax,Pmax,D); % Initialise husimi grid
Hus=get_husimi(D,phin,n_efn,z); % get the Husimi distribution

% Plot things

% figure
% plot(1:1:N,real(log(Es)),'k.')
% return

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
