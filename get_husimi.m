%==========================================================================
% Function to calculate the Husimi distribution of an eigenstate in the 
% harmonic oscillator basis. The function can also take in a set of
%  M eigenstates and return the husimi distribution of each individual
% eigenstate as a seperate dimensnion of an NxNxM array.
%  
%   Input 
%       * N (integer): Dimension of the underlying phase space grid
%       * phin (NxN):  array of quantum states
%       * n_efn (integer): Number of states in phin to calculate husimi for 
%       * z (NxN array): phase space coordinates z=q+ip
%   Output
%       * Hus_out (NxNxn_efn array): Husimi distribution of each eigenstate
%==========================================================================

function [Hus_out]=get_husimi(N,phin,n_efn,z)

% reverseStr = ''; % String for counter
% display(' ') % Formating in terminal

exa=exp(-0.5*abs(z).^2); % Preallocate for speed
Hus_out=zeros(N,N,length(phin));

for n=1:n_efn % The state |phi_n>

%   msg = sprintf('Constructing Husimi %d/%d', n, n_efn); % Tell you how many are done
%   fprintf([reverseStr, msg]);
%   reverseStr = repmat(sprintf('\b'), 1, length(msg));
   
for j=1:length(phin) % Calculate overlap <z|phi>
    k=j-1;
    Hus_out(:,:,n)=Hus_out(:,:,n)+exa.*conj(z).^k/sqrt(factorial(k))*phin(j,n);
end

Hus_out(:,:,n)=abs(Hus_out(:,:,n)).^2; % Husimi |<z|phi>|^2

end
% display(' ') 
end