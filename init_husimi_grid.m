%==========================================================================
%  Function to set up the Husimi grid
%   Input 
%       * D (integer): Grid spacing density 
%       * Qmax,Pmax (reals): Symmetric Limits of the q and p axis
%   Output
%       * q_out,p_out (1xN array): position and momentum intervals
%       * zmesh_out (NxN arrays): meshgrids
%       * dz_out (real): measure
%==========================================================================

function [q_out,p_out,zmesh_out,dz_out]=init_husimi_grid(Qmax,Pmax,D)

    q_out=linspace(-Qmax,Qmax,D);
    p_out=linspace(-Pmax,Pmax,D);
    dq=abs(q_out(1)-q_out(2));
    dp=abs(p_out(1)-p_out(2));
    dz_out=dq*dp/(2*pi);
    [qmesh,pmesh]=meshgrid(q_out,p_out);
    zmesh_out=sqrt(0.5)*(qmesh+1i*pmesh);

end