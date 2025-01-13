f_fiducial_running = f_fiducial;
fd_fiducial_running = fd_fiducial;
z = zs;
pn = [0];
for n=1:numel(z)
%   if draw
%       pcolor(fran,fdran,gamma(:,:,n)),shading interp
%       title(num2str(n))
%       drawnow
%   end
    f_fiducial_running = f_fiducial_running - fd_fiducial_running*z(n) + 0.5*fdd_fiducial*z(n)^2;
    fd_fiducial_running = fd_fiducial_running - fdd_fiducial*z(n);
    phase_fiducial = (z(n)*f_fiducial_running + 0.5*fd_fiducial_running*z(n)^2 + 1/6*fdd_fiducial*z(n)^3);


    phase = (z(n)*freqs(path(n, 2))+1/2*z(n).^2*fdots(path(n, 1))) + phase_fiducial;
    pn(n+1) = pn(n) + round(phase);
end
dlmwrite('pn.dat', pn, 'delimiter', '\n', 'precision', '%d');