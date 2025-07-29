function [A,B] = crank_nicholson_sparsematrix(drift,nd,ny,dy,dt)

    D = 0.5;% (half?) the variance of the momentary evidence; not the variance as indicated in kiani2009
    sigma = D*dt/(2*dy^2);

    auxA = nan(nd*ny,3);
    auxB = nan(nd*ny,3);
    for idrift = 1:nd
        a = -1*drift(idrift);
        rho = a*dt/(4*dy);
        inds = [1:ny]+ny*(idrift-1);
        auxA(inds,:) = repmat([-sigma+rho,(1+2*sigma),-(rho+sigma)],ny,1);
        auxB(inds,:) = repmat([sigma-rho,(1-2*sigma),sigma+rho],ny,1);
    end
    A = spdiags(auxA,-1:1,ny*nd,ny*nd);
    B = spdiags(auxB,-1:1,ny*nd,ny*nd);
    
end