function f = rayleighPDF(t,dt,fshift,sigma)
    % Reyleigh distribution for time-out.
    f = (t-fshift).*exp(-(t-fshift).^2./(2.*sigma.^2))./(sigma.^2).*dt;
    f(f(:)<0) = 0;
end