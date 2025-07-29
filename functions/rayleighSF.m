function f = rayleighSF(t,fshift,sigma)
    % Reyleigh survivor function for time-out.
    f = exp(-(t - fshift).^2/(2.*sigma.^2));
    pick = logical(t<=fshift);
    f(pick) = 1;
end