function [phi_res residual1 residual2]=phase_residual(phi1,phi2)
if phi1>phi2
    factor = -2*pi;
else 
    factor = 2*pi;  
end
% phi_res = min(abs(phi1-phi2),abs(phi1+factor-phi2));
residual1 = abs(phi1-phi2);
residual2 = abs(phi1+factor-phi2);
phi_res = min(residual1,residual2);

end