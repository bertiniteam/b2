function corrector_test()

%% Test parameters
corr_znp1 = [2.3+1i*0.2; 1.1+1i*1.87]; % Point in space at time t=t_n
current_time = 0.9; 
digits(33);  %Precision used
N = 2; %Number of newton iterations in correction step



%% Homotopy system
num_vars = 2;  % number of variables
z = sym('z',[num_vars,1]);
syms t

%%%%%%%%%%%%%%%%%%%% polynomials that make up the homotopy%%%%%%%%%%%%%%%%%%%
H(1) = t*(z(1)^2-1) + (1-t)*(z(1)^2+z(2)^2-4);
H(2) = t*(z(2)-1) + (1-t)*(2*z(1)+5*z(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\frac{dH}{dt}
dHt = diff(H,t);

% Jacobian of H(z,t) w.r.t z
for ii = 1:num_vars
    for jj = 1:num_vars
        JH(ii,jj) = diff(H(ii),z(jj));
    end
end
% Inverse of the Jacobian
JHinv = inv(JH);





%% Predictor Corrector step

% corr_znp1 = the corrected approximation of z_{n+1}


for ii = 1:N
    corr_znp1 = corr_znp1 - vpa(subs(JHinv,[z;t],[corr_znp1;current_time]))*vpa(subs(H,[z;t],[corr_znp1;current_time])).';
	display(corr_znp1);
end


residual = subs(H,[z;t],[corr_znp1;current_time]);

end
