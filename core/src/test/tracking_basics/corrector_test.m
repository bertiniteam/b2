function corrector_test()

%% Test parameters
%corr_znp1 = [2.3+1i*0.2; 1.1+1i*1.87]; % Point in space at time t=t_n
%corr_znp1 = [vpa(329.200131365025965858932555654814) - vpa(0.0000000000000000627817514906492713836728689129202)*1i]
corr_znp1 = [vpa(256185069753.408853236449454838735) - vpa(387520022558.051912233172695037026)*1i,
 vpa(-0.0212298348984663761753389579320967) - vpa(0.177814646531698303094367770239999)*1i];

current_time = vpa(.1); 
digits(16);  %Precision used
N = 10; %Number of newton iterations in correction step



%% Homotopy system
num_vars = 2;  % number of variables
z = sym('z',[num_vars,1]);
syms t

%%%%%%%%%%%%%%%%%%%% polynomials that make up the homotopy%%%%%%%%%%%%%%%%%%%
H(1) = vpa((29/16)*z(1)^3 - 2*z(1)*z(2)) + t;
H(2) = z(2) - z(1)^2;
% H(1) = t*(z(1)^2-1) + (1-t)*(z(1)^2+z(2)^2-4);
% H(2) = t*(z(2)-1) + (1-t)*(2*z(1)+5*z(2));
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
    residual = subs(H,[z;t],[corr_znp1;current_time])

end


residual = subs(H,[z;t],[corr_znp1;current_time]);

end
