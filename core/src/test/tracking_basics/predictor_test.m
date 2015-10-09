function predictor_test()

%% Test parameters
zn = [4.641588833612776e-1; 7.416198487095662e-1]; % Point in space at time t=t_n
tn = 0.7; dt = -.01; tnp1 = tn + dt; %Starting time t_n, and time step dt
digits(30);  %Precision used




%% Homotopy system
num_vars = 2;  % number of variables
z = sym('z',[num_vars,1]);
syms t

%%%%%%%%%%%%%%%%%%%% polynomials that make up the homotopy%%%%%%%%%%%%%%%%%%%

H(1) = t*(z(1)^3-1) + (1-t)*(z(1)^3+2);
H(2) = t*(z(2)^2-1) + (1-t)*(z(2)^2+.5);

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




%\frac{dH}{dt}(z_n, t_n)
z
zn
dHtn = vpa(subs(dHt,[z], [zn])).';
% Inverse of Jacobian at (zn,tn)
JHinvn = vpa(subs(JHinv, [z;t], [zn;tn]));




%% Predictor Corrector step

% pred_znp1 = the prediction step approximation of z_{n+1}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -1 = Constant(no predictor)
% 0 = Forward Euler
% 1 = Heun
% 2 = RK4
% 3 = Heun Prediction with 1st order Euler for error estimation
% 4 = Runge-Kutta-Norsett
% 5 = RK fehlburg
% 6 = RK-Cash-Karp
% 7 = RK-Dormand-Prince
% 8 = RK-Verner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ODEPredictor = 0;
switch(ODEPredictor)
    case 0
        % Forward Euler
        pred_znp1 = zn - dt*JHinvn*dHtn;
end

display(pred_znp1);


end
