//This file is part of Bertini 2.
//
//corrector_test.m is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//corrector_test.m is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with corrector_test.m.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// jeb collins, west texas A&M
// dani brake, university of wisconsin eau claire




function corrector_test()

%% Test parameters
%corr_znp1 = [2.3+1i*0.2; 1.1+1i*1.87]; % Point in space at time t=t_n
%corr_znp1 = [vpa(329.200131365025965858932555654814) - vpa(0.0000000000000000627817514906492713836728689129202)*1i]
corr_znp1 = [ 0.000000000000001,
 0.000000000000001];
%                 vpa(0)+vpa(0)*1i];

current_time = vpa(.9); 
digits(33);  %Precision used
N = 1; %Number of newton iterations in correction step



%% Homotopy system
num_vars = 2;  % number of variables
z = sym('z',[num_vars,1]);
syms t

%%%%%%%%%%%%%%%%%%%% polynomials that make up the homotopy%%%%%%%%%%%%%%%%%%%
H(1) = 0 + vpa(1e-15)*z(1);
H(2) = 0 + vpa(1e-15)*z(2);

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
