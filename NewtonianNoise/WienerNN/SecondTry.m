%% Seismometer array optimization
% from eq 51 in "Improving the sensitivity of future GW observatories in
% the 1-10Hz band: Newtonian and seismic noise"

cd /Users/jdrigger/LIGO/C1_Jenne/NN_3rd4thGen/SensorArrayOptimization/

load IterationNumber
iterationNumber = iterationNumber + 1;
SaveFolder = '';
SaveFile   = ['ArrayOptimum_IterationNum_' num2str(iterationNumber)];
save IterationNumber iterationNumber

% Picking some stuff
N = 20;  % Number of sensors to use in array.

Max_u = 3;  % Max distance for seismometers from test mass, in units of correlation length
Max_theta = 2*pi;  % In case we want to constrain sensors to be within a certain opening angle.

sigma = 0;  % Instrument noise of sensors. Assumed same for all.
Gamma = 1;  % sensor-to-ground correlation.  Assumed same for all.  If sigma=0, Gamma not used, just can't be 0 b/c (sigma/Gamma)^2

%% do one-by-one optimization

u_min     = zeros(1,N);
theta_min = zeros(1,N);

for ii = 1:N
    if ii == 1
        u_ii = 0;
        theta_ii = 0;
    else
        u_ii     = Max_u     * rand(1);
        theta_ii = Max_theta * rand(1);
    end
    
    SearchableParameterGuesses = [u_ii, theta_ii];
    
    %% The actual optimization

    FixedParameters(1) = N;
    FixedParameters(2) = sigma;
    FixedParameters(3) = Gamma;

    fminopts                     = optimset('TolFun',1e-5,'MaxIter',100000,'Display','off');
    [resout, err_val, exitflag]  = fminsearch(@(x) minimize_one_minus_epsilon_step(x, FixedParameters, u_min, theta_min), SearchableParameterGuesses, fminopts);

    u_min(ii)     = resout(1);
    theta_min(ii) = resout(2);
    
end

u_minimized     = u_min;
theta_minimized = theta_min;

OneMinusEpsilon = one_minus_epsilon(u_minimized, theta_minimized,FixedParameters);
Epsilon         = 1 - OneMinusEpsilon;


save([SaveFolder SaveFile '.mat'], 'u_minimized', 'theta_minimized', 'OneMinusEpsilon', 'Epsilon', 'FixedParameters')
%% This is how you can make the balloon plot:

xx = u_minimized.*cos(theta_minimized);
yy = u_minimized.*sin(theta_minimized);

figure(104)
clf
colormap([.4 .6 .1])
[X,Y,Z] = sphere(100);
surf(0.1*X+xx(1),0.1*Y+yy(1),0.1*Z)
hold on
for k = 2:N
    surf(0.1*X+xx(k),0.1*Y+yy(k),0.1*Z)
end
hold off
shading interp
grid
daspect([1 1 1])
view([1 1 1])
grid
light('Position',[0 1 1])
xlabel('X: Along Beam Axis')
ylabel('Y: Perpendicular to Beam')


















