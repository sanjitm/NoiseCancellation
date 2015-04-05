% This program calculates the coherence of ground displacement at two 
% points at positions x1 and x2. x2 is varied while x1 stays constant. 
% 
% A number of waves with random amplitudes, phases and directions are added
% to give a total displacement at two points. Each wave in the simulation 
% is given by its quadrature form: 
% 
%   C = A*cos(omega*x/c)+B*sin(omega*x/c) 
% 
% where c is the phase velocity along the direction of the separation vector
% between the two points which depends on the propagation direction of the 
% wave. 

set(0,'DefaultAxesFontSize',20);
set(0,'DefaultTextFontSize',20);

f = 10; %frequency of wave
speed = 3000; %speed of seismic wave

x1 = 4000; %position of first point of correlation pair

nspec = 100; %number of averages 
nwaves = 100; %number of waves
ndist = 100; %number of distance values

corr = zeros(1,ndist);
psd1 = zeros(1,ndist);
psd2 = zeros(1,ndist);
distances = linspace(0,speed/f*30,ndist);

%phase offset of waves
p0 = 2*pi*rand(1,nwaves);    
    
for k = 1:nspec   
    %Gaussian distribution of amplitudes
    A = 1+0.001*randn(1,nwaves);
    B = 1+0.001*randn(1,nwaves);

    %distribution of wave propagation directions; phi is the angle between
    %the direction of propagation and the separation vector of the two
    %points
    frac = 1;
    phi = -pi/frac+2*pi/frac*rand(1,nwaves);

    %The phase velocity depending on direction of propagation
    %relative to the vector x1 -> x2
    cx = speed./cos(phi);

    %convert quadrature representation of wave into amplitude, phase
    %representation and then make use of algebraic sum-of-sines formula
    Amp = sqrt(A.^2+B.^2);
    p1 = p0+atan2(A,B)+2*pi*f*x1./cx; %phase of wave at position x1
        
    for j = 1:length(distances)
        x2 = distances(j);
        
        p2 = p0+atan2(A,B)+2*pi*f*x2./cx; %phase of wave at position x2

        %these two formulas calculate the complex amplitudes of the
        %oscillation at x1 and x2 by summing complex amplitudes of all
        %waves
        S1 = sum(Amp.*exp(i*p1));
        S2 = sum(Amp.*exp(i*p2));

        
        corr(j) = corr(j) + S1*conj(S2)/nspec;
        psd1(j) = psd1(j) + abs(S1)^2/nspec;
        psd2(j) = psd2(j) + abs(S2)^2/nspec;
    end
end

coherence = abs(corr)./sqrt(psd1.*psd2);

plot(distances,coherence)
grid
axis tight
xlabel('Position of P2 [m]')
ylabel('Coherence')
set(gca,'YLim',[0 1])




