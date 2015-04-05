function [r, model] = makeGrid(model)

radius = model.radius;
ng = 2*floor(model.ng/2);
model.ng = ng+1;

x_minus = linspace(-radius(1),0,ng(1)/2+1);
x_plus = linspace(0,radius(1),ng(1)/2+1);
x = [x_minus x_plus(2:end)];
y_minus = linspace(-radius(2),0,ng(2)/2+1);
y_plus = linspace(0,radius(2),ng(2)/2+1);
y = [y_minus y_plus(2:end)];
z = linspace(-radius(3),min(radius(3),model.surface),ng(3)+1);

[r(:,:,:,1) r(:,:,:,2) r(:,:,:,3)] = meshgrid(x,y,z);
