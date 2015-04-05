function model = addWaveModes(model)

%determine S and Rayleigh wave speed
nu = model.Poisson;

model.c_S = model.c_P*sqrt((1-2*nu)/(2-2*nu));
ray_R = @(n)n^6-8*n^4+8*n^2*(2-nu)/(1-nu)-8/(1-nu);
nu_R = fsolve(ray_R,0.92);
model.c_R = nu_R*model.c_S;
model.lambda_P = model.c_P/model.freq;
model.lambda_S = model.c_S/model.freq;
model.lambda_R = model.c_R/model.freq;