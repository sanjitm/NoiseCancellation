function time = calculateNN(time,params,sgrid,array,seism,plots,cnt)

global plcnt

h = waitbar(0,strcat({'Run '},num2str(cnt),': calculating time series...'),...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

N = time.T*time.fs;
fs = time.fs;
T = time.T;
fc = time.fc;

X = sgrid.X;
Y = sgrid.Y;
R = sgrid.R;
dS = sgrid.dS;

tt = time.tt;

G = params.G;
rho0 = params.rho0;

AN = length(array);

acc_NN = zeros(N,1);
for ai = 1:AN
    time.array(ai).S = time.array(ai).noiseS;
end

for ti = 1:N
    if getappdata(h,'canceling')
        break
    end
    waitbar(ti/N);
    XI = zeros(size(X));
    % on grid
    XI = displace2D(XI,tt(ti),X,Y,seism.waves);
    % at seismometers
    for ai = 1:AN
        SP = array(ai).loc';
        time.array(ai).S(ti,:) = ...
            displace2D(time.array(ai).S(ti,:),tt(ti),SP(1,:),SP(2,:),seism.waves);
    end
    % gravity perturbation
    acc_NN_field = G*rho0*dS*XI./R.^2.*sgrid.cos;
    acc_NN(ti) = sum(sum(acc_NN_field));
    
    if plots.gravityfieldmovie
        figure(99)
        set(gcf, 'PaperSize',[9 6])
        set(gcf, 'PaperPosition', [0 0 9 6])
        clf
        contourf(X,Y,acc_NN_field,20)
        shading flat
        set(gca,'clim',[-3e-16 3e-16])
        xlabel('X [m]')
        ylabel('Y [m]')
        saveas(gcf,strcat('./movieNN/Snap_',num2str(ti),'.jpg'))
    end
    
    if plots.seismicfieldmovie
        figure(99)
        set(gcf, 'PaperSize',[9 6])
        set(gcf, 'PaperPosition', [0 0 9 6])
        clf
        contourf(X,Y,XI,20)
        shading flat
        set(gca,'clim',[-3e-7 3e-7])
        xlabel('X [m]')
        ylabel('Y [m]')
        saveas(gcf,strcat('./moviexi/Snap_',num2str(ti),'.jpg'))
    end
    
end
delete(h);

%evaluate spectrum of first seismometer (supposedly at origin)
spec = fft(time.array(ai).S(:,1))/fs;
time.seismicspectrum = sqrt(abs(spec(1:N/2)).^2/T);


time.NN = acc2disp(acc_NN,T,fs,fc);
time.TM = time.NN+time.noiseTM;

if plots.displacementfield
    figure(plcnt)
    set(gcf, 'PaperSize',[9 6])
    set(gcf, 'PaperPosition', [0 0 9 6])
    clf
    plcnt = plcnt+1;
    contourf(X,Y,XI,50)
    shading flat
    grid
    xlabel('X [m]')
    ylabel('Y [m]')
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String','Surface displacement [m]');
end

