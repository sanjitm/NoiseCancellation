function results = calculateWiener(results,time,array,hTM,plots)

global plcnt

tt = time.tt;
TM = time.TM;
NN = time.NN;

AN = length(array);

for ai = 1:AN
    S = prefilter(time.array(ai).S,time.T,time.fs,hTM);

    %covariance between all sensors and test mass
    cov_s = cov_ff([S TM],time.T,time.fs,time.noisemodel); 
    results.array(ai).cov_s = cov_s;
    
    results.array(ai).cWF = results.array(ai).cov_s(...
        1:(end-1),1:(end-1))\cov_s(1:(end-1),end); %filter coefficients
    
    results.array(ai).WF = (results.array(ai).cWF'*S')';
    results.array(ai).residual = results.array(ai).WF-TM;

    results.array(ai).g_WF = 1/(1/(sum(NN.*results.array(ai).WF)^2 ...
        /(sum(NN.^2)*sum(results.array(ai).WF.^2)))-1);
    results.array(ai).g_S0 = 1/(1/(sum(NN.*S(:,1))^2/(sum(NN.^2)*sum(S(:,1).^2)))-1);

    if plots.timeseries
        figure(plcnt)
        set(gcf, 'PaperSize',[8 6])
        set(gcf, 'PaperPosition', [0 0 8 6])
        clf
        plcnt = plcnt+1;
        [AX, H1, H2] = plotyy(tt,TM,tt,time.array(ai).S(:,1));
        hold all
        plot(tt,results.array(ai).WF,'r','LineWidth',2)
        hold off
        set(H1,'LineWidth',2)
        set(H2,'LineStyle','--','LineWidth',2)
        grid
        xlabel('Time [sec]')
        set(AX(1),'YColor','k') 
        set(AX(2),'YColor','k') 
        set(get(AX(1),'Ylabel'),'String','TM disp [m]') 
        set(get(AX(2),'Ylabel'),'String','Seismic displ at TM [m]') 
        legend('Test mass','Wiener filtered array data','Seismometer at test mass')
        title(array(ai).name)
    end
end