function results = calculateCoh(time,sgrid,array,plots)

global plcnt

TM = time.TM;

L = sgrid.L;

AN = length(array);
for ai = 1:AN
    S = time.array(ai).S;
    SN = length(S(1,:));
    SP = array(ai).loc;
    results.array(ai).coh_TM_S = ...
        sum(TM(:,ones(1,SN)).*S,1)./sqrt(sum(TM(:,ones(1,SN)).^2,1).*sum(S.^2,1));
    
    if plots.arraycoherence
        figure(plcnt)
        set(gcf, 'PaperSize',[9 6])
        set(gcf, 'PaperPosition', [0 0 9 6])
        clf
        plcnt = plcnt+1;
        scatter(SP(:,1),SP(:,2),50,results.array(ai).coh_TM_S,'filled') 
        axis([-L L -L L])
        set(gca,'clim',0.1*[-1 1])
        xlabel('X [m]')
        ylabel('Y [m]')
        grid
        title(array(ai).name)
        set(gca,'DataAspectRatio',[1 1 1])
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String','Coherence');
    end
end