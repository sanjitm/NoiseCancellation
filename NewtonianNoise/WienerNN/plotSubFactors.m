function plotSubFactors(sim,NN)

global plcnt

leg = strcat('<',NN.array(1).name,{'> = '},num2str(sim.g_WF_av(1)));
n = length(sim.run);
AN = length(sim.run(1).results.array);

y = zeros(AN,n);
x = zeros(1,n);
for k = 1:n
    x(k) = std(sim.run(k).NN);
    for ai = 1:AN
        y(ai,k) = sim.run(k).results.array(ai).g_WF;
    end
end

figure(plcnt)
set(gcf, 'PaperSize',[8 6])
set(gcf, 'PaperPosition', [0 0 8 6])
clf
plcnt = plcnt+1;
h = plot(1:n,y(1,:),'o','MarkerSize',9);
set(h, 'MarkerFaceColor', get(h, 'Color'),'MarkerEdgeColor', get(h, 'Color'));
hold all
for ai = 2:AN
    h = plot(1:n,y(ai,:),'o','MarkerSize',9);
    set(h, 'MarkerFaceColor', get(h, 'Color'),'MarkerEdgeColor', get(h, 'Color'));
    leg = [leg{:}, strcat('<',NN.array(ai).name,{'> = '},num2str(sim.g_WF_av(ai)))];
end
hold off
legend(leg,'Location','SouthWest')
set(gca,'yscale','log')
set(gca,'xlim',[0.8 n+0.2])
grid
set(gca,'Layer','top')
xlabel('Simulations')
ylabel('1/(1/\gamma^2-1)')

for ai = 1:AN
    figure(plcnt)
    set(gcf, 'PaperSize',[8 6])
    set(gcf, 'PaperPosition', [0 0 8 6])
    clf
    plcnt = plcnt+1;
    loglog(x,y(ai,:),'o','MarkerSize',9,...
        'MarkerEdgeColor','r','MarkerFaceColor','r')
    for k = 1:n
        text(0.95*x(k),1.1*y(ai,k),num2str(k),'FontSize',11); 
    end
    axis([0.8*min(x) 1.2*max(x) 0.8*min(y(ai,:)) 2*max(y(ai,:))])
    grid
    xlabel('NN amplitude')
    ylabel('Filter performance, 1/(1/\gamma^2-1)')
    title(strcat('<',NN.array(ai).name,{'> = '},num2str(sim.p_WF_av(ai),3)))
end
