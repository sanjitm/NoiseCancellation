function results = calculateMap(results,time,array,waves,plots,cnt)

global plcnt

N = time.T*time.fs;
T = time.T;
fs = time.fs;

tt = linspace(0,T,N)';
window = nuttallwin(N);

WN = length(waves);

nk = 101;
k0 = 1;
kx = linspace(-k0,k0,nk);
ky = linspace(-k0,k0,nk);
[kX, kY] = meshgrid(kx,ky);

ff = 10;

i = sqrt(-1);

AN = length(array);

for ai = 1:AN
    SP = array(ai).loc;
    SN = length(SP(:,1));
    results.array(ai).map_fk = zeros(length(ff),nk,nk);

    for fi = 1:length(ff)
        a = sum(time.array(ai).S.*window(:,ones(1,SN)).*exp(-2*pi*i*tt(:,ones(1,SN))*ff(fi)),1)/fs;
        csd = zeros(SN,SN);
        map_fk = zeros(nk,nk);
        for i1 = 1:SN
            for i2 = 1:SN
                csd(i1,i2) = a(i1).*conj(a(i2));
            end
        end

        for c1 = 1:SN
            for c2 = 1:SN
                 map_fk = map_fk + ...
                    +squeeze(csd(c1,c2)/sqrt(abs(csd(c1,c1)*csd(c2,c2))))...
                    *exp(-i*((SP(c2,1)-SP(c1,1))*kX+(SP(c2,2)-SP(c1,2))*kY));
            end
        end
        results.array(ai).map_fk(fi,:,:) = abs(map_fk);
    end

    if plots.arrayresponse
        figure(plcnt)
        set(gcf, 'PaperSize',[9 6])
        set(gcf, 'PaperPosition', [0 0 9 6])
        clf
        plcnt = plcnt+1;
        C = zeros(nk,nk);
        for k1 = 1:nk
            for k2 = 1:nk
                kv = [kx(k1); ky(k2)];
                C(k1,k2) = abs(mean(exp(i*SP*kv),1));
            end
        end
        contourf(kx,ky,C',30)
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String','Array Response'); 
        shading flat
        grid
        set(gca,'DataAspectRatio',[1 1 1])
        xlabel('Wavenumber k_x [1/m]')
        ylabel('Wavenumber k_y [1/m]')
        title(array(ai).name)
    end

    if plots.spatialspectra
        max_v = 0;
        for wi = 1:WN
            for wwi = 1:waves(wi).num
                thisone = abs(waves(wi).xi0(wwi));
                if thisone>max_v
                    max_v = thisone;
                end
            end
        end
        % spatial spectra
        for fi = 1:length(ff)
            figure(plcnt)
            set(gcf, 'PaperSize',[8 6])
            set(gcf, 'PaperPosition', [0 0 8 6])
            clf
            plcnt = plcnt+1;
            contourf(kx,ky,squeeze(results.array(ai).map_fk(fi,:,:)),30)
            hold on
            plot(2*pi*ff(fi)/1000*cos(linspace(0,2*pi,200)),...
                2*pi*ff(fi)/1000*sin(linspace(0,2*pi,200)),'w',...
                2*pi*ff(fi)/500*cos(linspace(0,2*pi,200)),...
                2*pi*ff(fi)/500*sin(linspace(0,2*pi,200)),'w',...
                2*pi*ff(fi)/250*cos(linspace(0,2*pi,200)),...
                2*pi*ff(fi)/250*sin(linspace(0,2*pi,200)),'w',...
                2*pi*ff(fi)/100*cos(linspace(0,2*pi,200)),...
                2*pi*ff(fi)/100*sin(linspace(0,2*pi,200)),'w')
            for wi = 1:WN
                if strcmp(waves(wi).type,'plane') && (waves(wi).freq==ff(fi))
                    for wwi = 1:waves(wi).num
                        plot(waves(wi).k(1,wwi),waves(wi).k(2,wwi),'w.',...
                            'MarkerSize',max(1,floor(25*abs(waves(wi).xi0(wwi))/max_v)));
                    end
                end
                if strcmp(waves(wi).type,'wavelet')
                    for wwi = 1:waves(wi).num
                        plot(waves(wi).k(1,wwi),waves(wi).k(2,wwi),'w.',...
                            'MarkerSize',max(1,floor(25*abs(waves(wi).xi0(wwi))/max_v)));
                    end
                end
                if strcmp(waves(wi).type,'spherical') && (waves(wi).freq==ff(fi))
                    for wwi = 1:waves(wi).num
                        plot(2*pi*ff(fi)/100*cos(waves(wi).az(wwi)),...
                            2*pi*ff(fi)/100*sin(waves(wi).az(wwi)),'w.',...
                            'MarkerSize',max(1,floor(25*abs(waves(wi).xi0(wwi))/max_v)));
                    end
                end
            end
            set(gca,'layer','top')
            hold off
            shading flat
            grid
            set(gca,'DataAspectRatio',[1 1 1])
            xlabel('Wavenumber k_x [1/m]')
            ylabel('Wavenumber k_y [1/m]')
            title(strcat(array(ai).name,{', f = '},num2str(ff(fi),2),'Hz'))
            % saveas(gcf,strcat('./plots/Map_',num2str(cnt),'.jpg'))
        end
    end
end

