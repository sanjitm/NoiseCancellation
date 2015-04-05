function sim = simProcessing(sim,time)

T = time.T;
fs = time.fs;
N = T*fs;
n = length(sim.run);
AN = length(sim.run(1).results.array);

sim.sub_ff = zeros(N/2,AN);
for ai = 1:AN
    y = zeros(1,n);
    res = zeros(1,n);    
    sub_ff  = zeros(N,n);
    for k = 1:n
        y(k) = sim.run(k).results.array(ai).g_WF;
        res(k) = std(sim.run(k).results.array(ai).residual);
        sub_ff(:,k) = fft(sim.run(k).results.array(ai).residual)/fs;
    end
    sim.sub_ff(:,ai) = sqrt(mean(abs(sub_ff(1:end/2,:)).^2,2)/T);
    sim.g_WF_av(ai) = prod(y)^(1/n);
    sim.p_WF_av(ai) = 1/sqrt(1/sim.g_WF_av(ai)+1);
    sim.res_av(ai) = mean(res);
end

NN_ff = zeros(N,n);
for k = 1:n
    NN_ff(:,k) = fft(sim.run(k).NN)/fs;
end
sim.NN_ff = sqrt(mean(abs(NN_ff(1:end/2,:)).^2,2)/T);