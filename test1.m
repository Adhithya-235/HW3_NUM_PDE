clear 
close all
clc

x = linspace(-1,1,101)';
W = zeros(length(x),1);

for i = 1:length(x)
    W(i) = wavetest(x(i));
end

plot(x,W,'k','linewidth',2)
grid on
xlim([-1 1])
ylim([-0.2,1.2])