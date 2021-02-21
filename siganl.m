signal=zeros(length(oneStimMultiRegResultsAno(2).LFP),1);
On=[100:20:400];
Off=[105:20:405];
for i = 1:length(On)
    signal(On(i):Off(i)) = 10000;
end
plot(signal)
title('Stimulation Signal')
xlabel('Time (ms)','FontSize', 14)
ylabel('Stimulation amplitude (mV)','FontSize', 14)

    