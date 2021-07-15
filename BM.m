%response on changing n in PLA
fc = 1e9;
lambda = physconst('LightSpeed')/fc;
array = phased.ULA('NumElements',100,'ElementSpacing',lambda/2); %no. of elements
array.Element.FrequencyRange = [8e8 1.2e9]; %frequency range of array

t = linspace(0,0.3,300)';  %sample pulse
testsig = zeros(size(t));
testsig(201:205) = 1;

angle_of_arrival = [30;0]; %elevation, azimuth of pulse and PLA
x = collectPlaneWave(array,testsig,angle_of_arrival,fc);
subplot(321)
plot(t,abs(x(:,1)))
title('Received Signal)')
axis tight
ylabel('Magnitude')

rng default
npower = 0.5;
x = x + sqrt(npower/2)*(randn(size(x)) + 1i*randn(size(x)));
subplot(322)
plot(t,abs(x(:,1)))
title('Element 1 (with Noise)')
axis tight
ylabel('Magnitude')
subplot(323)
plot(t,abs(x(:,2)))
title('Element 2 (with Noise)')
axis tight
ylabel('Magnitude')
subplot(324)
plot(t,abs(x(:,3)))
title('Element 3 (with Noise)')
axis tight
xlabel('Seconds')
ylabel('Magnitude')
subplot(325)
plot(t,abs(x(:,4)))
title('Element 4 (with Noise)')
axis tight
xlabel('Seconds')
ylabel('Magnitude')


beamformer = phased.PhaseShiftBeamformer('SensorArray',array,...
    'OperatingFrequency',1e9,'Direction',angle_of_arrival,...
    'WeightsOutputPort',true);

[y,w] = beamformer(x);

subplot(326)
plot(t,abs(y))
axis tight
title('Received Signal with Beamforming')
ylabel('Magnitude')
xlabel('Seconds')

%azang = -180:30:180;
%subplot(211)
%pattern(array,fc,[-180:180],0,'CoordinateSystem','rectangular',...
    %'Type','powerdb','PropagationSpeed',physconst('LightSpeed'))
%set(gca,'xtick',azang);
%title('Array Response without Beamforming Weights')
%subplot(212)
%pattern(array,fc,[-180:180],0,'CoordinateSystem','rectangular',...
    %'Type','powerdb','PropagationSpeed',physconst('LightSpeed'),...
    %'Weights',w)
%set(gca,'xtick',azang);
%title('Array Response with Beamforming Weights')