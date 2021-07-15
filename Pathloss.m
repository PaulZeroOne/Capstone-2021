close all;
clear;

%Wavelength (3 GHz)
lambda = 0.1;

%Propagation distances and angles for source and destination
d = 300;
eta = pi/3;
deltaValues = linspace(0.5,50,200);

%Area of each element
A = (lambda/5)^2;

%Number of RIS elements
Nsqrt = 20*5;
N = Nsqrt.^2;

%Side length of each element
a = sqrt(A);

%Set the antenna gain at the source (10 dBi)
antennaGainTx = 10;

%Set the extra propagation loss through the window (20 dB)
penetrationLoss = 100;

%Define location of the source
p_t = [d*sin(eta); 0; d*cos(eta)];


%Prepare to save simulation results
channelGain_RIS = zeros(length(deltaValues),1);
channelGain_RIS_approx = zeros(length(deltaValues),1);
channelGain_mirror = zeros(length(deltaValues),1);


%% Go through the different number elements
for j = 1:length(deltaValues)
    
    %Position of the destination
    p_r = [0; 0; deltaValues(j)];
    
    %Prepare to store channel gains for individual elements/antennas
    betaHn = zeros(N,1);
    betaGn = zeros(N,1);
    phaseHn = zeros(N,1);
    phaseGn = zeros(N,1);
    
    %Go through each element/antenna and compute pathlosses
    for n = 1:N
        
        %Compute location using Eqs. (22)-(23) in [11]
        x = -a*(sqrt(N)-1)/2 + a*mod(n-1,sqrt(N));
        y = a*(sqrt(N)-1)/2 - a*floor((n-1)/sqrt(N));
        
        %Compute channel gain for the n:th element
        betaHn(n) = channelgainGeneral(p_t,[x; y; 0],a);
        betaGn(n) = channelgainGeneral(p_r,[x; y; 0],a);
        
        %Compute phase-shift for the n:th element
        phaseHn(n) = mod(norm(p_t-[x; y; 0])/lambda,1)*2*pi;
        phaseGn(n) = mod(norm(p_r-[x; y; 0])/lambda,1)*2*pi;
        
    end
    
    
    %Compute the total channel gain with the RIS using Eq. (42) in [11],
    %by removing the P/sigma^2 term
    channelGain_RIS(j) = sum(sqrt(betaHn.*betaGn)).^2/penetrationLoss;
    
    %Compute the  total channel gain with the RIS using Eq. (21) in [11]
    %when mimicking a mirror using theta_n=0 and mu_n=1 (by removing the
    %P/sigma^2 term)
    channelGain_RIS_approx(j) = abs(sum(sqrt(betaHn.*betaGn).*exp(-1i*(phaseGn)))).^2/penetrationLoss;
    
    %Compute the mirror limit in Eq. (54) in [11]
    channelGain_mirror(j) = lambda^2/((4*pi)^2*(deltaValues(j)+d)^2)/penetrationLoss;
    
end

%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(deltaValues,pow2db(channelGain_RIS),'b-','LineWidth',2);
xlabel('Distance from surface [m]','Interpreter','Latex');
ylabel('Pathloss [dB]','Interpreter','Latex');
legend({'RIS '},'Interpreter','Latex','Location','NorthEast');
set(gca,'fontsize',18);
function zetaVal = channelgainGeneral(p_t,p_n,a)
% INPUT:
%p_t = Location (x_t,y_t,d) of the transmitter
%p_n = Location (x_n,y_n,0) of the receiver in the XY-plane
%a   = Area of the square-shaped receive antenna
%
% OUTPUT:
% zetaVal = Free-space channel gain


%Compute the X and Y sets defined in Lemma 1
X = [a/2+p_n(1)-p_t(1); a/2-p_n(1)+p_t(1)]';
Y = [a/2+p_n(2)-p_t(2); a/2-p_n(2)+p_t(2)]';

%Compute square distance in Z-dimension
d2 = (p_t(3) - p_n(3))^2;

%Compute the channel gain using Eq. (4)
zetaVal = 0;

for x = X
    
    x2 = x^2;
    
    for y = Y
        
        y2 = y^2;
        
        zetaVal = zetaVal + x*y/d2 / ( 3*(y2/d2+1)*sqrt(x2/d2+y2/d2+1)) + (2/3)*atan( x*y/d2 /  sqrt(x2/d2+y2/d2+1) );
        
    end
    
end

zetaVal = zetaVal/(4*pi);
end
