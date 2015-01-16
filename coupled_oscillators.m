clear;
load osc.txt
load order.txt
[~,NOSC] = size(osc);
NOSC=NOSC-1;

% get frequencies and phases
ws = osc(1,2:end);
phs = osc(2:end,2:end);

% time
t = osc(2:end,1);

% order parameter
r = order(:,2);
psi = order(:,3);

%% plot
figure(1); clf;
subplot(311);
plot(t,r);
subplot(312);
plot(t,sin(psi));
subplot(313);
plot(r.*cos(psi),r.*sin(psi))
axis square

%% raster plot
figure(2); clf;
imagesc(t,1:NOSC,sin(phs)');
xlabel('time')
ylabel('oscillator number')