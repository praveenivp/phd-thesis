GIRF_sens(3,30,40,200,512*10 );




function [sens_spectrum,freq_axis,cal_grad_nopad]=GIRF_sens(npul,start_time,increment,slew_rate,win_size,Grad_raster_time )

%npul, no of pulses
% start_amplitude, in us
% increment, in us
% slew_rate, Slew rate in  mT/m/ms
% win_size, resolution of spectrum
% Sampling_rate Gradient raster time or smapling time


if nargin<5, win_size=512; end
if nargin<6, Grad_raster_time=10; end



% slew_rate=-1*setting(j,4)/100;


cal_grad=zeros(npul, win_size);
% cal_grad_nopad=zeros(npul, win_size/2);
for i=1:npul
    
    ramp_samples= ceil((start_time + increment*(i-1))/Grad_raster_time);
    
    
    pulse=(0:ramp_samples)*slew_rate/100;
    
    pulse=[pulse  pulse(end)+1*(pulse(end)-pulse(end-1)) fliplr(pulse) ];
    
    pulse1=padarray(pulse,[0 ceil((win_size-length(pulse))/2)],0,'pre');
    
    cal_grad(i,:)=padarray(pulse1,[0 floor((win_size-length(pulse))/2)],0,'post');
    cal_grad_nopad(i,:)=padarray(pulse,[0 floor((win_size-length(pulse)))],0,'post');
end


if(abs(max(cal_grad(:)))>40)
    warning('GRadient amplitude exceeds 40 mT/m : %d',abs(max(cal_grad(:))))
end

%   figure,plot(cal_grad.','LineWidth',1)
freq_axis= linspace(-0.5*(1/Grad_raster_time),0.5*(1/Grad_raster_time),win_size)*1000; %in kHz


sens_spectrum=sum(abs(fftshift(fft(cal_grad,[],2),2)),1);

impulse=zeros(1,win_size);
impulse(win_size/2)=1;
% figure,subplot(211),stem(impulse(end,floor(win_size/2 -win_size/12):ceil(win_size/2+win_size/12)));
% xlabel('Time(us)'),ylabel('Amplitude(mT/m)'),title('Gradient pulses')
% subplot(212),plot(freq_axis,ones(1,length(freq_axis)));
%  xlabel('Frequency(kHz)'),ylabel('Amplitude(AU)'),title('Spectrum'),xlim([-30 30])
% 
% figure,subplot(211),plot(cal_grad(end,floor(win_size/2 -win_size/12):ceil(win_size/2+win_size/12)).','LineWidth',1)
% xlabel('Time(us)'),ylabel('Amplitude(mT/m)'),title('Gradient pulses')
% subplot(212),plot(freq_axis,abs(fftshift(fft(cal_grad(end,:),[],2),2)));
%  xlabel('Frequency(kHz)'),ylabel('Amplitude(AU)'),title('Spectrum'),xlim([-30 30])

figure,
tt=tiledlayout(1,2,"TileSpacing","tight",'Padding','compact');
nexttile(),
grad_plot=cal_grad(:,floor(win_size/2 -20):ceil(win_size/2+20)).';

plot(0:Grad_raster_time:Grad_raster_time*(size(grad_plot,1)-1),grad_plot,'LineWidth',1.5)
grid on,axis square
xlabel('time [\mus]'),ylabel('amplitude [mT/m]'),title('gradient pulses')
% xlim([0 60])

nexttile(),semilogy(freq_axis,abs(fftshift(fft(cal_grad,[],2),2)./ max(abs(sens_spectrum(:)))),'LineWidth',1.5);
 hold on,semilogy(freq_axis,   sens_spectrum./ max(abs(sens_spectrum(:))),'LineStyle','--','LineWidth',1.5)
 xlabel('frequency [kHz]'),ylabel('amplitude [n.u.]'),title('spectrum'),xlim([-1 1]*20)
grid on, axis square,ylim([10^-4,1])
set(gcf,'Color','w','Position', [-1625 485 1059 505])
fontsize(gcf,'scale',1.5)
annotation(gcf,'textbox',...
    [0.0273843248347505 0.924752475247526 0.0292728983703038 0.0920693069306944],...
    'String',{'A'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.487252124645896 0.924752475247527 0.0292728983703033 0.0920693069306938],...
    'String','B',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');




end