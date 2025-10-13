%%
formatplots={'linewidth',2};
actFormat={'Color',[235, 77, 77]./255};
restFormat={'Color',[63, 171, 99]./255};
wmFormat={'Color',[63, 171, 200]./255};

TR=[7]*1e-3;
TE=1e-3;
T1=2000e-3;
T2=35e-3;
shift=0;
flip=12;
phi=(0:1:359);
xin=Phase2Freq(deg2rad(phi),TR(1));
bSSFP=bSSFP_profile(TR,TE,deg2rad(flip),deg2rad(phi),shift,T1,T2);
figure(8),clf,

tt=tiledlayout(2,2,"TileSpacing","compact","Padding","compact");
nexttile()
plot(xin,abs(bSSFP),formatplots{:}),ylabel('amplitude [a.u.]')
yyaxis('right'),plot(xin,unwrap(angle(bSSFP)),'--',formatplots{:}),ylabel('Phase [rad]')
xlabel('frequency [Hz]')
title(sprintf('SSFP frequency response | TR = %d ms',TR*1e3))
legend('Magnitude','Phase','Location','northwest')

%% simulate spetra
% full width at half maximum is 2*Linwist
FWHM_hz=[18,18,18,180].*1;

T2Star=1000./(pi*FWHM_hz);  %ms
freq_Hz=[0, 20, 1380, 0];
amp_rest=[1,0.4,2,10];
amp_act=[1,0.,2,10];
% freq_Hz_wm=[0, 100, 1380, 0];
% amp_wm=[0.1,2,2,10];

% freq_Hz=[0,35,100];
% amp_rest=[0,1,0];
% amp_act=[1 0 0];


faxis=linspace(-2e3,2e3,1024*32);
model = @(freq_axis,amp1,freq,LW)amp1./((pi*LW/2)*(1 + ((freq_axis - freq)./(LW/2)).^2)) ;
spec_rest=zeros(size(faxis));
spec_act=zeros(size(faxis));
% spec_wm=zeros(size(faxis));
for ii=1:length(amp_rest)
spec_rest=spec_rest+model(faxis,amp_rest(ii),freq_Hz(ii),FWHM_hz(ii));
spec_act=spec_act+model(faxis,amp_act(ii),freq_Hz(ii),FWHM_hz(ii));
% spec_wm=spec_wm+model(faxis,amp_wm(ii),freq_Hz_wm(ii),FWHM_hz(ii));
if(ii==2)
spec_act=spec_act./sum(spec_act);
sclfac=max(spec_act);
spec_act=spec_act./sclfac;
spec_rest=spec_rest./sum(spec_rest);
 spec_rest=spec_rest./sclfac;

end
end


nexttile()
plot(faxis,spec_rest,formatplots{:},restFormat{:}),
hold on,plot(faxis,spec_act,formatplots{:},actFormat{:})
% plot(faxis,spec_wm,formatplots{:},wmFormat{:})
xlabel('frequency [Hz]')
ylabel('amplitude [a.u.]')
legend('resting','active')
title('Intravoxel spectrum')

% xlim([-20 20])

%% alias spectrum
faxis_bssfp=linspace(-0.5/TR,0.5/TR,length(phi));
center_freq=0;
max_f=max(faxis_bssfp);

[~,cIdx]=min(abs(faxis-center_freq));
[~,mIdx]=min(abs(faxis(cIdx:end)-max_f));
rIDX= [(cIdx-mIdx):(cIdx+mIdx)];
faxis_bssfp=linspace(-0.5/TR,0.5/TR,length(rIDX));
alias_spec_act=zeros(size(rIDX));
alias_spec_rest=zeros(size(rIDX));
for i= 1:floor(0.5*length(faxis)/length(faxis_bssfp))
    alias_spec_act=alias_spec_act+spec_act(rIDX+(i-1)*length(rIDX));
    alias_spec_act=alias_spec_act+spec_act(rIDX-(i-1)*length(rIDX));


        alias_spec_rest=alias_spec_rest+spec_rest(rIDX+(i-1)*length(rIDX));
    alias_spec_rest=alias_spec_rest+spec_rest(rIDX-(i-1)*length(rIDX));

end
alias_spec_rest=alias_spec_rest.*0.99; % additonl exchange faster relaxation
nexttile(),
plot(faxis_bssfp,alias_spec_rest,formatplots{:},restFormat{:}),hold on,
plot(faxis_bssfp,alias_spec_act,formatplots{:},actFormat{:})
% plot(faxis_bssfp,alias_spec_wm,formatplots{:},wmFormat{:})
xlabel('frequency [Hz]')
ylabel('amplitude [a.u.]')
legend('resting','active')
title('Aliased intravoxel spectrum | TR = 7 ms')




%% interp bSSFP to faxis grid;

xin=Phase2Freq(deg2rad(phi),TR(1));

bSSFP_interp=interp1(xin,bSSFP,faxis_bssfp-min(faxis_bssfp),'linear','extrap');


%% own convolution

conv_prof_act=zeros(size(faxis_bssfp));
conv_prof_rest=zeros(size(faxis_bssfp));
conv_prof_wm=zeros(size(faxis_bssfp));
for i=1:length(faxis_bssfp)
    shift=round(faxis_bssfp(i)/diff(faxis_bssfp(1:2)));
    conv_prof_act=conv_prof_act+alias_spec_act(i).*circshift(bSSFP_interp,-1*shift);
    conv_prof_rest=conv_prof_rest+alias_spec_rest(i).*circshift(bSSFP_interp,-1*shift);

end
nexttile()
plot(faxis_bssfp-min(faxis_bssfp),abs(conv_prof_rest),formatplots{:},restFormat{:})

hold on

plot(faxis_bssfp-min(faxis_bssfp),abs(conv_prof_act),formatplots{:},actFormat{:})
% plot(faxis_bssfp-min(faxis_bssfp),abs(conv_prof_wm),formatplots{:},actFormat{:})
,ylabel('amplitude [a.u.]')

yyaxis('right'),plot(faxis_bssfp-min(faxis_bssfp),angle(conv_prof_rest),'--',formatplots{:},restFormat{:}),ylabel('Phase [rad]')
yyaxis('right'),plot(faxis_bssfp-min(faxis_bssfp),angle(conv_prof_act),'--',formatplots{:},actFormat{:})
xlabel('frequency [Hz]')
title(sprintf(' Resulting SSFP profile| TR = %d ms',TR*1e3))

legend('resting [mag]','active [mag]','resting [phase]','active [phase]','Location','northwest')

fontsize(gcf,'scale',1.5),
set(gcf,'color','w')
%%
% figure(5),clf,
% % plot(xin,bSSFP),hold on,
% plot(faxis,bSSFP_interp)

function freq_Hz=Phase2Freq(phase_rad,TR_s)
% freq_Hz=Phase2Freq(phase_rad,TR_s)
% all these things are equivalent
% -1/TR to +1/TR  <====> -360 to 360 Hz  <===> -2*pi to 2*pi

if(phase_rad<0)
    phase_rad= -1*mod(phase_rad,2*pi);
else
    phase_rad= mod(phase_rad,2*pi);
end

freq_Hz=interp1(linspace(-2*pi,2*pi,1e3),linspace(-1/TR_s,1/TR_s,1e3),phase_rad(:));
end


function bSSFP=bSSFP_profile(TR,TE,flip,phi,shift,T1,T2)
%function bSSFP=bSSFP_profile(TR,TE,flip,phi,T1,T2)
%
%Calculation of the bSSFP profile (series of measurements with varying RF phase increment)
%
%INPUT
%
%repetition time in ms (TR)
%echo time in ms (TE)
%flip angle in rad (flip)
%vector of RF phase increments in rad (phi)
%off-resonance shift in rad (shift)
%spin-lattice relaxation time in ms (T1)
%spin-spin relaxation time in ms (T2)
%
%OUTPUT
%
%complex bSSFP signal for the given range of RF phase increments phi (bSSFP)

M0    = 1;
E1    = exp(-TR./T1);
E2    = exp(-TR./T2);

C     = E2.*(E1-1).*(1+cos(flip));
D     = (1-E1.*cos(flip))-(E1-cos(flip)).*E2.^2;
bSSFP = (1-E1).*sin(flip).*(1-E2.*exp(-complex(0,1)*(-phi+shift)))./(C.*cos(-phi+shift)+D);

%% Read out signal at t=TE

bSSFP = M0*bSSFP.*exp(-TE./T2).*exp(complex(0,-1)*phi*(TE/TR));
% bSSFP = M0*bSSFP.*exp(-TE./T2);
end
