%
%
%	Simulation of bSSFP using EPG.  Yes, you can do this!!
%   supporting fucntions from EPG package from Brian Hargreaves: https://web.stanford.edu/~bah/software/epg/
%
%%
clearvars,clc
Ntr = 80;	% Number of TRs
Nf = 144;
Nstates = 20;	% Number of states to simulate (easier to keep constant
		% if showing state evolution.


TR = 10e-3;	% s
T1 = 2000e-3;	% s
T2 = 35e-3;	% s
phaseinc=180;	% Degreees

E1=exp(-TR/T1);E2=exp(-TR/T2);
flipang = rad2deg(acos((E1-E2)/(1-E1*E2)));	% Degrees

P = zeros(3,Nstates);	% State matrix
P(3,1)=1;		% Equilibrium magnetization.

Pstore = zeros(2*Nstates,Ntr);	% Matrix to store for image of evolution.
Zstore = zeros(Nstates,Ntr);	% Matrix to store for image of evolution.

% -- excitation

flipphase = 0;			% Start with flip angle along x (arbitrary)
flipphaseincr = 0;		% Phase increment 0 to start.
flipphaseincrincr = 0;	% Phase increment-increment of 117 deg.


for k = 1:Ntr

  flipphase = flipphase + pi/180*(flipphaseincr);	% Incr. phase
  flipphaseincr = flipphaseincr + flipphaseincrincr;	% Increment increment!

  fp(k)=flipphase;
  if(Ntr>-1) %alpha/2 pulse
  P = epg_rf(P,pi/180*flipang,flipphase);		% RF pulse
  else
        P = epg_rf(P,pi/180*flipang*1/2,flipphase);		% RF pulse
  end


  s(k) = P(1,1);      			% Signal is F0 state.
  sd(k) = s(k)*exp(-i*flipphase);	% Phase-Demodulated signal.

  temp=epg_FZ2spins(P);	% There are many redundant (0) states here that could	
  s_all(k)=complex(temp(1,end),temp(2,end))*size(temp,2);
  plotMax=5;
  maxF=min(k,plotMax-1);
  
sd1(plotMax+(-maxF+1:maxF-1),k) = [fliplr(conj(P(2,2:maxF))) P(1,1:maxF)]; 		% Stretched out


  % -- Keep states for evolution diagram
  Pstore(Nstates:2*Nstates-1,k) = P(2,:).'; % Put in negative states
  Pstore(1:Nstates,k) = flipud(P(1,:).');  % Put in positive, overwrite center.
  Zstore(:,k) = P(3,:).';

  if(Ntr>-1)
  % -- Simulate relaxation and spoiler gradient
  P = epg_grelax(P,T1,T2,TR,0,0,1,1);   % spoiler gradient, relaxation.
  else
        P = epg_grelax(P,T1,T2,TR/2,0,0,1,1);   % spoiler gradient, relaxation.
  end
  


end;

figure(1);
plotstate = cat(1,Pstore,Zstore);
dispim(plotstate,0,0.2);
xlabel('Time'); ylabel('F(top) and Z(bottom) States');
colormap jet,clim([0 20])
%%
figure(2),clf
colororder(lines(6))
% magphase(sd);
sd1(sd1==0)=nan;
magphase(sd1(2:3,:)');
magphase(s_all);

legend(strsplit(sprintf('F_{%d} ',-plotMax:plotMax),' '),'NumColumns',2)

figure(3),clf
Np = 200;
if(flipphaseincrincr~=0)
for q=1%k-Nf:k
  Fp = [flipud(Pstore(1:Nstates,q)).'];
  Fn = [Pstore(Nstates:2*Nstates-1,q).']; Fn(1)=conj(Fn(1));
  Zn = [Zstore(:,q).'];
  M = epg_FZ2spins([Fp;Fn;Zn],Np)*Np;
  S = exp(-i*fp(q))*(M(1,:)+i*M(2,:));
  magphase(S);
  subplot(2,1,1);
  hold on;
%   h = plot([1,200],[1 1]*exrecsignal(T1,T2,0,TR,flipang),'y--');
  plot([1,200],[1 1]*abs(mean(S)),'w--');
  hold off;
%   setprops
  subplot(2,1,1);

  tt = sprintf('Signal Across Voxel, %g^\\circ Increment',phaseinc);
  h = title(tt);
  set(h,'FontSize',24);
  set(h,'Color',[1 1 1]*.8);
  drawnow;


  f = getframe(gcf);
  [fim,map] = frame2im(f);
  fn = sprintf('%s_%04d.tif','/Users/brian/tmp/ims/spgr',q+Nf-Ntr);
%   imwrite(fim,fn,'TIFF');

end;
end


%%
figure(5),clf

taxis=permute(linspace(0,TR,100),[1 3 2]);
T2star=20e-3;

sd_decay=cat(1,sd1(1:4,:).*flip(exp(-taxis./T2star),3),...
    sd1(plotMax:end,:).*exp(-taxis./T2star));

sd_decay=permute(sd_decay,[1 3 2]);

sig_decay=s_all.*exp(-taxis./T2star);
sig_decay=permute(sig_decay,[1 3 2]);

col=@(x) x(:);


cl=lines(4);
cl=cl([1 2 4],:);
for cM=plotMax+(-2:2)
    if(cM-plotMax>=0)
        format_str={'color',cl(abs(cM-plotMax)+1,:),'LineStyle','-','LineWidth',1.5};
    else
        format_str={'color',cl(abs(cM-plotMax)+1,:),'LineStyle','--','LineWidth',1.5};
    end

subplot(3,2,1:2),hold on
taxis=linspace(0,TR*10,100*10);
plot(taxis*1e2,col(abs(sd_decay(cM,:,1:10))),format_str{:})


subplot(3,2,[5,6]),hold on
taxis=linspace(0,TR*10,100*10)+70*TR;
plot(taxis*1e2,col(abs(sd_decay(cM,:,end-9:end))),format_str{:})
ylim([0,0.3])
grid on
end

subplot(3,2,1:2),hold on
taxis=linspace(0,TR*10,100*10);
plot(taxis*1e2,col(abs(sig_decay(1,:,1:10))),'LineWidth',4)
rf_pulses=0.3*(mod(taxis*1e3,TR*1e3)<1e-1);
rf_pulses(rf_pulses==0)=nan;
stem(taxis*1e2,rf_pulses,'Marker','none','LineWidth',5)
ylim([-1e-3,0.3])
grid on,grid minor
ylabel('M_0 [normalized]'),xlabel('TR [#]')

subplot(3,2,5:6),hold on
taxis=linspace(0,TR*10,100*10)+70*TR;
plot(taxis*1e2,col(abs(sig_decay(1,:,end-9:end))),'LineWidth',4)
stem(taxis*1e2,rf_pulses,'Marker','none','LineWidth',5)
ylim([-1e-3,0.2])
grid on,grid minor


legend({'F_{-2}','F_{-1}','F_{0}','F_{1}','F_{2}','Signal'},'NumColumns',2)
xlabel('TR [#]')
ylabel('M_0 [normalized]')


subplot(3,2,3)
imagesc(1:80,-20:20,abs(Pstore)),clim([0,0.10])
ylabel('F states'),xlabel('TR [#]')

subplot(3,2,4)
imagesc(1:80,0:20,abs(Zstore)),clim([0,0.3])
ylabel('Z states'),xlabel('TR [#]')
sgtitle('Time evolution of EPG states in bSSFP')
fontsize(gcf,'scale',1.5)
set(gcf,'color','w')
axis tight
