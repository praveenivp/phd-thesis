%% FA
T1=500e-3;
T2=T1;

TR=10e-3;
TE=TR/2;
dfreq=linspace(-1/TR,1/TR,1e4);
FA=deg2rad([10 40 90])
phi=0;
bSSFP=bSSFP_profile_Ganter2(FA(:),T1,T2,TE, TR,phi,dfreq*TR*2*pi);

figure(21)
set(gcf,'Position',[-1643 465 1060 468],'Color','w')
tt=tiledlayout(1,2,"TileSpacing","compact",'Padding','compact');
nexttile(tt)
l_h1=plot(dfreq*TR,abs(squeeze(bSSFP.')),'LineWidth',1.5);
ylim([0 0.52]),ylabel('M_{xy} magnitude [n.u]')
yyaxis('right'),
set(gca,'ColorOrder',lines(3),'LineStyleOrder','--')
l_h2=plot(dfreq*TR,angle(squeeze(bSSFP(1,:).')),'LineWidth',1.5,'LineStyle','--','Color','black'); hold on,
plot(dfreq*TR,angle(squeeze(bSSFP(1,:).')),'LineWidth',1.5,'LineStyle','--','Color',lines(1));

title('flip angle dependence')
xlabel('off-resonance [Hz]'),ylabel('M_{xy} phase [rad]'),axis square,grid minor
xticklabels({'-1/TR','-0.5/TR','0','0.5/TR','1/TR'})
yticks([-1 -0.5 0 0.5 1]*pi),yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
ax=gca;ax.YAxis(2).Color=[0 0 0];
legend([l_h1;l_h2],'10^\circ','40^\circ','90^\circ','phase','Position',[0.30 0.404365086058775 0.0825147338079797 0.186904756724835])


% T2/T1
T2=T1.*[0.1 0.5 1];

dfreq=linspace(-1/TR,1/TR,2000);
get_opt=@(t1,t2,tr) acos((exp(-tr./t1)-exp(-tr./t2))./(1- exp(-tr./t1).*exp(-tr./t2)));
FA=get_opt(T1,T2,TR);
phi=0;
bSSFP=bSSFP_profile_Ganter2(FA(:),T1(:),T2(:),TE, TR,phi,dfreq*TR*2*pi);


nexttile(tt)
l_h1=plot(dfreq*TR,abs(squeeze(bSSFP.')),'LineWidth',1.5);
ylim([0 0.52]),ylabel('M_{xy} magnitude [n.u]')
yyaxis('right'),
ax=gca;ax.YAxis(2).Color=[0 0 0];
l_h2=plot(dfreq*TR,angle(squeeze(bSSFP(1,:).')),'LineWidth',1.5,'LineStyle','--','Color','black'); hold on,
plot(dfreq*TR,angle(squeeze(bSSFP(:,:).')),'LineWidth',1.5,'LineStyle','--');
title('T2/T1 dependence')
set(gca,'ColorOrder',lines(3),'LineStyleOrder','--')
xlabel('off-resonance [Hz]'),ylabel('M_{xy} phase [rad]'),axis square,grid minor
xticklabels({'-1/TR','-0.5/TR','0','0.5/TR','1/TR'})
yticks([-1 -0.5 0 0.5 1]*pi),yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
legh=legend([l_h1;l_h2],'T2/T1=0.1','T2/T1=0.5','T2/T1=1','phase',...
    'Position',[0.782995255016859 0.132936513991583 0.105108053370403 0.165476185934884]);
fontsize(gcf,"scale",1.15)
%%

plotMinMaxgradients(1:2000, min(abs(bSSFP),[],1), max(abs(bSSFP),[],1));


%% only T2 


T1=[2000e-3,2500e-3];
T2=[20e-3];
TR=10e-3;
TE=5e-3;

dfreq=linspace(-1/TR,1/TR,2000);
get_opt=@(t1,t2,tr) acos((exp(-tr./t1)-exp(-tr./t2))./(1- exp(-tr./t1).*exp(-tr./t2)));
FA=get_opt(T1,T2,TR);
phi=0;
bSSFP=bSSFP_profile_Ganter2(FA(:),T1(:),T2(:),TE, TR,phi,dfreq*TR*2*pi);

figure(25),clf
l_h1=plot(dfreq*TR,abs(squeeze(bSSFP.')),'LineWidth',1.5);
% ylim([0 0.52]),ylabel('M_{xy} magnitude [n.u]')
yyaxis('right'),
ax=gca;ax.YAxis(2).Color=[0 0 0];
% l_h2=plot(dfreq*TR,angle(squeeze(bSSFP(1,:).')),'LineWidth',1.5,'LineStyle','--','Color','black'); hold on,
% plot(dfreq*TR,angle(squeeze(bSSFP(:,:).')),'LineWidth',1.5,'LineStyle','--');
title('T2/T1 dependence')
% set(gca,'ColorOrder',lines(3),'LineStyleOrder','--')
xlabel('off-resonance [Hz]'),ylabel('M_{xy} phase [rad]'),axis square,grid minor
 xticklabels({'-1/TR','-0.5/TR','0','0.5/TR','1/TR'})
% yticks([-1 -0.5 0 0.5 1]*pi),yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% legh=legend([l_h1;l_h2],'T2/T1=0.1','T2/T1=0.5','T2/T1=1','phase',...
%     'Position',[0.782995255016859 0.132936513991583 0.105108053370403 0.165476185934884]);
fontsize(gcf,"scale",1.15)


%%

function bSSFP=bSSFP_profile_Ganter2(flip,T1,T2,TE, TR,phi,off_resonance)
%function bSSFP=bSSFP_profile_Ganter(TR,TE,flip,phi,T1,T2)
%
%Calculation of the bSSFP profile (series of measurements with varying RF phase increment)
%Following Ganter, MRM, 2006: 
%Steady State of Gradient Echo Sequences with Radiofrequency Phase Cycling: Analytical Solution, Contrast Enhancement with Partial Spoiling
%
%INPUT
%
%repetition time in s (TR)
%echo time in s (TE)
%flip angle in rad (flip)
%vector of RF phase increments in rad (phi)
%off-resonance-related shift in rad (off_resonance)
%spin-lattice relaxation 
%spin-spin relaxation time in s (T2)
%
%OUTPUT
%
%complex bSSFP signal for the given range of RF phase increments phi (bSSFP)

M0 = 1;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);

theta = off_resonance - phi;
% Equation 44
D = (1-E1.*cos(flip)).*(1-E2.*cos(theta))-(E1-cos(flip)).*(E2-cos(theta)).*E2;
% equation 42
bSSFP = -(1i./D).*(1-E1).*sin(flip).*(1-E2.*exp(-1i*theta));

%% Read out signal at t=TE and 1i is probably RF-ADC phase offset

bSSFP = 1i*M0*bSSFP.*exp(-TE./T2).*exp(1i*off_resonance*((TE)/TR));
end



function plotMinMaxgradients(xaxis, ymin, ymax, yaxis)
    % PLOTMINMAXGRADIENTS Plots min-max range with monotonic gradient
    
    % Make sure inputs are column vectors
    xaxis = xaxis(:);
    ymin = ymin(:);
    ymax = ymax(:);
    
    % Create figure
    figure;
    hold on;
    
    % Number of steps for smooth gradient
    n_steps = 100;
    
    % Create x coordinates for the patch
    x_patch = [xaxis; flipud(xaxis)];
    
    % Initialize y and color matrices
    y_patch = zeros(n_steps, 2*length(xaxis));
    c_patch = zeros(n_steps, 2*length(xaxis));
    
    % Fill in y and color values
    for i = 1:length(xaxis)
        % Y values from min to max
        y_vals = linspace(ymin(i), ymax(i), n_steps);
        
        % Store forward and backward
        y_patch(:, i) = y_vals;
        y_patch(:, 2*length(xaxis) - i + 1) = flipud(y_vals);
        
        % Color values (monotonic from 0 to 1)
        c_vals = linspace(0, 1, n_steps);
        c_patch(:, i) = c_vals;
        c_patch(:, 2*length(xaxis) - i + 1) = flipud(c_vals);
    end
    
    % Create X matrix by repeating xaxis
    x_matrix = repmat(x_patch, n_steps, 1);
    
    % Plot the patch with interpolated colors
    patch(x_matrix, y_patch, c_patch, 'EdgeColor', 'none', 'FaceColor', 'interp');
    
    % Set colormap
    colormap jet;
    colorbar;
    
    % Plot mean line
    mean_vals = (ymin + ymax)/2;
    plot(xaxis, mean_vals, 'k', 'LineWidth', 2);
    
    % Labels
    xlabel('X');
    ylabel('Y');
    title('Min-Max Range with Monotonic Gradient');
    
    hold off;
end