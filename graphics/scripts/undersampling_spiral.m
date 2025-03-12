%%
sn='Y:\MRDATA_ARCHIVE\mri_etl_s94t\studies\42\experiments\DZZN-Z3Q3\TWIX';
pn='C:\Users\pvalsala\Documents\ISMRM2024\skope_stuff';
%sync pulse medium
 fn= 'allData#S94Tuebingen#F67352#M80#D061123#T165214#peSpiral_R1_p8_GRE_med.dat'; 
%sync pulse high
% fn='allData#S94Tuebingen#F67354#M82#D061123#T165342#peSpiral_R1_p8_GRE_high3.dat';

r=SpiralReco(fullfile(sn,fn),'doCoilCombine','sos','CompMode','CPU2DHybrid', 'RepSel',5);

%%
figure(6),clf
tt=tiledlayout(2,2,'TileSpacing','tight','Padding','compact');
nexttile()
im=ndflip(squeeze(r.img(1,:,:,1,1,5)),[1,2]);
im=im./max(im(:));
imagesc(im)
axis image,colormap gray,xticks([]),yticks([])
clim([0 0.6])
title('Nyquist sampled')
nexttile()
imf=myfft(im,[1,2]);
imf(1:2:end,:)=0;
imr=myifft(imf,[1,2]);
imr=imr./max(imr(:));
imagesc(abs(imr)),axis image,colormap gray,xticks([]),yticks([])
clim([0 0.6])
title('Cartesian 2x undersampling ')
nexttile()
s=r.sig;
s(:,:,1:2:end)=0;
imr2=sos(r.NUFFT_obj'*permute(s,[2 3 4 1]),4);
imr2=imr2./max(imr2(:));
 imagesc(ndflip(abs(imr2),[1,2])),axis image,colormap gray,xticks([]),yticks([])
clim([0 0.6])
title('Spiral 2x undersampling')
nexttile()
s=r.sig;
s(:,:,2:4:end)=0;
s(:,:,3:4:end)=0;
s(:,:,4:4:end)=0;
imr2=sos(r.NUFFT_obj'*permute(s,[2 3 4 1]),4);
imr2=imr2./max(imr2(:));
% plotk(r.KTraj(:,1:2:end,1:1:end))
% axis square,xticks([]),yticks([])
 imagesc(ndflip(abs(imr2),[1,2])),axis image,colormap gray,xticks([]),yticks([])
clim([0 0.6])
title('Spiral 4x undersampling')

set(gcf,'color','w')
fontsize(gcf,"scale",1.3)

%  Position: [-1240 397 712 701]