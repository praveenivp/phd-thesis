 figure(8),clf,
 
tt=tiledlayout(2,6);
% rout=SpiralReco('allData#S94Tuebingen#F26564#M1262#D260822#T172245#peSpiral_R1_TR10_p6iso_test.dat');
% im=squeeze(rout.img(1,:,:,9,2));
k=rout.KTraj.*0.9;
im1=2*ndflip(im,[1,2])./max(im(:));

nexttile(1,[2 2])
 imagesc(linspace(-5e3,5e3,340),linspace(-5e3,5e3,340),im1),axis image,colormap('gray'),clim([0,1])
% hold on,scatter(col(k(1,1:10:end,1:10:end,1)),col(k(2,1:10:end,1:10:end,1)),'r.')
xticks([]),yticks([]),xlabel('K_x [rad/m]'),ylabel('K_y [rad/m]'),title('Image')

nexttile(3,[2 2])
  a=colormap('jet');
    a=interp1(linspace(0,1,length(a)),a, linspace(0,1,size(k,3)),'linear');
 imagesc(linspace(-5e3,5e3,340),linspace(-5e3,5e3,340),abs(log(myfft(im,[1,2])))),axis image,colormap('gray'),clim([3 5])
hold on,
for cintlv=1:10:size(k,3)
scatter(col(k(1,1:10:end,cintlv,1)),col(k(2,1:10:end,cintlv,1)),20,"Marker",".",'LineWidth',1,'MarkerFaceColor',a(cintlv,:),'MarkerEdgeColor',a(cintlv*1,:))
end
xticks([]),yticks([]),xlabel('K_x [rad/m]'),ylabel('K_y [rad/m]'),title('k-space')

[xx,yy]=ndgrid(linspace(-5e3,5e3,340));

nexttile(5,[1,2])
kb=real(exp(1i*(2e-3*xx+0*yy)));

imagesc([kb,(im1.*kb)]),colormap(gca,'parula'),axis image
xticks([]),yticks([]),title('Fourier basis | weighted image')

nexttile(11,[1,2])
kb=real(exp(1i*(3e-3*xx+2e-3*yy)));

imagesc([kb,(im1.*kb)]),colormap(gca,'parula'),axis image
xticks([]),yticks([]),title('Fourier basis | weighted image')

fontsize(gcf,"scale",1.4)

annotation(gcf,'arrow',[0.678932178932179 0.498556998557],...
    [0.766899766899767 0.564102564102564],...
    'Color',[1 0.411764705882353 0.16078431372549]);

annotation(gcf,'arrow',[0.676767676767677 0.51010101010101],...
    [0.27039627039627 0.461538461538462],...
    'Color',[1 0.411764705882353 0.16078431372549]);

%    Position: [-1531 467 1386 429]
