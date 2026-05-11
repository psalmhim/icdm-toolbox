function show_coronal_tiles(R,T,ratio)
if nargin<3, ratio=0.0; end
K1=size(R,4);
dim=size(R);hdim=round(dim/2);
Rs=squeeze(R(:,hdim(2),:,:));
Ts=squeeze(T(:,hdim(2),:));
mR=max(abs(Rs(:)))*0.9; 
thr=mR*ratio;
t = tiledlayout(ceil(sqrt(K1)), ceil(sqrt(K1)));
t.TileSpacing = 'none';
t.Padding = 'none';
for d = 1:K1
    nexttile;
    rR = corient(Rs(:,:,d));
    rT = corient(Ts); 
    monet_overlay(rT,gray(256),[],rR,[-thr thr],jet(256),[-mR mR],0.2);
    axis image off;
end
colormap coolwarm;
end

function show_multi_slices(R,T,region,cmap,ratio)
if nargin<4,cmap=[]; end
if nargin<5, ratio=0.1; end
R = squeeze(R(:,:,:,region));
[m,mid]=max(R(:)); [x,y,z]=ind2sub(size(R),mid);
mR=max(abs(R(:)))*0.9; 
thr=mR*ratio;
if isempty(cmap)
    if min(R(:))<0
     cmap=coolwarm(256);
    else
        cmap=hot(256);
    end
end
subplot(1,3,1);
rR = corient(R(:,:,z));
rT = corient(T(:,:,z));
monet_overlay(rT,gray(256),[],rR,[-thr thr],cmap,[-mR mR],0.2);
axis image off;
title(sprintf('Region %d: Axial', region));

subplot(1,3,2);
rR = corient(R(:,y,:));
rT = corient(T(:,y,:));
monet_overlay(rT,gray(256),[],rR,[-thr thr],cmap,[-mR mR],0.2);
axis image off;
title('coronal');

subplot(1,3,3);
rR = corient(R(x,:,:));
rT = corient(T(x,:,:));
monet_overlay(rT,gray(256),[],rR,[-thr thr],cmap,[-mR mR],0.2);
axis image off;
title('Sagittal');

end


function img=corient(img)
    img = rot90(squeeze(img), -1);
end