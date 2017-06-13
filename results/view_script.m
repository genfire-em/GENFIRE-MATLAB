%% DFT : 146s
%% OS25: 45s
%% OS10: 6.4s
%% OS3 : 1.7s


markersize = 20;
load('GENFIRE_rec180_gridOS3.mat')
rec3 = GENFIRE_parameters.reconstruction;
load('GENFIRE_rec180_gridOS10.mat')
rec10 = GENFIRE_parameters.reconstruction;
load('GENFIRE_rec180_gridOS25.mat')
rec25 = GENFIRE_parameters.reconstruction;
load('GENFIRE_rec180_gridDft.mat')
recDft = GENFIRE_parameters.reconstruction;
load ../data/model
model = padarray(model,[0 0 60]);
[cc3, ir] = FourierShellCorrelate(model,rec3,15,.5);
[cc10, ir] = FourierShellCorrelate(model,rec10,15,.5);
[cc25, ir] = FourierShellCorrelate(model,rec25,15,.5);
[ccDft, ir] = FourierShellCorrelate(model,recDft,15,.5);


figure
plot(ir,cc3,'g','LineWidth',3)
hold on

plot(ir,cc10,'b:','LineWidth',3)
plot(ir,cc25,'k--','LineWidth',3)
plot(ir,ccDft,'r.','LineWidth',3)
plot(ir,cc3,'g+','MarkerSize',markersize)
plot(ir,cc10,'bo','MarkerSize',markersize)
plot(ir,cc25,'k^','MarkerSize',markersize)
plot(ir,ccDft,'r*','MarkerSize',markersize)
legend('OS3','OS10','OS25','DFT')
ht = title('FSC for varying 2D OS ratios');
set(ht,'FontSize',20)

rec3   = rec3(:,:,91-30:91+29);
rec10  = rec10(:,:,91-30:91+29);
rec25  = rec25(:,:,91-30:91+29);
recDft = recDft(:,:,91-30:91+29);
model  = model(:,:,91-30:91+29);

smodel1 = squeeze(sum(model(29:31,:,:),1));
smodel2 = squeeze(sum(model(:, 29:31,:),2));
smodel3 = squeeze(sum(model(:, :, 29:31),3));
figure
subplot(2,3,1), imagesc(smodel1),axis image,title('model')
subplot(2,3,2), imagesc(squeeze(sum(rec3(29:31,:,:),1))),axis image,title('OS 3')
subplot(2,3,3), imagesc(squeeze(sum(rec10(29:31,:,:),1))),axis image,title('OS 10')
subplot(2,3,4), imagesc(squeeze(sum(rec25(29:31,:,:),1))),axis image,title('OS 25')
subplot(2,3,5), imagesc(squeeze(sum(recDft(29:31,:,:),1))),axis image,title('DFT')

figure
subplot(2,3,1), imagesc(smodel2),axis image,title('model')
subplot(2,3,2), imagesc(squeeze(sum(rec3(:, 29:31,:),2))),axis image,title('OS 3')
subplot(2,3,3), imagesc(squeeze(sum(rec10(:, 29:31,:),2))),axis image,title('OS 10')
subplot(2,3,4), imagesc(squeeze(sum(rec25(:, 29:31,:),2))),axis image,title('OS 25')
subplot(2,3,5), imagesc(squeeze(sum(recDft(:, 29:31,:),2))),axis image,title('DFT')

figure
subplot(2,3,1), imagesc(smodel3),axis image,title('model')
subplot(2,3,2), imagesc(squeeze(sum(rec3(:, :, 29:31),3))),axis image,title('OS 3')
subplot(2,3,3), imagesc(squeeze(sum(rec10(:, :, 29:31),3))),axis image,title('OS 10')
subplot(2,3,4), imagesc(squeeze(sum(rec25(:, :, 29:31),3))),axis image,title('OS 25')
subplot(2,3,5), imagesc(squeeze(sum(recDft(:, :, 29:31),3))),axis image,title('DFT')

figure
subplot(2,3,1), imagesc(smodel1),axis image,title('model')
subplot(2,3,2), imagesc(squeeze(sum(rec3(29:31,:,:),1)) - smodel1),axis image,title('OS 3')
subplot(2,3,3), imagesc(squeeze(sum(rec10(29:31,:,:),1)) - smodel1),axis image,title('OS 10')
subplot(2,3,4), imagesc(squeeze(sum(rec25(29:31,:,:),1)) - smodel1),axis image,title('OS 25')
subplot(2,3,5), imagesc(squeeze(sum(recDft(29:31,:,:),1)) - smodel1),axis image,title('DFT')

figure
subplot(2,3,1), imagesc(smodel2),axis image,title('model')
subplot(2,3,2), imagesc(squeeze(sum(rec3(:, 29:31,:),2)) - smodel2),axis image,title('OS 3')
subplot(2,3,3), imagesc(squeeze(sum(rec10(:, 29:31,:),2)) - smodel2),axis image,title('OS 10')
subplot(2,3,4), imagesc(squeeze(sum(rec25(:, 29:31,:),2)) - smodel2),axis image,title('OS 25')
subplot(2,3,5), imagesc(squeeze(sum(recDft(:, 29:31,:),2)) - smodel2),axis image,title('DFT')

figure
subplot(2,3,1), imagesc(smodel3),axis image,title('model')
subplot(2,3,2), imagesc(squeeze(sum(rec3(:, :, 29:31),3)) - smodel1),axis image,title('OS 3')
subplot(2,3,3), imagesc(squeeze(sum(rec10(:, :, 29:31),3)) - smodel1),axis image,title('OS 10')
subplot(2,3,4), imagesc(squeeze(sum(rec25(:, :, 29:31),3)) - smodel1),axis image,title('OS 25')
subplot(2,3,5), imagesc(squeeze(sum(recDft(:, :, 29:31),3)) - smodel1),axis image,title('DFT')
