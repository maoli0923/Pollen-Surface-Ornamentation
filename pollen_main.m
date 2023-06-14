%This is code for "Grass pollen surface ornamentation is diverse across 
%the phylogeny: evidence from northern South America and the global literature" 
list=dir('*_3.3');
s=1;
fid=fopen('name.txt','wt');
for k=1:length(list)
    cd(list(k).name)
    list2=dir('*_*');
    for kk=1:length(list2)
        cd(list2(kk).name)
        file=dir('*.tif');
        for j=1:length(file)
            tmp=imread(file(j).name);
            if size(tmp,1)~=1855
                I(:,:,s)=imresize(tmp(:,:,1),[1855 1855]);
                s
            else
                I(:,:,s)=tmp(:,:,1);
            end
            fprintf(fid,'%s\n',file(j).name(1:end-4));
            s=s+1;
        end
        cd ..
    end
    cd ..
end
fclose(fid);

for j=1:size(I,3)
I_adj(:,:,j)=adapthisteq(I(:,:,j),'ClipLimit',0.005,'Distribution','rayleigh');
end

cropindex=randi([1 856],size(I,3)*5,2);
save cropindex.mat cropindex
s=1;
se = strel('disk',3);
for j=1:size(I,3)
    %figure;
    for i=1:5
        h=I_adj(cropindex(s,1):cropindex(s,1)+999,cropindex(s,2):cropindex(s,2)+999,j);
        BW = edge(h,'Sobel');
        hc=imclose(BW,se);
        CC = bwconncomp(hc);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [m1 m2] = sort(numPixels,'descend');
        for ii=1:10
            hc2=ones(size(hc));
            hc2(CC.PixelIdxList{m2(ii)}) = 0;
            [x y]=find(hc2==0);
            [COEFF SCORE latent]=pca([x y]);
            tmp_size_sub(ii)=latent(1);
        end
        tmp_size(i)=mean(tmp_size_sub);
        h_rgb = cat(3,h,h,h);
        [X_no_dither,map] = rgb2ind(h_rgb,4,'nodither');
        [m1 m2]=max(map(:,1));
        a=find(X_no_dither==m2-1);
        BW2=zeros(size(h));
        BW2(a) = 1;
        h2c=imclose(BW2,se);
        CC = bwconncomp(h2c);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        tmp_density(i)=length(find(numPixels>37));
        s=s+1;
        %subplot(2,3,i)
        %imshow(hc)
    end
    Feature_adj(j,1)=mean(tmp_size);
    Feature_adj(j,2)=mean(tmp_density);
end


%This is code for the SC features developed for "Micro-ornamentation of the Poaceae pollen wall as a tool for identifying grasses in the fossil record" 
s=1;
for j=1:size(I_adj,3)
    for i=1:5
        h=I_adj(cropindex(s,1):cropindex(s,1)+999,cropindex(s,2):cropindex(s,2)+999,j);
        h=imresize(h,0.2);
        h=h(41:160,41:160);
        h_rgb = cat(3,h,h,h);
        [X_no_dither,map] = rgb2ind(h_rgb,4,'nodither');
        [m1 m2]=sort(map(:,1),'descend');
        a=find(X_no_dither==m2(1)-1|X_no_dither==m2(2)-1);
        BW=zeros(size(h));
        BW(a) = 1;      
        hh=bwareaopen(BW,5);
        hh=bwareaopen(~hh,5);
        hh=~hh;
        g(i,:)=pollen_newSC(hh);
        s=s+1;
    end
    SC2(j,1:20)=mean(g);
end 
save SC2.mat SC2

s=1;
for j=1:size(I_adj,3)
    for i=1:5
        h=I_adj(cropindex(s,1):cropindex(s,1)+999,cropindex(s,2):cropindex(s,2)+999,j);
        h=imresize(h,0.2);
        h=h(41:160,41:160);
        h_rgb = cat(3,h,h,h);
        [X_no_dither,map] = rgb2ind(h_rgb,4,'nodither');
        [m1 m2]=sort(map(:,1),'descend');
        a=find(X_no_dither==m2(1)-1|X_no_dither==m2(2)-1|X_no_dither==m2(3)-1);
        BW=zeros(size(h));
        BW(a) = 1;      
        hh=bwareaopen(BW,5);
        hh=bwareaopen(~hh,5);
        hh=~hh;
        g(i,:)=pollen_newSC(hh);
        s=s+1;
    end
    SC3(j,1:20)=mean(g);
end 
save SC3.mat SC3