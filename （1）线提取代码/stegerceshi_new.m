
close all;
clear all;
clc;

%%
%sigma��low high������ֵ�������С�����ã�����һ���̶���Ӱ��ȡ�����������������ú�֮�󣬿���ͨ���ⲿ����ȡ��������Ż�
%ʵ�ڲ��У�������ͼ�л����������⼸�������ڣ��������е㣨������Թ���β����ͷ����

% һ�����������߿�����������ȡ���� �� �����ܸ�ǿ ��ƪ����ԭ���ĺ�ϸ�¡����
% sigma��ѡȡ���߿��й�
%%
num=1;
for ii=1:num

I=imread([num2str(ii),'y.bmp']);
% I=imread(['E:\���ⱨ��+������\���������ϻ���\������ش���\��ƴ�Ӵ���\����ȡ����\',num2str(ii),'y.bmp']);
% figure,imshow(255-I);
% I=im2bw(I,0.99);
% figure,imshow(I);
I=I*2;

figure,imshow(255-I);
gp=ginput();
zs=gp(1,:);
zx=gp(2,:);
yx=gp(3,:);
ys=gp(4,:);
save('hy0114','zs','zx','yx','ys');    %���г��򣬵��ĸ��ǵ�ȷ����ȡ����

% load('hy0114');  %hy0114��������ʾ����

numpix=10;
thresoldd=80;
% c=[zs(1) zx(1) yx(1) ys(1) ];
% r=[zs(2) zx(2) yx(2) ys(2) ];
c=[zs(1) zx(1) yx(1) ys(1) ];
r=[zs(2) zx(2) yx(2) ys(2) ];
BW=roipoly(I,c,r);
BW2=uint8(BW);
I=I.*BW2;
% figure,imshow(I);
% I=I(zs(2):yx(2),zs(1):yx(1));
% figure,imshow(I);
% [tranXY,img]=imgcrop(I); %��һ�������roi�������С������ֵ���ڶ�������ǻ�����roi���򣬰����������ȡ������
img=double(I);
tranXY=[1,1];
% img=I;
% t2=clock;
% time1=etime(t2,t1);%�ü�ͼ���ʱ�䣻

%%%%%%%%%%%%%%%%%%%�ϰ汾������%%%%%%%%%
tranx(1)=tranXY(1,1);
trany(1)=tranXY(1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(img);
sigma=2.3; %�����С����ʹȡ�������


%sigma=3.96;
%���׵�������ֵΪ���ܵ����ĵ�
%�����������ֵ��Ӧ����������
%���߷���Ϊnx ny
%��Ϊpx py
[a,eigVec,eigVal,n1,n2,p1,p2,dx,dy,ismax]=compute_line_points(img,sigma,1,0.34,0.85); %�����С����ʹȡ�������
% figure,imshow(img,[]);
% hold on
% plot(p2,p1,'.');
% hold off
% [a,eigVec,eigVal,n1,n2,p1,p2,dx,dy,ismax]=compute_line_points(img,sigma,1,0.4,0.9);
% t21=clock;
% time21=etime(t21,t2);%����Hessian����Ͷ��׵�������ֵ���ʱ��
%��������ֵƽ��ֵ
% ���������������һֱ��posx��posy����û���õģ�û���õ����溯�������
% �����������Ӧ���������õĶ����ģ������Ҵ����
% contours.row��.col ����������������һ�Խ���᲻��������

mean=0;
for i=1:height
    for j=1:width
        if(eigVal(i,j)>0)
        mean=mean+eigVal(i,j);
        end
    end
end
mean_eig=mean/(width*height);
% [num_junc,num_cont,cont,junction,cross]=get_contours(img,ismax,eigVal,n1,n2,p1,p2,sigma,100,1,thresoldd,1,dx,dy,width,height,mean_eig);
% ������һ�飬ʵ���Ͼ��ǰ���ɢ��ȥ����ֻ����һ���ߵĹ���
% t22=clock;
% time22=etime(t22,t21);%�������ӳ��ߵ�ʱ��
% [~,num_contours]=size(cont);
% [contours,correction,asymm]=compute_line_width(dx,dy,width,height,sigma,1,cont,num_contours);
% 
% qqq=cont.row;
% www=cont.col;
% figure,imshow(img,[]);
% hold on
% plot(www,qqq,'.');
% hold off
% t3=clock;
% time23=etime(t3,t22);%�����߿��ò������ĵ�ʱ��
% time2=etime(t3,t2);%Hessian����ʱ�䣻
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ͼָ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t4=clock;
% time3=etime(t4,t3);%Hessian����ʱ�䣻
[k,h]=size(p1);
p=1;
for i=1:k
    for j=1:h
        if((p1(i,j)~=0)||(p2(i,j)~=0))
        posx(p)=p1(i,j);
        posy(p)=p2(i,j);
        end
        p=p+1;
    end   
end
% figure,imshow(I);
% hold on
% plot(posy,posx,'.');
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',@NewCallback) %���������α�ȡѡȡ�ĵ㣬����ȡ6λ��Ч����
% hold off
% figure(2);
% imshow(I);
% hold on;
% plot(posy,posx,'.');
% save(gcf,
% hold on;
% [~,n]=size(cutcont);
% for i=1:n
%     plot(cutcont(i).col,cutcont(i).row,'-g');
%     hold on;
% end
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ݴ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ɾ��ȡ����0,0��
j=1;
px=[];
py=[];
for i=1:size(posy,2)
    if posx(1,i)~=0;
        px(j,1)=posx(1,i);
        py(j,1)=posy(1,i);
        j=j+1;
    end
end
%% ѡ����ĻҶ�ֵҪ����֮ǰ����ֵ
j=1;
pxcho=[];
pycho=[];
for i=1:size(px,1)
    a1=round(px(i,1));
    a2=round(py(i,1));
    if a1==0
        a1=1;
    end
    if a1>size(img,1)
        a1=size(img,1);
    end
    if a2==0
        a2=1;
    end
    if a2>size(img,2)
        a2=size(img,2);
    end
        
    if I(a1,a2)>thresoldd
        pycho(j,1)=py(i,1);
        pxcho(j,1)=px(i,1);
        j=j+1;
    end
end
% hold on
% plot(pycho,pxcho,'.');
% hold off
%% ϸ������ʴ��С�ı߽磬�ж�֮ǰѡ���ĵ��Ƿ�Ϊ1��ȥ���߽��ϵ���ɢ�㣩
% I2=I;
% I2=im2bw(I2,100/255);
% % figure,imshow(I2);
% % se=strel('disk',2);
% % I2ero=imerode(I2,se);
% % figure,imshow(I2ero);
% I2thin=bwmorph(I2,'thin',1); %ϸ�� %����������Ϊϸ������
% % figure,imshow(I2thin);
% % se=strel('disk',1);
% % I2ero=imerode(I2thin,se); %��ʴ
% % figure,imshow(I2ero);
% I2p=I2thin;
% j=1;
% pxcho2=[];
% pycho2=[];
% for i=1:size(pycho,1)
%     if I2p(round(pxcho(i,1)),round(pycho(i,1)))==1
%         pycho2(j,1)=pycho(i,1);
%         pxcho2(j,1)=pxcho(i,1);
%         j=j+1;
%     end
% end
% % hold on
% % plot(pycho2,pxcho2,'.');
% % hold off
%% ��������ȫ��Ϊ1����������㣬����ɾ��
% %ϸ�����࣬�ᵼ�°�����ȡ��ʱ�򣬰ѹ�ϸλ�ô��ĵ�����
% I3=I2thin;
% nei=8;
% for i=1:nei/8
%     for j=1:size(pxcho2,1)
%         xnei=round(pxcho2(j,1));
%         ynei=round(pycho2(j,1));
%         neipo(1,j)=I3(xnei+i,ynei);
%         neipo(2,j)=I3(xnei-i,ynei);
%         neipo(3,j)=I3(xnei,ynei+i);
%         neipo(4,j)=I3(xnei,ynei-i);
%         neipo(5,j)=I3(xnei+i,ynei+i);
%         neipo(6,j)=I3(xnei+i,ynei-i);
%         neipo(7,j)=I3(xnei-i,ynei+i);
%         neipo(8,j)=I3(xnei-i,ynei-i);
%     end
% end
%      neiposum=sum(neipo);
%      j=1;
%      pxcho3=[];
% pycho3=[];
%      for i=1:size(neiposum,2)
%          if neiposum(1,i)==nei
%             pxcho3(j,1)=pxcho2(i,:);
%             pycho3(j,1)=pycho2(i,:);
%              j=j+1;
%          end
%      end
%    hold on
% plot(pycho3,pxcho3,'.');
%% ���е�ĵط�Ϊ1�������ط�Ϊ0��������ͨ����ͨ��С�ģ�ɾ������ɾ�����ܱߵ���ɢ��
Ifilt=zeros(size(I));
for i=1:size(pxcho,1)
    Ifilt(round(pxcho(i)),round(pycho(i)))=1;
end
% figure,imshow(Ifilt);
Ifilt2=bwareaopen(Ifilt,numpix);  %С��ָ��������Ŀ�Ļᱻɾ��
% figure,imshow(Ifilt2);
pxcho3=[];
pycho3=[];
j=1;
for i=1:size(pycho,1)
    if Ifilt2(round(pxcho(i,1)),round(pycho(i,1)))==1
        pycho3(j,1)=pycho(i,1);
        pxcho3(j,1)=pxcho(i,1);
        j=j+1;
    end
end
% hold on
% plot(pycho3,pxcho3,'.');
% hold off
%% �����ı����ڣ�ɾ����ȥ����β��Щ�㣩
% xvs= [961 963 991 992 961]; %x���� ͷ���ı���
% yvs= [188 222 222 191 188];%y����
% xvw= [994 945 1032 1034 994]; %x���� β���ı���
% yvw= [566 629 635 566 566];%y����
xvs= [0 0 0 0 0]; %x���� ͷ���ı���
yvs= [0 0 0 0 0];%y����
xvw= [0 0 0 0 0]; %x���� β���ı���
yvw= [0 0 0 0 0];%y����

% xv3=[923 959 960 919 923];
% yv3=[786 786 835 837 786];



pxcho4=[];
pycho4=[];
j=1;
for i=1:size(pycho3,1)
    inparaw(i,1)=inpolygon(pxcho3(i),pycho3(i),yvw,xvw); %1Ϊ�ڲ��ͱ߽磬0Ϊ�ⲿ
    inparas(i,1)=inpolygon(pxcho3(i),pycho3(i),yvs,xvs);
%     inpara3(i,1)=inpolygon(pxcho3(i),pycho3(i),yv3,xv3);
%     if inparaw(i,1)==0&&inparas(i,1)==0&&inpara3(i,1)==0
    if inparaw(i,1)==0&&inparas(i,1)==0
        pycho4(j,1)=pycho3(i,1);
        pxcho4(j,1)=pxcho3(i,1);
        j=j+1;
    end
end
% figure,imshow(I2thin);
% hold on
% plot(pycho4,pxcho4,'.');
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',@NewCallback) %���������α�ȡѡȡ�ĵ㣬����ȡ6λ��Ч����
% hold off
%% ���ù����ĵ� �ҵ�,������ȡ��ֵ
[pxcho4so,indexso] = sort(pxcho4);
pxcho5=pxcho4;
pycho5=pycho4;
% for i=1:size(pxcho4,1)
%     numso=size(pxcho4,1)-i;
%     xdif=[];
%     j=1;
%     if numso>=3 %��������3���㣬����3���㣬����ʣ�µĵ�
%         for j=1:3
%             xdif(j)=abs(pxcho4so(i,1)-pxcho4so(i+j,1));
%             if xdif(j)<=0.3 %���x�Ĳ�ֵ��0.3����
%                 locso1=indexso(i,1);
%                 locso2=indexso(i+j,1);
%                 dify=abs(pycho4(locso1,1)-pycho4(locso2,1));
%                 if dify<0.3 %���y�Ĳ�ֵ��0.3����
%                     pxcho5(locso1,1)=(pxcho5(locso1,1)+pxcho5(locso2,1))/2;
%                     pxcho5(locso2,1)=0;
%                     pycho5(locso1,1)=(pycho5(locso1,1)+pycho5(locso2,1))/2;
%                     pycho5(locso2,1)=0;
%                 end
%             end
%         end
%     else
%         for j=1:numso
%             xdif(j)=abs(pxcho4so(i,1)-pxcho4so(i+j,1));
%             if xdif(j)<=0.3 %���x�Ĳ�ֵ��0.3����
%                 locso1=indexso(i,1);
%                 locso2=indexso(i+j,1);
%                 dify=abs(pycho4(locso1,1)-pycho4(locso2,1));
%                 if dify<0.3 %���y�Ĳ�ֵ��0.3����
%                     pxcho5(locso1,1)=(pxcho5(locso1,1)+pxcho5(locso2,1))/2;
%                     pxcho5(locso2,1)=0;
%                     pycho5(locso1,1)=(pycho5(locso1,1)+pycho5(locso2,1))/2;
%                     pycho5(locso2,1)=0;
%                 end
%             end
%         end
%     end
%     
% end

for i=1:size(pxcho4,1)-1
    xdif=abs(pxcho4so(i,1)-pxcho4so(i+1,1));
    if xdif<=0.4 %���x�Ĳ�ֵ��0.3����
        locso1=indexso(i,1);
        locso2=indexso(i+1,1);
        dify=abs(pycho4(locso1,1)-pycho4(locso2,1));
        if dify<0.4 %���� ���y�Ĳ�ֵ��0.3����
            pxcho5(locso1,1)=(pxcho5(locso1,1)+pxcho5(locso2,1))/2;
            pxcho5(locso2,1)=0;
            pycho5(locso1,1)=(pycho5(locso1,1)+pycho5(locso2,1))/2;
            pycho5(locso2,1)=0;
        end
    end
end
pxcho5(find(sum(abs(pxcho5),2)==0),:)=[];
pycho5(find(sum(abs(pycho5),2)==0),:)=[];
k=120;
figure,imshow(I);
hold on
plot(pycho5,pxcho5,'.');
% gfframe=getframe(gcf);
% imwrite(gfframe.cdata,strcat('F:\��άƴ��\��άƴ��ʵ������\0103����������\˫Ŀ�߱궨\�߼���ͼ\11\',num2str(ii),'y.jpg'));
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',@NewCallback) %���������α�ȡѡȡ�ĵ㣬����ȡ6λ��Ч����
hold off
% close 
end
%           j=1;          
%      for i=1:297
%          for k=1:605
%              if pxcho5(i)==pxcho5s(k)&&pycho5(i)==pycho5s(k)
%                  j=j+1;
%              end
%          end
%      end
 
% ppy=floor(pycho5);
% ppx=floor(pxcho5);
% 
% ii=zeros(size(I,1),size(I,2));
% for i=1:size(ppy,1)
%     ii(ppx(i,1),ppy(i,1))=255;
% end
% figure,imshow(ii);
    