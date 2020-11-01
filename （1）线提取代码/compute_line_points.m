function [a,eigVec1,ev,nx,ny,px,py,k0,k1,ismax]=compute_line_points(img,sigma,mode,low,high)
[m,n]=size(img);
%计算grdx grdy grdxx grdxy grdyy
% [k0]=convol(img,sigma,1);
% [k1]=convol(img,sigma,2);
% [k2]=convol(img,sigma,3); 
% [k3]=convol(img,sigma,4);
% [k4]=convol(img,sigma,5);
[k0,k1,k2,k3,k4] = Hessian2D(img,sigma); %高斯卷积 得到hessian矩阵，不止是hessian矩阵，hessian只得到Dxx，Dyy，Dxy

H1=cell(m,n);
eigVec1=cell(m,n);
eigVec2=cell(m,n); %这两个应该是特征向量
eigVal1=zeros(m,n);
eigVal2=zeros(m,n);%这两个应该是特征值
for i=1:m
   for j=1:n 
     H1{i,j}=[k2(i,j) k3(i,j);k3(i,j) k4(i,j)];    %这个是hessian矩阵 roi区域中每个点的
   end
end
% Compute the eigenvalues and eigenvectors of the Hessian matrix given by
 %  H{i,j} and sort them in descending order according to
 %  their absolute values. 
for i=1:m
   for j=1:n 
    [X,B]=eig(H1{i,j}); %这个点的hessian矩阵的特征值和特征向量
    % Hessian矩阵的特征值就是形容其在该点附近特征向量方向的凹凸性，特征值越大，凸性越强，特征向量是方向
    % 每一行或每一列，灰度是由小变大再变小的
  %  C=eig(H1);
    if(abs(B(1,1))>abs(B(2,2)))  %特征值绝对值大的，放在1中，小的放在2中，特征值和特征向量
        eigVec1{i,j}=X(:,1);
        eigVec2{i,j}=X(:,2);
        eigVal1(i,j)=B(1,1);
        eigVal2(i,j)=B(2,2);
    elseif(abs(B(1,1))<abs(B(2,2)))
        eigVec1{i,j}=X(:,2);
        eigVec2{i,j}=X(:,1);
        eigVal1(i,j)=B(2,2);
        eigVal2(i,j)=B(1,1);
    else 
        if(B(1,1)<B(2,2))    %特征值绝对值相同，负的放在1中，正的放在2中
        eigVec1{i,j}=X(:,1);
        eigVec2{i,j}=X(:,2);
        eigVal1(i,j)=B(1,1);
        eigVal2(i,j)=B(2,2);
        else
        eigVec1{i,j}=X(:,2);
        eigVec2{i,j}=X(:,1);
        eigVal1(i,j)=B(2,2);
        eigVal2(i,j)=B(1,1);    
        end
    end   
   end
end
a=ones(m,n);
b=ones(m,n);
nx=zeros(m,n);
ny=zeros(m,n);
px=zeros(m,n);
py=zeros(m,n);
ismax=zeros(m,n);
ev=zeros(m,n);
for i=1:m
   for j=1:n 
       if(mode==1)
        val=-eigVal1(i,j);
       else
           val=eigVa1(i,j); %这句话原始代码写错了，没用到 
           %其实这句话的意思是，mode==1就是暗背景亮条纹，否则的话就是亮背景暗条纹
       end
       if(-eigVal1(i,j)>0)  % Hessian矩阵是实对称的，一个特征向量，它的这个方向的二阶导数就是对应的特征值，小于零，为凸
          ev(i,j)= val;   
           nd=eigVec1{i,j};
            n1=nd(1); 
            n2=nd(2); %下面这一块就是二阶泰勒展开，具体原理看那个论文
            a(i,j)=k2(i,j)*n1*n1+2.0*k3(i,j)*n1*n2+k4(i,j)*n2*n2;
            b(i,j)=k0(i,j)*n1+k1(i,j)*n2;
%             if(abs(a)==0)
%                 num=0;
%                 break;
%             else
%                 num=1;
%                 t=-b(i,j)/a(i,j);
%             end
            [t,num]=solve_linear(a(i,j),b(i,j));
            if(num~=0)
                p1=t*n1; %保证二阶泰勒展开后，求得的一阶对x和y的导数均为0的位置，在原位置附近
                p2=t*n2;
                if((abs(p1)<0.6)&&(abs(p2)<0.6))
                    if(val>=low)
                        if(val<=high)
                            ismax(i,j)=2; %特征值在预设范围内的点，ismax为2 否则为1
                        else
                            ismax(i,j)=1;
                        end
                    end
                    nx(i,j)=n1; %法线方向
                    ny(i,j)=n2;
                    px(i,j)=i+p1; %求得的极值点
                    py(i,j)=j+p2;
                end
            end
       end
   end
end
end