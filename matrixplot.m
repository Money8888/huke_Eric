function imh=matrixplot(matrix)
%将矩阵转化为彩色图像，返回句柄
%矩阵中元素范围为0～6
B(:,:,1)=ones(size(matrix))*255;
B(:,:,2)=ones(size(matrix))*255;
B(:,:,3)=ones(size(matrix))*255;      %新建彩色矩阵B，颜色值为(255,255,255)白色
[m0,n0]=find(matrix==0);              %寻找元胞0的坐标
[m1,n1]=find(matrix==1);              %寻找元胞1的坐标
[m2,n2]=find(matrix==2);              %寻找元胞2的坐标
[m3,n3]=find(matrix==3);              %寻找元胞3的坐标
[m4,n4]=find(matrix==4);              %寻找元胞4的坐标
[m5,n5]=find(matrix==5);              %寻找元胞5的坐标
[m6,n6]=find(matrix==6);              %寻找元胞6的坐标
for k=1:numel(m0)                     %numel表示向量m的大小
    B(m0(k),n0(k),1)=0;
    B(m0(k),n0(k),2)=255;
    B(m0(k),n0(k),3)=255;             %膜孔为暖灰色
end
for k=1:numel(m1)                     
    B(m1(k),n1(k),1)=255;
    B(m1(k),n1(k),2)=0;
    B(m1(k),n1(k),3)=0;               %哌嗪为红色
end
for k=1:numel(m2)                     
    B(m2(k),n2(k),1)=0;
    B(m2(k),n2(k),2)=0;
    B(m2(k),n2(k),3)=255;             %水为蓝色
end
for k=1:numel(m3)                     
    B(m3(k),n3(k),1)=255;
    B(m3(k),n3(k),2)=255;
    B(m3(k),n3(k),3)=255;             %TMC为白色
end
for k=1:numel(m4)                     
    B(m4(k),n4(k),1)=255;
    B(m4(k),n4(k),2)=255;
    B(m4(k),n4(k),3)=0;               %油相溶剂为黄色
end
for k=1:numel(m5)                     
    B(m5(k),n5(k),1)=255;
    B(m5(k),n5(k),2)=0;
    B(m5(k),n5(k),3)=255;             %基膜为蛋壳色
end
for k=1:numel(m6)                     
    B(m6(k),n6(k),1)=0;
    B(m6(k),n6(k),2)=255;
    B(m6(k),n6(k),3)=0;               %PA为绿色
end
imh=image(B);