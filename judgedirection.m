function [swapi,swapj]=judgedirection(index,i,j,n)    %判断邻居向量中选择交换的元胞具体属于哪个方位
%index为方向数组的下标，i,j为中心元胞的行号和列号,n为列数
switch index
    case 1
        swapi=i-1;swapj=j;
    case 2
        swapi=i+1;swapj=j;
    case 3
        swapi=i;swapj=j-1;
    case 4
        swapi=i;swapj=j+1;
end
if swapj==0
    swapj=swapj+n;
end
if swapj==n+1
    swapj=swapj-n;
end
end
