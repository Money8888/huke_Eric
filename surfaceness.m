function R_a=surfaceness(statusmat)
%计算表面粗糙度
%statusmat为状态矩阵
%R_a为表面粗糙度
statusmat(find(statusmat==1))=0;
statusmat(find(statusmat==2))=0;
statusmat(find(statusmat==3))=0;
statusmat(find(statusmat==4))=0;
statusmat(find(statusmat==-1))=0;%将基膜和PA以及膜孔之外的元素全部赋值为0
[m,n]=size(statusmat);
meanval=zeros(1,n);
meancol=zeros(1,n);
for i=1:n
    index=find(statusmat(:,i)~=0);
    if length(index)~=0
        meanval(i)=mean(index);
    else
        meanval(i)=0;
    end
end
meankey=mean(meanval);%y平均
for i=1:n
    index=find(statusmat(:,i)~=0);
    if length(index)==0
        index=zeros(m,1);
    end
    meancol(i)=mean(index-meankey*ones(length(index),1));
end
R_a=mean(abs(meancol));


