function nextlist=neighbor(rindex,cindex,matrix)
%采用Von Neumann邻居规则
%nextlist为邻居四维数组，按顺序为上下左右
%rindex，cindex为元胞在矩阵matrix的行号和列号
%求元胞的邻居，-1的元胞(边界元胞)邻居全定义为inf
nextlist=zeros(4,1);
[m,n]=size(matrix);
if rindex==1 || rindex==m || cindex==1 || cindex==n
    if rindex==1 || rindex==m
        nextlist=inf*ones(4,1);
    end
    if cindex==1 && rindex~=1 && rindex~=m
       nextlist(1)=matrix(rindex-1,cindex);
       nextlist(2)=matrix(rindex+1,cindex);
       nextlist(3)=matrix(rindex,cindex+n-1);
       nextlist(4)=matrix(rindex,cindex+1);
    end
    if cindex==n && rindex~=1 && rindex~=m
       nextlist(1)=matrix(rindex-1,cindex);
       nextlist(2)=matrix(rindex+1,cindex);
       nextlist(3)=matrix(rindex,cindex-1);
       nextlist(4)=matrix(rindex,cindex-n+1);
    end
else
    nextlist(1)=matrix(rindex-1,cindex);
    nextlist(2)=matrix(rindex+1,cindex);
    nextlist(3)=matrix(rindex,cindex-1);
    nextlist(4)=matrix(rindex,cindex+1);
end
end