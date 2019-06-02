function [f1,f2]=intercept(statusmat)
%求函数
%f1为基膜截留
%f2为纳滤膜截留
%statusmat为状态矩阵
%m1,m2,m3为各物相矩阵的深度
%此时状态4背景状态，不是油相
m1=22;m2=5;m3=53;
[~,n]=size(statusmat);
core=statusmat(2:m1+1,2:n-1);
core(ceil(rand(ceil(0.2*(n-2)*m1),1)*m1*(n-2)))=7;
member=statusmat(1:m1+m2+1,:);%member基膜
member(2:m1+1,2:n-1)=core;
member=[member;4*ones(m3,n);-1*ones(1,n)];
member(find(member==1))=2;%将哌嗪替换为水分子
member(find(member==6))=2;%将PA替换为水分子

poly=statusmat;%聚合之后的矩阵
poly(find(poly==3))=4;%将TMC替换成4
poly(find(poly==1))=2;%将哌嗪替换为水分子
water=find(poly==2);
randindex=randperm(length(water));
randindex=randindex(1:ceil(0.2*length(water)));
poly(water(randindex))=7;
[rm,rn]=size(member);
nextcellm=cell(rm,rn);
nextcellp=cell(rm,rn);

for i=1:rm
    for j=1:rn
        list=neighbor(i,j,member);%获取邻居
        nextcellm{i,j}=list;
    end
end
for i=1:rm
    for j=1:rn
        list=neighbor(i,j,poly);
        nextcellp{i,j}=list;
    end
end
num=0;%迭代次数
number=uicontrol('style','text','string','1','fontsize',12,'position',[100,400,50,20]);
imh=matrixplot(member(2:m1+m2+m3+1,1:n));
%imh=matrixplot(poly(2:m1+m2+m3+1,1:n));
while num<=70
    for i=1:rm
        for j=1:rn
            switch member(i,j)
                case 7
                    times=0;
                    while(times<3)
                        if nextcellm{i,j}(2)==4 || nextcellm{i,j}(2)==2 && i<=m1+m2+1
                            [swapi,swapj]=judgedirection(2,i,j,rn);
                            temp=member(swapi,swapj);
                            member(swapi,swapj)=member(i,j);
                            member(i,j)=temp;
                            times=times+1;
                        end
                    end
            end
            list=neighbor(i,j,member);
            nextcellm{i,j}=list;
            %             imh=matrixplot(member(2:m1+m2+m3+1,1:n));
            %             pause(0.002);
            switch poly(i,j)
                case 7
                    if (nextcellp{i,j}(2)==4 || nextcellp{i,j}(2)==2 || nextcellp{i,j}(2)==6) && i<=m1+m2+1
                        if nextcellp{i,j}(2)==6
                            if  i>2 && i<m1+m2+1
                                [swapi,swapj]=judgedirection(2,i,j,rn);
                                temp=poly(2*i-swapi,swapj);                            %交换元胞状态
                                poly(2*i-swapi,j)=poly(i,j);
                                poly(i,j)=temp;
                            else
                                continue;
                            end
                        else
                            times=0;
                            while(times<3)
                                if nextcellm{i,j}(2)==4 || nextcellm{i,j}(2)==2 && i<=m1+m2+1
                                    [swapi,swapj]=judgedirection(2,i,j,rn);
                                    temp=member(swapi,swapj);
                                    member(swapi,swapj)=member(i,j);
                                    member(i,j)=temp;
                                    times=times+1;
                                end
                            end
                        end
                    end
            end
        end
    end
end
imh=matrixplot(member(2:m1+m2+m3+1,1:n));
pause(0.05);
num=num+1;
stepnumber = 1 + str2num(get(number,'string'));
set(number,'string',num2str(stepnumber));

rest1=member(m1+m2+2:end,:);
count1=length(find(rest1==7));
f1=count1/size(find(member==7),1);
rest2=poly(m1+m2+2:end,:);
count2=length(find(rest2==7));
f2=count2/size(find(poly==7),1);