%% --------------------------------------------------------
clc;clear;
%主函数
% writerObj=VideoWriter('cellautoma.avi');
% writerObj.FrameRate=2;%帧速率
% open(writerObj);
number=uicontrol('style','text','string','1','fontsize',12,'position',[100,400,50,20]);
%% 构成元胞
%0为膜孔空隙，1为哌嗪分子，2为水分子，3为TMC分子，4为油相溶剂分子，5为基膜，6为聚酰胺
m1=22;m2=5;m3=53;   %各区域深度
n=50;               %模拟时的宽度
water=2*ones(m1*n,1);  %水相
membr=zeros(m2*n,1);  %界面层
oil=4*ones(m3*n,1);    %油相
active=zeros(m1+m2+m3,n); %PA活性设置
c1=0.2;c2=0.50;c3=0.02; %各个区域单体所占的百分比
membr(ceil(rand(ceil(m2*n*c2),1)*m2*n))=5;
membr=reshape(membr,m2,n); %含有基膜的界面矩阵
water(ceil(rand(ceil((m1*n+m2*n*(1-c2))*c1),1)*m1*n))=1;
water=reshape(water,m1,n); %含有单体的水相矩阵
oil(ceil(rand(ceil(m3*n*c3),1)*m3*n))=3;
oil=reshape(oil,m3,n);    %含有单体的油相矩阵
allmatrix=[water;membr;oil];%组合
%% 处理边界
simulation=zeros(m1+m2+m3+2,n+2); %多加了一层外层边界
simulation([1 m1+m2+m3+2],:)=-1*ones(2,n+2);%上下层边界元胞赋值为-1
simulation(:,1)=[-1;allmatrix(:,end);-1];
simulation(:,n+2)=[-1;allmatrix(:,1);-1];   %左右层采取周期性原则，模拟无限宽流体
simulation(2:m1+m2+m3+1,2:n+1)=allmatrix;
w_pip=0.8; %pip不形成团簇的概率
w_TMC=0.6; %TMC不形成团簇的概率
%% 元胞邻居
%调用neighbor函数，将邻居存储到每个元胞对应的元胞数组中
[rm,rn]=size(simulation);
nextcell=cell(rm,rn);
for i=1:rm
    for j=1:rn
        list=neighbor(i,j,simulation);
        nextcell{i,j}=list;
    end
end
%% 团簇元胞
grouppip=cell(rm,rn);  %%定义哌嗪团簇元胞
for i=1:rm
    for j=1:rn
        grouppip{i,j}=zeros(1,2);
    end
end
groupTMC=cell(rm,rn);  %%定义TMC团簇元胞
for i=1:rm
    for j=1:rn
        groupTMC{i,j}=zeros(1,2);
    end
end
%% 制定元胞规则
num=0;%迭代次数
imh=matrixplot(simulation(2:m1+m2+m3+1,1:n+2));
while num<=98
    %tempmatrix = simulation;                  %定义新的矩阵变量暂时保存当前画面
    for i=1:rm
        for j=1:rn
            switch simulation(i,j)
                case -1                       %边界元胞，不作处理
                case 0                        %膜孔元胞，不作处理
                case 1                        %哌嗪元胞
                    if ismember(simulation(i,j),grouppip{i,j})
                        continue;
                    else
                        if ismember(3,nextcell{i,j})        %若元胞上面邻居含有元胞3，则结合生成6
                            direction=find(nextcell{i,j}==3);
                            randindex=randperm(length(direction));
                            resultindex=direction(randindex(1));             %假设水分子运动为布朗运动
                            [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                            simulation(i,j)=6;                               %哌嗪位被PA替代
                            simulation(swapi,swapj)=4;                          %TMC位被油相溶剂替代
                            active(i,j)=2;  %假设pa只反应两次
                        elseif ismember(0,nextcell{i,j}) || ismember(2,nextcell{i,j}) || ismember(4,nextcell{i,j})
                            direction=union(find(nextcell{i,j}==0),union(find(nextcell{i,j}==2),find(nextcell{i,j}==4))); %寻找四个方向中含有元胞0,1,2,4的方向
                            p_max=1;p_min=0.5;                              %概率变化范围
                            if i<=m1    %若哌嗪在水相
                                p_i=p_max-(i-1)*(p_max-p_min)/m1;
                                resultindex=choosedirection(direction,2,p_i);
                                %resultindex=2;
                            elseif i>m1+m2 %若哌嗪在油相
                                p_i=p_min+(i-m1-m2)*(p_max-p_min)/m3;
                                resultindex=choosedirection(direction,1,p_i);
                            else           %若哌嗪在界面中，等概率扩散
                                %p_i=p_max-(i-1)*(p_max-p_min)/m1;
                                resultindex=choosedirection(direction,2,1);
                                %resultindex=2;
                            end
                            %由于界面聚合，水相哌嗪在界面处移动概率最小，两边概率最大，概率随着行数逐渐变化,依据概率选择最佳方向
                            [swapi,swapj]=judgedirection(resultindex,i,j,rn);   %获得选定方向元胞的行列号
                            if simulation(swapi,swapj)==0                    %如果扩散到膜孔,膜孔元胞置为哌嗪元胞，哌嗪元胞变成水分子元胞
                                simulation(swapi,swapj)=1;
                                simulation(i,j)=2;
                            elseif simulation(swapi,swapj)==4
                                if swapi<=m1+m2
                                    continue;
                                elseif swapi>m1+m2+min(size(find(simulation==1)),size(find(simulation==3)))/n+1   %哌嗪只能活跃于油相的表层,该层为期望膜厚度
                                    continue;
                                end
                            else
                                temp=simulation(i,j);                            %交换元胞状态
                                simulation(i,j)=simulation(swapi,swapj);
                                simulation(swapi,swapj)=temp;
                            end
                            
                        elseif ismember(1,nextcell{i,j})   %若邻居有自身，则依玻尔兹曼因子可能形成团簇
                            if rand>=w_pip && i>2 && i<rm-1 && j>2 && j<rn-1 %弹性碰撞
                                direction=find(nextcell{i,j}==1);
                                randindex=randperm(length(direction));
                                resultindex=direction(randindex(1));             %假设水分子运动为布朗运动
                                [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                                temp=simulation(2*i-swapi,2*j-swapj);                            %交换元胞状态
                                simulation(2*i-swapi,2*j-swapj)=simulation(i,j);
                                simulation(i,j)=temp;
                            else             %形成团簇
                                grouppip{i,j}=[simulation(i,j),simulation(swapi,swapj)];
                            end
                        end
                    end
                case 2                        %水分子元胞
                    if ismember(0,nextcell{i,j})  %若水分子扩散到膜孔
                        direction=find(nextcell{i,j}==0);
                        randindex=randperm(length(direction));
                        resultindex=direction(randindex(1));             %假设水分子运动为布朗运动
                        [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                        if simulation(swapi,swapj)==0
                            simulation(swapi,swapj)=2;
                        end
                    end
                case 3                        %TMC元胞
                    if ismember(simulation(i,j),groupTMC{i,j})
                        continue;
                    else
                        if ismember(1,nextcell{i,j})        %若元胞邻居含有元胞1，则结合生成6
                            direction=find(nextcell{i,j}==1);
                            randindex=randperm(length(direction));
                            resultindex=direction(randindex(1));
                            [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                            simulation(i,j)=4;
                            simulation(swapi,swapj)=6; %哌嗪位被PA替代
                            active(swapi,swapj)=2;
                        elseif ismember(4,nextcell{i,j})
                            direction=find(nextcell{i,j}==4);
                            p_max=1;p_min=0.5;
                            p_i=p_min+(i-m1-m2)*(p_max-p_min)/m3;%由于界面聚合，油相哌嗪要往界面去的概率更大，假设概率为0.7
                            resultindex=choosedirection(direction,1,p_i);  %1表示向上，依据概率选择最佳方向
                            [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                            if (simulation(swapi,swapj)==4 && swapi>m1+m2) || simulation(swapi,swapj)==3
                                temp=simulation(i,j);
                                simulation(i,j)=simulation(swapi,swapj);
                                simulation(swapi,swapj)=temp;
                            end
                        elseif ismember(3,nextcell{i,j})
                            if rand>=w_TMC && i>2 && i<rm-1 && j>2 && j<rn-1 %弹性碰撞
                                direction=find(nextcell{i,j}==3);
                                randindex=randperm(length(direction));
                                resultindex=direction(randindex(1));             %假设水分子运动为布朗运动
                                [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                                temp=simulation(2*i-swapi,2*j-swapj);                            %交换元胞状态
                                simulation(2*i-swapi,2*j-swapj)=simulation(i,j);
                                simulation(i,j)=temp;
                            else             %形成团簇
                                groupTMC{i,j}=[simulation(i,j),simulation(swapi,swapj)];
                            end
                        end
                    end
                case 4                        %HEX油相溶剂元胞
                case 5                        %基膜元胞
                case 6                        %PA元胞
                    p1=0.4286;p2=0.5714;           %经计算，p1表示6和1反应再生成6的概率，p2表示6和3反应再生成6的概率
                    if active(i,j)>0
                        if ismember(1,nextcell{i,j}) && ismember(3,nextcell{i,j}) %如果周围仍然有单体
                            flag=(p1*rand>p2*rand);
                            if flag==1  %和1发生反应
                                direction=find(nextcell{i,j}==1);
                                randindex=randperm(length(direction));
                                resultindex=randindex(1);
                                [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                                simulation(swapi,swapj)=6;
                                active(swapi,swapj)=1;
                                active(i,j)=active(i,j)-1;
                            else        %和3发生反应
                                direction=find(nextcell{i,j}==3);
                                resultindex=direction(1);
                                [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                                simulation(swapi,swapj)=6;
                                active(swapi,swapj)=1;
                                active(i,j)=active(i,j)-1;
                            end
                        elseif ismember(1,nextcell{i,j}) || ismember(3,nextcell{i,j}) %若6既和1又和3相邻
                            direction=union(find(nextcell{i,j}==1),find(nextcell{i,j}==3));
                            resultindex=direction(1);
                            [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                            simulation(swapi,swapj)=6;
                            active(swapi,swapj)=1;
                            active(i,j)=active(i,j)-1;
                        end
                    else
                        continue;
                    end
                    if i>m1+m2       %若生成的PA处于油相中，则应该以一定的规律附着在界面
                        direction=[1 2 3 4];
                        p_max=1;p_min=0.5;
                        p_i=p_min+(i-m1-m2)*(p_max-p_min)/m3;
                        resultindex=choosedirection(direction,1,p_i);
                        [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                        temp=simulation(i,j);
                        if simulation(swapi,swapj) ==0 || simulation(swapi,swapj)==5
                            continue;
                        elseif swapi<=m1+m2 && simulation(swapi,swapj)==3
                            continue;
                        elseif swapi>=m1-1 && simulation(swapi,swapj)==2
                            continue;
                        else
                            simulation(i,j)=simulation(swapi,swapj);
                            simulation(swapi,swapj)=temp;
                        end
                        %矩阵元素聚类
                    elseif i<m1
                        direction=[1 2 3 4];
                        p_max=1;p_min=0.5;
                        p_i=p_min-(i-m1)*(p_max-p_min)/m1;
                        resultindex=choosedirection(direction,2,p_i);
                        [swapi,swapj]=judgedirection(resultindex,i,j,rn);
                        temp=simulation(i,j);
                        if simulation(swapi,swapj) ==0 || simulation(swapi,swapj)==5
                            continue;
                        elseif swapi<=m1+m2 && simulation(swapi,swapj)==3
                            continue;
                        elseif swapi>=m1-1 && simulation(swapi,swapj)==2
                            continue;
                        else
                            simulation(i,j)=simulation(swapi,swapj);
                            simulation(swapi,swapj)=temp;
                        end
                    end
            end
            list=neighbor(i,j,simulation);
            nextcell{i,j}=list;
        end
    end
    [row,col]=find(simulation==6); %PA的电荷效应
    len=length(row);
    for i=1:len
        if row(i)>m1+m2
            active(row(i),col(i))=active(row(i),col(i))+1;
        elseif row(i)<m1
            active(row(i),col(i))=active(row(i),col(i))-1;
        end
    end

    
    imh=matrixplot(simulation(2:m1+m2+m3+1,1:n+2));
    pause(0.05);
%     frame=getframe(gcf);
%     writeVideo(writerObj,frame);
    num=num+1;
    stepnumber = 1 + str2num(get(number,'string'));
    set(number,'string',num2str(stepnumber));
    %legend('膜孔空隙','哌嗪分子','水分子','TMC分子','油相溶剂分子','基膜','聚酰胺');
end
% close(writerObj)
%% 可视化创建图形界面
% imh=imagesc(simulation(2:m1+m2+m3+1,2:n+1));
% plotbutton=uicontrol('style','pushbutton','string','Run','fontsize',12,'position',[150,400,50,20],'callback','run=1;');%运行
% erasebutton=uicontrol('style','pushbutton','string','Stop','fontsize',12,'position',[300,400,50,20],'callback','freeze=1;');%终止
% quitbutton=uicontrol('style','pushbutton','string','Quit','fontsize',12,'position',[450,400,50,20],'callback','freeze=1;');%退出
% number=uicontrol('style','text','string','1','fontsize',12,'position',[50,400,50,20]);

