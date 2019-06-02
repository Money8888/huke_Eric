function resultindex=choosedirection(direction,loc,lr_refer)
%最佳选择方向函数
%direction为可选的方向向量
%loc为待赋权的方向
%lr_refer指定方向上的权值
if ismember(loc,direction)
    if lr_refer>0.5
        resultindex=loc;
    else
        direction(find(direction==loc))=[];
        if (sum(direction)~=0)
            randindex=randperm(length(direction));
            resultindex=direction(randindex(1));
        else
            resultindex=loc;
        end
    end
else
    randindex=randperm(length(direction));
    resultindex=direction(randindex(1));
end
    