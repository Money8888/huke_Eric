function yita=porosity(statusmat)
statusmat(find(statusmat==1))=-1;
statusmat(find(statusmat==3))=-1;
statusmat(find(statusmat==4))=-1;
tempmat=statusmat(23:27,:);
coll=union(find(tempmat==2),find(tempmat==0));

all=find(tempmat~=-1);
yita=length(coll)/length(all);
