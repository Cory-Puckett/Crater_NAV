%{
Data reduction on crater database
%}
function craterDatabase =  DataReduce(maxDia,minDia, latTol, lonTol)
%{
maxDia = maximum diameter of crater to keep
minDia = minimum diameter of crater to keep
latTol = tolerance in latitude from about large crater center degrees
lonTol = tolerance in longitude from about large crater center degrees
%}
load craterDatabaseAll.mat craterDatabase
% maxDia = 55; %Maximum crater diameter to fit within frame
% minDia = 5; %minimum crater diameter to look for
ids = (craterDatabase.DIAM_CIRC_IMG <= maxDia)&(craterDatabase.DIAM_CIRC_IMG >= minDia);
craterDatabase = craterDatabase(ids,:);

%Renaming columns to shorter names
colnames = craterDatabase.Properties.VariableNames;
craterDatabase = renamevars(craterDatabase, colnames, ["ID","Lat","Lon","Diam"]);

%Limiting to area near large crater in image
bigCraterCord = [63.5, 93]; %center of big crater
%if in latitude range
ids = (craterDatabase.Lat >= bigCraterCord(1)-latTol);
ids = ids&(craterDatabase.Lat <= bigCraterCord(1)+latTol);
%if in long range
ids = ids & (craterDatabase.Lon >= bigCraterCord(2)-lonTol);
ids = ids & (craterDatabase.Lon <= bigCraterCord(2)+lonTol);

craterDatabase = craterDatabase(ids,:);

%Saving new table
save craterDatabaseSmall.mat craterDatabase
end