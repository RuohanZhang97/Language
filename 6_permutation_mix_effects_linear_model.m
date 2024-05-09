clear
clc

datapath = '/home1/zhangruohan/Documents/LanguageComprehension/';
load(fullfile(datapath,'results','sorted_dLI_temporal_properties.mat'));
load(fullfile(datapath,'data_supply','ArticleClassification.mat'));
load(fullfile(datapath,'data_supply','SubSexInfo.mat'));
savefolder = fullfile(datapath,'results','permutaion_mix_effects_results');
if ~exist(savefolder,'dir')
    mkdir(savefolder);
end

ind_sensible_article = find([ArticleClassification.ArticleClass1] == 1);
NumSensibleArticle = length(ind_sensible_article);
ind_rational_article = find([ArticleClassification.ArticleClass1] == 2);
NumRationalArticle = length(ind_rational_article);
State1_OccurrenceRate_Sensible = State1_OccurrenceRateResults(ind_sensible_article);
State1_OccurrenceRate_Rational = State1_OccurrenceRateResults(ind_rational_article);

SubAge = [SexInfo.Age];
ind_male = find([SexInfo.Sex] == 0);
Age_male = SubAge(ind_male);
ind_female = find([SexInfo.Sex] == 1);
Age_female = SubAge(ind_female);

NumMale = length(ind_male);
NumFemale = length(ind_female);
NumSub = NumMale + NumFemale;
NoSub = [1:NumSub];

occurrence_rate_sensible_male = [];
occurrence_rate_sensible_female = [];
for i = 1:length(State1_OccurrenceRate_Sensible)
    occurrence_rate_sensible = [State1_OccurrenceRate_Sensible(i).Result.occurrence_rate];
    occurrence_rate_sensible_male = [occurrence_rate_sensible_male occurrence_rate_sensible(1,ind_male)];   
    occurrence_rate_sensible_female = [occurrence_rate_sensible_female occurrence_rate_sensible(1,ind_female)];   
end
State1_OccurrenceRate_sensible_male = occurrence_rate_sensible_male;
State1_OccurrenceRate_sensible_female = occurrence_rate_sensible_female;

occurrence_rate_rational_male = [];
occurrence_rate_rational_female = [];
for i = 1:length(State1_OccurrenceRate_Rational)
    occurrence_rate_rational = [State1_OccurrenceRate_Rational(i).Result.occurrence_rate];
    occurrence_rate_rational_male = [occurrence_rate_rational_male occurrence_rate_rational(1,ind_male)];   
    occurrence_rate_rational_female = [occurrence_rate_rational_female occurrence_rate_rational(1,ind_female)];   
end
State1_OccurrenceRate_rational_male = occurrence_rate_rational_male;
State1_OccurrenceRate_rational_female = occurrence_rate_rational_female;

NumWin_sensible_male = [];
NumWin_sensible_female = [];
for i = 1:length(State1_OccurrenceRate_Sensible)
    NumWin_sensible_male = [NumWin_sensible_male repelem(State1_OccurrenceRate_Sensible(i).NumWin,NumMale)];
    NumWin_sensible_female = [NumWin_sensible_female repelem(State1_OccurrenceRate_Sensible(i).NumWin,NumFemale)];
end

NumWin_rational_male = [];
NumWin_rational_female = [];
for i = 1:length(State1_OccurrenceRate_Rational)
    NumWin_rational_male = [NumWin_rational_male repelem(State1_OccurrenceRate_Rational(i).NumWin,NumMale)];
    NumWin_rational_female = [NumWin_rational_female repelem(State1_OccurrenceRate_Rational(i).NumWin,NumFemale)];
end

Age_sensible_male = repmat(Age_male,1,length(State1_OccurrenceRate_Sensible));
Age_sensible_female = repmat(Age_female,1,length(State1_OccurrenceRate_Sensible));
Age_rational_male = repmat(Age_male,1,length(State1_OccurrenceRate_Rational));
Age_rational_female = repmat(Age_female,1,length(State1_OccurrenceRate_Rational));

%% t-test on Occurrence Rate
NumPerm = 10000;
% Male - Female all articles
Subject = [repmat(NoSub(1,ind_male),1,NumSensibleArticle) repmat(NoSub(1,ind_male),1,NumRationalArticle)...
    repmat(NoSub(1,ind_female),1,NumSensibleArticle) repmat(NoSub(1,ind_female),1,NumRationalArticle)];
y1 = [State1_OccurrenceRate_sensible_male State1_OccurrenceRate_rational_male];
NumWin1 = [NumWin_sensible_male NumWin_rational_male];
Age1 = [Age_sensible_male Age_rational_male];
ArticleClass1 = [zeros(1,length(State1_OccurrenceRate_sensible_male)) ones(1,length(State1_OccurrenceRate_rational_male))];

y2 = [State1_OccurrenceRate_sensible_female State1_OccurrenceRate_rational_female];
NumWin2 = [NumWin_sensible_female NumWin_rational_female];
Age2 = [Age_sensible_female Age_rational_female];
ArticleClass2 = [zeros(1,length(State1_OccurrenceRate_sensible_female)) ones(1,length(State1_OccurrenceRate_rational_female))];

DependentVariable = [y2'; y1'];
Covariate = [NumWin2' Age2'; NumWin1' Age1'];
[~,~,r] = regress(DependentVariable,[ones(size(Covariate,1),1) Covariate]);
GroupLabel = [zeros(size(y2',1),1); ones(size(y1',1),1)];
ArticleClass = [ArticleClass2 ArticleClass1];

data = table;
data.y = r;
data.group = GroupLabel;
data.subject = Subject';
data.article_class = ArticleClass';
data = rmmissing(data);
glme = fitglme(data,'y ~ 1 + group + (1|subject) + (1|article_class)');
State1_OccurrenceRate_beta_value1 = glme.Coefficients.Estimate(2);
State1_OccurrenceRate_P_value1 = glme.Coefficients.pValue(2);
State1_OccurrenceRate_T_value1 = glme.Coefficients.tStat(2);

for i = 1:NumPerm
    data = table;
    data.y = r(randperm(length(r))');
    data.group = GroupLabel;
    data.subject = Subject';
    data.article_class = ArticleClass';
    data = rmmissing(data);
    glme = fitglme(data,'y ~ 1 + group + (1|subject) + (1|article_class)');
    perm_State1_OccurrenceRate_beta_value1(i) = glme.Coefficients.Estimate(2);
    perm_State1_OccurrenceRate_P_value1(i) = glme.Coefficients.pValue(2);
    perm_State1_OccurrenceRate_T_value1(i) = glme.Coefficients.tStat(2);
end

if State1_OccurrenceRate_beta_value1 > 0
    perm_pVal1 = (1+length(find(perm_State1_OccurrenceRate_beta_value1 >= State1_OccurrenceRate_beta_value1)))/(NumPerm+1);
else
    perm_pVal1 = (1+length(find(perm_State1_OccurrenceRate_beta_value1 <= State1_OccurrenceRate_beta_value1)))/(NumPerm+1);
end    

% Rational - Sensible all subjects
Subject = [repmat(NoSub(1,ind_male),1,NumSensibleArticle) repmat(NoSub(1,ind_female),1,NumSensibleArticle)...
    repmat(NoSub(1,ind_male),1,NumRationalArticle) repmat(NoSub(1,ind_female),1,NumRationalArticle)];
y3 = [State1_OccurrenceRate_sensible_male State1_OccurrenceRate_sensible_female];
NumWin3 = [NumWin_sensible_male NumWin_sensible_female];
Age3 = [Age_sensible_male Age_sensible_female];
Sex3 = [ones(1,length(State1_OccurrenceRate_sensible_male)) zeros(1,length(State1_OccurrenceRate_sensible_female))];

y4 = [State1_OccurrenceRate_rational_male State1_OccurrenceRate_rational_female];
NumWin4 = [NumWin_rational_male NumWin_rational_female];
Age4 = [Age_rational_male Age_rational_female];
Sex4 = [ones(1,length(State1_OccurrenceRate_rational_male)) zeros(1,length(State1_OccurrenceRate_rational_female))];

DependentVariable = [y3'; y4'];
Covariate = [NumWin3' Age3'; NumWin4' Age4'];
[~,~,r] = regress(DependentVariable,[ones(size(Covariate,1),1) Covariate]);
GroupLabel = [zeros(size(y3',1),1); ones(size(y4',1),1)];
Sex = [Sex3 Sex4];

data = table;
data.y = r;
data.group = GroupLabel;
data.subject = Subject';
data.sex = Sex';
data = rmmissing(data);
glme = fitglme(data,'y ~ 1 + group + (1|subject) + (1|sex)');
State1_OccurrenceRate_beta_value2 = glme.Coefficients.Estimate(2);
State1_OccurrenceRate_P_value2 = glme.Coefficients.pValue(2);
State1_OccurrenceRate_T_value2 = glme.Coefficients.tStat(2);    

for i = 1:NumPerm
    data = table;
    data.y = r(randperm(length(r))');
    data.group = GroupLabel;
    data.subject = Subject';
    data.sex = Sex';
    data = rmmissing(data);
    glme = fitglme(data,'y ~ 1 + group + (1|subject) + (1|sex)');
    perm_State1_OccurrenceRate_beta_value2(i) = glme.Coefficients.Estimate(2);
    perm_State1_OccurrenceRate_P_value2(i) = glme.Coefficients.pValue(2);
    perm_State1_OccurrenceRate_T_value2(i) = glme.Coefficients.tStat(2);  
end

if State1_OccurrenceRate_beta_value2 > 0
    perm_pVal2 = (1+length(find(perm_State1_OccurrenceRate_beta_value2 >= State1_OccurrenceRate_beta_value2)))/(NumPerm+1);
else
    perm_pVal2 = (1+length(find(perm_State1_OccurrenceRate_beta_value2 <= State1_OccurrenceRate_beta_value2)))/(NumPerm+1);
end 

% Male - Female Sensible
Subject = [repmat(NoSub(1,ind_male),1,NumSensibleArticle) repmat(NoSub(1,ind_female),1,NumSensibleArticle)];
y5 = State1_OccurrenceRate_sensible_male;
NumWin5 = NumWin_sensible_male;
Age5 = Age_sensible_male;
y6 = State1_OccurrenceRate_sensible_female;
NumWin6 = NumWin_sensible_female;
Age6 = Age_sensible_female;

DependentVariable = [y6'; y5'];
Covariate = [NumWin6' Age6'; NumWin5' Age5'];
[~,~,r] = regress(DependentVariable,[ones(size(Covariate,1),1) Covariate]);
GroupLabel = [zeros(size(y6',1),1); ones(size(y5',1),1)];

data = table;
data.y = r;
data.group = GroupLabel;
data.subject = Subject';
data = rmmissing(data);
glme = fitglme(data,'y ~ 1 + group + (1|subject)');
State1_OccurrenceRate_beta_value3 = glme.Coefficients.Estimate(2);
State1_OccurrenceRate_P_value3 = glme.Coefficients.pValue(2);
State1_OccurrenceRate_T_value3 = glme.Coefficients.tStat(2);    

for i = 1:NumPerm
    data = table;
    data.y = r(randperm(length(r))');
    data.group = GroupLabel;
    data.subject = Subject';
    data = rmmissing(data);
    glme = fitglme(data,'y ~ 1 + group + (1|subject)');
    perm_State1_OccurrenceRate_beta_value3(i) = glme.Coefficients.Estimate(2);
    perm_State1_OccurrenceRate_P_value3(i) = glme.Coefficients.pValue(2);
    perm_State1_OccurrenceRate_T_value3(i) = glme.Coefficients.tStat(2); 
end

if State1_OccurrenceRate_beta_value3 > 0
    perm_pVal3 = (1+length(find(perm_State1_OccurrenceRate_beta_value3 >= State1_OccurrenceRate_beta_value3)))/(NumPerm+1);
else
    perm_pVal3 = (1+length(find(perm_State1_OccurrenceRate_beta_value3 <= State1_OccurrenceRate_beta_value3)))/(NumPerm+1);
end 

% Male - Female Rational
Subject = [repmat(NoSub(1,ind_male),1,NumRationalArticle) repmat(NoSub(1,ind_female),1,NumRationalArticle)];
y7 = State1_OccurrenceRate_rational_male;
NumWin7 = NumWin_rational_male;
Age7 = Age_rational_male;
y8 = State1_OccurrenceRate_rational_female;
NumWin8 = NumWin_rational_female;
Age8 = Age_rational_female;

DependentVariable = [y8'; y7'];
Covariate = [NumWin8' Age8'; NumWin7' Age7'];
[~,~,r] = regress(DependentVariable,[ones(size(Covariate,1),1) Covariate]);
GroupLabel = [zeros(size(y8',1),1); ones(size(y7',1),1)];

data = table;
data.y = r;
data.group = GroupLabel;
data.subject = Subject';
data = rmmissing(data);
glme = fitglme(data,'y ~ 1 + group + (1|subject)');
State1_OccurrenceRate_beta_value4 = glme.Coefficients.Estimate(2);
State1_OccurrenceRate_P_value4 = glme.Coefficients.pValue(2);
State1_OccurrenceRate_T_value4 = glme.Coefficients.tStat(2);        

for i = 1:NumPerm
    data = table;
    data.y = r(randperm(length(r))');
    data.group = GroupLabel;
    data.subject = Subject';
    data = rmmissing(data);
    glme = fitglme(data,'y ~ 1 + group + (1|subject)');
    perm_State1_OccurrenceRate_beta_value4(i) = glme.Coefficients.Estimate(2);
    perm_State1_OccurrenceRate_P_value4(i) = glme.Coefficients.pValue(2);
    perm_State1_OccurrenceRate_T_value4(i) = glme.Coefficients.tStat(2);
end

if State1_OccurrenceRate_beta_value4 > 0
    perm_pVal4 = (1+length(find(perm_State1_OccurrenceRate_beta_value4 >= State1_OccurrenceRate_beta_value4)))/(NumPerm+1);
else
    perm_pVal4 = (1+length(find(perm_State1_OccurrenceRate_beta_value4 <= State1_OccurrenceRate_beta_value4)))/(NumPerm+1);
end 

% Rational - Sensible Male
Subject = [repmat(NoSub(1,ind_male),1,NumSensibleArticle) repmat(NoSub(1,ind_male),1,NumRationalArticle)];
y9 = State1_OccurrenceRate_sensible_male;
NumWin9 = NumWin_sensible_male;
Age9 = Age_sensible_male;
y10 = State1_OccurrenceRate_rational_male;
NumWin10 = NumWin_rational_male;
Age10 = Age_rational_male;

DependentVariable = [y9'; y10']; 
Covariate = [NumWin9' Age9'; NumWin10' Age10'];
[~,~,r] = regress(DependentVariable,[ones(size(Covariate,1),1) Covariate]);
GroupLabel = [zeros(size(y9',1),1); ones(size(y10',1),1)];

data = table;
data.y = r;
data.group = GroupLabel;
data.subject = Subject';
data = rmmissing(data);
glme = fitglme(data,'y ~ 1 + group + (1|subject)');
State1_OccurrenceRate_beta_value5 = glme.Coefficients.Estimate(2);
State1_OccurrenceRate_P_value5 = glme.Coefficients.pValue(2);
State1_OccurrenceRate_T_value5 = glme.Coefficients.tStat(2); 

for i = 1:NumPerm
    data = table;
    data.y = r(randperm(length(r))');
    data.group = GroupLabel;
    data.subject = Subject';
    data = rmmissing(data);
    glme = fitglme(data,'y ~ 1 + group + (1|subject)');
    perm_State1_OccurrenceRate_beta_value5(i) = glme.Coefficients.Estimate(2);
    perm_State1_OccurrenceRate_P_value5(i) = glme.Coefficients.pValue(2);
    perm_State1_OccurrenceRate_T_value5(i) = glme.Coefficients.tStat(2);
end

if State1_OccurrenceRate_beta_value5 > 0
    perm_pVal5 = (1+length(find(perm_State1_OccurrenceRate_beta_value5 >= State1_OccurrenceRate_beta_value5)))/(NumPerm+1);
else
    perm_pVal5 = (1+length(find(perm_State1_OccurrenceRate_beta_value5 <= State1_OccurrenceRate_beta_value5)))/(NumPerm+1);
end 

% Rational - Sensible Female
Subject = [repmat(NoSub(1,ind_female),1,NumSensibleArticle) repmat(NoSub(1,ind_female),1,NumRationalArticle)];
y11 = State1_OccurrenceRate_sensible_female;
NumWin11 = NumWin_sensible_female;
Age11 = Age_sensible_female;
y12 = State1_OccurrenceRate_rational_female;
NumWin12 = NumWin_rational_female;
Age12 = Age_rational_female;

DependentVariable = [y11'; y12']; 
Covariate = [NumWin11' Age11'; NumWin12' Age12'];
[~,~,r] = regress(DependentVariable,[ones(size(Covariate,1),1) Covariate]);
GroupLabel = [zeros(size(y11',1),1); ones(size(y12',1),1)];

data = table;
data.y = r;
data.group = GroupLabel;
data.subject = Subject';
data = rmmissing(data);
glme = fitglme(data,'y ~ 1 + group + (1|subject)');
State1_OccurrenceRate_beta_value6 = glme.Coefficients.Estimate(2);
State1_OccurrenceRate_P_value6 = glme.Coefficients.pValue(2);
State1_OccurrenceRate_T_value6 = glme.Coefficients.tStat(2); 

for i = 1:NumPerm
    data = table;
    data.y = r(randperm(length(r))');
    data.group = GroupLabel;
    data.subject = Subject';
    data = rmmissing(data);
    glme = fitglme(data,'y ~ 1 + group + (1|subject)');
    perm_State1_OccurrenceRate_beta_value6(i) = glme.Coefficients.Estimate(2);
    perm_State1_OccurrenceRate_P_value6(i) = glme.Coefficients.pValue(2);
    perm_State1_OccurrenceRate_T_value6(i) = glme.Coefficients.tStat(2);
end

if State1_OccurrenceRate_beta_value6 > 0
    perm_pVal6 = (1+length(find(perm_State1_OccurrenceRate_beta_value6 >= State1_OccurrenceRate_beta_value6)))/(NumPerm+1);
else
    perm_pVal6 = (1+length(find(perm_State1_OccurrenceRate_beta_value6 <= State1_OccurrenceRate_beta_value6)))/(NumPerm+1);
end 

savefile = fullfile(savefolder,'State1_OccurrenceRate_between_conditions.mat');
save(savefile);
