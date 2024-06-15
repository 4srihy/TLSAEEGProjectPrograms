%% fig7_v2_without iterations

clear;
load('stimFeaturesHP.mat')
numIter = 1; % no. of iterations of auc


% centre freq 7,31,55 with freqwidth 12 for hfd spec power 
features = [powerAlpha_median powerSG_median powerFG_median exponent_LFR_median exponent_HFR_median  HFDAlpha HFDSG HFDFG hfd_LFR hfd_HFR];

features(isnan(features)) = 0.01;
featureNames = [{'P_{\alpha}'},{'P_{SG}'},{'P_{FG}'},{'Slope_{LFR}'},{'Slope_{HFR}'},{'HFD_{\alpha}'},{'HFD_{SG}'},{'HFD_{FG}'},{'HFD_{LFR}'},{'HFD_{HFR}'}];

% zn = zscore(features,0,1);
new_features = features;

figure;

%% auc
    % 5-fold cross validation sets
    k = 5; % "k" of kfold cross-validation
    subjParts =cvpartition(group,'kfold',k);    % stratified by default

    % Single feature AUC
    nfeats = size(new_features,2);
    auc = zeros(k,nfeats+3);

    for sel_feat = 1:nfeats+3
        classSubj = cell(1,k);
        cp = cell(1,k); ratio_self = zeros(1,k);

        post_p = cell(1,k);
        
        if sel_feat<=nfeats
            featIndex = sel_feat;
        elseif sel_feat==nfeats+1
            featIndex = 1:5;
        elseif sel_feat==nfeats+2
            featIndex = 6:10;
        else
            featIndex = 1:10;
        end

        features_considered = new_features(:,featIndex);
        numFeatures = size(features_considered,2);
        Mdl = fitcdiscr(features_considered,group);

        %auc = zeros(1,k);
        for i=1:k
            testIDs = test(subjParts,i);
            labels = group(testIDs);
            [~,scores,~] = predict(Mdl,features_considered(testIDs,:));
            [~,~,~,auc(i,sel_feat)] = perfcurve(labels,scores(:,2),1);

        end

        % m_acc(sel_feat) = mean(acc);  %accuracy
        % m_ratio_self(sel_feat) = (mean(ratio_self));    %sensitivity
        % m_ratio_cross(sel_feat) = mean(ratio_cross);    %specificity

    end
SingleFeatAUC= mean(auc,1);            %auc
SingleFeatAUCStd= std(auc,1); 

% Visualizing result
subplot(1,1,1)
x = 1:nfeats+3;
y = SingleFeatAUC;
bar(x,y); hold on;
err = SingleFeatAUCStd/sqrt(k);
%er = errorbar(x,y,err,err,'LineStyle','none','LineWidth',1,'Color',[0 0 0]);
xticks(1:nfeats+3);

featureNames = [featureNames, {'Spec Features'},{'HFD Features'},{'All Features'}];
xticklabels(featureNames);
ylim([0.5 0.75]);
ylabel('AUC','Fontweight','bold');
for i=1:length(y)
text(i-0.2,y(i)+0.012,num2str(round(y(i),2)),'FontSize',13,'FontWeight','bold');
end
set(gca,'FontWeight','bold','FontSize',12,'Box','off');





