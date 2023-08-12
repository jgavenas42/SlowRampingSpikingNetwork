%% Quantifying the Charactieristics of Spontaneous Fluctuations

% By quantifying the presence of fluctuations using CV, and the presence
% of bistability using the ST measure (Schaub et al., 2015), we will
% compare the effects of clustering on fluctuations in the presence of fast
% and slow synapses

%% Testing set-up
Ree_list = 1:0.1:6;
connMat_list = 1:50;
Ksyn_list = 0:0.1:0.5;

N = 400;
Nclusters = 4;
noisescale = 0.5;
Jee = 0.2;
Jclust = 1.9;
tsteps = 5000;

table_length = length(Ree_list) *length(connMat_list) * length(Ksyn_list);
tab_Ree = zeros(table_length,1);
tab_rng = zeros(table_length,1);
tab_ksyn = zeros(table_length,1);
tab_cv = zeros(table_length,1);
tab_st = zeros(table_length,1);
tab_ss = zeros(table_length,1);
tab_ff = zeros(table_length,1);

counter = 1;
tic
for ii =1:length(Ree_list)
    for jj = 1:length(connMat_list)
        for kk = 1:length(Ksyn_list)
            
            Ree = Ree_list(ii);
            connMat_rng = connMat_list(jj);
            Ksyn = Ksyn_list(kk);




            % Create the connectivity matrix
            W = connMatrix('N',N,...
                'Nclusters',Nclusters,...
                'Ree',Ree,...
                'Jee',Jee,...
                'Jclust',Jclust,...
                'connMat_rng',connMat_rng);


            [spks,Rs] = snap_singletrial(W,'tsteps',tsteps,...
                            'noisescale',noisescale,...
                            'Ksyn',Ksyn); 

            ST = calculateST(spks(1:320,501:end),Nclusters);
            CV = calculateCV(spks(1:320,501:end));
%             SS = calculateSilhouetteScore(spks(1:320,501:end),2);
%             FF = calculateFanoFactor(spks(1:320,501:end));


            tab_Ree(counter) = Ree;
            tab_rng(counter) = connMat_rng;
            tab_ksyn(counter) = Ksyn;
            tab_cv(counter) = CV;
            tab_st(counter) = ST;
            
            
            counter = counter + 1;
            toc
            
        end
    end
end


fluxQuant = table;

fluxQuant.Ree = tab_Ree;
fluxQuant.Network = tab_rng;
fluxQuant.Ksyn = tab_ksyn;
fluxQuant.CV = tab_cv;
fluxQuant.ST = tab_st;

writetable(fluxQuant,'./fluxQuant/fluxquant.csv')