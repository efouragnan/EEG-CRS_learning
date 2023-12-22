function [ npar ] = get_npar( modelID)
% Lookup table to get number of free parameters per model
% MKW 2018


%%%%%
if       strcmp(modelID,'ms_RL'),            npar = 3;
elseif   strcmp(modelID,'ms_RL_cl'),         npar = 5;
elseif   strcmp(modelID,'ms_RL_cs'),         npar = 5;
elseif   strcmp(modelID,'ms_RL_rt'),         npar = 5;
elseif   strcmp(modelID,'ms_RL_cl_cs'),      npar = 7;
elseif   strcmp(modelID,'ms_RL_cl_rt'),      npar = 7;
elseif   strcmp(modelID,'ms_RL_cs_rt'),      npar = 7;   
elseif   strcmp(modelID,'ms_RL_cl_cs_rt'),   npar = 9;
   
   
   
   
%%%%%%

elseif   strcmp(modelID,'ms_RL_cl_cs_rt_unboundrt'),   npar = 9;

elseif   strcmp(modelID,'RL'),               npar = 2;    
elseif   strcmp(modelID,'RL_fo'),            npar = 3; 
elseif   strcmp(modelID,'RL_foAs'),          npar = 4;   
elseif   strcmp(modelID,'RL_re'),            npar = 3;
elseif   strcmp(modelID,'RL_ch'),            npar = 3;
elseif   strcmp(modelID,'RL_re_ch'),         npar = 4;
elseif   strcmp(modelID,'RL_re_ch_foAs'),    npar = 6;
elseif   strcmp(modelID,'RL_cl0'),           npar = 3;
elseif   strcmp(modelID,'RL_sbLr'),          npar = 3;
elseif   strcmp(modelID,'RL_cl0_sbLr'),      npar = 4;
elseif   strcmp(modelID,'RL_cl1_sbLr'),      npar = 5;
elseif   strcmp(modelID,'RL_cl1_sbLr_ma2rt'),      npar = 6;   
elseif   strcmp(modelID,'RL_sbLr_2a'),        npar = 4;
elseif   strcmp(modelID,'RL_cl1_sbLr_2a'),   npar = 6;
elseif   strcmp(modelID,'RL_cl1_sbLr_maloss'),   npar = 8;   
elseif   strcmp(modelID,'RL_cl1_sbLr_maloss2'),   npar = 7;   
elseif   strcmp(modelID,'RL_cl1_sbLr_2art'),   npar = 8;   
elseif   strcmp(modelID,'RL_cl1_sbLr_avg2a'),npar = 7;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs0'),  npar = 6;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1'),  npar = 7;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_ymaxQ'),  npar = 8;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_ma2rt'),  npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_ma2rtpn'),  npar = 10;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_ma3rt'),  npar = 9;
elseif   strcmp(modelID,'RL_cl1_sbLr_arr'),  npar = 7;
elseif   strcmp(modelID,'RL_cl1_sbLr_ma2rtpn'),  npar = 8;   
   
elseif   strcmp(modelID,'RL_cl1_sbLr_arr2'), npar = 7;
elseif   strcmp(modelID,'RL_sbLr_cs1'),      npar = 5;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs2'),  npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1ch'),            npar = 7;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1ch1'),           npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1ch12p'),         npar = 9;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1ch1av'),         npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_re'),           npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_dyna'),         npar = 8;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_reFix'),        npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_avg'),          npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_davgV1'),       npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_davgV2'),       npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_davgV3'),       npar = 7;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_davgV4'),       npar = 8;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_davgV5'),       npar = 9;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_davgV6'),       npar = 8;  
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_davgV5_ch1'),   npar = 10;    
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_davgV5_ch1a'),  npar = 10;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr'),          npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arrSeq'),          npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arrSeq2'),          npar = 9;    
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arrSeqm'),          npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arrSeqs'),          npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arrs'),          npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arrymaxQ'),          npar = 10; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr_iss'),       npar = 10;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr_is'),       npar = 11;    
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arrB1'),          npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr2'),         npar = 9;     
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr2sum'),      npar = 8;  
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr2sumR'),     npar = 7;  
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr3'),         npar = 9;     
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr4'),         npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr5'),         npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arr_ucs'),   npar = 10;     
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_arrN'),         npar = 9;     
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_2a'),           npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_2a_avg'),       npar = 9;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_2a_avg2a'),     npar = 10;   
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_avg2a'),        npar = 9;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_avg2a_ch1'),    npar = 10;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_reM'),          npar = 8;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_reXscs'),       npar = 9;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_reXccs'),       npar = 9;
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_reXscsMod'),    npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_reXccsMod'),    npar = 8; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_reXcsc2'),      npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_re2p'),         npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_re3p'),         npar = 10; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_rede'),         npar = 8; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_rede2'),        npar = 9; 
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_rede2p'),       npar = 9;    
elseif   strcmp(modelID,'RL_cl1_sbLr_cs1_rede2pM'),      npar = 9;    
elseif   strcmp(modelID,'RL_cl1'),                       npar = 4;






end




end

