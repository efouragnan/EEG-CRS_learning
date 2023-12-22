function [ qout,qbounds ] = norm2par(modelID,qin)
% Lookup table to transform gaussian parameters to real parameter values for all models;
%  
% columns are different free parameters, have to be in the order defined by the model
% MKW 2017
%
%
% INPUT:    - qin: input parameters, observation x parameter
% OUTPUT:   - qout: qin, transformed
%           - qbounds: sensible bounds for the parameter, rows are min/max, cols are parameters

%% in case input is a vector

if size(qin,2)==1, qin = qin'; end                                           % transpose input if it has wrong dimension, but only if it is a n x 1 vector


% do checking
if       strcmp(modelID,'ms_RL')   
   if size(qin,2)~=3, disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2sblr(qin(:,3))];  
   qbounds  = [0 1; 0 15; -.7 .7]';
   
elseif       strcmp(modelID,'ms_RL_cl')
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5))];
   qbounds = [0 1; 0 15; -.7 .7; -.7 .7; 0 1]';
      
elseif       strcmp(modelID,'ms_RL_cs')
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2sblr(qin(:,3)) (qin(:,4)) norm2alpha(qin(:,5))];
   qbounds = [0 1; 0 15; -.7 .7; -1 1; 0 1]';   
   
elseif       strcmp(modelID,'ms_RL_rt')
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2sblr(qin(:,3)) norm2sblr(qin(:,4)) norm2alpha(qin(:,5))];
   qbounds = [0 1; 0 15; -.7 .7; -2 2; 0 1]';        
   
elseif       strcmp(modelID,'ms_RL_cl_cs')
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7))];
   qbounds = [0 1; 0 15; -.7 .7; -.7 .7; 0 1; -1 1;0 1]';
   
elseif       strcmp(modelID,'ms_RL_cl_rt')
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) norm2sblr(qin(:,6)) norm2alpha(qin(:,7))];
   qbounds = [0 1; 0 15; -.7 .7; -.7 .7; 0 1; -1 1;0 1]';   
   
 elseif       strcmp(modelID,'ms_RL_cs_rt')
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2sblr(qin(:,3)) (qin(:,4)) norm2alpha(qin(:,5)) norm2sblr(qin(:,6)) norm2alpha(qin(:,7))];
   qbounds = [0 1; 0 15; -.7 .7; -1 1;0 1; -2 2; 0 1]';  
   
   
elseif       strcmp(modelID,'ms_RL_cl_cs_rt')
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2sblr(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 15; -1 1; -1 1; 0 1; -1.5 1.5;0 1; -1 1; 0 1]';



   
   
   
   
   
   
   
   
   
%% all other models below

elseif       strcmp(modelID,'ms_RL_cl_cs_rt_unboundrt')
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) qin(:,8) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -1 1;0 1; -2 2; 0 1]';




elseif       strcmp(modelID,'RL'),     
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2))];
   qbounds  = [0 1; 0 20]';
   
elseif       strcmp(modelID,'RL_cl0'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3)];
   qbounds  = [0 1; 0 20; -.7 .7]';  
   
elseif       strcmp(modelID,'RL_cl0_sbLr'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7]';

elseif       strcmp(modelID,'RL_cl1_sbLr'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1]';   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_maloss'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) norm2alpha(qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1;0 1; 0 1; 0 1]';   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_maloss2'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) norm2alpha(qin(:,6)) norm2alpha(qin(:,7)) ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1;0 1; 0 1]';    
   
elseif       strcmp(modelID,'RL_cl1_sbLr_ma2rt'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) norm2alpha(qin(:,6))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; 0 1]';  
   
elseif       strcmp(modelID,'RL_cl1_sbLr_ma2rtpn'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) norm2alpha(qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; 0 1; 0 1; 0 1]';     
   
elseif       strcmp(modelID,'RL_sbLr_2a'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2))  norm2sblr(qin(:,3))  norm2alpha(qin(:,4))];
   qbounds = [0 1; 0 20;  -.7 .7;0 1]';       
   
elseif       strcmp(modelID,'RL_cl1_sbLr_2a'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) norm2alpha(qin(:,6))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1;0 1]';   
      
elseif       strcmp(modelID,'RL_cl1_sbLr_2art'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) norm2alpha(qin(:,6)) norm2alpha(qin(:,7))   qin(:,8)];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1;0 1;0 1;-1 1]';     
   
elseif       strcmp(modelID,'RL_cl1_sbLr_avg2a'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) norm2alpha(qin(:,6)) norm2alpha(qin(:,7))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; 0 1; 0 1]';     
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs0'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7]';    
      
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1]'; 
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_dyna'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2beta(qin(:,8))+.00001];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 20]';   
   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_ymaxQ'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1;0 1]';    
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_ma2rt'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';  
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_ma2rtpn'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9)) norm2alpha(qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;0 1; 0 1]';    
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_ma3rt'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]';     
   
elseif       strcmp(modelID,'RL_cl1_sbLr_arr'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -1.7 1.7;0 1]';  
     
   
elseif       strcmp(modelID,'RL_cl1_sbLr_arr2'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -1.7 1.7;0 1]';     
   
   
elseif       strcmp(modelID,'RL_sbLr_cs1'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1]';  

   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs2'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';       
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1ch'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1]';   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1ch1av'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1;0 1]';      

elseif       strcmp(modelID,'RL_cl1_sbLr_cs1ch1'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1ch12p'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]';   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_re'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';  
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_davgV1'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';  
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_davgV2'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';   
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_davgV3'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1]';  
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_davgV4'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) abs(qin(:,8))+6 ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 200]';   
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_davgV5'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1]';  
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_davgV6'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7))  norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1;  0 1]';     
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_davgV5_ch1'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9)) norm2alpha(qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1;0 1]';      
 
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_davgV5_ch1a'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9)) norm2alpha(qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1;0 1]';
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1]'; 
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arrymaxQ'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9)) norm2alpha(qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1;0 1]';    
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr4'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2sblr(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1]';  
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr5'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2sblr(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1]';     
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arrSeq'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2sblr(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]';    
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arrSeq2'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2sblr(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]';    
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arrSeqm'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2sblr(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]';  

elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arrSeqs'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2sblr(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]'; 
 
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arrs'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1]'; 
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr_iss'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9))  (qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1; -1 1]'; 
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr_is'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9)) norm2alpha(qin(:,10)) (qin(:,11))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1; 0 1; -1 1]'; 
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arrB1'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1]';   
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr2'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1]';

elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr2sum'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2]';
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr2sumR'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7))  ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1]';
   
      
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr3'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]';    
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arr_ucs'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9)) (qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1;-.3 .3]';      
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_arrN'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) (qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; -2 2; 0 1]';    
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_reFix'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';    
 
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_avg'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';    
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_avg2a'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))  norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;0 1]'; 
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_2a'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';    
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_2a_avg'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]';    
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_2a_avg2a'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9)) norm2alpha(qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1; 0 1]';     
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_avg2a_ch1'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8))  norm2alpha(qin(:,9)) norm2alpha(qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;0 1;0 1]';    
   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_reM'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2beta(qin(:,8))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';  
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_rede2pM'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2beta(qin(:,8)) norm2beta(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;0 1]';     
   
         
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_re2p'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1]';  
      
 
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_re3p'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9)) norm2alpha(qin(:,10))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1; 0 1; 0 1]'; 
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_rede'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_rede2'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) exp(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;-2 2]';    
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_rede2p'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) ((qin(:,9)))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;-2 6]';      
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_reXscs'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) qin(:,9)];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;1 2]'; 
   
 elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_reXscsMod'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) qin(:,9)];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;1 2]';   
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_reXccs'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) qin(:,9)];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;1 2]';  
   
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_reXcsc2'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) norm2alpha(qin(:,9))];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1;0 3]';   
                           
elseif       strcmp(modelID,'RL_cl1_sbLr_cs1_reXccsMod'),    
   if size(qin,2)~=get_npar(modelID), disp('ERROR'); keyboard; end; 
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) qin(:,3) norm2sblr(qin(:,4)) norm2alpha(qin(:,5)) (qin(:,6)) norm2alpha(qin(:,7)) norm2alpha(qin(:,8)) ];
   qbounds = [0 1; 0 20; -.7 .7; -.7 .7; 0 1; -.7 .7;0 1; 0 1]';      
   
elseif       strcmp(modelID,'RL_sbLr'),   
   if size(qin,2)~=3, disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2sblr(qin(:,3))];  
   qbounds  = [0 1; 0 20; -.7 .7]';
               
elseif   strcmp(modelID,'RL_fo'),     
   if size(qin,2)~=3, disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3))];
   qbounds  = [0 1; 0 40; 0 1]';  
   
elseif   strcmp(modelID,'RL_foAs'),   
   if size(qin,2)~=4, disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3)) norm2alpha(qin(:,4))];
   qbounds  = [0 1; 0 40; 0 1; 0 1]';
   
elseif   strcmp(modelID,'RL_re'),  
   if size(qin,2)~=3, disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3)) ];
   qbounds  = [0 1; 0 40; 0 1]';

elseif   strcmp(modelID,'RL_ch'),  
   if size(qin,2)~=3, disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3)) ];
   qbounds  = [0 1; 0 40; 0 1]';

elseif   strcmp(modelID,'RL_re_ch'),  
   if size(qin,2)~=4, disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3)) norm2alpha(qin(:,4))];
   qbounds  = [0 1; 0 40; 0 1; 0 1]';
      
elseif   strcmp(modelID,'RL_re_ch_foAs'), 
   if size(qin,2)~=6, disp('ERROR'); keyboard; end;
   qout     = [norm2alpha(qin(:,1)) norm2beta(qin(:,2)) norm2alpha(qin(:,3))  norm2alpha(qin(:,4))  norm2alpha(qin(:,5))  norm2alpha(qin(:,6)) ];
   qbounds  = [0 1; 0 40; 0 1; 0 1; 0 1; 0 1]';
         
end;






end

