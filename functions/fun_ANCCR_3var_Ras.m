function error = fun_ANCCR_3var_Ras(xpar,input_set)
eventlog = input_set{1};

% parameters
R_pn = xpar(1);
eventlog(eventlog(:,3)<0,3) = R_pn;

T = input_set{3};
samplingperiod = xpar(3);
alpha_anccr = xpar(2);
alpha_r = 0.2;
w = 0.5;
k = 1;
minimumrate = 10^(-3);
maximumjitter = 0.1;
beta = [0,0,0,1,1,0,0,0];
threshold = 0.6;
Tratio = 1.2;
omidx = [6 4; 7 5];

% estimate value like variable
[~,~,~,SRC,~,Rs] = calculateANCCR(eventlog, T*Tratio, alpha_anccr, k,...
    samplingperiod,w,threshold,minimumrate,beta,alpha_r,maximumjitter,nan,omidx);

SRC_R_data = cell(1,8);
SRC_R = SRC.*Rs;
SRC_R_data{1} = SRC_R(1,4,eventlog(:,1)>3); % CSrw-pn & Rw % outcome
SRC_R_data{2} = SRC_R(1,5,eventlog(:,1)>3); % CSrw-pn & Pn
SRC_R_data{3} = SRC_R(2,4,eventlog(:,1)>3); % CSpn-rw & Rw
SRC_R_data{4} = SRC_R(2,5,eventlog(:,1)>3); % CSpn-rw & Pn
SRC_R_data{5} = SRC_R(1,4,eventlog(:,1)<4); % CSrw-pn & Rw % cue
SRC_R_data{6} = SRC_R(1,5,eventlog(:,1)<4); % CSrw-pn & Pn
SRC_R_data{7} = SRC_R(2,4,eventlog(:,1)<4); % CSpn-rw & Rw
SRC_R_data{8} = SRC_R(2,5,eventlog(:,1)<4); % CSpn-rw & Pn
SRC_R_data = cellfun(@squeeze,SRC_R_data,'uniformoutput',false);

temp_data_csrwpn_rw = SRC_R_data{5};
temp_data_csrwpn_pn = SRC_R_data{6};
temp_data_cspnrw_rw = SRC_R_data{7};
temp_data_cspnrw_pn = SRC_R_data{8};

if R_pn<0
    temp_data_csrwpn = temp_data_csrwpn_rw+temp_data_csrwpn_pn;
    temp_data_cspnrw = temp_data_cspnrw_rw+temp_data_cspnrw_pn;
else
    temp_data_csrwpn = temp_data_csrwpn_rw-temp_data_csrwpn_pn;
    temp_data_cspnrw = temp_data_cspnrw_rw-temp_data_cspnrw_pn;
end

temp_data_cs_set = {temp_data_csrwpn temp_data_cspnrw};

% compare with lick
delay_licknum = input_set{2};

cs_idx = eventlog(eventlog(:,1)<4,1);
r_set = zeros(1,2);
for idata = 1:2
    temp_lick = delay_licknum;
    temp_lick(cs_idx~=idata)=nan;
    temp_lick = movmean(temp_lick,40,'omitnan');
    temp_lick = temp_lick(cs_idx==idata);
    
    temp_srcr = temp_data_cs_set{idata};
    temp_srcr(cs_idx~=idata)=nan;
    temp_srcr = movmean(temp_srcr,40,'omitnan');
    temp_srcr = temp_srcr(cs_idx==idata);
    
    [r,~] = corrcoef(temp_lick,temp_srcr);
    r_set(idata) = r(1,2);
end

error = -mean(r_set);