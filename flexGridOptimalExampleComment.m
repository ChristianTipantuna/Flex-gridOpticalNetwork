% -------------------------------------------------------------------------
%Universidad Politecnica de Catalunya (UPC)
%Doctoral Program in Computer Networks
%Author: Christian Tipantuna
%Advisor: Dr. Xavier Hesselbach
%Date: April 2020
%Description: Optimal Spectrum Allocation in Flex-grid optical
%communications considering time shifting capabilities
%This wavelength scheduling provides the optimal result and it has been
%proven be NP-Hard.
% -------------------------------------------------------------------------

%--- System parameters:
%   Central frequency of the system         freqCentral = 193100 GHz
%                                           Recommendation ITU-T G.694.1
%   Power reference associate with lambdas  PwRefLevel>0 (PwRefLevel = 1)
%   Number of services (demands or jobs):   n, k = 1,...,n
%   Service identifier:                     sk
%   Priority identifier:                    lk = 1
%   Power demanded by lambdas:              uk = 1
%   Central frequency of services:          fcInitk
%   Channel bandwidth:                      bk = 12.5, 25, 37.5, 50, 100; or wider (integer multiples of 100 GHz)
%   Lifetime duration of demands:           dk = 1,2,3,4
%   Frecuency shifting:                     fsBackward, fsForward = 0,1,2,3,4


%Services: Service identifier, priority, power demand, central frequency, bandwidth channel,
%frequency shifting
%Fields: sk, lk, uk, fcInitk, bk, fsBackward, fsForward

%--- Performance metrics:
%Acceptance ratio:                      AR
%Bandwidth utilization:                 BWUse
%Number of changes of frequency:        NumFreqChange
%Cumulative frequency shifting value:   CumFreShifting

% -------------------------------------------------------------------------

clc, clear, close all;

tic;

%--- Initialization of arrays to store variables
caseInfo = []; %information of case under analysis
frecShiftingInfo = []; %value of frequency shifting of services
acceptanceRatioInfo = []; %value of acceptance ratio of processed services
useAvaBWInfo = []; %used bandwidth in GHz
useAvaBWpercentInfo = []; %utilization percentage of available bandwidth
numFreqShiftinInfo = []; %number of times frequency shifting has been used
cumFreqShiftinInfo = []; %cummulative value of the frequency shifting

%--------------------------------------------------------------------------
% Parameters of the system
%--------------------------------------------------------------------------

for iterCase = 1:1 %iterations for the case under analysis
    
    for iterFs = 0:3 %iterations for Fs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Fs: frequecy shifting
        
        %Setting seed values to perform a comparison with the sub-optimal
        %solution
        rng default;
        auxSeed = 100 + iterCase;
        rng(auxSeed);
        
        %Option for computing the acceptance ratio
        infoOptionAR = 1;
        
        %--- Central frequency of the system
        freqCentral = 193100; %GHZ, based on Recommendation ITU-T G.694.1 for dense DWDM networks
        
        %--- Frequency slot width
        fsWidth = 6.25; %GHz, based on Recommendation ITU-T G.694.1
        
        %--- Power level associate with the transmission of lambdas
        PowerReference = 1; %reference power level
        
        %--- Number of lambdas/jobs
        n = 8;
        
        %--------------------------------------------------------------------------
        % Parameters of services (lambdas)
        %--------------------------------------------------------------------------
        
        %Services: Service identifier, priority, power demand, central frequency, bandwidth channel,
        %frequency shifting
        %Fields: sk, lk, uk, fcInitk, bk, fsBackward, fsForward
        services = zeros(n,7); %7 number of parameters
        
        
        %--- Number of services: sk
        for i = 1:n
            services(i,1) = i;
        end
        
        %--- Priority of the service: lk
        %In case of the system uses prioritization of lambdas
        for i = 1:n
            services(i,2) = 1;
        end
        
        %--- Power demand by services, reference power: uk
        %The power demanded can be different for each service, however, the
        %idea in this work is optimize the use of available bandwidth, so this
        %power level is only used as a reference
        for i = 1:n
            services(i,3) = 1; %reference power level
        end
        
        
        %--- Bandwidth channel: bk
        % Bandwidth channels = 12.5 × m, where m is a positive integer and 12.5 is the slot width granularity in GHz
        serviceChannel = [50 12.5 12.5 25 12.5 25 37.5 50]; %Channels bandwidth example
        for i = 1:n
            %channels from 12.5 Ghz till 100 GHz, wider channels change 8
            %by a larger number
            %services(i,5) = 12.5 * randi([1 8],1,1);
            services(i,5) = serviceChannel(i);
        end
%         auxChannelBW = sum(services(:,5),1); %auxiliar variable for computing the central frequency
%         auxChannelBWMaxMin = auxChannelBW/2; 
%         auxChannelBWFs = auxChannelBW/fsWidth;
%         auxmaxNumChannel = auxChannelBWFs/2; %maximum number of frequency slots below or above the cetral frequency
%         auxMinChannel = min(min(services(:,5))); %information of the minimum BW of services
%         auxMinChannelFs = auxMinChannel/fsWidth;
%         numChannels = auxmaxNumChannel - auxMinChannelFs; %maximum frequency slot for central frequencies
%         auxFmin = freqCentral - auxChannelBWMaxMin; %minimum range in the frequency horizon of analysis
%         auxFmax = freqCentral + auxChannelBWMaxMin; %maximum range in the frequency horizon of analysis
        
        %--- Central frequency of services (lambdas): FcInitk
        %Central frequency services: 193100 + n * 6.25,
        servicesFc = [193018.75 193087.5 193118.75 193200 192993.75... %Central frequencies example
            193043.75 193106.25 193181.25];
        for i = 1:n
            %where n is a positive or negative integer including 0
            services(i,4) = servicesFc(i);
            %services(i,4) = freqCentral + randi([-n n],1,1) * fsWidth;
            %services(i,4) = freqCentral + randi([-numChannels numChannels],1,1) * fsWidth; %serviceFc(i);
            
%             temCentralFrew = 1;
%             
%             while temCentralFrew == 1
%                 
%                 auxCentralFreq = freqCentral + randi([-numChannels numChannels],1,1) * fsWidth;
%                 
%                 if (auxCentralFreq - (services(i,5)/2)) >= auxFmin && (auxCentralFreq + (services(i,5)/2)) <= auxFmax
%                     disp('Valid central frequency');
%                     break;
%                 else
%                     disp('Not valid central frequency');
%                 end
%             end
%      
%             services(i,4) = auxCentralFreq;
            
        end
        
        
        %--- Frequency Shifting Forward: fFs (in frequency slots, 1 slot = 6.25 GHz)
        for i = 1:n
            services(i,6) = iterFs; %0 means no frequency shifting performed
        end
        
        %--- Frequency Shifting Backward: bFs (in frequency slots, 1 slot = 6.25 GHz)
        for i = 1:n
            services(i,7) = iterFs; %0 means no frequency shifting performed
        end
        %disp(services); %parameters in terms of GHz
        
        
        %--- Theoretical total bandwidth (Available spectrum)
        BW =  sum(services(:,5)); %measurement given in GHz
        
        %-Badwidth to the left of the reference point (193100)
        bwLow = freqCentral - (BW/2); %lower bound
        
        %-Badwidth to the right of the reference point (193100)
        bwHigh = freqCentral + (BW/2); %higher bound
        
        
        %--------------------------------------------------------------------------
        % Information expressed as Frequency Slots
        %--------------------------------------------------------------------------
        %- Central frequency of the system
        freqCentralFs = 193100/fsWidth; %central frequency of the system
        
        %Total available bandwidth (available spectrum)
        BWFs = BW/fsWidth;
        
        %-Lower and higher frequencies of the system
        
        %Lower badwidth value
        bwLowFs = bwLow/fsWidth;
        
        %Higher badwidth value
        bwHighFs = bwHigh/fsWidth;
        
        
        %--- Services information considering the minimum granularity of frequency
        %slots
        servicesFs = zeros(n,7); %7 number of parameter
        
        servicesFs(:,1) = services(:,1); %sk: Service identifier
        servicesFs(:,2) = services(:,2); %lk: Service priority
        servicesFs(:,3) = services(:,3); %uk: Power demanded by services
        servicesFs(:,4) = services(:,4)/fsWidth; %fcInitk: Central frequency of services
        servicesFs(:,5) = services(:,5)/fsWidth; %bk: Channel bandwidth
        servicesFs(:,6) = services(:,6); %fsBackward: Backward frequency shifting
        servicesFs(:,7) = services(:,7); %fsForward: Forward frequency shifting
        
        
        %--------------------------------------------------------------------------
        %Description of the input parameters of the system
        %--------------------------------------------------------------------------
%         disp('------ Information of the system -------');
%         disp(['Central frequencies DWDM systems: ',num2str(freqCentral),' GHz']);
%         disp(['Frequency slot granularity: ',num2str(fsWidth),' GHz']);
%         disp (' ');
%         disp('---- Parameters of services (wavelengths):');
%         disp(['Number of services: ',num2str(n)]);
%         disp(['Maximum available bandwidth: ',num2str(BW),' GHz']);
%         disp(['Lower band frequency: ',num2str(bwLow),' GHz']);
%         disp(['Higher band frequency: ',num2str(bwHigh),' GHz']);
%         disp (' ');
%         disp('---- Parameters of services in frequency slots:');
%         disp(['Maximum available bandwidth in freq. slots: ',num2str(BWFs),' slots']);
%         disp(['Lower band frequency: ',num2str(bwLowFs),' slots']);
%         disp(['Higher band frequency: ',num2str(bwHighFs),' slots']);
%         disp (' ');
%         disp('- Information of services:');
%         disp('Service(Sk), Priority(lk), Power demanded(uk), Central frequency channel (fck), Channel bandwidth (bk)');
%         disp('Backward frequency shifting (fsBackward), Forward frequency shifting fsForward');
%         disp(' ');
%         disp('TABLE I: Description of services');
%         disp('--------------------------------');
%         disp('           Sk         lk           uk        fck          bk     fsBackward  fsForward ');
%         disp(servicesFs);
        
        %--- Initial and final frequencies of services
        %- Initial frequencies
        initFreqServ = servicesFs(:,4) - servicesFs(:,5)/2;
        
        %- Final frequencies
        finalFreqServ = servicesFs(:,4) + servicesFs(:,5)/2;
        
        
        %--- Information about different frequencies of services
        %Fields: Sk, initFk, fcInitk, finalFk, bk
        freqServInfo = zeros(n,9);
        
        freqServInfo(:,1) = servicesFs(:,1); %Sk
        freqServInfo(:,2) = servicesFs(:,2); %lk
        freqServInfo(:,3) = servicesFs(:,3); %uk
        freqServInfo(:,4) = initFreqServ; %initFk
        freqServInfo(:,5) = servicesFs(:,4); %fcInitk
        freqServInfo(:,6) = finalFreqServ; %finalFk
        freqServInfo(:,7) = servicesFs(:,5); %bk
        freqServInfo(:,8) = servicesFs(:,6); %bFs
        freqServInfo(:,9) = servicesFs(:,7); %fFs
        
        
%         disp('TABLE II: Description of services (Frequency slots)');
%         disp('--------------------------------');
%         disp('           Sk         lk           uk      initFk      fck         finalFk        bk   fsBackward  fsForward');
%         disp(freqServInfo);
        
        %-Bandiwth of services
        bkInfo = freqServInfo(:,7);
        
        %-Lower and higher frequencies of services (wavelengths)
        
        %Lower frequency value
        minFreqServ = min(min(freqServInfo(:,4)));
        
        %Higher frequency value
        maxFreqServ = max(max(freqServInfo(:,6)));
        
        
        %--- Minimum and maximum frequencies in the system
        %Reference values to perform the spectrum allocation
        
        %Minimum frequency of reference in the system
        minFreqRef = min(bwLowFs,minFreqServ);
        
        %Maximum frequency of reference in the system
        maxFreqRef = max(bwHighFs,maxFreqServ);
        
        
        %-Frequency horizon
        %Maximum frequency range for analyzing the system
        maxBW = maxFreqRef - minFreqRef;
        
        
        %--------------------------------------------------------------------------
        % Initial distribution of services (wavelengths)
        %--------------------------------------------------------------------------
        
        %---Information of available bandwidth (BW) within the maximum
        %frequency horizon
        
        %BW refereced at frequency slot 0
        availableBWFs = zeros(1,maxBW);
        bwInitFreq = (bwLowFs - minFreqRef) + 1; %1 for avoiding staring on 0 value
        availableBWFs(1,bwInitFreq:(bwInitFreq+BWFs)-1) = 1;
        PwRefLevel = availableBWFs;
        
        %---Information of services (wavelengthd) within the maximum
        %frequency horizon
        initDistributionServ = zeros(n,maxBW);
        
        refInitFreqServ = zeros(1,n); %initial frequencies of services referenced at frequency slot 0
        for i = 1:n
            %Reference frequency slot according the minimum reference in
            %the system
            minFreServ = (initFreqServ(i)- minFreqRef) + 1;
            refInitFreqServ(1,i) = minFreServ;
            %Initial service distribution referenced at frequency slot 0
            initDistributionServ(i,minFreServ:(minFreServ+servicesFs(i,5))-1) = i;
        end
        
        
        %- Figures: Available bandwidth (BW) and initial service distribution
        %(wavelengths)
%         figure('name','Available bandwidth (BW) and initial distribution of services');
%         
%         subplot(2,1,1) %Available bandwidth (BW)
%         imagesc(availableBWFs);
%         %imagesc([minFreqRef maxFreqRef],1,availableBWFs)
%         title(['Available bandwidth (BW): ',num2str(BWFs), ' frequency slots']);
%         xlabel('Frequency slots');
%         ylabel('Available bandwidth (BW)');
%         grid on;
%         
%         subplot(2,1,2) %Service distribution (wavelengths)
%         %imagesc([minFreqRef maxFreqRef],[1 8],availableBWFs)
%         imagesc(initDistributionServ);
%         title('Initial service distribution (wavelengths)');
%         xlabel('Frequency slots');
%         ylabel('Services');
%         grid on;
        
        
        
        %--------------------------------------------------------------------------
        %Information of services (wavelengths) considering the frequency
        %shifting
        %--------------------------------------------------------------------------
        
        %--- Information about the services and specifically about the variations
        %of services considering the frequency shifting interval
        %numVar: Number of services
        %idVar: Identification of variation within all variations
        %varInfo: Information of the service to which the variation belongs to
        %priorVarInfo: Priority information of variation
        %mFsVarInfo: Maximum frequency shifting of variation
        %powVarServ: Power demanded per variation and per frequency slot
        
        [numVar,idVar,varInfo,priorVarInfo,mFsVarInfo,...
            powVarServ] = funcVarServ(n,maxBW,freqServInfo,refInitFreqServ);
        
%         disp(' ');
%         disp('-----------------------------------------------');
        %disp('-- Valid combinations per service and valid combinations among services --');
%         disp(['Number of valid variations of services: ', num2str(numVar)]);
%         disp(' ');
%         disp('TABLE II: Valid Variations per Service');
%         disp('---------------- Variations per service ---------------------');
%         disp('- Number - Sk - lk - mts - Power demanded per time slot)');
        %disp([idVar', varInfo', priorVarInfo', mFsVarInfo', powVarServ]);
%         disp([idVar', varInfo', priorVarInfo', mFsVarInfo']);
        %imagesc(powVarServ)
        %grid on;
        
        
        %--------------------------------------------------------------------------
        %Partial costs of variations per services
        %--------------------------------------------------------------------------
        %Cost function: Linear cost function ck = car + cl + cmts
        
        %--- Cost related to the priority level: cl
        costLk = priorVarInfo'; %the corresponds to priority level
        
        %--- Cost related to the time shifting performed: cmts
        costMk = mFsVarInfo'; %the cost corresponds to the time shifting value
        
        
        %--------------------------------------------------------------------------
        % Computation of combinations between variations
        %--------------------------------------------------------------------------
        
        %--- Combinations of variations of services
        [numRowCombVar,numColCombVar,infCombVar,combVar] = funcCombVarSer(varInfo);
        
%         disp(' ');
%         disp('------------------------------------------------------');
%         if length(unique(varInfo)) < n
%             disp('There is no valid combination ');
%         else
%             disp(['Number of combinations of variations of services: ', num2str(numRowCombVar)]);
%         end
%         
%         disp('TABLE III: Combinations of variations of services');
%         disp(' ------ Combinations of variations of services -------');
%         disp([infCombVar',combVar]);
        
        
        %--------------------------------------------------------------------------
        % Computation of power consumed per combinations of variations of services
        %--------------------------------------------------------------------------
        
        %--- Consumed power of variations of services
        consumedPowerComb = zeros(numRowCombVar,maxBW); %
        
        %--- Cost related to the priority level (cl) of combinations of variations
        costLkCombVar = zeros(numRowCombVar,1);
        
        %--- Cost related to the time shifting performed (cmts) of combinations of variations
        costMkCombVar = zeros(numRowCombVar,1);
        
        for i = 1:numRowCombVar % i number of combinations
            %consumedPower(i,:) = demandPower(combinations(i,1),:) + ...
            %    demandPower(combinations(i,2),:) + demandPower(combinations(i,3),:);
            
            aux = zeros(1,maxBW); %auxiliar variable to store the current power of each variation
            auxCostLk = 0; %auxiliar variable to store the priority information of a variation
            auxCostMts = 0; %auxiliar variable to store the time shifting of a variation
            consumedPowerComb(i,:) = zeros(1,maxBW); %power consumer per combination of variations
            costLkCombVar(i,1) = 0; %cost of priority per combination of variations
            costMkCombVar(i,1) = 0; %cost of time shifting per combination of variations
            
            for j = 1:numColCombVar %i number of columns of combinations
                aux = powVarServ(combVar(i,j),:); %current value of a variation
                auxCostLk = costLk(combVar(i,j),1); %current value of priority
                auxCostMts = costMk(combVar(i,j),1); %current value of time shifting
                consumedPowerComb(i,:) = consumedPowerComb(i,:) + aux;
                costLkCombVar(i,1) = costLkCombVar(i,1) + auxCostLk;
                costMkCombVar(i,1) = costMkCombVar(i,1) + auxCostMts;
            end
            
        end
        
%         disp(' ');
%         disp('-------------------------------------------------------------');
%         disp('Power demanded by combinations and cost of lk and cost of mTs');
%         disp(' ');
%         disp('TABLE IV: Power Consumed of Combinations of variations of services');
%         disp(' ---------- Power consumed of combination of variations ----------');
%         disp('- NumComb - Consumed Power per Frequency Slot -------------------------');
        %disp([infCombVar',consumedPowerComb]);
%         disp(infCombVar');
%         
%         disp(' ');
%         disp('TABLE V: Cost related to priority of Combinations of variations of services');
%         disp(' -------------- Cost of priority of combination of variations -------------');
%         disp('- NumComb - costLk --------------------------------------------------------');
%         disp([infCombVar',costLkCombVar]);
        
%         disp(' ');
%         disp('TABLE VI: Cost related to time shifting of Combinations of variations of services');
%         disp(' -------------- Cost of time shifting of combination of variations --------------');
%         disp('- NumComb - costMk --------------------------------------------------------------');
%         disp([infCombVar',costMkCombVar]);
        
        
        %--------------------------------------------------------------------------
        % Computation of residual power of combinations of variations of services
        %--------------------------------------------------------------------------
        
        resPowCombVar = zeros(numRowCombVar,maxBW);
        
        for i = 1:numRowCombVar % i number of combinations
            resPowCombVar(i,:) = PwRefLevel - consumedPowerComb(i,:);
        end
        
%         disp(' ');
%         disp('TABLE VII: Residual Power of Combinations of Variations of Services');
%         disp('--- Residual power of Combinations of variations ------------------');
%         disp('Residual power = PES - Consumed Power Combination of Variations');
%         disp('- NumComb - Residual Power per Time Slot -----------------------------------');
        %disp(residualPower);
        %disp([infCombVar',resPowCombVar]);
        
        
        %--------------------------------------------------------------------------
        % Computation of acceptance ratio and the related cost of combinations
        %--------------------------------------------------------------------------
        
        %For the computation of this metric, we have three option
        %1) AR computation that firstly maximizes the use of available spectrum and
        %secondly ensure the highest possible AR value
        %2) AR that considers the discarding of services with higher bandwidth
        %3) AR that considers the discarding of services with lower bandwidth
        %Option: 1,2,3
        %optionAR = 1;
        optionAR = infoOptionAR;
        %--- AR, cost related and power demanded by the processed services
        [arCombVar,costARCombVar,powCombVar,comVarServ] = funcARCostCombVarSer(numRowCombVar,...
            resPowCombVar,combVar,n,maxBW,PwRefLevel,BWFs,priorVarInfo,...
            powVarServ,consumedPowerComb,bkInfo,optionAR);
        
        
%         disp(' ');
%         disp('TABLE VIII: Acceptance Ratio, Cost related to rejection of service, and Real power consumed combinations ');
%         disp(' --------------------------------- Cost of acceptance ratio --------------------------------------------');
%         disp('- NumComb -     AR     -    costAR -     Power demanded per Combination of variations per time slot');
        
        %fprintf('%10d      %s     \t %s  \t \t %8.3f            \t %s\n',numValComb(i),num2str(validInfoComb(i,:)),...
        %num2str(sigmaRes(i),'% .3f'),acceptanceRatio(i), num2str(sigmaShift(i),'% .3f'));
        
%         for i = 1:numRowCombVar
%             fprintf('%10d   %8.3f  \t  %8.3f \t     %s\n',infCombVar(i),arCombVar(i),...
%                 costARCombVar(i), num2str(powCombVar(i,:)));
%         end
        
        
        %--------------------------------------------------------------------------
        %Total cost of Combination of variations of services
        %--------------------------------------------------------------------------
        %Cost function: Linear cost function ck = car + cl + cmts
        
        totCostComVar = costLkCombVar + costMkCombVar + costARCombVar;
        
        
        %--------------------------------------------------------------------------
        % Summarized information of combinations of variations of services
        %--------------------------------------------------------------------------
        
%         disp(' ');
%         disp('TABLE IX: Summarized information of Combinations of Variations');
%         disp(' --------------------------------- Cost of acceptance ratio -----------------------------------------------------------');
%         disp('- NumComb -  Comb. Var.  -  costLk  -  costMk  - CostAR    -  totCost -   AR -     Power demanded per Combination of variations');
        
%         for i = 1:numRowCombVar
%             fprintf('%10d    %s %10d %10d   %8.3f     %8.3f    %8.3f \t %s\n',infCombVar(i),num2str(combVar(i,:)),...
%                 costLkCombVar(i), costMkCombVar(i), costARCombVar(i), totCostComVar(i), arCombVar(i),num2str(powCombVar(i,:)));
%         end
        
        
        %--------------------------------------------------------------------------
        %Selection of the best Combination of variations of services
        %--------------------------------------------------------------------------
        
        %--- The criteria for The best combnation is the minimum total cost
        [sortComVar,indSortCombVar] = sort(totCostComVar);
        
        
        %--- Information of the best Combination of variations of services
        idBestCombVar = indSortCombVar(1); %id of the best combinations of variations of services
        infBestCombVar = comVarServ(indSortCombVar(1),:); % information of variations that compose the best combination
        
        costLkBestCombVar = costLkCombVar(indSortCombVar(1)); %information about the cost of best combination
        costMkBestCombVar = costMkCombVar(indSortCombVar(1));
        costARBestCombVar = costARCombVar(indSortCombVar(1));
        totCostBestComVar = totCostComVar(indSortCombVar(1));
        
        arBestCombVar = arCombVar(indSortCombVar(1)); %AR of the best combination
        powBestCombVar = powCombVar(indSortCombVar(1),:); %Real power consumed of the best combination
        resPowBestCombVar = PwRefLevel - powBestCombVar; %residual power of the best combination
        
        
%         disp(' ');
%         disp(' ---------------------------------   Information of the best Combination of variations of services --------------------------------');
%         disp('- NumComb -  Comb. Var.  -  costLk  - costMk  - CostAR  - totCost -   AR -    Power demanded   - Residual power');
        %fprintf('%10d    %s %10d %10d   %8.3f  %8.3f  %8.3f     %s \t %s \t   %s\n',idBestCombVar,num2str(infBestCombVar),...
        %    costLkBestCombVar, costMkBestCombVar, costARBestCombVar, totCostBestComVar, arBestCombVar,...
        %    num2str(powBestCombVar),num2str(resPowBestCombVar));
%         fprintf('%10d    %s %10d %10d   %8.3f  %8.3f  %8.3f     %s \n',idBestCombVar,num2str(infBestCombVar),...
%             costLkBestCombVar, costMkBestCombVar, costARBestCombVar, totCostBestComVar, arBestCombVar);
        
        
        %--------------------------------------------------------------------------
        %Information about the optimal allocation of services with the aim
        %of maximizing the utilization of available spectrum
        %--------------------------------------------------------------------------
        
        %Information of variations of services processed
        infCombVarSerProc = infBestCombVar;
        infCombVarSerProc(infCombVarSerProc == 0) = []; %delete the columns with zeros
        
        %-Information of services processed
        infServProc = varInfo(infCombVarSerProc);
        
        %-Information of power of processed services
        tempPowServProc = powVarServ(infCombVarSerProc,:);
        
        %Service ID per frequency slot
        for i = 1:length(infServProc)
            infPowServProc(i,:) = tempPowServProc(i,:)*infServProc(i);
        end
        
        %--- Frequency shfting of the processed services and cummulative
        %value of the total frequency shifting performed
        
        %-Number of frequency shifting performed
        freqShifVarSerProc = mFsVarInfo(infCombVarSerProc);
        numFreqShitServ = (length(find(freqShifVarSerProc > 0)))/n; %n to present a normalized value
        
        %-Cummulative value of total frequency shifting
        totalFreqShifVal = (sum(freqShifVarSerProc,2))/n; %n to present a normalized value
        
        
        %--- Utilization of available spectrum
        
        %-Use of spectrum in GHz
        useBWServ = sum(services(infServProc,5));
        
        %-Percentage of utilization of available spectrum considering only
        %the processed services
        useBWServFreSlots = sum(servicesFs(infServProc,5));
        UBW = (useBWServFreSlots /BWFs) * 100;
        
        
        %--------------------------------------------------------------------------
        % Graphs of allocation of services within the available spectrum
        %--------------------------------------------------------------------------
        
        %- Figures: Summary available bandwidth (BW), initial service distribution
        %(wavelengths) and optimal spectrum allocation
        figure('name','Available bandwidth (BW) and initial distribution of services');
        
        subplot(3,1,1) %Available bandwidth (BW)
        imagesc(availableBWFs);
        %imagesc([minFreqRef maxFreqRef],1,availableBWFs)
        title(['Available bandwidth (BW): ',num2str(BWFs), ' frequency slots']);
        xlabel('Frequency slots');
        ylabel('Available bandwidth (BW)');
        grid on;
        
        subplot(3,1,2) %Service distribution (wavelengths)
        %imagesc([minFreqRef maxFreqRef],[1 8],availableBWFs)
        imagesc(initDistributionServ);
        title('Initial service distribution (wavelengths)');
        xlabel('Frequency slots');
        ylabel('Services');
        grid on;
        
        subplot(3,1,3) %Service distribution (wavelengths)
        %imagesc([minFreqRef maxFreqRef],[1 8],availableBWFs)
        imagesc(infPowServProc);
        title(['Optimal spectrum allocation (wavelengths)',' AR: ',num2str(arBestCombVar),' %']);
        xlabel('Frequency slots');
        ylabel('Services');
        grid on;
        
        
        %--- Information of each case
        idCase = strcat('Case',num2str(iterCase),'_Fs',num2str(iterFs));
        parsave(idCase,iterCase,iterFs,arBestCombVar,useBWServ,UBW,numFreqShitServ,totalFreqShifVal);
        
        %--- Store of information parameters
        caseInfo = [caseInfo;iterCase]; %information of case under analysis
        frecShiftingInfo = [frecShiftingInfo;iterFs]; %value of frequency shifting of services
        acceptanceRatioInfo = [acceptanceRatioInfo;arBestCombVar]; %value of acceptance ratio of processed services
        useAvaBWInfo = [useAvaBWInfo;useBWServ]; %used bandwidth in GHz
        useAvaBWpercentInfo = [useAvaBWpercentInfo;UBW]; %utilization percentage of available bandwidth
        numFreqShiftinInfo = [numFreqShiftinInfo;numFreqShitServ]; %number of times frequency shifting has been used
        cumFreqShiftinInfo = [cumFreqShiftinInfo;totalFreqShifVal]; %cummulative value of the frequency shifting
        
    end %iterations for Ts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end %iterations for case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infElapsedTime = toc;

% ---- Saving variables for getting the pictorial results
solOptFlexGridExample = [caseInfo,frecShiftingInfo,acceptanceRatioInfo,useAvaBWInfo,...
    useAvaBWpercentInfo,numFreqShiftinInfo,cumFreqShiftinInfo];
save('solOptFlexGridExample.mat','solOptFlexGridExample',...
    'caseInfo','frecShiftingInfo','acceptanceRatioInfo','useAvaBWInfo',...
    'useAvaBWpercentInfo','numFreqShiftinInfo','cumFreqShiftinInfo','infElapsedTime');




%--------------------------------------------------------------------------
% Functions
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%------ Function for computing the power demanded and all information of
%variation per services considering the time shifting interval
%--- Input parameters:
%inN: Number of services
%inMaxBW: Maximum frequency horizon
%inServices: Information of services
%inRefInitFreqServ: Inital frequencies of services referenced at frequency
%slot 0
%--- Out parameters:
%- outNumVariat: Total number of valid variations
%- outInfIdVariat: Identification of all possible variation per services (row
%vector)
%- outServInfo: Information of the service accoring to its variation. Each
%service can be originated one or more than one variations.
%- outPriorityInfo: Information about the priority level of a variation of
%service.
%outMtsInfo: Information about the maximum time shifting of a variation of
%service.
%outPowerVar: Power demanded per variation and per time slot
%Example: [numVar,idVar,varInfo,priorVarInfo,mtsVarInfo,powVarServ] = funcInfoDemands(n,tmax,services);

function [outNumVariat,outInfIdVariat,outServInfo,...
    outPriorityInfo,outMtsInfo,outPowerVar] = funcVarServ(inN,inMaxBW,inServices,inRefInitFreqServ)


%Power demanded for services
powerDemanded = zeros(inN,inMaxBW);
backwardDemand = zeros(1,inMaxBW);
forwardDemand = zeros(1,inMaxBW);
cont = 0; % variable to know the number of possible combinations

demandPower = []; %variable array for the information of all demands (services)
serviceInfo = []; %information of the processed services
outMtsInfo = []; %number of shifting per service in terms of time slots
outPriorityInfo = []; %information of the priority of each service
auxFShif = 0;

%auxiliar variables for time shifting
auxBackshifted = zeros(inN,inMaxBW);
auxForshifted = zeros(inN,inMaxBW);

for i = 1:inN % i number of services
    
    finit = inRefInitFreqServ(i);
    bk = inServices(i,7) - 1;
    powerDemanded(i,finit : finit + bk) = 1 * inServices(i,3);
    cont = cont + 1;
    serviceInfo = [serviceInfo,i];
    demandPower =[demandPower; powerDemanded(i,:)];
    auxFShif = 0; %indicates no time shifting(numer slots)
    outMtsInfo = [outMtsInfo,auxFShif];
    outPriorityInfo = [outPriorityInfo,inServices(i,2)];
    %disp(['Demand ',num2str(i), ' without time shifting']);
    %disp(powerDemanded(i,:));
    
    %auxBackshifted = zeros(n,m);
    auxBackshifted = powerDemanded(i,:);
    %auxForshifted = zeros(n,m);
    auxForshifted = powerDemanded(i,:);
    
    if inServices(i,9) > 0 %if there exist time shifting
        
        for j = 1 : inServices(i,9)
            
            %Forward frequency shifting
            if ((finit + bk + j) > inMaxBW)
                powerDemanded(i,:)= Inf;
                %disp('No valid demand, demand out of limit');
                %disp(' ');
                demandPower =[demandPower; powerDemanded(i,:)];
            else
                backwardDemand(:,j+1:end) = auxBackshifted(:,1:end-j);
                %disp(['Demand ',num2str(i), ' backward time shifting ', num2str(j), ' time slots']);
                %disp(backwardDemand);
                demandPower =[demandPower; backwardDemand];
            end
            backwardDemand = zeros(1,inMaxBW); %n for erase de previous result
            
            cont = cont + 1;
            serviceInfo = [serviceInfo,i];
            outMtsInfo = [outMtsInfo,j];
            outPriorityInfo = [outPriorityInfo,inServices(i,2)];
            
        end
    end
    
    if inServices(i,8) > 0 %if there exist time shifting
        
        for j = 1: inServices(i,8)
            %Backward frequency shifting
            if ((finit - j) <= 0)
                powerDemanded(i,:)= Inf;
                demandPower =[demandPower; powerDemanded(i,:)];
                
            else
                forwardDemand(:,1:end-j) = auxForshifted(:,j+1:end);
                %disp(['Demand ',num2str(i), ' forward time shifting ', num2str(j), ' time slots']);
                %disp(forwardDemand);
                demandPower =[demandPower; forwardDemand];
            end
            forwardDemand = zeros(1,inMaxBW); %n for erase de previous result
            
            cont = cont + 1;
            serviceInfo = [serviceInfo,i];
            outMtsInfo = [outMtsInfo,j];
            outPriorityInfo = [outPriorityInfo,inServices(i,2)];
            
        end
        
    end
    
end

%disp('--------------------------------------------------------------------------------');
%disp('-- Possible variations per service and possible combinations among services --');
%disp(' ');
%disp(['Number of possible combinations per service: ', num2str(cont)]);
%disp(' ');
%disp('Order of service combinations: without Ts, Ts backward, Ts forward');

infoDemands = demandPower; %only to show the results of demand demandPower
infoDemands(infoDemands == Inf) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('TABLE II: Possible Combinations per service');
%disp(' ------------- Possible Combinations per service ------------------');
%disp(' "---" means not valid demand, invalid time shifting');
%disp('In this table in column 3 the time slot used is given in terms of power consumed by demand');
%disp('- Number - Sk - TimeSlot used - Absolute Time shifting Performed');
%numberServices = 1:cont;
%[rInfoD,cInfoD] = size(infoDemands);
%disp([numberValidServices', serviceInfo', demandPower]);
%disp('- Number - Sk - Type of demand (in terms of power consumed by demand)');
%disp([numberServices', serviceInfo', infoDemands]);

% for i = 1:cont
%     if infoDemands(i,:) == 0
%         disp(['     ',num2str(numberServices(i)), '     ',...
%             num2str(serviceInfo(i)) 9, '  ', '-------' blanks(10),...
%             '-']); %num2str(shiftingServ(i))
%     else
%         disp(['     ',num2str(numberServices(i)), '     ',...
%             num2str(serviceInfo(i)) 9, '  ',num2str(infoDemands(i,:)) blanks(10), ...
%             num2str(shiftingServ(i))]);
%
%     end
% end
% disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Valid demands avoiding the invalid type of demands (time shifting)

%--- Valid type of demands
infoDemandsTrans = infoDemands';
outPowerVar = infoDemands(find(sum(infoDemandsTrans)>0),:);
[outNumVariat, valdc] = size(outPowerVar);
orderValidDemands = find(sum(infoDemandsTrans)>0);
noValidDemands = find(sum(infoDemandsTrans)==0);
[noValr, noValc] = size(noValidDemands);

outServInfo = serviceInfo(find(sum(infoDemandsTrans)>0));

outPriorityInfo = outPriorityInfo(find(sum(infoDemandsTrans)>0));

%--- Power demanded
demandPower = demandPower(find(sum(infoDemandsTrans)>0),:);

%--- Time shifting information
outMtsInfo = outMtsInfo(find(sum(infoDemandsTrans)>0));

%disp(' ');
%disp('--------------------------------------------------------------------------');
%disp('-- Valid combinations per service and valid combinations among services --');
%disp(' ');
%disp(['Number of valid combinations per service: ', num2str(outNumVariat)]);
%disp(' ');
%disp('TABLE III: Valid Combination per Service');
%disp('---------------- Valid combination per service ---------------------');
%disp('- Number - Sk - lk - mts - Power demanded per time slot)');
outInfIdVariat = 1:outNumVariat;
%disp([numberValidServices', serviceInfo', demandPower]);
%disp([orderValidDemands', valserviceInfo', valinfoDemads]);
%disp([outInfIdVariat', outServInfo', outPriorityInfo', outMtsInfo', outPowerVar]);


%Type of demands according to the time slot where the service is processed
%This information is needed later on, in the computation of power
%inforDemPerSlot = infoDemands;
inforDemPerSlot = outPowerVar;
for i = 1:outNumVariat
    for j = 1:valdc
        if (inforDemPerSlot(i,j) > 0)
            inforDemPerSlot(i,j) = j;
        end
    end
end
%inforDemPerSlot
%disp(' ');
%disp('Types of demands according to the time slot where the service is processed');
%disp('- Number - Sk - Type of demand -----');
%disp([numberServices', serviceInfo', inforDemPerSlot]);


%Type of demands according to the service to which they belong
%inforDemPerSer = infoDemands;
inforDemPerSer = outPowerVar;
inforDemPerSer(inforDemPerSer>0) = 1; %to produce a matrix with unity elements
[rdps cdps] = size(inforDemPerSer);

for i = 1:rdps
    %inforDemPerSer(i,:) = inforDemPerSer(i,:)* serviceInfo(:,i);
    inforDemPerSer(i,:) = inforDemPerSer(i,:)* outServInfo(:,i);
end
%disp(' ');
%disp('Types of demands according to the service to which they belong');
%disp('- Number - Sk - Type of demand -----');
%disp([numberServices', serviceInfo', inforDemPerSer]);

end


%--------------------------------------------------------------------------
%------ Function for computing the combinations among variations of
%services
%--- Input parameters:
%inVarInfo: Information of all variations per service, ID of each variation
%--- Out parameters:
%outNumRowCombVar: Number of combinations of variations of services
%outNumColCombVar: Number of processed services at the same time (n)
%outInfCombVar: Identification of combination
%outCombVar: Information of combinations of variations per services
%Example: [numRowCombVar,numColCombVar,infCombVar,combVar] =
%funcCombVarSer(varInfo;

function [outNumRowCombVar,outNumColCombVar,outInfCombVar,outCombVar] = funcCombVarSer(inVarInfo)

%Values of services
valServ = unique(inVarInfo);
[rvalSer,cvalSer] = size(valServ);

%Temporal information of service indexes
tmpInfIndSer = {};

for i = 1:cvalSer
    auxPosServ = find(inVarInfo == i);
    tmpInfIndSer{end+1} = auxPosServ; %to append the indexes of all services
end

numSer = cellfun(@numel,tmpInfIndSer); %number of demands per service
cumSumSer = cumsum(numSer); %cummulative sum of number of services
combSerPos = fullfact(numSer); % different combinations of elements among services
%Expansion function to obtain combinations
outCombVar = bsxfun(@plus,combSerPos,[0, cumSumSer(1:end-1)]);
[outNumRowCombVar,outNumColCombVar] = size(outCombVar);
outInfCombVar = 1:outNumRowCombVar; %information id combination according variations

end


%--------------------------------------------------------------------------
%------ Function for computing the acceptance ratio and the cost related
%per combination of variations of services considering priority
%--- Input parameters:
%inNumRowCombVar: Number of variations of services
%inResPowCombVar: Residual information of combination of variations
%inCombVar: Information of Combination of variations
%inN: Number of variations per combination
%inTmax: Maximum time horizon
%inPes: Information of Pes
%inTpes: Time window related to Pes
%inPriorVarInfo: Priority of variations of service
%inPowVarServ: Power of variations of service
%inConsumedPowerComb: Consumed power per combination
%inBkInfo: Bandwidth of services
%inOptionAR: Option for the computation of the AR

%--- Out parameters:
%arCombVar: Acceptance ratio combination of variations of services
%outCostCombVar: Cost related to rejected service (acceptance ratio)
%outPowCombVar: Power consumed by only the processed services
%outComVarServ: Combination of variations processed

%Example: [arCombVar,costCombVar,powCombVar] = funcARCostCombVarSer(numRowCombVar,resPowCombVar,combVar,n,tmax,Pes,tpes,priorVarInfo,powVarServ,consumedPowerComb);

function [outARCombVar,outCostCombVar,outPowCombVar,outComVarServ] = funcARCostCombVarSer(inNumRowCombVar,...
    inResPowCombVar,inCombVar,inN,inTmax,inPes,inTpes,inPriorVarInfo,...
    inPowVarServ,inConsumedPowerComb,inBkInfo,inOptionAR)

outARCombVar = zeros(inNumRowCombVar,1); %Acceptance ration of combinations of variations
outCostCombVar = zeros(inNumRowCombVar,1); %Cost of combinations of variations
outPowCombVar = zeros(inNumRowCombVar,inTmax); %Power of combinations of variations
outComVarServ = zeros(inNumRowCombVar,inN); %Information of services processed

[increChBW, indIncreChBW] = sort(inBkInfo); %services bandwidth ordered in increasing order
[decreChBW, indDecreChBW] = sort(inBkInfo,'descend'); %services bandwidth ordered in decreasing order

for i = 1:inNumRowCombVar %exploration of all combinations
    
    auxInPes = inPes; %auxiliar varible to compute the remainder available spectrum
    
    %----------------------------------------------------------------------
    %--- Combinations of variations with negative residual power
    %The cost of each service that cannot be processed will be 1000(big M)
    if (find(inResPowCombVar(i,:)<0) > 0)
        %disp('Negative combination of variations of services');
        %outARCombVar(i,1) = 0;
        %outCostCombVar(i,1) = inN*1000;
        
        %--- Variations that compose the negative combination of variations
        auxNegServCombVar = inCombVar(i,:);
        
        %--- Priority info of the services that compose the negative
        %combination of variations
        auxPriorNegSerCombVar = inPriorVarInfo(1,auxNegServCombVar);
        
        
        switch inOptionAR
            
            case 1 %--- AR option 1
                %disp('AR option 1');
                
                %------------------------------------------------------------------
                %--- Processing of negative combination of variatiosn with equal
                %priorities in services
                
                %Combination of possible rejectetable services
                allRejSerCombVar = [];
                
                %Maximun number of colums among combinations
                auxLarRejServComVar = zeros(1,inN);
                
                %--- Analysis of all possible rejectable services
                for j=1:inN
                    combRejSer = combnk(1:inN,j); %all combinations of possible rejected services
                    [rCombRejSer,cCombRejSer] = size(combRejSer);
                    
                    if cCombRejSer < inN
                        combRejSer(:,numel(auxLarRejServComVar)) = 0;
                    end
                    
                    %All possible combinations of rejectable services
                    allRejSerCombVar = [allRejSerCombVar; combRejSer];
                end
                %disp(allRejSerCombVar);
                [rAllRejSerCombVar,cAllRejSerCombVar] = size(allRejSerCombVar);
                
                %--- Analysis of power demanded by the non-rejectable services
                %and the cost related to the rejection of the services
                
                %After the computation of the residual power with the information of
                %non-rejectable services if the resulting Pres is negative cost of 10000
                %Cost of Pres non negative cero
                %Therefore the cost of ar is equal to CostAR= costNegPres + cost for each service rejected
                %Example: CostAR = 10000 + 2*1000 (2 services rejected)
                %If after the rejection, and fater the computation of the power consumed by
                %the non-rejectable services the Pres is positive the cost related is the
                %the Pres normalized
                
                tempCostRejServ = zeros(rAllRejSerCombVar,1); %temporal information about the cost if a service is rejected
                tempCostResPow = zeros(rAllRejSerCombVar,1); %temporal information about the cost if the residual power after
                %the assigned cost will be a
                %big value (1000), if the
                %residual power is positive
                %this value can be se a metric
                %(shigma Pres)
                tempAcceptanceRatio = zeros(rAllRejSerCombVar,1); %temporal information of acceptance ratio
                tempPowerDemandNonRejSev = zeros(rAllRejSerCombVar,inTmax); %temporal information of power consumed
                %by non-rejectable
                %services
                tempCombVarSer = zeros(rAllRejSerCombVar,inN); %temporal information of variations processed
                
                for j = 1:rAllRejSerCombVar
                    
                    %Information about the number of rejected services
                    numRejServ = length(find(allRejSerCombVar(j,:) > 0));
                    
                    %--- Cost of rejection of all possible rejectable services
                    tempCostRejServ(j,1) = numRejServ * 1000;
                    
                    %--- Computation of power of nonrejectable services
                    %Variable to save all possible services of combination
                    varAuxRejecServ = auxNegServCombVar;
                    
                    %Variable to be deleted according the combination of rejected services
                    tempAuxRejecServ = varAuxRejecServ;
                    
                    %Variable to store the temporal information of rejectable services
                    tempRejServ = allRejSerCombVar(j,:);
                    tempRejServ(tempRejServ == 0) = []; %delete the columns with zeros
                    
                    %Information of non-rejectable services per combination of variations
                    tempAuxRejecServ(:,tempRejServ) = [];
                    
                    
                    %Information of combination of variations processed
                    auxCombVarser = zeros(1,inN);
                    auxCombVarser(1:length(tempAuxRejecServ)) = tempAuxRejecServ;
                    
                    
                    %Power consumed by non-rejectable services
                    powNonRejServ = sum(inPowVarServ(tempAuxRejecServ,:),1);
                    
                    %Residual power of non rejectable services
                    resPowNonRejSer = inPes - powNonRejServ;
                    
                    if (find(resPowNonRejSer<0) > 0)
                        %disp('The combination of non-rejectable services produce a negative residual power');
                        %--- Cost related to the negative residual power of all
                        %combinations
                        tempCostResPow(j,1) = 10000; %cost asigned if the combination of non-rejectable services
                        %still produce a negative residual
                        %power
                        tempAcceptanceRatio(j,1) = 0; %if Pres still negative no processing of demands
                        tempPowerDemandNonRejSev(j,:) = zeros(1,inTmax); %no real consumption
                        
                        tempCombVarSer(j,:) = zeros(1,inN); %no combinations of services processed
                    else
                        %disp('The combination of non-rejectable services produce a positive residual power');
                        %--- Cost related to the negative residual power of all
                        %combinations
                        auxSigResRejServ = ((sum(resPowNonRejSer.^2))/inTpes)^(1/2); %metric for resiudal power of rejected services
                        tempCostResPow(j,1) = auxSigResRejServ; %cost related to the remaining Pres
                        tempAcceptanceRatio(j,1) = ((inN - numRejServ)/inN)*100; %computation of AR
                        tempPowerDemandNonRejSev(j,:) = powNonRejServ; %power consumed according the services processed
                        
                        tempCombVarSer(j,:) = auxCombVarser; %information of variations of services processed
                    end
                end
                
                %--- Total cost related to rejection per Combination of Variations of
                %services
                tempTotalCostRejServ = tempCostRejServ + tempCostResPow; %temporal total cost related to rejection
                %tempTotalCostRejServ =
                %tempCostRejServ + tempCostResPow
                
                %Sorted information of all analysis of non-rejectable services
                [tempCostRej,indTempCosRej] = sort(tempTotalCostRejServ); %The best value is the first in the list
                
                %Information of AR for services of equal priority
                infARNonRejSer = tempAcceptanceRatio(indTempCosRej(1,1),1);
                outARCombVar(i,1) = infARNonRejSer;
                
                %Information of cost for services of equal priority
                infCosNonRejSer = tempTotalCostRejServ(indTempCosRej(1,1),1);
                outCostCombVar(i,1) = infCosNonRejSer;
                
                %Information of power consumed for services of equal priority
                infConsumPowNonRejSer = tempPowerDemandNonRejSev(indTempCosRej(1,1),:);
                outPowCombVar(i,:) = infConsumPowNonRejSer;
                
                %Information of variation of services processed within the
                %combination
                infCombVarSerProcess = tempCombVarSer(indTempCosRej(1,1),:);
                outComVarServ(i,:) = infCombVarSerProcess;
                
            case 2 %--- AR option 2
                %disp('AR option 2');
                
                tempVarSer = []; %temporal information of processed services
                auxCombVarser = zeros(1,inN); %variation of services porcessed
                tempPowerDemandNonRejSev = zeros(1,inTmax); %variable to store the power of processed services
                numRejServ = 0; %storage of rejected services
                
                for j = 1:inN
                    auxVarServ = auxNegServCombVar(indIncreChBW(j));
                    auxPowVarServ = inPowVarServ(auxVarServ,:);
                    auxResPow = auxInPes - auxPowVarServ;
                    
                    if (find(auxResPow<0) > 0)
                        numRejServ = numRejServ + 1;
                        continue;
                    else
                        %Variations of the service processed
                        tempVarSer = [tempVarSer,auxVarServ];
                        %Power of variations of services
                        tempPowerDemandNonRejSev = tempPowerDemandNonRejSev + auxPowVarServ;
                        %Update of power for the nex variation of service
                        auxInPes = auxResPow;
                        auxInPes(auxInPes<0)=0; %to ensure only possitive values
                    end
                    
                end
                
                tempAcceptanceRatio = ((inN - numRejServ)/inN)*100; %computation of AR
                tempCostRejServ = numRejServ * 1000;
                auxCombVarser(1:length(tempVarSer)) = tempVarSer;
                
                %Information of AR for services of equal priority
                outARCombVar(i,1) = tempAcceptanceRatio;
                
                %Information of cost for services of equal priority
                outCostCombVar(i,1) = tempCostRejServ;
                
                %Information of power consumed for services of equal priority
                outPowCombVar(i,:) = tempPowerDemandNonRejSev;
                
                %Information of variation of services processed within the
                %combination
                outComVarServ(i,:) = auxCombVarser;
                
                
            case 3 %--- AR option 3
                %disp('AR option 3');
                
                tempVarSer = []; %temporal information of processed services
                auxCombVarser = zeros(1,inN); %variation of services porcessed
                tempPowerDemandNonRejSev = zeros(1,inTmax); %variable to store the power of processed services
                numRejServ = 0; %storage of rejected services
                
                for j = 1:inN
                    auxVarServ = auxNegServCombVar(indDecreChBW(j));
                    auxPowVarServ = inPowVarServ(auxVarServ,:);
                    auxResPow = auxInPes - auxPowVarServ;
                    
                    if (find(auxResPow<0) > 0)
                        numRejServ = numRejServ + 1;
                        continue;
                    else
                        %Variations of the service processed
                        tempVarSer = [tempVarSer,auxVarServ];
                        %Power of variations of services
                        tempPowerDemandNonRejSev = tempPowerDemandNonRejSev + auxPowVarServ;
                        %Update of power for the nex variation of service
                        auxInPes = auxResPow;
                        auxInPes(auxInPes<0)=0; %to ensure only possitive values
                    end
                    
                end
                
                tempAcceptanceRatio = ((inN - numRejServ)/inN)*100; %computation of AR
                tempCostRejServ = numRejServ * 1000;
                auxCombVarser(1:length(tempVarSer)) = tempVarSer;
                
                %Information of AR for services of equal priority
                outARCombVar(i,1) = tempAcceptanceRatio;
                
                %Information of cost for services of equal priority
                outCostCombVar(i,1) = tempCostRejServ;
                
                %Information of power consumed for services of equal priority
                outPowCombVar(i,:) = tempPowerDemandNonRejSev;
                
                %Information of variation of services processed within the
                %combination
                outComVarServ(i,:) = auxCombVarser;
                
                
        end
        
        
        %----------------------------------------------------------------------
        %--- Combinations of variations with positive residual power
    else
        %disp('Possitive combination of variations of services');
        outARCombVar(i,1) = 100; %Information of Acceptance Ratio
        outCostCombVar(i,1) = 0; %Information of cost related to acceptance ratio
        outPowCombVar(i,:) = inConsumedPowerComb(i,:); %Final power consumed per combination
        outComVarServ(i,:) = inCombVar(i,:); %Variation of services processed
    end
    
end

end


% --- Function to store data, compatible with parfor instruction
%par_Case = [idCase,iterCase,iterFs,arBestCombVar,useBWServ,...
%UBW,numFreqShitServ,totalFreqShifVal];
function parsave(fname,iterCase,iterFs,arBestCombVar,useBWServ,UBW,...
    numFreqShitServ,totalFreqShifVal)
par_Case = [iterCase,iterFs,arBestCombVar,useBWServ,UBW,numFreqShitServ,totalFreqShifVal];
save(fname, 'par_Case')
end
