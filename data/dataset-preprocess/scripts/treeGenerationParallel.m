function treeGenerationParallel(rnd_seed, Nsim_input)
%-------------------------------------------------------------------------
% DESCRIPTION
% simulate the deletions using parameters from the experimental data
% and the Drosophila cell division times
% the experimental data is based on the "FAST" target in the embryo
% for explanation see the notebook "Fast_sequence_Embryo_Reader.ipynb"
%-------------------------------------------------------------------------
% set the random seed to 123 so it can be reproducible
rng(rnd_seed);

% Define the MAIN parameters
%-------------------------------------------------------------------------
% cells divisions
% Ndiv = 16;
% Ndiv = 10;
Ndiv = 7;

%Matrix of cell division times
CellDiv_Time = zeros(Ndiv,1);
CellDiv_Time(1:13) = 0:60:720;  % First 13 cell divisions are every 10 mins
% CellDiv_Time(1:13) = 0:10:120;  % First 13 cell divisions are every 10 mins
CellDiv_Time(13:16) = 120:160:600; %Last 3 divisions are every 160 mins
CellDiv_Time(17) = 1440;  % end of Embryonic development

% Number of targets
targets = 200;
% Mutation rate per minute
Lambda = 0.0014;
%NUmber of sims
% Nsims = 1000;
Nsims = Nsim_input;

%% to account for the "missing" mutations (~5%)
states = 30;
rvs = gamrnd(0.1,2, [1 states]);
rvs = rvs/sum(rvs);
Nmer_prob = sort(rvs);

    for r = 1:Nsims                      % TO HAVE MULTIPLE SAMPLES

        %the array is a 2D array with the targets (within the cell) as columns and
        % each row is a daughter cell (matrix changes in each iteration)
        % ----------Initialize some variables ------------------------------------
        sequences = struct('Header',{},'Sequence',{});
        cell_name = ones(1);
        A1 = zeros(1,targets);

        %-----------------------------------------------------------------
        % -------------- MAIN LOOP ---------------------------------------
        %------------------------------------------------------------------
        for n_div = 1:Ndiv                           %for each division
            %%% ---Determine the time interval in each cell division
            t1 = CellDiv_Time(n_div+1);        % end   time in min to check for events wit POISSON
            t0 = CellDiv_Time(n_div);      % start time in min to check for events wit POISSON

            %%% -------
            Ndaug_1 = 2^(n_div-1);
            Ndaugh  = 2^n_div;                            % calculate daughters number
            ncells_1 = size(A1,2);

            % create the matrix of daughters by copying the mother
            A2 = vertcat(A1,A1);
            A2_old = A2;
            A3 = char(zeros(size(A2)));
            % keep track of the cell names
            cell_name = horzcat(cell_name,cell_name);
            for cn = 1:Ndaug_1
                cell_name(cn)   = (cell_name(cn) *2)-1;
                cell_name(cn+Ndaug_1)   = cell_name(cn+Ndaug_1) *2;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%
            %events = zeros(size(Aind));
            events = zeros(Ndaugh*targets,1);

            for i = 1:size(events) 
                a = poisson_fixed_events(0.0014,1);
                events(i,1) = a(2);
            end

            Nevents = size(events(events > t0 & events <= t1),1);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % for each copied target, mutate randomly with a fixed probability
            % the site can mutate only once, into one of 30 possible states
            Aind = find(A2 == 0);
            AR = 0;
            if size(Aind,1) > 0 %%%% only if there are unmutated targets!
                if Nevents >= size(Aind,1)
                Eind = Aind;
                end

                if Nevents < size(Aind,1)
                    Eind = randsample(Aind,Nevents);
                end

            AR = rand(size(Eind));
            AR = arrayfun(@(z)sum(z <= Nmer_prob), AR);

            A2(Eind) = AR;     %%% MODIFY HEREEE!! Aind should be the ref
            end

            %%%%%%%%% MAKE THE INTER-TARGETS HAPPEN  %%%%%%%%%%%%%%%%%%%%%%%%%%

            for j = 1:Ndaugh
                C_ne = find(A2(j,:)~=A2_old(j,:));
                 if length(C_ne) > 1 %% if more than 1 events happened
                    for k = 1:length(C_ne) %% Parse all the mutation events
                        C_ne_dist = abs(C_ne - C_ne(k)); 
                        idx = find((0< C_ne_dist) & (C_ne_dist < 20)); 
                        if ~isempty(idx) %% the distance of mutation is less than 20
                            for l = 1:length(idx)
                                if idx(l) > k
                                    % assign to the targets inbetween the state 200 
                                    A2(j,C_ne(k)+1:C_ne(idx(l))-1) = 200;
                                end
                            end
                        end
                    end
                 end
            end       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make the daughter matrix the mother matrix for the next iteration
            A1 = A2;

            % ------------- IO --------------------------------------------
            % name of the simulation
            sim = [sprintf('%02d',n_div),'div_',...
                sprintf('%02d',targets),'_targets_'...
                sprintf('%02d',states),'_states_'...
                ];

            %define the output file name for the alignment (fasta)
            out = [sim, sprintf('%04d', r),'rep_FAST_', sprintf('%04d', rnd_seed), 'seed','.fas'];
            freq_file = [sim, sprintf('%04d', r),'rep_FAST_', sprintf('%04d', rnd_seed), 'seed','_freqs','.txt'];

            %---------- produce the output for each division!!
            %%%% Converts to ASCII characters, including 0-9, A-Z, a-z
            %%%% Avoids troublesome characters as ?, \, etc..

            %%% most common, last is most common (A-Z)
            A3(A2 > 0 & A2 < 27) = char(A2(A2 > 0 & A2 < 27) + 64);

            % 11th below most common (a-d)
            A3(A2 >= 27 & A2 <= 30) = char(A2(A2 >= 27 & A2 <= 30 ) + 69);


            %DELETED (MISSING) targets
            A3(A2==200)= 45;

            %unmutated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            A3(A2==0)= char(A2(A2==0) + 48);

            %%%%%% EXPORT THE SIMULATED FREQS TO A FILE %%%%%%%%%%%%

            unqA = unique(A3(1:numel(A3)));
            Nmers = cellstr(unqA');
            countNmers = histc(A3(1:numel(A3)),unqA); 
            relFreq=countNmers/numel(A3);

            relFreqCell = num2cell(relFreq)';
            %transp(relFreq);
            %num2cell(relFreq,dim)
            %T = table(Nmers, transp(relFreq));
            T = table(Nmers, relFreqCell);
            if n_div == Ndiv
                % create the file
                writetable(T,['../simulations/',freq_file], ...
                'WriteVariableNames',false);            
            end
            clear A2 AR;
            %%%%%%%%%%%%%
            % define the "true" tree with the number of final daughter cells
            %Totdaugh = 2^n_div;


            if n_div == Ndiv
                Nbranch = Ndaugh -1;
                canonic = 1:(2*Nbranch);
                canonic = reshape(canonic,[2,Nbranch]);
                canonic = canonic';
                % create the name of the cells
                myfasta = fopen(['../simulations/',out],'w');
                fclose(myfasta);

                myfasta = fopen(['../simulations/',out],'a');

                for ii = 1:Ndaugh
                    sequences(ii).Header = ['c_',sprintf('%07d',cell_name(ii))];
                    sequences(ii).Sequence = strjoin(cellstr(A3(ii,:)));
                end
                %%%% SORT THE CELLS IN THE OUTPUT SO IT IS EASY TO READ
                [~,index] = sortrows({sequences.Header}.');
                sequences = sequences(index); clear index

                for ii = 1:Ndaugh
                    % append to a file
                    fprintf(myfasta,'%s\t',sequences(ii).Header);
                    fprintf(myfasta,'%s\n',sequences(ii).Sequence);
                end
                fclose(myfasta);
                if r == 1
                    %define the output file name for the reference tree (newick)
                    RefOut = [sim, 'REF.nw'];
                    % construct a nw file to compare it with the simulated ones
                    RefTree = phytree(canonic);
                    d = get(RefTree,'distances');
                    p = get(RefTree,'pointers');
                    n = sort({sequences.Header});
                    RefTree2 = phytree(p,d,n);
                    %view(RefTree2)
                    % write newick files to compare with ete software externally
                    phytreewrite(['../Ref_Trees/',RefOut], RefTree2, ...
                        'Distances', 'false','BranchNames','false');
                end
            end
        end
    end

end


function [ t, w ] = poisson_fixed_events ( lambda, event_num )

%*****************************************************************************
%
%% POISSON_FIXED_EVENTS simulates the occurrence of a given number of Poisson events.
%  Parameters:
%
%    Input, real LAMBDA, the average number of events per unit time.
%
%    Input, integer EVENT_NUM, the number of events to wait for.
%
%    Output, real T(EVENT_NUM+1), the time at which a total of 0, 1, 2, ...
%    and EVENT_NUM events were observed.
%
%    Output, real W(EVENT_NUM+1), the waiting time until the I-th event
%    occurred.
%
%
%  Poisson waiting times follow an exponential distribution.
%
  w(1) = 0.0;
  w(2:event_num+1) = - log ( rand ( event_num, 1 ) ) / lambda;
%
%  The time til event I is the sum of the waiting times 0 through I.
%
  t(1:event_num+1) = cumsum ( w(1:event_num+1) );

  return
end