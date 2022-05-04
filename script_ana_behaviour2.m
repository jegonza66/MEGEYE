%% Quick script to run some preliminary analysis on behaviour
%% 4/7/2018
%% 
clear all
close all
clc
platform=computer;
switch platform
    case{'PCWIN64'}
        addpath(genpath('H:/MEGEYE/analysis'))
        DIR.ET_data = 'H:/MEGEYE/ET_DATA/';
        DIR.out     = 'H:/MEGEYE/ET_DATA/out';
        DIR.codes   = 'D:/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/';
        addpath(genpath(DIR.codes))
        if(~exist(DIR.out))
            mkdir(DIR.out)
        end
    case{'GLNXA64'}
%         addpath(genpath('H:/MEGEYE/analysis'))
        DIR.ET_data = '/home/juank/Dropbox/big_files_shared/ET_DATA/';
        DIR.out     = '/home/juank/Dropbox/big_files_shared/ET_DATA/out';
        DIR.codes   = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/';
        addpath(genpath(DIR.codes))
        if(~exist(DIR.out))
            mkdir(DIR.out)
        end
    otherwise
        display('Check OS and add paths')
end
%%
subj_ids = {'KS1','KS2','KS3','KS4','KS5','KS6','KS7','KS8'};
sess_ids = {'KS10622';
            'KS20622';
            'KS30622';
            'KS40629';
            {'KS50629';'S530629'}; %'S520629' only contains one trial (disk full)
            'KS60629';
            'KS70629';
            'KS80702'};
%%
for isub=8:8 %length(subj_ids) %loop over subjects
    Nsubj_sessions=size(sess_ids{isub},1); %loop over sessions
    for isess=1:Nsubj_sessions
        if(Nsubj_sessions==1)
            session_name=sess_ids{isub};
        else
            session_name=sess_ids{isub}{isess};
        end
        beh.global_beh_file = fullfile(DIR.ET_data,[session_name '.mat']);
        beh.trials_dir      = fullfile(DIR.ET_data,[session_name '_Trials']);
        listfile_raw=dir(fullfile(beh.trials_dir,[session_name '_Trial*']))  ;
        cd(beh.trials_dir)
        button_key=-1*ones(length(listfile_raw),1);
        button_correct=-1*ones(length(listfile_raw),1);
        true_key=-1*ones(length(listfile_raw),1);
        stim_type=-1*ones(length(listfile_raw),1);
        clear trial_type;clear button_time;clear button_keyCode;clear trial_acc;clear button_trial_timing;
        %trial_acc=[];
        beh_ground_truth = struct;
        for ikk=1:length(listfile_raw)
            [tik tok] = regexp(listfile_raw(ikk).name,'[Trial]+\d+[.mat]');
            itrial    = str2double(listfile_raw(ikk).name(tik+5:tok-1)); %looks for number in the middle
            clear SingleTrial;
            load(listfile_raw(ikk).name);
            trial_acc{ikk}                = itrial;
            %button_key(itrial)            = SingleTrial.resp.key-1; %0=A, 1=P
            button_key(itrial)            = mod(SingleTrial.resp.key+1,2); %0=A, 1=P
            [beh_ground_truth]            = calc_beh_ground_truth(beh_ground_truth,SingleTrial,itrial);   
           if (strcmp(cfg.recond{SingleTrial.info.resptype},'A')) 
               true_key(itrial)           = 0;                      %0=A, 1=P
           elseif(strcmp(cfg.recond{SingleTrial.info.resptype},'P'))
               true_key(itrial)           = 1;
           else
               disp('cannot find true response')
           end            
            if(isfield(SingleTrial.resp,'correct'))
                button_correct(itrial)  = SingleTrial.resp.correct;
            end
            if(isfield(SingleTrial.resp,'keyCode'))
                button_keyCode{itrial}  = SingleTrial.resp.keyCode;
            end
            button_time(itrial)          = SingleTrial.resp.time;
            button_trial_timing(itrial)  = SingleTrial.trial_timing(7)-SingleTrial.trial_timing(5);
            trial_type{itrial}           = SingleTrial.info.trialtype;
            stim_type(itrial)            = SingleTrial.info.stimtype; %1=F, 2=O
        end
        %
        true_key = beh_ground_truth.true_key(1:length(true_key))'; %OVERWRITE true_key because of previous implementation errors
        % Checks if CFS buttons file exists
        buttons_CFS_file=fullfile(DIR.ET_data,['buttons_CFS_' session_name(1:3) '.mat']);
        if(exist(buttons_CFS_file))
            load(buttons_CFS_file);
            buttons_only_CFS=buttons_CFS(:,2);
            right_button_key=mode(buttons_only_CFS);
            if(mod(right_button_key+1,2)==1)
                buttons_only_CFS=mod(buttons_only_CFS+1,2); %0=A, 1=P
            else
                buttons_only_CFS=mod(buttons_only_CFS,2); %0=A, 1=P
            end
            trial_min_size=min(length(buttons_only_CFS),length(true_key));
            button_key=buttons_only_CFS(1:trial_min_size);
        end
        trial_min_size=min(length(button_key),length(true_key));
        button_key=button_key(1:trial_min_size);
        true_key=true_key(1:trial_min_size);
        trial_type=trial_type(1:trial_min_size);
        stim_type=stim_type(1:trial_min_size);
        perf_button1=sum(button_key==true_key & button_key~=-1)/length(button_key); %BUTTONS OPTION 1
        true_key2=ones(length(true_key),1)-true_key;
        perf_button2=sum(button_key==true_key2 & button_key~=-1)/length(button_key); %BUTTONS OPTION 2
        threshold_fixed=0.55;
        if(perf_button1>threshold_fixed)
            disp(['Session: ' session_name '. Normal button box used'])
        elseif(perf_button1<1-threshold_fixed)
            true_key=true_key2;
            disp(['Session: ' session_name '. Inverted button box used'])
        else
           disp(['Session: ' session_name '. Warning: button box used to be investigated!'])
        end 

        button_correct=true_key; %TO DO: check if the original button_correct 
        %doesn't work because of inverted keys or if there's another reason
        %calculate performance with both button options, report and keep
        %the highest value. We can later verify that this is the true
        %answer from the EM data.
        %quick and dirty plot with button responses
        %figure(102);plot(1:length(button_correct),beh_ground_truth.true_key,'ro');hold on;plot(1:length(button_correct),button_key,'bx')
        figure(102);clf;plot(1:length(button_correct),true_key,'ro');hold on;plot(1:length(button_correct),button_key,'bx')
        xlabel('Trial Number','FontSize',18,'FontWeight','bold') 
        ylabel('Response','FontSize',18,'FontWeight','bold')
        set(gca, 'Box', 'off');
        left_button_pressed=sum(button_key==0)/length(button_key);
        right_button_pressed=sum(button_key==1)/length(button_key);
        title(['Session: ' session_name '     %LeftPressed:  ' num2str(left_button_pressed)])
        outfile=fullfile(DIR.out,['buttons_' session_name]); print(gcf,'-dpng',outfile);
        %
        figure;clf
        figure;plot(button_trial_timing,'bo')
        hold on;plot(button_time,'rx')
        legend('button trial timing','button time')
        legend boxoff
        outfile=fullfile(DIR.out,['buttons_time' session_name]); print(gcf,'-dpng',outfile);
        %% eCDF
        clear acc_button_key;
        acc_button_key=zeros(length(button_key ),1);
        acc_VS=0;acc_EX=0;acc_ME=0;
        acc_button_VS=0;acc_button_EX=0;acc_button_ME=0;
        clear perf_trial;
        if(button_key(1)==button_correct(1) && button_key(1)~=-1)
            acc_button_key(1)=1;
        else
            acc_button_key(1)=0;
        end
        switch trial_type{1}
            case('VS')
                acc_VS=acc_VS+1;
                if(button_key(1)==button_correct(1) && button_key(1)~=-1 )
                    acc_button_VS(1)=1;
                else
                    acc_button_VS(1)=0;
                end
                perf_trial(1)           =acc_button_VS(1);
            case('EX')
                acc_EX=acc_EX+1;
                if(button_key(1)==button_correct(1) && button_key(1)~=-1 )
                    acc_button_EX(1)=1;
                else
                    acc_button_EX(1)=0;
                end
                perf_trial(1)           =acc_button_EX(1);
            case('ME')
                acc_ME=acc_ME+1;
                if(button_key(1)==button_correct(1) && button_key(1)~=-1 )
                    acc_button_ME(1)=1;
                else
                    acc_button_ME(1)=0;
                end
                perf_trial(1)           =acc_button_ME(1);
        end
        for ii=2:length(button_key)
            if(button_key(ii)==button_correct(ii) && button_key(ii)~=-1 )
                acc_button_key(ii)=1+acc_button_key(ii-1);
            else
                acc_button_key(ii)=acc_button_key(ii-1);
            end
            switch trial_type{ii}
                case('VS')
                    acc_VS=acc_VS+1;
                    if(button_key(ii)==button_correct(ii) && button_key(ii)~=-1 )
                        acc_button_VS(ii)=1+acc_button_VS(ii-1);
                    else
                        acc_button_VS(ii)=acc_button_VS(ii-1);
                    end
                    perf_trial(ii)           =acc_button_VS(ii)/acc_VS;
                    acc_button_EX(ii)    =acc_button_EX(ii-1);
                    acc_button_ME(ii)    =acc_button_ME(ii-1);
                case('EX')
                    acc_EX=acc_EX+1;
                    if(button_key(ii)==button_correct(ii) && button_key(ii)~=-1 )
                        acc_button_EX(ii)=1+acc_button_EX(ii-1);
                    else
                        acc_button_EX(ii)=acc_button_EX(ii-1);
                    end
                    perf_trial(ii)           =acc_button_EX(ii)/acc_EX;
                    acc_button_VS(ii)    =acc_button_VS(ii-1);
                    acc_button_ME(ii)    =acc_button_ME(ii-1);
                case('ME')
                    acc_ME=acc_ME+1;
                    if(button_key(ii)==button_correct(ii) && button_key(ii)~=-1 )
                        acc_button_ME(ii)=1+acc_button_ME(ii-1);
                    else
                        acc_button_ME(ii)=acc_button_ME(ii-1);
                    end
                    perf_trial(ii)           =acc_button_ME(ii)/acc_ME;
                    acc_button_EX(ii)    =acc_button_EX(ii-1);
                    acc_button_VS(ii)    =acc_button_VS(ii-1);
            end
        end
        acc_button_key=acc_button_key/length(button_key);
        %find specific task trials
        if(subj.cntbalance==1)
            block_borders{1}=[1 20 61 80 121 140]; %VS
            block_borders{2}=[41 60 101 120]; %EX
            block_borders{3}=[21 40 81 100 141 160]; %ME
        else
            block_borders{1}=[21 40 81 100 141 160]; %VS
            block_borders{2}=[41 60 101 120]; %EX
            block_borders{3}=[1 20 61 80 121 140]; %ME
        end
        
        %% [ind_VS] = arrayfun(@(x)find(strmatch(x,'VS')==1),trial_type,'uniformoutput',false)
        figure;clf;
        h1=stairs(1:length(button_key),acc_button_key,'k','LineWidth',4);
        set(gca,'YLim',[0 1]);
        hold on;
        alpha_value=0.6;
        %VS
        x = [block_borders{1}(1) block_borders{1}(1) block_borders{1}(2) block_borders{1}(2)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{1},'FaceAlpha',alpha_value);
        x = [block_borders{1}(3) block_borders{1}(3) block_borders{1}(4) block_borders{1}(4)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{1},'FaceAlpha',alpha_value)
        x = [block_borders{1}(5) block_borders{1}(5) block_borders{1}(6) block_borders{1}(6)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{1},'FaceAlpha',alpha_value);
        %EX
        x = [block_borders{2}(1) block_borders{2}(1) block_borders{2}(2) block_borders{2}(2)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{2},'FaceAlpha',alpha_value);
        x = [block_borders{2}(3) block_borders{2}(3) block_borders{2}(4) block_borders{2}(4)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{2},'FaceAlpha',alpha_value);
        %ME
        x = [block_borders{3}(1) block_borders{3}(1) block_borders{3}(2) block_borders{3}(2)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{3},'FaceAlpha',alpha_value);
        x = [block_borders{3}(3) block_borders{3}(3) block_borders{3}(4) block_borders{3}(4)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{3},'FaceAlpha',alpha_value);
        x = [block_borders{3}(5) block_borders{3}(5) block_borders{3}(6) block_borders{3}(6)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{3},'FaceAlpha',alpha_value);
        ylabel(['Norm Acc Correct, Session: ' session_name]);
        xlabel(['Trial']);
        outfile=fullfile(DIR.out,['cdf_correct_' session_name]); print(gcf,'-dpng',outfile);
        %%         %% [ind_VS] = arrayfun(@(x)find(strmatch(x,'VS')==1),trial_type,'uniformoutput',false)
        figure;clf;
        h1=stairs(1:length(button_key),perf_trial,'k','LineWidth',2);
        hold on;plot(1:length(button_correct),button_key==true_key,'kx')
        set(gca,'YLim',[0 1]);
        hold on;
        alpha_value=0.6;
        %VS
        x = [block_borders{1}(1) block_borders{1}(1) block_borders{1}(2) block_borders{1}(2)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{1},'FaceAlpha',alpha_value);
        x = [block_borders{1}(3) block_borders{1}(3) block_borders{1}(4) block_borders{1}(4)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{1},'FaceAlpha',alpha_value)
        x = [block_borders{1}(5) block_borders{1}(5) block_borders{1}(6) block_borders{1}(6)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{1},'FaceAlpha',alpha_value);
        %EX
        x = [block_borders{2}(1) block_borders{2}(1) block_borders{2}(2) block_borders{2}(2)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{2},'FaceAlpha',alpha_value);
        x = [block_borders{2}(3) block_borders{2}(3) block_borders{2}(4) block_borders{2}(4)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{2},'FaceAlpha',alpha_value);
        %ME
        x = [block_borders{3}(1) block_borders{3}(1) block_borders{3}(2) block_borders{3}(2)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{3},'FaceAlpha',alpha_value);
        x = [block_borders{3}(3) block_borders{3}(3) block_borders{3}(4) block_borders{3}(4)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{3},'FaceAlpha',alpha_value);
        x = [block_borders{3}(5) block_borders{3}(5) block_borders{3}(6) block_borders{3}(6)];
        y = [0 1 1 0];
        fill(x,y,cfg.colortasks{3},'FaceAlpha',alpha_value);
        ylabel(['Performance, Session: ' session_name]);
        xlabel(['Trial']);
        outfile=fullfile(DIR.out,['perfo_trial_' session_name]); print(gcf,'-dpng',outfile);  
        
        %% performance in each task
        hits_ME=0;FA_ME=0;nobutt_ME=0;
        hits_EX=0;FA_EX=0;nobutt_EX=0;
        hits_VS=0;FA_VS=0;nobutt_VS=0;
        for i=1:trial_min_size
            switch trial_type{i}
                case('ME')
                    if(button_key(i)==button_correct(i) && button_key(i)~=-1) %valid and correct response
                        hits_ME=hits_ME+1;
                    elseif(button_key(i)~=button_correct(i) && button_key(i)~=-1) %valid and incorrect
                        FA_ME=FA_ME+1;
                    else
                        nobutt_ME=nobutt_ME+1;
                    end
                case('EX')
                    if(button_key(i)==button_correct(i) && button_key(i)~=-1) %valid and correct response
                        hits_EX=hits_EX+1;
                    elseif(button_key(i)~=button_correct(i) && button_key(i)~=-1) %valid and incorrect
                        FA_EX=FA_EX+1;
                    else
                        nobutt_EX=nobutt_EX+1;
                    end
                case('VS')
                    if(button_key(i)==button_correct(i) && button_key(i)~=-1) %valid and correct response
                        hits_VS=hits_VS+1;
                    elseif(button_key(i)~=button_correct(i) && button_key(i)~=-1) %valid and incorrect
                        FA_VS=FA_VS+1;
                    else
                        nobutt_VS=nobutt_VS+1;
                    end
            end
        end
        perf_valid_ME=hits_ME/(hits_ME+FA_ME);
        perf_all_ME=hits_ME/(hits_ME+FA_ME+nobutt_ME);
        perf_valid_EX=hits_EX/(hits_EX+FA_EX);
        perf_all_EX=hits_EX/(hits_EX+FA_EX+nobutt_EX);
        perf_valid_VS=hits_VS/(hits_VS+FA_VS);
        perf_all_VS=hits_VS/(hits_VS+FA_VS+nobutt_VS);
        %% performance for Faces (stim_type=1) and Objects (stim_type=2)
        %VS
        VS_trials = strmatch('VS', trial_type);
        %F
        hits_VS_F=sum(button_key(VS_trials)==button_correct(VS_trials) & button_key(VS_trials)~=-1 &stim_type(VS_trials)==1);
        FA_VS_F=sum(button_key(VS_trials)~=button_correct(VS_trials) & button_key(VS_trials)~=-1 &stim_type(VS_trials)==1);
        perf_valid_VS_F=hits_VS_F/(hits_VS_F+FA_VS_F);
        %O
        hits_VS_O=sum(button_key(VS_trials)==button_correct(VS_trials) & button_key(VS_trials)~=-1 &stim_type(VS_trials)==2);
        FA_VS_O=sum(button_key(VS_trials)~=button_correct(VS_trials) & button_key(VS_trials)~=-1 &stim_type(VS_trials)==2);
        perf_valid_VS_O=hits_VS_O/(hits_VS_O+FA_VS_O);
        %EX
        EX_trials = strmatch('EX', trial_type);
        %F
        hits_EX_F=sum(button_key(EX_trials)==button_correct(EX_trials) & button_key(EX_trials)~=-1 &stim_type(EX_trials)==1);
        FA_EX_F=sum(button_key(EX_trials)~=button_correct(EX_trials) & button_key(EX_trials)~=-1 &stim_type(EX_trials)==1);
        perf_valid_EX_F=hits_EX_F/(hits_EX_F+FA_EX_F);
        %O
        hits_EX_O=sum(button_key(EX_trials)==button_correct(EX_trials) & button_key(EX_trials)~=-1 &stim_type(EX_trials)==2);
        FA_EX_O=sum(button_key(EX_trials)~=button_correct(EX_trials) & button_key(EX_trials)~=-1 &stim_type(EX_trials)==2);
        perf_valid_EX_O=hits_EX_O/(hits_EX_O+FA_EX_O);
        %ME
        ME_trials = strmatch('ME', trial_type);
        %F
        hits_ME_F=sum(button_key(ME_trials)==button_correct(ME_trials) & button_key(ME_trials)~=-1 &stim_type(ME_trials)==1);
        FA_ME_F=sum(button_key(ME_trials)~=button_correct(ME_trials) & button_key(ME_trials)~=-1 &stim_type(ME_trials)==1);
        perf_valid_ME_F=hits_ME_F/(hits_ME_F+FA_ME_F);
        %O
        hits_ME_O=sum(button_key(ME_trials)==button_correct(ME_trials) & button_key(ME_trials)~=-1 &stim_type(ME_trials)==2);
        FA_ME_O=sum(button_key(ME_trials)~=button_correct(ME_trials) & button_key(ME_trials)~=-1 &stim_type(ME_trials)==2);
        perf_valid_ME_O=hits_ME_O/(hits_ME_O+FA_ME_O);
        %% TEMPORARY CALCULATIONS
%         F_trials = find(stim_type==1);
%         O_trials = find(stim_type==2);
%         EX_F_trials = intersect(EX_trials,F_trials);
%         EX_O_trials = intersect(EX_trials,O_trials);
%         beh_ground_truth.check_resptype(EX_O_trials)
%         beh_ground_truth.gender_size(EX_trials)
%         
%         gender_size = beh_ground_truth.gender_size;
%         for ifoo=1:4
%             footrials = cell2mat(gender_size)==ifoo & isEX_trial;  
%             sum(beh_ground_truth.check_resptype(footrials))
%         end
%         
%         isVS_trial = cell2mat(cellfun(@(x) strcmp(x,'VS'), trial_type, 'UniformOutput',false));
%         isEX_trial = cell2mat(cellfun(@(x) strcmp(x,'EX'), trial_type, 'UniformOutput',false));
%         isME_trial = cell2mat(cellfun(@(x) strcmp(x,'ME'), trial_type, 'UniformOutput',false));


        
        
        %% figure with performance for valid and all trials
        figure(103);clf;
        y1 = [perf_valid_VS perf_valid_EX perf_valid_ME];
        y2 = [perf_all_VS perf_all_EX perf_all_ME]; 
        subplot(1,2,1)
        hold on
        for itask=1:length(cfg.tasks)
            bar(itask,y1(itask),0.3,'FaceColor',cfg.colortasks{itask},'EdgeColor',[0 0 0],'LineWidth',1.5);
        end
        hold off
        title('Valid Trials')
        ylabel(['Performance Session: ' session_name])
        set(gca,'XLim',[0 4],'visible','off')
        set(gca,'YLim',[0 1],'visible','on')
        set(gca,'XTickLabels',[])
        yticks([0 0.2 0.4 0.6 0.8 1.0])
        set(gca,'YTickLabels',[0 0.2 0.4 0.6 0.8 1.0])
        legend(cfg.tasks)
        legend boxoff
        subplot(1,2,2)
        hold on
        for itask=1:length(cfg.tasks)
            bar(itask,y2(itask),0.3,'FaceColor',[1 1 1],'EdgeColor',cfg.colortasks{itask},'LineWidth',1.5);
        end
        hold off
        title('All Trials')
        set(gca,'XLim',[0 4],'visible','off')
        set(gca,'YLim',[0 1],'visible','on')
        set(gca,'XTickLabels',[])
        yticks([0 0.2 0.4 0.6 0.8 1.0])
        set(gca,'YTickLabels',[0 0.2 0.4 0.6 0.8 1.0])
        legend(cfg.tasks)
        legend boxoff
        outfile=fullfile(DIR.out,['performance_' session_name]); print(gcf,'-dpng',outfile);
         %% figure with performance for valid trials faces and objects
        figure(104);clf;
        y1 = [perf_valid_VS_F perf_valid_EX_F perf_valid_ME_F];
        y2 = [perf_valid_VS_O perf_valid_EX_O perf_valid_ME_O]; 
        subplot(1,2,1)
        hold on
        for itask=1:length(cfg.tasks)
            %bar(itask,y1(itask),0.3,'FaceColor',cfg.colortasks{itask},'EdgeColor',[0 0 0],'LineWidth',1.5);
            bar(itask,y1(itask),0.4,'FaceColor',[1 0 1],'EdgeColor',cfg.colortasks{itask},'LineWidth',3); %magenta F
        end
        hold off
        title('Valid FACE Trials')
        ylabel(['Performance Session: ' session_name])
        set(gca,'XLim',[0 4],'visible','off')
        set(gca,'YLim',[0 1],'visible','on')
        set(gca,'XTickLabels',[])
        yticks([0 0.2 0.4 0.6 0.8 1.0])
        set(gca,'YTickLabels',[0 0.2 0.4 0.6 0.8 1.0])
        legend(cfg.tasks)
        legend boxoff
        subplot(1,2,2)
        hold on
        for itask=1:length(cfg.tasks)
            %bar(itask,y2(itask),0.3,'FaceColor',[1 1 1],'EdgeColor',cfg.colortasks{itask},'LineWidth',1.5);
            bar(itask,y2(itask),0.4,'FaceColor',[0 1 1],'EdgeColor',cfg.colortasks{itask},'LineWidth',3); %cyan O
        end
        hold off
        title('Valid OBJECT Trials')
        set(gca,'XLim',[0 4],'visible','off')
        set(gca,'YLim',[0 1],'visible','on')
        set(gca,'XTickLabels',[])
        yticks([0 0.2 0.4 0.6 0.8 1.0])
        set(gca,'YTickLabels',[0 0.2 0.4 0.6 0.8 1.0])
        legend(cfg.tasks)
        legend boxoff
        outfile=fullfile(DIR.out,['perf_F_OBJ_' session_name]); print(gcf,'-dpng',outfile);
    end
end
