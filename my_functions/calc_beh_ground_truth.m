function [beh_ground_truth]       = calc_beh_ground_truth(beh_ground_truth,SingleTrial,itrial)
%function [beh_ground_truth]       = calc_beh_ground_truth(SingleTrial)
%this function calculates the ground truth for the generated responses
%based on the actual stimuli displayed (taken from the file SingleTrial)
%stimtype:          1=F, 2=O
%resptype:          1=A, 2=P
%gendertype:        1=M, 2=F
%sizetype:          1=S, 2=B
%gender_size:       1=M, 2=F, 3=S, 4=B
%MJI, 11/7/2-18

trialtype               = SingleTrial.info.trialtype;
stimtype                = SingleTrial.info.stimtype;
resptype                = SingleTrial.info.resptype;
target_A                = SingleTrial.info.indtarget==0;

[tik tok] = regexp(SingleTrial.info.filenames_target,'\w*\w*\w*','match');
%the first element of tik is 'experiment', the second one 'faces' or
%'objects', the third one starts with male_ female_ small_ big_
if(strcmp(tik{3},'refined')) %for 2019 images, it detects the extra 'refined' field
    itc=4;
else
    itc=3;
end
if(strcmp(tik{2},'faces'))
    if(stimtype==1)
        beh_ground_truth.check_stimtype(itrial) = 1;
    else
        beh_ground_truth.check_stimtype(itrial) = 0;
    end
    gendertik = regexp(tik{itc},'female_','match');
    if(~isempty(gendertik))
        beh_ground_truth.gendertype{itrial} = 2; %'female' string found
        beh_ground_truth.gender_size{itrial}= 2; 
    else
        beh_ground_truth.gendertype{itrial} = 1; %'female' string not found
        beh_ground_truth.gender_size{itrial}= 1;
    end
elseif(strcmp(tik{2},'objects'))
    if(stimtype==2)
        beh_ground_truth.check_stimtype(itrial) = 1;
    else
        beh_ground_truth.check_stimtype(itrial) = 0;
    end
    gendertik = regexp(tik{itc},'big_','match');
    if(~isempty(gendertik))
        beh_ground_truth.sizetype{itrial} = 2; %'big' string found
        beh_ground_truth.gender_size{itrial}= 4;
    else
        beh_ground_truth.sizetype{itrial} = 1; %'big' string not found
        beh_ground_truth.gender_size{itrial}= 3;
    end
else
    beh_ground_truth.check_stimtype(itrial) = 0;
end
%check A label (resptype=1) or P (resptype=2) label
if(resptype==1 && target_A==1)
    beh_ground_truth.check_resptype(itrial) = 1;
elseif(resptype==2 && target_A==0)
    beh_ground_truth.check_resptype(itrial) = 1;
else
    beh_ground_truth.check_resptype(itrial) = 0;
end
%Estimate current responses in EX trials
if(strcmp(trialtype,'EX'))
    %beh_ground_truth.true_key(itrial) = mod(beh_ground_truth.gender_size{itrial}+1,2);
    beh_ground_truth.true_key(itrial) = mod(beh_ground_truth.gender_size{itrial}+1,2);
else
    beh_ground_truth.true_key(itrial) = mod(resptype+1,2);
end
end