function [runpath,code_path,session_path,whoisrunning]=add_paths_matlab_EEG(print_info)
%add_paths_matlab(print_info)
%  
if(nargin == 0)
    print_info = 1;
end
[~,whoisrunning] = system('whoami');

switch strtrim(whoisrunning)
    case 'duip66033\lpzmi_local'
        dropbox_path='D:/Dropbox';
        local_path  ='H:/EEGEYE';
        runpath=fullfile(local_path,'data');
        code_path.fieldtrippath=fullfile(dropbox_path,'big_files_shared','fieldtrip-20170827');
        %code_path.fieldtrippath=fullfile('D:\toolbox','fieldtrip-20180711');
        code_path.eeglabpath   =fullfile(dropbox_path,'toolbox','eeglab13_6_5b');
        code_path.my_functions=fullfile(local_path,'my_functions');
        session_path.raw=fullfile(local_path,'eeg_data');
        %         session_path.directories = {'';
        %                '';
        %                };
        session_path.rawet      = fullfile(local_path,'data');
        session_path.matfiles   = fullfile(local_path,'data');
        session_path.subjname = {'E01','E03','E04','E05','E06','E07','E08','E09',...
            'E10','E11','E12','E13','E14'};
        session_path.sessionfilenames = {'E010620';'E030622';'E040625';'E050625';'E060626';...
            'E070627';'E080628';'E090702';'E100703';...
            'E110703';'E120703';'E130705';'E140709'};
        session_path.Trialfilenames = {'S10620_Trials';'E010622_Trials';'E040625_Trials';'E050625_Trials';'E060626_Trials';...
            'E070627_Trials';'E080628_Trials';'E080702_Trials';'E100703_Trials';...
            'E110703_Trials';'E120703_Trials';'E13missing';'E140709'};
        
    case 'root' %'juank-Desktop'
        dropbox_path='/home/juank/Dropbox';
        local_path  ='/home/juank/Dropbox/big_files_shared/EEGEYE';
        runpath=fullfile(local_path,'data');
        code_path.fieldtrippath='/home/juank/toolbox/fieldtrip-20170827';
        code_path.eeglabpath   ='/home/juank/toolbox/eeglab14_1_1b';
        code_path.my_functions='my_functions';
        session_path.raw=fullfile(local_path,'eeg_data');
        %         session_path.directories = {'';
        %                '';
        %                };
        session_path.rawet      = fullfile(local_path,'data');
        session_path.matfiles   = fullfile(local_path,'data');
        session_path.subjname = {'E01','E03','E04','E05','E06','E07','E08','E09',...
            'E10','E11','E12','E13','E14'};
        session_path.sessionfilenames = {'E010620';'E030622';'E040625';'E050625';'E060626';...
            'E070627';'E080628';'E090702';'E100703';...
            'E110703';'E120703';'E130705';'E140709'};
        session_path.Trialfilenames = {'S10620_Trials';'E010622_Trials';'E040625_Trials';'E050625_Trials';'E060626_Trials';...
            'E070627_Trials';'E080628_Trials';'E080702_Trials';'E100703_Trials';...
            'E110703_Trials';'E120703_Trials';'E13missing';'E140709'};
        
    otherwise
        warning(['User unknown: ' whoisrunning ' , check paths']);
end

if(print_info)
    disp(whoisrunning)
    disp(code_path)
    disp(runpath)
end

end

