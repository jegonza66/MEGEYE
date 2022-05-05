function [runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG(print_info)
%add_paths_matlab(print_info)
%  
if(nargin == 0)
    print_info = 1;
end
[~,whoisrunning] = system('whoami')

switch strtrim(whoisrunning)
    case 'duip66033\lpzmi_local'
        main_path                       = 'D:/Dropbox';
        runpath                         = fullfile(main_path,'fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/');
        code_path.fieldtrippath         = fullfile(main_path,'big_files_shared','fieldtrip-20170827');
        %code_path.fieldtrippath=fullfile('D:\toolbox','fieldtrip-20180711');
        code_path.my_functions          = fullfile(runpath,'my_functions');
        %session_path.raw=fullfile('H:\MEGEYE\CTF_DATA');
        
        session_path.raw        = fullfile(main_path,'big_files_shared/CTF_DATA');
        session_path.rawet      = fullfile(main_path,'big_files_shared/ET_DATA');
        session_path.matfiles   = fullfile(main_path,'big_files_shared/MAT_DATA');
        
    case 'ad\lpajg1' % Joac Notts
        main_path                       = 'H:\';
        runpath                         = fullfile('H:\','MEGEYE\MEG_analysis\');
        code_path.fieldtrippath         = fullfile('H:\fieldtrip-20220104');
        code_path.my_functions          = fullfile(runpath,'my_functions');
        
        session_path.raw        = fullfile(main_path,'MEGEYE\CTF_DATA');
        session_path.rawet      = fullfile(main_path,'MEGEYE\ET_DATA');
        session_path.matfiles   = fullfile(main_path,'MEGEYE\MAT_DATA');
    
    case 'laptop-5i5qsv76\joaco' % Joac Asus
        data_path                       = 'E:\Doc\MEGEYE';
        runpath                         = fullfile('H:\','MEGEYE\MEG_analysis\');
        code_path.fieldtrippath         = fullfile('H:\fieldtrip-20220104');
        code_path.my_functions          = fullfile(runpath,'my_functions');
        
        session_path.raw        = fullfile(data_path,'\CTF_DATA');
        session_path.rawet      = fullfile(data_path,'\ET_DATA');
        session_path.matfiles   = fullfile(data_path,'\MAT_DATA');
        
    case 'root' %'juank-Desktop'
        main_path                       = '~/Dropbox';
        runpath                         = fullfile(main_path,'fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/');
        code_path.fieldtrippath         = fullfile('/home/juank/toolbox/fieldtrip-20170827');
        code_path.my_functions          = fullfile(runpath,'my_functions');
        
        session_path.raw        = fullfile(main_path,'big_files_shared/CTF_DATA');
        session_path.rawet      = fullfile(main_path,'big_files_shared/ET_DATA');
        session_path.matfiles   = fullfile(main_path,'big_files_shared/MAT_DATA');

    otherwise
        warning(['User unknown: ' whoisrunning ' , check paths']);
end

session_path.directories = {''; ''};
% session_path.sessionfilenames = {'13434001_MarkusBauer_20180629_01.ds';
%             '13434001_MarkusBauer_20180629_02.ds';
%             };

session_path.subjname = {'11214-008','13229-001',...
                        'KS1','KS2','KS3','KS4','KS5','KS6','KS7','KS8', ...
                        '11804','13130','13416','13555','13556','13557',...
                        '13575','15028_01','10925_01','13444_01','11535_01','14343_01',...
                        '13891_01','13993_01','14359_01'};
session_path.subjcode = {'11214-008','13229-001','13891-002','13993-004'...
                        '10845-017','13397-001','13392-001','13433-001',...
                        '11081-018','13434-001','10335-003','13129-002',...
                        '11804_01','13130_01','13416_01','13555_04','13556_03','13557_01',...
                        '13575_002','15028_01','10925_01','13444_01','11535_01','14343_01',...
                        '13891_01','13993_01','14359_01'};
session_path.sessionfilenames = {{''};{''}; ...
    {'10845-17_MarkusBauer_20180622_01.ds',  '10845-17_MarkusBauer_20180622_02.ds'};...
    {'13392001_MarkusBauer_20180622_02.ds'};...
    {'13392001_MarkusBauer_20180622_01.ds',  '13392001_MarkusBauer_20180622_02.ds',  '13392001_MarkusBauer_20180622_03.ds',  '13392001_MarkusBauer_20180622_04.ds'};...
    {'1343301_MarkusBauer_20180629_01.ds',   '1343301_MarkusBauer_20180629_02.ds'};...
    {'11081018_MarkusBauer_20180629_01.ds',  '11081018_MarkusBauer_20180629_02.ds'};...
    {'13434001_MarkusBauer_20180629_01.ds',  '13434001_MarkusBauer_20180629_02.ds'};...
    {'10335003_MarkusBauer_20180629_01.ds',  '10335003_MarkusBauer_20180629_02.ds',  '10335003_MarkusBauer_20180629_03.ds'};...
    {'13129002_MarkusBauer_20180702_01.ds',  '13129002_MarkusBauer_20180702_02.ds'};...
    {'11804007_MarkusBauer_20180718_01.ds'};...
    {'13130002_MarkusBauer_20180718_01.ds',  '13130002_MarkusBauer_20180718_02.ds',  '13130002_MarkusBauer_20180718_03.ds',  '13130002_MarkusBauer_20180718_04.ds', '13130002_MarkusBauer_20180718_05.ds',  '13130002_MarkusBauer_20180718_06.ds'};...
    {'13416003_MarkusBauer_20180718_01.ds',  '13416003_MarkusBauer_20180718_02.ds',  '13416003_MarkusBauer_20180718_03.ds'};...
    {'13555001_MarkusBauer_20180730_01.ds','13555001_MarkusBauer_20180730_02.ds','13555001_MarkusBauer_20180730_03.ds'};...
    {'13556001_MarkusBauer_20180730_01.ds','13556001_MarkusBauer_20180730_02.ds','13556001_MarkusBauer_20180730_03.ds'};...
    {'13556001_MarkusBauer_20180730_01.ds','13556001_MarkusBauer_20180730_02.ds','13556001_MarkusBauer_20180730_03.ds'};...
    {'13575002_MarkusBauer_20180820_01.ds','13575002_MarkusBauer_20180820_02.ds','13575002_MarkusBauer_20180820_03.ds','13575002_MarkusBauer_20180820_04.ds','13575002_MarkusBauer_20180820_05.ds'};...
    {'15028005_MarkusBauer_20180807_01.ds'};...
    {'10925041_MarkusBauer_20180913_01.ds'};...
    {'13444001_MarkusBauer_20190710_01.ds'};...
    {'13444001_MarkusBauer_20190710_02.ds'};...
    {'13444001_MarkusBauer_20190710_03.ds'};...
    {'11535003_MarkusBauer_20190711_01.ds'};...
    {'11535003_MarkusBauer_20190711_02.ds'};...
    {'11535003_MarkusBauer_20190711_03.ds'};...
    {'14343001_MarkusBauer_20190711_01.ds'};...
    {'14343001_MarkusBauer_20190711_02.ds'};...
    {'14343001_MarkusBauer_20190711_03.ds'};...
    {'13891002_MarkusBauer_20190717_01.ds'};...
    {'13891002_MarkusBauer_20190717_02.ds'};...  
    {'13891002_MarkusBauer_20190717_03.ds'};...  
    {'13993004_MarkusBauer_20190717_01.ds'};...
    {'13993004_MarkusBauer_20190717_02.ds'};...
    {'14359001_MarkusBauer_20190718_01.ds'};...          
    };
session_path.behavfilenames = {{''}; {''}; ...
    {'KS10622'};...
    {'KS20622'};...
    {'KS30622'};...
    {'KS40629'};...
    {'KS50629','S530629'};...
    {'KS60629'};...
    {'KS70629'};...
    {'KS80702'};...
    {'11804_01'};...
    {'13130_01'};...
    {'13416_01'};...
    {'13555_04'};...
    {'13556_03'};...
    {'13557_01'};
    {'13575_002'};...
    {'15028_01'};...
    {'10925_01'};...
    {'13444_01'};...
    {'11535_01'};...
    {'14343_01'};...
    {'13891_01'};...
    {'13993_01'};...
    {'14359_01'};...
    };     
if(print_info)
    disp(whoisrunning)
    disp(code_path)
    disp(runpath)
end

end

