%% Script to check button presses

%% ---------------------- Keyboard -----------------------------
KbName('UnifyKeyNames');
cfg.escapeKey   = KbName('ESCAPE');
cfg.pauseKey    = KbName('p');
%cfg.leftKey     = [KbName('1') KbName('2') KbName('3') KbName('4') KbName('LeftArrow')]; %KbName('LeftArrow');
%cfg.rightKey    = [KbName('6') KbName('7') KbName('8') KbName('9') KbName('RightArrow')]; %KbName('RightArrow');
cfg.leftKey     = KbName('6^'); %left hand
cfg.rightKey     = KbName('1!'); %right hand
%%
maxtime=10;
exitcond = 1;
tStart  = GetSecs;
t       = GetSecs - tStart;
response= 0;
exitflag= 0;
if (exitcond==1);       cond = (t < maxtime) && (response==0);
elseif (exitcond==0);   cond = (t < maxtime);
end
while ( cond )
    if (exitcond==1);       cond = (t < maxtime) && (response==0);
    elseif (exitcond==0);   cond = (t < maxtime);
    end
    % Check the keyboard. The person should press the
    [~,~, keyCode] = KbCheck;
    if keyCode(cfg.escapeKey)
        resp.time           = nan;
        resp.key            = nan;
        ShowCursor;
        sca;
        exitflag = 1;
        break
    elseif sum(keyCode(cfg.rightKey))
        response = 1;
        exitflag = 2;
    elseif sum(keyCode(cfg.leftKey))
        response = 2;
        exitflag = 2;
    end
    t       = GetSecs - tStart;
    WaitSecs(0.001);
end
tEnd = GetSecs;
resp.time           = tEnd - tStart
resp.key            = response
resp.keyCode        = find(keyCode)

