function [xc_xy,mxcf,xcf,tlags,t1,z1,t2,z2,z2i]=fun_xcorr_lag(tr,i,data,indch,eyedata)
    % Analog Channel MEG data
%     t1 = data.time{tr};
%     z1 = data.trial{tr}(indch(i),:);
    t1 = data.time{tr}(data.time{tr}>=0);
    z1 = data.trial{tr}(indch(i),data.time{tr}>=0);
    % Eye Tracking data
    t2 = eyedata(tr).samples(:,1)/1000;
    z2 = eyedata(tr).samples(:,i+1);
    % Interpolate the Eye Tracking data into the MEG timepoints
    z2i = interp1(t2,z2,t1);

    % I'm going to interpolate the NaN values with a rect only for
    % the lag estimation
    if any(isnan(z2i))
        indini = find(diff(isnan(z2i))==1);
        indend = find(diff(isnan(z2i))==-1);
        if ( isempty(indini) );  indini = 1;  end
        if ( isempty(indend) );  indend = length(z2i)-1;  end

        if ( indini(1)  > indend(1) );  indini = [1 indini];            end
        if ( indini(end)> indend(end) );indend = [indend length(z2i)-1];end

        for j = 1:length(indini)
            ti = t1(indini(j));
            tf = z2i(indini(j));
            yi = t1(indend(j)+1);
            yf = z2i(indend(j)+1);

            if isnan(yi); yi=yf; end
            if isnan(yf); yf=yi; end

            a = (yf-yi)/(tf-ti);
            b = yf - a*tf;

            z2i(indini(j):(indend(j)+1)) = a*t1(indini(j):(indend(j)+1)) + b;
        end
    end

    % Calculate the maximum of the cross-correlation 
    [xcf,lags] = xcorr( (z2i-mean(z2i))/std(z2i) , (z1-mean(z1))/std(z1));
    tlags = 1000*lags/data.fsample;
    [mxcf,imxcf]=max(xcf);
%     xc_xy(tr,i) = tlags(imxcf);
    xc_xy = tlags(imxcf);
end