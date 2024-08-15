function output = process_terasmart(file,time_axis,td_bg,w_mode,w_center,w_length,filter_bd,filter_width,padding,padding_dir,denoise)

% This code is for processing averaged or non-averaged data by Terasmart
% LabVIEW program.

% The processing includs loading, interpolation, low/high pass filter,
% window, signal padding and fourier transmform

% [time_axis]: give a defined time axis (will be auto shortened), or input
% 0 to be automatically generated

% [td_bg] is an array of time-domain background. The length of each row
% should be identical to the time axis. The number of rows should be equal
% to the number of data read from the file.

% [w_mode]=0: no window function applied; =1:normal chebwin function;
% =2:sharpened chebwin function; =3:flat pass chebwin

% [w_center]: specify the data point (of the interpolated data) as the
% center for the window function. =0 to be automatically generated

% [w_length]: window length in ps unit: ~6-15 ps. Input 0 to be
% autoamtically determined according to the available length (not recommended)

% [filter_bd]=[low_bd, high_bd], the passing bandwidth for the low/high
% pass filter (transmission reaches 1)

% [filter_width]=[low_width, high_width], the bandwidth the let the filter
% falls to e^-1

% [paading]: the total number of data points desired for the padded TD data

% [padding_dir]: padding direction. =0: padding 0 in front, =1: 0 after

% [denoise]: =1 to use wavelet denoising

%%%% template
% time_axis = 0;
% w_mode = 0; w_center = 0; w_atten = 0;
% filter_bd = [0,0]; filter_width = [0,0];
% padding = 0; padding_dir = 0;
% denoise = 0;

%% check input
if nargin < 11
    denoise = 0;
    if nargin < 9
        padding = 0;
        padding_dir = 0;
        if nargin < 7
            filter_bd = [0,0];
            filter_width = [0,0]; 
            if nargin <4
                w_mode = 0;
                w_center = 0;
                w_length = 0;
                if nargin <3
                    td_bg = 0;
                    if nargin <2
                        time_axis = 0;
                    end
                end
            end
        end
    end
end

%% load data
data_load_itp = read_terasmart(time_axis,file); 
num_rows = numel(data_load_itp.data_itp(:,1));
if numel(w_center) == 1
    w_center = linspace(w_center, w_center, num_rows);
end
if numel(w_length) == 1
    w_length = linspace(w_length, w_length, num_rows);
end
%% build filter for later use
time_axis=data_load_itp.data_itp(1,:);
td_length = numel(time_axis);
maxF=1/(time_axis(2)-time_axis(1));
freq=linspace(0,maxF/2,numel(rfft(time_axis)));

if filter_bd(1)~=0
    LP_filter=exp(-((freq-filter_bd(1)).^2)./filter_width(1)^2);
else
    LP_filter=linspace(1,1,numel(freq));
end

if filter_bd(2)~=0
    HP_filter=exp(-((freq-filter_bd(2)).^2)./filter_width(2)^2);
else
    HP_filter=linspace(1,1,numel(freq));
end

filter_applied=LP_filter+HP_filter;
filter_applied(filter_applied>1)=1;
filter_applied(freq>filter_bd(1) & freq<filter_bd(2))=1;

%% process each pulse

data_td_filtered = zeros(num_rows, td_length);
signal_win = zeros(num_rows, td_length);
if padding<td_length || w_mode==0
    signal_padding = zeros(num_rows, td_length);
    data_fd = zeros(num_rows, numel(rfft(time_axis)));
else
    signal_padding = zeros(num_rows, padding);
    data_fd = zeros(num_rows, numel(rfft(linspace(1,1,padding))));
end

w_center_out = linspace(0,0,num_rows-1);

for k=2:num_rows % k starts from 2!!
    %% wavelet denoising
    if denoise==1
        data_load_itp.data_itp(k,:)=wden(data_load_itp.data_itp(k,:),'modwtsqtwolog','s','mln',4,'sym4');
    end

    %% background subtration
    if numel(td_bg)~=1
        data_load_itp.data_itp(k,:) = data_load_itp.data_itp(k,:)-td_bg(k-1,:);
    end

    %% LP / HP filter
    data_fd_filtered=rfft(data_load_itp.data_itp(k,:)).*filter_applied;
    data_td_filtered(k,:)=irfft(data_fd_filtered,numel(time_axis));

    %% TD window
    w_atten=130; % enough to make window edges decayed to zero
    
    % find the window center
    if w_center(k-1)<4 && w_mode~=0 % a window center needs to be found
        w_center_use = win_center_hilbert(data_td_filtered(k,:), w_center(k-1));
    else % window center has been defined
        w_center_use = w_center(min([numel(w_center),k-1]));
    end

    % define the window length
    if w_mode~=0
        if w_length(k-1)==0 % window length not defined
            end_length = numel(time_axis)-w_center_use-1;
            if w_mode~=3 % symmetric window without a flat top
                front_length = w_center_use-1;
                length_15ps = round(15/(time_axis(2)-time_axis(1)));
                w_length(k-1) = min([end_length,front_length,length_15ps]);
            else % asymetric winow with a flat passing top
                w_length(k-1) = end_length; % This length correspond to the flat length
            end
        else % window length already defined but in the unit of ps
            w_length(k-1) = fix(w_length(k-1)/(time_axis(2)-time_axis(1))); % convert time length to index
        end
        
        % define the window function
        if w_mode ==1
            win_ftn=chebwin(2*w_length(k-1),w_atten);
            w_peak_pos = w_length(k-1)+1;
        elseif w_mode==2
            win_ftn0=chebwin(2*w_length(k-1),w_atten);
            win_ftn=0.5*sin(win_ftn0*pi-pi/2)+0.5; % modify to make it sharpper but flatter on top
            w_peak_pos = w_length(k-1)+1;
        elseif w_mode==3
            length_3ps = round(3/(time_axis(2)-time_axis(1)));
            win_ftn0=chebwin(2*length_3ps, w_atten);
            win_ftn=[win_ftn0(1:length_3ps); ones(w_length(k-1),1); win_ftn0(length_3ps+1:end)];
            w_peak_pos = 2*length_3ps; % i.e. there will be 3ps flat window before w_center
        end
        front_zeros = w_center_use-w_peak_pos; % if this is negative it will produce an error
        end_zeros = numel(data_td_filtered(k,:))-numel(win_ftn)-front_zeros;
        if front_zeros<1 || end_zeros<1
            msg = 'signal too short to apply a window'
        end
        win_ftn_extend = [zeros(front_zeros, 1); win_ftn; zeros(end_zeros, 1)];
        signal_win(k,:) = data_td_filtered(k,:).*win_ftn_extend';
        td_changed = 1;
    else
        signal_win(k,:)=data_td_filtered(k,:);
        td_changed = 0;
    end
    w_center_out(k-1) = w_center_use;

    %% zero padding and FFT
    time_step=time_axis(2)-time_axis(1);
    if padding<td_length || w_mode==0
        signal_padding(k,:)=signal_win(k,:);
        time_padding=time_axis;
    elseif padding_dir==0
        signal_padding(k,:)=[linspace(0,0,padding-td_length),signal_win(k,:)];
        time_padding=0:time_step:(padding-1)*time_step;
    else
        signal_padding(k,:)=[signal_win(k,:),linspace(0,0,padding-td_length)];
        time_padding=0:time_step:(padding-1)*time_step;
    end
    
    if td_changed == 1
        data_fd(k,:) = rfft(signal_padding(k,:));
    else
        data_fd(k,:) = data_fd_filtered;
    end
    
end

data_td_filtered(1,:) = time_axis;
signal_win(1,:) = time_axis;
data_fd(1,:)=linspace(0,maxF/2,numel(data_fd(2,:)));
signal_padding(1,:) = time_padding;

%% data output
output.data_read=data_load_itp.data_read;
output.data_itp=data_load_itp.data_itp;
output.data_filter=data_td_filtered;
output.data_win=signal_win;
output.data_padding=signal_padding;
output.data_fd=data_fd;
output.win_center=w_center_out;
output.win_length=w_length;
output.filter=filter_applied;