function output=read_terasmart(time_axis_define,file_name)

format_flag = isempty(intersect(file_name(end-3:end),'x')); % =0 for txt file
if format_flag == 0
    data_read=load(file_name);
else  % .mat file
    data_read=load(file_name);
    data_read = data_read.signal;
end

% define the time axis
if numel(time_axis_define)>1
    time_axis=time_axis_define;
    if time_axis(end)>data_read(1,end)+0.1 % +0.1 to avoid the two values are extremely close but still larger
        cut_ind=find(abs(time_axis-(round(data_read(1,end))-1))==min(abs(time_axis-(round(data_read(1,end))-1))));
        time_axis=time_axis(1:cut_ind);
    end
    % time-domain interpolation
    data_itp(1,:)=time_axis;
    for k=2:numel(data_read(:,1))
        data_itp(k,:)=interp1(data_read(1,:),data_read(k,:),data_itp(1,:),'spline', 'extrap');
    end
else
    data_itp=data_read;
end

output.data_read=data_read;
output.data_itp=data_itp;