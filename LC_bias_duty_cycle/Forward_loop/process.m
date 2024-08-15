clear all

file_num1 = [0:7,9];
file_num2 = [11:4:27,37,47,67,87,999];

data_num = numel(file_num1)+numel(file_num2);
colorset = ColorGradient(data_num, 'full');

for k = 1:data_num
    if k<=numel(file_num1)
        sam_file = ['0.0',num2str(file_num1(k)), '_LC.mat'];
        ref_file = ['0.0',num2str(file_num1(k)), '_air.mat'];
    else
        sam_file = ['0.',num2str(file_num2(k-numel(file_num1))), '_LC.mat'];
        ref_file = ['0.',num2str(file_num2(k-numel(file_num1))), '_air.mat'];
    end
    sam_load = process_terasmart(sam_file, 0,0,2,0,4);
    ref_load = process_terasmart(ref_file, 0,0,2,0,4);
    
    freq = sam_load.data_fd(1,:);
    sam_ps = sam_load.data_fd(2,:) ./ sam_load.data_fd(3,:);
    ref_ps = ref_load.data_fd(2,:) ./ ref_load.data_fd(3,:);
    ps_norm(k,:) = sam_ps./ref_ps;
    ps_phase(k) = unwrap_TDS(freq, ps_norm(k,:), [0.4,1.8], 3);
end

figure
plot([file_num1,file_num2(1:end-1),100], ps_phase, 'o-', 'linewidth', 1)

%%
forward_dc = [file_num1,file_num2(1:end-1),100];
forward_phase = ps_phase;
save('forward_bias.mat', 'forward_dc', 'forward_phase');