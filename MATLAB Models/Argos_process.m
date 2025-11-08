%need to download Argos CSI data from their website
close all
clear
path = "ArgosCSI-96x7-2016-11-04-00-57-12_5GHz_track_left_to_right_outdoor_clocksync.hdf5";
info = h5info(path);
disp(info);

FrameCompleteTime = h5read(path,"/FrameCompleteTime");
Pilot_Samples = h5read(path,"/Pilot_Samples");
RSSI = h5read(path,"/RSSI");
% save("ArgosCSI-96x7-2016-11-04-00-57-12_5GHz_track_left_to_right_outdoor_clocksync_parameter.mat","-v7.3")
% load("ArgosCSI-96x7-2016-11-04-00-57-12_5GHz_track_left_to_right_outdoor_clocksync_parameter.mat")
% 
% lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
% lts_t = ifft(lts_f, 64);
% 
% Ant_index = 1;
% frame_index = 70;
% mobile_index = 4;
% samples_per_user = 224;
% cp_length = 16;
% nfft = 64;
% thresh = 0.9;
% peak_spacing = 64;
% interpolated_fft_size = 1024;
% 
% 
% Ant_start = 1;
% Ant_end = 96;
% frame_start = 101;
% frame_end = 1100;
% mobile_start = 1;
% mobile_end = 8;
% 
% Channel_result_f = zeros(Ant_end - Ant_start + 1, mobile_end - mobile_start + 1, frame_end - frame_start + 1, interpolated_fft_size);
% Channel_result_t = zeros(Ant_end - Ant_start + 1, mobile_end - mobile_start + 1, frame_end - frame_start + 1, interpolated_fft_size);
% 
% for Ant_index = Ant_start : Ant_end
%     l = Ant_index - Ant_start + 1;
%     for mobile_index = mobile_start : mobile_end
%         m = mobile_index - mobile_start + 1; 
%         for frame_index = frame_start : frame_end
%             n = frame_index - frame_start + 1;
% 
%             Pilot_Sample = squeeze(Pilot_Samples(:,samples_per_user*(mobile_index - 1)+1:samples_per_user*(mobile_index),Ant_index,frame_index));
%             Pilot_Sample_complex = (double(Pilot_Sample(1,:)) + 1i*double(Pilot_Sample(2,:))) * 2^(-15);
% 
%             lts_flip = flip(lts_t);
%             lts_flip_conj = conj(lts_flip);
% 
%             sign_fct = Pilot_Sample_complex ./ abs(Pilot_Sample_complex);
%             sign_fct(isnan(sign_fct)) = 0;      
% 
%             lts_corr = abs(conv(lts_flip_conj, sign_fct));
%             lts_pks = find(lts_corr(1:length(Pilot_Sample_complex)) > (thresh * max(lts_corr)));
% 
%             [x_vec, y_vec] = meshgrid(lts_pks, lts_pks);
% 
%             [second_peak_idx, y] = find((y_vec - x_vec) == peak_spacing);
% 
%             if not(isempty(y))
% 
%                 max_corr_idx = lts_pks(y(1));
% 
%                 recv_lts = Pilot_Sample_complex(max_corr_idx:max_corr_idx + nfft - 1);
% 
%                 h = fftshift(fft(recv_lts))./fftshift(lts_f);
%                 used_h = [h(7:32), h(34:59)];
%                 interpolated_h = interp1(1:length(used_h),used_h,linspace(1,52,interpolated_fft_size),'cubic');
%                 P = polyfit(1:interpolated_fft_size,unwrap(angle(interpolated_h)),1);
%                 slope = P(1);
%                 interpolated_h = interpolated_h .* exp(-1i*(0:interpolated_fft_size-1)*slope);
% 
%                 Channel_result_f(l,m,n,:) = interpolated_h;
%                 Channel_result_t(l,m,n,:) = ifft(fftshift(interpolated_h));
%             end
%         end
%     end
% end
% 
% figure
% mobile_idx = 5;
% frame_num = 1000;
% phase_diff = zeros(96,frame_num);
% for frame_idx = 1:frame_num
%     for ant_idx = 1:96
%         % phase_diff(ant_idx + 1,frame_idx) = mean(unwrap(angle(squeeze(Channel_result_f(ant_idx + 1,mobile_idx,frame_idx,:))) - angle(squeeze(Channel_result_f(ant_idx,mobile_idx,frame_idx,:)))));
%         phase_diff(ant_idx,frame_idx) = mean(unwrap(angle(squeeze(Channel_result_f(ant_idx,mobile_idx,frame_idx,:)))));
%     end
% end
% 
% figure
% 
% plot(unwrap(phase_diff(2,:) - phase_diff(1,:)))
% hold on
% plot(unwrap(phase_diff(3,:) - phase_diff(1,:)))
% hold on
% plot(unwrap(phase_diff(4,:) - phase_diff(1,:)))
% hold on
% hold off