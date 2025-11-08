
% 加载无 predistortion 数据
data_no_predistortion1 = load('results_noPredistortion_CQI1-3.mat');
throughput_no_predistortion1 = data_no_predistortion1.throughputResults(:, 1:3);

data_no_predistortion2 = load('results_noPredistortion_CQI4-6.mat');
throughput_no_predistortion2 = data_no_predistortion2.throughputResults(:, 1:3);

data_no_predistortion3 = load('results_noPredistortion_CQI7-9.mat');
throughput_no_predistortion3 = data_no_predistortion3.throughputResults(:, 1:3);



% 合并所有无 predistortion 数据
throughput_no_predistortion = [throughput_no_predistortion1, throughput_no_predistortion2, throughput_no_predistortion3];

% 定义 SNR 值
snrValues = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

% 定义颜色（避免使用红色）
colors = [0 0.4470 0.7410; ... % 蓝色
          0.8500 0.3250 0.0980; ... % 橙色
          0.9290 0.6940 0.1250; ... % 黄色
          0.4940 0.1840 0.5560; ... % 紫色
          0.4660 0.6740 0.1880; ... % 绿色
          0.3010 0.7450 0.9330; ... % 青色
          0.6350 0.0780 0.1840; ... % 深棕色
          0.25 0.25 0.25; ... % 灰色
          0.75 0.75 0];       % 橄榄色

% 绘制图形
figure;
hold on;

for cqi = 1:9
    % 绘制无 predistortion 的数据
    plot(snrValues, throughput_no_predistortion(:, cqi), '-o', 'Color', colors(cqi, :), ...
        'DisplayName', sprintf('CQI%d', cqi), 'LineWidth', 3, 'MarkerSize', 8);
end
% 添加标签和图例，统一字体设置
xlabel('SNR (dB)', 'FontSize', 28);
ylabel('Throughput (Mbps)', 'FontSize', 28);
set(gca, 'FontSize', 20);  % 加粗坐标轴字体
set(gcf, 'Color', 'w');
%title('Throughput vs SNR for CQI Levels with Waterfilling - Argos Dataset', 'FontSize', 22, 'FontWeight', 'bold');
% % 添加标签和图例
% xlabel('SNR (dB)');
% ylabel('Throughput (Mbps)');
% title('Throughput vs SNR for CQI Levels with Waterfilling');
legend('Location',  'best','FontSize', 18);
grid on;
hold off;