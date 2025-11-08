% 加载各个方法下的结果（第3列是 CQI3）
data_predistortion = load('results_withPredistortion_CQI7-9.mat');
cqi3_predistortion = data_predistortion.throughputResults(:, 3);  % CQI3

data_without = load('results_noPredistortion_CQI7-9.mat');
cqi3_without = data_without.throughputResults(:, 3);  % CQI3

data_water = load('results_withwater_CQI7-9.mat');
cqi3_water = data_water.throughputResults(:, 3);  % CQI3

% SNR 横轴
snrValues = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
cqi3_predistortion = cqi3_predistortion(2:end);
cqi3_without = cqi3_without(2:end);
cqi3_water = cqi3_water(2:end);
% 定义颜色
color_predistortion = [0.4940 0.1840 0.5560]; % 紫色
color_without       = [0 0.4470 0.7410];      % 蓝色
color_water         = [0.9290 0.6940 0.1250]; % 黄色

% 绘图
figure;
hold on;

plot(snrValues, cqi3_predistortion, '-o', 'Color', color_predistortion, ...
    'DisplayName', 'Power-Equalization', 'LineWidth', 5, 'MarkerSize', 8);

plot(snrValues, cqi3_without, '-s', 'Color', color_without, ...
    'DisplayName', 'Default setup', 'LineWidth', 5, 'MarkerSize', 8);

plot(snrValues, cqi3_water, '-^', 'Color', color_water, ...
    'DisplayName', 'Waterfilling', 'LineWidth', 5, 'MarkerSize', 8);

% 图形美化
xlabel('SNR (dB)', 'FontSize', 35);
ylabel('Throughput (Mbps)', 'FontSize', 35);
set(gca, 'FontSize', 24);
legend('Location', 'best', 'FontSize', 30);
set(gcf, 'Color', 'w');
grid on;
hold off;
