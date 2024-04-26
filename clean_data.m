%% Plot line profiles for each .mat file to check for poorly aligned cells
data = profiles{1};
figure();
cells = size(data, 1);
for i = 1:cells
p = cell2mat(data(i, :));
p = reshape(p, size(p, 1), [], 3);
p = mean(p, 3);
subplot(2, ceil(cells/2), i);
plot(p);
title(sprintf("Cell %d", i));
end
sgtitle("11-4")

%% Clean sheets 2-end after manually cleaning sheet 1
dT = 13;
T_keys = readtable('/Users/allisonlam/Downloads/osmotic_shock_rep1_temp.xls','Sheet',1);
T_keys(idx, :) = [];
writetable(T_keys, '/Users/allisonlam/Downloads/osmotic_shock_rep1_cleaned.xls','Sheet',1);
T_keys = T_keys(:, [1 2]);
for i = 2:dT
T2 = readtable('/Users/allisonlam/Downloads/osmotic_shock_rep1_temp.xls','Sheet',i);
T2 = innerjoin(T_keys, T2);
writetable(T2, '/Users/allisonlam/Downloads/osmotic_shock_rep1_cleaned.xls','Sheet',i);
end

%% Plot radius over time to check for cells w/ weird dynamics
dT = 13;
T = readtable('/Users/allisonlam/Downloads/osmotic_shock_rep1_cleaned.xls','Sheet',1);
T = T(:, [1 2 3]);
for i = 2:dT
T2 = readtable('/Users/allisonlam/Downloads/osmotic_shock_rep1_cleaned.xls','Sheet',i);
T2 = T2(:, [1 2 3]);
T2 = renamevars(T2, "Radius" , sprintf("Radius%d", i));
[T,ileft,iright] = innerjoin(T, T2);
end
Tplot = T(:, 3:15);

testtable = table2array(Tplot);
[~, columns] = size(testtable);
for column = 2:columns
testtable(:, column - 1) = testtable(:, column) - testtable(:, column - 1);
end
testtable = array2table(testtable);
testtable(:, 13) = [];
testidx = any(testtable{:,:}>30,2);
idx = find(testidx.');
Tplot(idx, :) = [];

plot(Tplot{:,:}.')
xlim([1, 13])
title("Radius by timepoint per cell")
xlabel("Timepoint")
ylabel("Radius (pixels)")

%% Small vs Large
Tsorted = sortrows(Tplot, 1);
ind = Tsorted(:, 1);
ind = table2array(ind);
smallmax = min(find(ind > 120));
largemin = min(find(ind > 150));
Tsmall = table2array(Tsorted(1:smallmax, 2:end));
Tlarge = table2array(Tsorted(largemin:end, 2:end));

% % Uncomment me for percent change
% Tlarge = Tlarge ./ Tlarge(:, 1);
% Tsmall = Tsmall ./ Tsmall(:, 1);

% Calculate mean, std
Tlarge_mean = mean(Tlarge, 1);
Tlarge_std = std(Tlarge, 1);
Tsmall_mean = mean(Tsmall, 1);
Tsmall_std = std(Tsmall, 1);
maxstd_s = Tsmall_mean + Tsmall_std;
minstd_s = Tsmall_mean - Tsmall_std;
maxstd_l = Tlarge_mean + Tlarge_std;
minstd_l = Tlarge_mean - Tlarge_std;

% Plotting
figure();
x = linspace(1, 12, 12);
plot(x, Tsmall_mean)
hold on;
plot(x, Tlarge_mean)
fill([x fliplr(x)],[minstd_s fliplr(maxstd_s)], [0.2 0.5 0.9], "FaceAlpha", 0.2, "EdgeColor", "none")
fill([x fliplr(x)],[minstd_l fliplr(maxstd_l)], [0.6350 0.0780 0.1840], "FaceAlpha", 0.2, "EdgeColor", "none")
legend("Small", "Large", '', '')
xlim([1 12])
xlabel("Timepoint"); ylabel("Radius (px)"); title("Radius by time (percent change)");

