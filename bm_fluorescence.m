T = readtable("/Users/allisonlam/Downloads/20211209_l1210 control rep1.txt");
% T = T(string(T.ellipse_gate) == 'pass', :);
%samples = unique(T.sample_name);
%samples = samples(4:13, :);


%for i = 1:size(samples, 1)

    figure();
    set(gcf,'Position',[100 100 500 1000])
    %name = string(samples{i});
    name = "20211209_l1210_control_rep1";
    filepath = strcat("/Users/allisonlam/Downloads/", name);
    %data = T(string(T.sample_name) == name, :);
    data = T(T.pmt3_mV ./ T.buoyant_mass_pg < 5*mean(T.pmt3_mV ./ T.buoyant_mass_pg), :);
    data = data(data.buoyant_mass_pg > (mean(data.buoyant_mass_pg) - 1.5 * std(data.buoyant_mass_pg)) & data.buoyant_mass_pg < mean(data.buoyant_mass_pg) + 2 * std(data.buoyant_mass_pg), :);
    subplot(3, 1, 1);
    hold on;
    scatter(data.buoyant_mass_pg, data.pmt3_mV ./ data.buoyant_mass_pg);
    xlabel("Buoyant Mass");
    ylabel("Amine Intensity (normalised)")
    title(sprintf("%s Amine vs Buoyant Mass", name), 'Interpreter', 'none');
    f1 = subplot(3, 1, 2);
    hold on;
    scatter(data.buoyant_mass_pg, data.pmt3_mV);
    xlabel("Buoyant Mass");
    ylabel("Amine Intensity")
    title(sprintf("%s Amine vs Buoyant Mass", name), 'Interpreter', 'none');
    f3 = subplot(3, 1, 3);
    f2 = scatter(data.buoyant_mass_pg, data.pmt2_mV);
    hold on;
    brush on;
    set(gca,'yscale','log')
    ylim([0.01, 300]);
    xlabel("Buoyant Mass");
    ylabel("FUCCI Intensity")
    title(sprintf("%s FUCCI vs Buoyant Mass", name), 'Interpreter', 'none');
    pause;
    xd = get(f2, 'XData');
    yd = get(f2, 'YData');
    b = logical(get(f2, 'BrushData'));
    brushed_x = xd(b);
    brushed_y = yd(b);
    brushed_locs = find(get(f2, 'BrushData'));
    other_x = xd;
    other_x(brushed_locs) = [];
    other_y = yd;
    other_y(brushed_locs) = [];
    other_x = other_x(other_y > 5 * mean(brushed_y));
    other_x = other_x(other_x>min(brushed_x) & other_x<max(brushed_x));
    
    other_x = other_x.';
    brushed_x = brushed_x.';
    differences = pdist2(brushed_x, other_x);
    [M, idx] = min(differences, [], 2);
    u=unique(idx);
    other_x = other_x(u);
    brush off;

    idxc = ismember(data.buoyant_mass_pg, brushed_x);
    cytokinetic = data(idxc, :);
    idxg = ismember(data.buoyant_mass_pg, other_x);
    g2 = data(idxg, :);
    scatter(f1, cytokinetic.buoyant_mass_pg, cytokinetic.pmt3_mV);
    scatter(f3, cytokinetic.buoyant_mass_pg, cytokinetic.pmt2_mV);
    scatter(f1, g2.buoyant_mass_pg, g2.pmt3_mV);
    scatter(f3, g2.buoyant_mass_pg, g2.pmt2_mV);
    hold off;
    print(filepath,'-dpdf','-fillpage')
    writetable(cytokinetic, strcat(filepath, ".xls"),'Sheet',1);
    writetable(g2, strcat(filepath, ".xls"),'Sheet',2);
%end
close all;