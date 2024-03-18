%% Run me ONCE per experiment/folder!!
% Change these paths to wherever you have bfmatlab and images
addpath /Users/allisonlam/Downloads/bfmatlab 
addpath(genpath('/Volumes/My Passport/Light microscopy/normal_vs_polyploidy')) 

% Parameters
warning('off','all');
dZ = 15; % Number of Z-slices per timepoint
dT = 13; % Number of timepoints
output_file = '/Users/allisonlam/Downloads/osmotic_shock_rep2_20_29.xls'; % Spreadsheet name

% Setup
[listOfFolderNames, listOfFileNames, ~] = find_files("_R3D_D3D.dv");
timepoints = cell(1, dT);
tables = cell(1,dT);
profiles = cell(1, 3);
all_radii = cell(1, 3);
all_centers = cell(1, 3);
empty_table = table('Size',[0 8],'VariableTypes',{'string','int32','double', 'double', 'double', 'double', 'double', 'double'});
slice_table = table('Size',[0 3],'VariableTypes',{'string','string','string'});
slice_table.Properties.VariableNames = ["Filename", "Early-timepoint Z-slices", "Later-timepoint Z-slices"];
empty_table.Properties.VariableNames = ["File Name", "Cell ID", "Radius", "Inner Baseline", "Outer Baseline", "Max Peak Height", "FWHM", "AUC"];
[tables{:}] = deal(empty_table);

%% Run me for every image
% Parameters
for imageNumber = 23:29 % starting image:end image in folder
%imageNumber = 17; % Image in selected folder
sensitivity = 0.97; % Circle finding sensitivity; decrease if too many circles, increase if not enough
frame = floor(dZ/2); % Starting frame to find cells
checkcircles = 1; % Setting whether to check found circles; 0 to display found circles, 1 to skip
checkpair = 1; % Setting whether to check paired cells; 0 to check, 1 to skip
usemanual = 0;

% Read images
filepath = strcat(listOfFolderNames{1}, '/', listOfFileNames{imageNumber});
imageFiles = dir(filepath);
dv_file = bfopen(fullfile(imageFiles.folder, imageFiles.name));
images = dv_file{1};

for i = 1:dT 
    timepoints{i} = images((dZ * (i - 1) + 1):dZ * i, :);
end

for i = 1:size(timepoints, 2)
    img = [];
    for j = 1:size(timepoints{i}, 1)
        img = cat(3, img, timepoints{i}{j});
    end
    timepoints{i} = img;
end
%%
% Main analysis
for time = 1:dT
img = timepoints{time};
slice_image = img(:, :, frame);
slice_image = rescale(slice_image,0,1);

if usemanual == 0
[centers, radii] = imfindcircles(imfill(slice_image, "holes"),[90 350], ...
    'Sensitivity', sensitivity,'Method','TwoStage', "ObjectPolarity","bright");

% Remove circles, radii cut off by the edge of the image
remove_idx = find(sum((centers - radii - 100) < 0, 2) + sum((centers + radii + 100) > size(slice_image, 2), 2));
centers(remove_idx, :) = [];
radii(remove_idx) = [];
end

% Match cells
if (time ~= 1 && usemanual == 0)
    dup = [];
    differences = pdist2(centers, centers_1);
    [M, idx] = min(differences, [], 2);
    u=unique(idx);
    [n,bin]=histc(idx,u);
    for v=find(n>1).'
        dup=find(bin==v);
    end

    if size(dup, 1) > 0
        max_dup = M(dup(size(dup, 1)));
        dup_id = dup(size(dup, 1));
        for i = 1:size(dup, 1)-1
            next_dup = M(dup(size(dup, 1)-i));
            if next_dup > max_dup
                idx(dup(size(dup, 1)-i)) = [];
                centers(dup(size(dup, 1)-i), :) = [];
                radii(dup(size(dup, 1)-i), :) = [];
            else
                idx(dup_id) = [];
                centers(dup_id, :) = [];
                radii(dup_id, :) = [];
                min_dup = dup;
                dup_id = dup(size(dup, 1)-i);
            end
        end
    end
end

if (time == 1)
    start_slice = slice_image;
    centers_1 = centers;
    radii_1 = radii;
    checkcircles = 0;
elseif (size(centers,1) ~= last_size)
    checkcircles = 0;
    checkpair = 0;
elseif (time == 4)
    checkpair = 0;
end

while checkcircles == 0
figure();
subplot(1, 2, 1)
imshow(start_slice);
viscircles(centers_1, radii_1,'EdgeColor','b');
title(sprintf('Identified cells in frame 1'));
subplot(1, 2, 2)
imshow(slice_image);
viscircles(centers, radii,'EdgeColor','b');
title(sprintf('Identified cells in frame %d', time));
checkcircles = input("Did I find the right number of circles? Enter 1 if yes, 0 if no.");
if (checkcircles == 0)
    select = CROIEditor(slice_image);
    pause;
    roi = select.roi;
    stats = regionprops("table",roi,"Centroid", ...
    "MajorAxisLength","MinorAxisLength");
    centers = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2;
    usemanual = input("Should I use these circles for later images? Enter 1 if yes, 0 if no.");
    % Match cells
    if (time ~= 1)
    dup = [];
    differences = pdist2(centers, centers_1);
    [M, idx] = min(differences, [], 2);
    u=unique(idx);
    [n,bin]=histc(idx,u);
    for v=find(n>1).'
        dup=find(bin==v);
    end

    if size(dup, 1) > 0
        max_dup = M(dup(size(dup, 1)));
        dup_id = dup(size(dup, 1));
        for i = 1:size(dup, 1)-1
            next_dup = M(dup(size(dup, 1)-i));
            if next_dup > max_dup
                idx(dup(size(dup, 1)-i)) = [];
                centers(dup(size(dup, 1)-i), :) = [];
                radii(dup(size(dup, 1)-i), :) = [];
            else
                idx(dup_id) = [];
                centers(dup_id, :) = [];
                radii(dup_id, :) = [];
                min_dup = dup;
                dup_id = dup(size(dup, 1)-i);
            end
        end
    end
    end
end
end

if (time == 1)
    idx = [linspace(1, size(centers, 1), size(centers, 1))];
    idx = idx.';
    imshow(slice_image);
    viscircles(centers, radii,'EdgeColor','b');
    
    % Cut subimages based on circle, radii coordinates
    subimages = cell(size(centers, 1), 1);
    for j = 1:size(subimages, 1)
        buffer = radii(j) + 100;
        subimages{j} = img(centers(j, 2) - buffer:centers(j, 2) + buffer, centers(j, 1) - buffer:centers(j, 1) + buffer, :);
    end
    
    figure();
    for k = 1:size(subimages, 1)
        subplot(size(subimages, 1), 1, k);
        montage(rescale(subimages{k},0,1), "Size", [1 size(subimages{k}, 3)]);
        hold on;
    end
    starting_slices = input("Enter a list of Z-slices corresponding to each early-timepoint cell:");
    slices = starting_slices;
elseif (time == 4)
    imshow(slice_image);
    viscircles(centers, radii,'EdgeColor','b');
    
    % Cut subimages based on circle, radii coordinates
    subimages = cell(size(centers, 1), 1);
    for j = 1:size(subimages, 1)
        buffer = radii(j) + 100;
        subimages{j} = img(centers(j, 2) - buffer:centers(j, 2) + buffer, centers(j, 1) - buffer:centers(j, 1) + buffer, :);
    end
    
    figure();
    for k = 1:size(subimages, 1)
        subplot(size(subimages, 1), 1, k);
        montage(rescale(subimages{k},0,1), "Size", [1 size(subimages{k}, 3)]);
        hold on;
    end
    later_slices = input("Enter a list of Z-slices corresponding to each later-timepoint cell:");
    slices = later_slices;
end
close all;

% Check that cells were paired correctly
if (time > 1)
while checkpair == 0
subplot(1, 2, 1)
imshow(start_slice);
hold on
for k = 1:size(centers_1, 1)
    text(centers_1(k, 1), centers_1(k, 2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color','white');
end
title("Labelled cells in frame 1")
subplot(1, 2, 2)
imshow(slice_image)
hold on
for k = 1:size(centers, 1)
    text(centers(k, 1), centers(k, 2), sprintf('%d', idx(k)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color','white');
end
title(sprintf('Labelled cells in frame %d', time));
checkpair = input("Are cells paired correctly? Enter 1 if yes, 0 if no.");
if (checkpair ~= 1) % Manually change cell assignment
    idx = input("Enter corrected cell IDs as a list:"); 
end
end
end

% Take correct image slice and black out all other cells
masked_images = cell(size(centers, 1), 1);
for k = 1:size(centers, 1)
   data = rescale(img(:, :, slices(k)-1:slices(k)+1), 0, 1);
   mask_centers = centers;
   mask_radii = radii;
   mask_centers(k, :) = [];
   mask_radii(k, :) = [];
   mask_radii = mask_radii + 70;
   mask = createCirclesMask(data(:, :, 1),mask_centers,mask_radii);
   data = bsxfun(@times, data, cast(~mask, 'like', data));
   buffer = radii(k) + 100;
   masked_images{k} = data;
end

% Line profiles
aligned = cell(size(centers, 1), 3, 3);
export_data = zeros(size(centers, 1), 6, 3);
for i = 1:size(masked_images, 1)
xcenter = centers(i, 1);
ycenter = centers(i, 2);
[X, Y] = make_endpoints(centers(i, 1), centers(i, 2), radii(i)*2);

for slice = 1:3
figure('visible','off');
mask_img = masked_images{i};
imshow(mask_img(:, :, slice))
hold on
    for j = 1:360
        plot(X(j,:),Y(j, :),'Color','r','LineWidth', 0.0002)
        hold on
    end
hold off;

% Get pixel intensity along each line
radius = radii(i);
c = double.empty(0, 360);
for k = 1:360
    ci = improfile(mask_img(:, :, slice), X(k,:),Y(k, :),radii(i)*2);
    c = [c ci];
end

 % make all true zeros into NaNs
    indices2 = find(abs(c)==0.0);
    c(indices2) = NaN;
    
    profcontainsnan = zeros(360, 1);
    for ii = 1:360
       if sum(isnan(c(1:radius, ii)), "all") ~= 0
            c(:, ii) = NaN(size(c(:, ii)));
            profcontainsnan(ii) = 1;
        end
    end
    
    
    cleaned_profiles = c;
    zero_sec = NaN(round(radius-radius*0.30), 360);
    cleaned_profiles(1:round(radius-radius*0.30), 1:360) = zero_sec;
     
    [peaks, pos] = max(cleaned_profiles);
    loc = pos';
    xrelativecir = zeros(360, 1);
    yrelativecir = zeros(360, 1);
    
    for jj = 1:360
        xrelativecir(jj, 1) = cos(jj*2*pi/360)*loc(jj) +xcenter;
        yrelativecir(jj, 1) = sin(jj*2*pi/360)*loc(jj) +ycenter;
    end

    [prof_i, prof_j] = size(cleaned_profiles);
    aligned_profiles = NaN(prof_i, prof_j);
    amount_moved = NaN(360, 1);
    
 
    for kk = 1:prof_j

        if profcontainsnan(kk) == 0
            col = cleaned_profiles(:, kk);
            amount_moved(kk, 1) = -(radius - int32(loc(kk)));
            colprime = circshift(col, radius - int32(loc(kk)));
            if loc(kk) - radius < 0
                colprime = vertcat(NaN(abs(int32(loc(kk))-radius), 1), colprime(abs(int32(loc(kk))-radius-1):size(colprime, 1)));
            elseif loc(kk) - radius > 0
                colprime((length(colprime) - abs(int32(loc(kk))-radius) + 1):length(colprime)) = NaN(abs(int32(loc(kk))-radius), 1);
            end
            aligned_profiles(:, kk) = colprime;
        end
        
    end

meanprofile = mean(aligned_profiles, 2, 'omitnan');
aligned{i, slice} = aligned_profiles; % Export profiles
plot(meanprofile);

% Find inner, outer baselines
innerbaselinepos = 0.9;
outerbaselinepos = 1.20;
innerBaseline = meanprofile(int32(radius*innerbaselinepos)); %inner baseline
outerBaseline = meanprofile(int32(radius*outerbaselinepos)); %outer baseline

% Find AUC
ind = find(meanprofile(round(radius*0.8):end)>=innerBaseline);	% Find indices where profile>innerbaseline
leftindex = min(ind) + round(radius*0.8);
rightindex = min(max(ind) + round(radius*0.8), size(meanprofile, 1));
intwithbase = trapz(meanprofile(leftindex:rightindex, 1));
base = (rightindex - leftindex)*innerBaseline;
AUC = intwithbase - base;

% Export metrics
export_data(i, 1, slice) = radius + mean(amount_moved, "omitnan"); % true radius
export_data(i, 2, slice) = innerBaseline; %inner baseline
export_data(i, 3, slice) = outerBaseline; %outer baseline
export_data(i, 4, slice) = max(meanprofile); % Maximum peak height
export_data(i, 5, slice) = findFWHM(meanprofile, innerBaseline, radius); % Full width half maximum
export_data(i, 6, slice) = AUC; % AUC
end
end
close all;
export_data = mean(export_data, 3);
last_size = size(centers, 1);


if (time == 1)
    centers_1 = centers;
    radii_1 = radii;
    meanAligned_1 = aligned;
    start_export = export_data;
    frame1_export = array2table([linspace(1, size(centers_1, 1), size(centers_1, 1)).' start_export]);
    frame1_export = addvars(frame1_export,repmat(imageFiles.name, size(centers_1, 1), 1),'Before',"Var1");
    frame1_export.Properties.VariableNames = ["File Name", "Cell ID", "Radius", "Inner Baseline", "Outer Baseline", "Max Peak Height", "FWHM", "AUC"];
else
    framen_export = array2table([idx export_data]);
    framen_export = addvars(framen_export,repmat(imageFiles.name, size(idx, 1), 1),'Before',"Var1");
    framen_export.Properties.VariableNames = ["File Name", "Cell ID", "Radius", "Inner Baseline", "Outer Baseline", "Max Peak Height", "FWHM", "AUC"];
end

% Add to big spreadsheet
if (time == 1)
    profiles{1} = aligned(:, :, 1);
    all_radii{1} = radii;
    all_centers{1} = centers;
    tables{1} = [tables{1}; frame1_export];
    writetable(tables{1}, output_file,'Sheet',1);
elseif (time == 4)
    profiles{2} = aligned(:, :, 2);
    all_radii{2} = radii;
    all_centers{2} = centers;
    tables{time} = [tables{time}; framen_export];
    writetable(tables{time}, output_file,'Sheet',time);
elseif (time == dT)
    profiles{3} = aligned(:, :, 3);
    all_radii{3} = radii;
    all_centers{3} = centers;
    tables{time} = [tables{time}; framen_export];
    writetable(tables{time}, output_file,'Sheet',time);
    slice_data = table(string(imageFiles.name), string(num2str(starting_slices)), string(num2str(later_slices)));
    slice_data.Properties.VariableNames = ["Filename", "Early-timepoint Z-slices", "Later-timepoint Z-slices"];
    slice_table = [slice_table; slice_data];
    writetable(slice_table, output_file, 'Sheet', (dT + 1));
    outfile = string(strcat(imageFiles.name, "_profiles.mat"));
    save(outfile, "profiles", "all_radii", "all_centers");
else
    tables{time} = [tables{time}; framen_export];
    writetable(tables{time}, output_file,'Sheet',time);
end
%end

% %% Plot
% for i = 1:size(centers, 1)
%     subplot(2, size(centers, 1), i);
%     plot(meanAligned_1{i});
%     xlim([0 350])
%     subplot(2, size(centers, 1), i+size(centers, 1));
%     plot(aligned{i});
%     xlim([0 350])
end
end

%%
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
sgtitle("17-5")

%% Functions
function [X, Y] = make_endpoints(xcenter, ycenter, radius)
    
    X1 = xcenter * ones(360,1);
    Y1 = ycenter * ones(360,1);
    
    temp = transpose(0:(2*pi/360):(2*pi - 2*pi/360));
    X2 = radius*cos(temp) + X1;
    Y2 = radius*sin(temp) + Y1;

    X = [X1 X2];
    Y = [Y1 Y2];
end

function [FWHM, mover2, leftind, rightind] = findFWHM(fx, baseline, radius)
    %findFWHM: Finds the full width half max (FWHM) of a function
    [m, n] = max(fx);		%	Find maximum value and index
    mover2 = (m-baseline)/2 + baseline;
    ind = find(fx(round(radius*0.8):end)>=((m-baseline)/2 + baseline));	%	Find indices where I>=max(I)/2
    leftind = min(ind) + round(radius*0.8);			%	Leftmost index
    rightind = max(ind) + round(radius*0.8);		%	Rightmost index
    %	Get FWHM
    FWHM = abs(rightind-leftind);
    if isempty(FWHM)
        FWHM = 0;
    end
end