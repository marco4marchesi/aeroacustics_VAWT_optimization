%{
Set graphical values for better looking plots
%}

%% set some values
LW = 2.5;
FS = 16;
MS = 6;

%% Interpreter:

% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%% Font
% This script changes all FontNames to Palatino Linotype. 
list_factory = fieldnames(get(groot,'factory'));
index_fontName = find(contains(list_factory,'FontName'));
for i = 1:length(index_fontName)
    default_name = strrep(list_factory{index_fontName(i)},'factory','default');
    set(groot, default_name,'Palatino Linotype');
end

% This script changes all FontSizes to FS value defined at the beginning. 
list_factory = fieldnames(get(groot,'factory'));
index_fontSize = find(contains(list_factory,'FontSize'));
index_fontSize_remove = [find(contains(list_factory,'FontSizeMode')); find(contains(list_factory,'FontSizeMultiplier'))];
index_fontSize_save = [];
for j= 1: length(index_fontSize)
    if sum(index_fontSize_remove == index_fontSize(j)) == 0
        index_fontSize_save = [index_fontSize_save, index_fontSize(j)];
    end
end
for i = 1:length(index_fontSize_save)
    default_name = strrep(list_factory{index_fontSize_save(i)},'factory','default');
    set(groot, default_name,FS);
end

%% Figure properties:

set(0, 'defaultFigureColormap',turbo(256));
set(0, 'defaultFigureColor', [1; 1; 1]);

%% Surfaces:
% transparency
set(0, 'defaultSurfaceEdgeAlpha', 0.3);

%% Lines:

% This script changes all Linewidth to value LW set at the beginning. 
list_factory = fieldnames(get(groot,'factory'));
index_lineWidth = find(contains(list_factory,'LineWidth'));
for i = 1:length(index_lineWidth)
    default_name = strrep(list_factory{index_lineWidth(i)},'factory','default');
    set(groot, default_name,LW);
end

%% Markers:

% This script changes all Linewidth to value LW set at the beginning. 
list_factory = fieldnames(get(groot,'factory'));
index_markerSize = find(contains(list_factory,'MarkerSize'));
for i = 1:length(index_markerSize)
    default_name = strrep(list_factory{index_markerSize(i)},'factory','default');
    set(groot, default_name,MS);
end
%% legend:
set(0, 'defaultLegendLocation','best');


%% grid 
set(0, 'defaultAxesXMinorGrid', 'off');
set(0, 'defaultAxesYMinorGrid', 'off');


%% color
set(0, 'defaultAxesColor', 'none');
%% fontSize


index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

clearvars i FS j LW MS index_fontName index_fontSize index_fontSize_remove index_fontSize_save index_interpreter index_lineWidth index_markerSize list_factory default_name


