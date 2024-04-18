function ClearAllFigures
%% Find all open figures
figureHandles=findall(0,'Type','figure');
for ii=1:1:length(figureHandles)
clf(figureHandles(ii));
end
end