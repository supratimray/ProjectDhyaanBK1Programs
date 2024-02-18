function [r,p,N]=getCorrelation(data1,subjects1,data2,subjects2,hPlot,hColor,showSignificanceFlag,linearFitFlag)

if ~exist('hPlot','var');                   hPlot=gca;                  end
if ~exist('hColor','var');                  hColor=[0 0 0];             end
if ~exist('showSignificanceFlag','var');    showSignificanceFlag=0;     end
if ~exist('linearFitFlag','var');           linearFitFlag=0;            end

d1 = []; d2 = [];
for i=1:length(subjects1)
    if ~isnan(data1(i))
        posS2 = find(strcmp(subjects1{i},subjects2));
        if ~isempty(posS2)
            if ~isnan(data2(posS2))
                d1 = cat(2,d1,data1(i));
                d2 = cat(2,d2,data2(posS2));
            end
        end
    end
end

% plot the data
scatter(hPlot,d1,d2,[],hColor,'filled','o');
hold(hPlot,'on');
[r,p]= corrcoef(d1,d2);
if showSignificanceFlag
    text(min(d1),max(d2),['r=' num2str(round(r(1,2),2)) ', p=' num2str(round(p(1,2),2)) ', N=' num2str(length(d1))],'Color',hColor,'parent',hPlot);
end

% linear fit of the line
if linearFitFlag
    linearCoefficients = polyfit(d1,d2, 1);
    xFit = linspace(min(d1), max(d1), 100);
    yFit = polyval(linearCoefficients, xFit);
    plot(hPlot,xFit, yFit,'LineStyle',":",'Color',hColor,'MarkerSize', 15,'LineWidth', 2);
end

% set the output
r=round(r(1,2),2); % correlation coefficient
p=round(p(1,2),2); % p-value
N=length(d1);      % no of subjects

end