function getCorrelation(data1,subjects1,data2,subjects2)

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

plot(d1,d2,'ko');
[r,p]= corrcoef(d1,d2);
text(0,0,['r=' num2str(r(1,2)) ',p=' num2str(p(1,2)) ',N=' num2str(length(d1))]);
end