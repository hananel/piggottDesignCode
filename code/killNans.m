function X = killNans(X,flag)
if nargin == 1
    for i=1:length(X)
        if isnan(X(i)) X(i)=0;
        end
    end
else
    if flag ==1
        for i=1:length(X)
            if isnan(X(i)) X(i)=X(i-1);
            end
        end
    end
end