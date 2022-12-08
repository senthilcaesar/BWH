function V_x = prctile2(Varray,p)
% from https://www.mathworks.com/matlabcentral/answers/184895-prctile-function-vs-excel-percentile
p=p/100;

% if ~exist('C')
%     C=1;
% end
%
% if C==1
% The linear interpolation between closest ranks method in the case of `C=1`
% [Percentile - Wikipedia](https://en.wikipedia.org/wiki/Percentile)
V_x = nan(length(p),size(Varray,2));

for j=1:size(Varray,2)
    
    V = Varray(:,j);
    for i=1:length(p)
        
        if ~isvector(p(i)) || numel(p(i)) == 0 || any(p(i) < 0 | p(i) > 1) || ~isreal(p(i))
            error('Make sure the Second digit within the [0,1] interval');
        end
        
        V = sort(V,'ascend');
        V(isnan(V))=[];
        N = length(V);

        x = p(i)*(N-1)+1; % position x
        if floor(x) < N
            V_x(i,j) = V(floor(x)) + mod(x,1)*(V(floor(x)+1) - V(floor(x))); % value
        else
            V_x(i,j) = V(N); % position N
        end
    end
    
end
%
% elseif C==0.5
%     % The linear interpolation between closest ranks method in the case of `C=0.5`
%     % [Percentile - Wikipedia](https://en.wikipedia.org/wiki/Percentile)
%     if ~isvector(p) || numel(p) == 0 || any(p < 0 | p > 1) || ~isreal(p)
%         error('Make sure the Second digit within the [0,1] interval');
%     end
%     V = sort(V,'ascend');
%     N = length(V);
%     p_1 = 1/(2*N); % position 1
%     p_N = 1 - 1/(2*N); % position N
%     if p_1<=p && p<=p_N
%         x = N*p + 0.5; % position x
%         V_x = V(floor(x)) + mod(x,1)*(V(floor(x)+1) - V(floor(x))); % value
%     else if 0<=p && p<p_1
%             V_x = V(1); % value 1
%         else if p_N<p && p<=1
%                 V_x = V(N); % value N
%             end
%         end
%     end
% end
