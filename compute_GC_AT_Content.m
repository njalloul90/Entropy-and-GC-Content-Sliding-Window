function [GC,AT] = compute_GC_AT_Content(sequence,window_length,shift)
%GC content is usually calculated as a percentage value and sometimes 
% called G+C ratio or GC-ratio. GC-content percentage is calculated as:
    %   Count(G + C)/Count(A + T + G + C) * 100%
% AT computation is similar

% input: 
    % sequence: character array of DNA sequence (include only A,T,C,G and N)
    % window_length: size of window to compute the GC content for
    % shift: overlap difference between two consecutive windows
% outout:
    % GC: G+C ratio
    % AT: A+T ratio

GC = [];
AT = [];
sequence = upper(sequence);                                                % uppercase
if isempty(sequence)
    error('sequence is empty');
end

l = length(sequence);
sprintf("analyzing sequence of length %d...",l)
bases = unique(sequence);
sprintf("summary:\n")
for k = 1 : length(unique(sequence))
    freq_bases(k) = length(find(ismember(sequence,bases(k))==1))/l*100;
    sprintf("%s: %0.2f %%\n",bases(k),freq_bases(k))
end


for i1 = 1 : shift : l
    if i1 + window_length - 1 <= l
        wind = sequence(i1:i1 + window_length - 1);                        % temp window
        alphab = unique(wind);                                             % bases present
        count_alphab = zeros(length(alphab),1);                            % count of each present base
        for a = 1 : length(count_alphab)
            count_alphab(a) = length(find(ismember(wind,alphab(a))==1));
            gc = count_alphab(find(alphab=='C')) + count_alphab(find(alphab=='G'));
            at = count_alphab(find(alphab=='A')) + count_alphab(find(alphab=='T'));
        end
        GC = [GC (gc/(gc+at))*100];
        AT = [AT (at/(gc+at))*100];
    end
end

end

