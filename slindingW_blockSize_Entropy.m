function [H] = slindingW_blockSize_Entropy(sequence,window_length,block_size,shift)

% a function to compute shannon entropy of a dna sequence
% input: 
    % sequence: chraacter array
    % window_length: length of sliding window
    % shift: overlap difference between two consecutive windows
    % block_size: number of adjacent nucleotides to compute entropy for
% output:
    % H: vector containing shannon entropy calculation for each window
    % based on the selected block size
H = [];
    
if isempty(sequence)
    error('sequence is empty');
end

l = length(sequence);
sprintf("analyzing sequence of length %d...",l)


% % initialize parameters
% i1 = 1;                                                                    % start index
% i2 = i1 + window_length - 1;                                               % end index


for i1 = 1 : shift : l
    if i1 + window_length - 1 <= l
    wind = sequence(i1:i1 + window_length - 1);                            % temp window
    alphab = unique(wind);                                                 % bases present
    freq_alphab = zeros(length(alphab),1);                                 % frequency of each present base
    prob_alphab = zeros(length(alphab),1);                                 % probability of each present base
    for a = 1 : length(freq_alphab)
        freq_alphab(a) = length(find(ismember(wind,alphab(a))==1));
        prob_alphab(a) = freq_alphab(a)/length(wind);
    end
    
    if block_size == 1
        h_x = -sum(prob_alphab.*log2(prob_alphab));                        % entropy for block size of 1
    else
        % count all possible words of a given block size within the window
        if rem(length(wind),block_size) ~= 0                               % can't divide window into equal blocks
            error('window should be divided into equal-sized blocks; adjust either window length or block size accordingly...');
        else
            words = {};
            n_words = length(wind)/block_size;
            for kk = 1 : length(wind) - (block_size - 1)
                block = wind(kk:(kk+length(wind)/n_words -1));
                words = [words;block];                                     % cell array containing all words of given block length where word consists of adjacent bases
            end
            
            freq_word = zeros(length(unique(words)),1);                    % frequency of each unique block
            prob_word = zeros(length(unique(words)),1);                    % probability of each unique block
            for b = 1 : length(freq_word)
                freq_word(b) = count(string(wind),string(words{b}));
                prob_word(b) = freq_word(b)/n_words;
            end
            h_x = -sum(prob_word.*log2(prob_word));
        end
    end
    else
        return;
    end
H = [H ; h_x];
end

end

