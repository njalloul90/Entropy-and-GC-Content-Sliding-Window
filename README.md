# Sequence Complexity Measures

## slindingW_blockSize_Entropy(sequence,window_length,block_size,shift):

computation of entropy as a measure of complexity for a DNA sequence.

function to compute shannon entropy of a dna sequence;

input:

    % sequence: character array
    
    % window_length: length of sliding window
    
    % shift: overlap difference between two consecutive windows
    
    % block_size: number of adjacent nucleotides to compute entropy for
    
    
output:


    % H: vector containing shannon entropy calculation for each window
    
    % based on the selected block size
    
## compute_GC_AT_Content(sequence,window_length,shift):

GC and AT content calculation; calculated as : Count(G + C)/Count(A + T + G + C) * 100%

input:

    % sequence: character array of DNA sequence (include only A,T,C,G and N)
    
    % window_length: size of window to compute the GC content for
    
    % shift: overlap difference between two consecutive windows
    
output:

    % GC: G+C ratio
    
    % AT: A+T ratio
    

    
    
