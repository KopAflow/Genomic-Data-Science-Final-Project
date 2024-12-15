from tqdm import tqdm

def find_max_overlap(s1, s2):
    # Finds the maximum overlap between two strings.
    # Returns the length of the overlap and the merged string.
    max_overlap = 0
    merged_string = ""

    # Check overlap where s1 suffix matches s2 prefix
    for i in range(1, len(s1)):
        if s2.startswith(s1[-i:]):
            if i > max_overlap:
                max_overlap = i
                merged_string = s1 + s2[i:]
    
    # Check overlap where s2 suffix matches s1 prefix
    for i in range(1, len(s2)):
        if s1.startswith(s2[-i:]):
            if i > max_overlap:
                max_overlap = i
                merged_string = s2 + s1[i:]
    
    return max_overlap, merged_string

def merge_scaffolds(scaffolds, min_overlap=10):
    
    # Iteratively merges scaffolds based on maximum overlap.
    max_iter = len(scaffolds) - 1
    for _ in tqdm(range(max_iter)):
        max_overlap = 0
        best_pair = (None, None)
        best_merged = None

        # Find the pair of scaffolds with the maximum overlap
        for i in range(len(scaffolds)):
            for j in range(i + 1, len(scaffolds)):
                overlap, merged = find_max_overlap(scaffolds[i], scaffolds[j])
                if overlap > max_overlap:
                    max_overlap = overlap
                    best_pair = (i, j)
                    best_merged = merged
        
        i, j = best_pair
        scaffolds.pop(j)
        scaffolds.append(best_merged) 
    return scaffolds