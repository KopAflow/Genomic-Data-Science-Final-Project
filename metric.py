from tqdm import tqdm
def edit_dist(s1, s2):
    # Minimal edit dist algorithm
    m, n = len(s1), len(s2)
    prev = list(range(n + 1))

    # Rest of the rows
    for i in (range(1, m + 1)):
        curr = [i]  # j = 0
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                curr.append(prev[j - 1])
            else:
                curr.append(1 + min(curr[-1], prev[j], prev[j - 1]))
        prev = curr
    return prev[n]



def completeness(assembled_seq):
    cnt = 0
    for g in assembled_seq:
        if g == 'N':
            cnt += 1
    return cnt/len(assembled_seq)