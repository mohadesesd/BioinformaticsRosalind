from Bio import SeqIO

def semiglobal_alignment(v, w, sigma):
    S = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    backtrack = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            scores = [S[i-1][j] - sigma, S[i][j-1] - sigma, S[i-1][j-1] + [-1, 1][v[i-1] == w[j-1]]]
            S[i][j] = max(scores)
            backtrack[i][j] = scores.index(S[i][j])

    last_row_index = max(range(len(w)+1), key=lambda x: S[len(v)][x])
    last_column_index = max(range(len(v)+1), key=lambda x: S[x][len(w)])

    if S[len(v)][last_row_index] >= S[last_column_index][len(w)]:
        i = len(v)
        j = last_row_index
    else:
        i = last_column_index
        j = len(w)
    max_score = S[i][j]

    def insert_indel( word, i):
        return word[:i] + '-' + word[i:]

    v_aligned, w_aligned = v, w

    for _ in range(len(v) - i):
        w_aligned += '-'
    for _ in range(len(w) - j):
        v_aligned += '-'

    while i*j != 0:
        if backtrack[i][j] == 0:
            i -= 1
            w_aligned = insert_indel(w_aligned, j)
        elif backtrack[i][j] == 1:
            j -= 1
            v_aligned = insert_indel(v_aligned, i)
        else:
            i -= 1
            j -= 1

    for _ in range(i):
        w_aligned = insert_indel(w_aligned, 0)

    for _ in range(j):
        v_aligned = insert_indel(v_aligned, 0)

    with open('SMGB1.txt', 'w') as output_data:
        output_data.write(str(max_score))
        output_data.write('\n')
        output_data.write(str(v_aligned.seq))
        output_data.write('\n')
        output_data.write(str(w_aligned.seq))


input_handle = "/home/msd/Downloads/rosalind_smgb.txt"
data = list(SeqIO.parse(input_handle, 'fasta'))
word1, word2 = data[0], data[1]
semiglobal_alignment(word1, word2, 1)
