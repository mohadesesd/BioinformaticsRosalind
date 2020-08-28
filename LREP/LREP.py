def Finding_the_Longest_Multiple_Repeat(edges,s,k,NodeInformation):
    previous_node = {'node1':'ROOT'}
    for edge in edges:
        previous_node[edge[1]] = edge[0]

    node_str = {'node1':''}
    for info in NodeInformation:
        node_str[info[0][1]] = s[info[1]-1:info[1]+info[2]-1].strip('$')

    heads = set([edge[0] for edge in edges])
    tails = set([edge[1] for edge in edges])
    leaves = tails - heads

    num_nodes = max([int(node[4:]) for node in leaves])
    descendants = [0]*num_nodes
    for leaf in leaves:
        descendants[int(leaf[4:])-1] += 1
        temp_node = previous_node[leaf]
        while temp_node != 'ROOT':
            descendants[int(temp_node[4:])-1] += 1
            temp_node = previous_node[temp_node]

    candidate_nodes = []
    for i, num in enumerate(descendants[1:]):
        if num >= k:
            candidate_nodes.append('node'+str(i+2))

    candidate_strings = []
    for node in candidate_nodes:
        temp_str = ''
        temp_node = node
        while temp_node != 'ROOT':
            temp_str = node_str[temp_node] + temp_str
            temp_node = previous_node[temp_node]
        candidate_strings.append(temp_str)

    lrep = max(candidate_strings, key=len)
    return lrep


with open('/home/msd/Downloads/rosalind_lrep.txt') as input_data:
    data = input_data.readlines()

s = data[0].strip()
k = int(data[1].strip())
NodeInformation = [line.strip().split() for line in data[2:]]
NodeInformation = [[[line[0],line[1]], int(line[2]),int(line[3])]  for line in NodeInformation]
edges = [info[0] for info in NodeInformation]
lrep= Finding_the_Longest_Multiple_Repeat(edges,s,k,NodeInformation)
with open('LREP.txt', 'w') as output_file:
    output_file.write(lrep)