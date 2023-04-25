from multiprocessing.connection import Listener
import sys
recombination_file_name = "filtering/data/combinedCatOnlyBestWithPVals.txt"
recombination_file_sock = "filtering/data/combinedCatOnlyBestWithPVals.sock"

with Listener(recombination_file_sock, family='AF_UNIX',backlog=int(sys.argv[1])) as listener:
    F = open(recombination_file_name, 'r')
    lines = F.readlines()
    F.close()

    labels = lines[0].split('\t')
    for i in range(len(labels)):
        labels[i] = labels[i].strip()

    del lines[0]
    parsimony_change = []

    for i in range(len(lines)):
        lines[i] = lines[i].split('\t')
        parsimony_change.append(int(lines[i][10]) - int(lines[i][11]))

    lines_sorted_by_parsimony_change = [x for _, x in sorted(zip(parsimony_change, lines) , reverse=True)]

    while (True):
        with listener.accept() as conn:
            conn.send(labels)
            idx=conn.recv()
            try:
                conn.send(lines_sorted_by_parsimony_change[idx])
            except:
                conn.send(([],[]))
