import sys
import random
import math

def gen_bfs_queries(inputFile, outputFile, q_num):
    v_num = 0
    with open(inputFile, 'r') as fin:
        for line in fin:
            u,v = map(int, line.strip().split())
            v_num = max(v_num, max(u, v))

    with open(outputFile, 'w') as fout:
        for i in xrange(q_num):
            u = random.randint(0, v_num)
            fout.write(str(u) + '\n')

def gen_shortest_distance_queries(inputFile, outputFile, q_num):
    v_num = 0
    with open(inputFile, 'r') as fin:
        for line in fin:
            u,v = map(int, line.strip().split())
            v_num = max(v_num, max(u, v))

    with open(outputFile, 'w') as fout:
        for i in xrange(q_num):
            u = random.randint(0, v_num)
            v = random.randint(0, v_num)
            fout.write(str(u) + ' ' + str(v) + '\n')

def trans_newid_queries(queryFile, idFile, outputFile):
    org2newid = {}
    with open(idFile, 'r') as fin:
        for line in fin:
            org_id, new_id = map(int, line.rstrip().split())
            org2newid[org_id] = new_id

    with open(queryFile, 'r') as fin, open(outputFile, 'w') as fout:
        for line in fin:
            org_ids = map(int, line.strip().split())
            for i in xrange(len(org_ids)):
                new_id = org2newid[org_ids[i]]
                fout.write(str(new_id))
                if (i + 1 == len(org_ids)): fout.write('\n')
                else: fout.write(' ')

if __name__ == "__main__":
    if (sys.argv[1] == '-b'):
        gen_bfs_queries(sys.argv[2], sys.argv[3], 1000)
    elif (sys.argv[1] == '-s'):
        gen_shortest_distance_queries(sys.argv[2], sys.argv[3], 50000)
    elif (sys.argv[1] == '-t'):
        trans_newid_queries(sys.argv[2], sys.argv[3], sys.argv[4])