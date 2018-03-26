import sys
import random

def gen_labels(inputFile, outputFile, l_num):
    v_num = 0
    with open(inputFile, 'r') as fin:
        for line in fin:
            u,v = map(int, line.strip().split())
            v_num = max(v_num, max(u, v))
    v_num += 1
    with open(outputFile, 'w') as fout:
        for u in xrange(v_num):
            l = random.randint(0, l_num - 1)
            fout.write(str(u) + ' ' + str(l) + '\n')

def trans_newid_labels(labelFile, idFile, outputFile):
    v_num = 0
    org2newid = {}
    with open(idFile, 'r') as fin:
        for line in fin:
            org_id, new_id = map(int, line.rstrip().split())
            org2newid[org_id] = new_id
            v_num = max(v_num, org_id)
    v_num += 1
    
    newid2label = {}
    with open(labelFile, 'r') as fin:
        for line in fin:
            org_id,label = map(int, line.strip().split())
            newid2label[org2newid[org_id]] = label

    with open(outputFile, 'w') as fout:
        for u in xrange(v_num):
            new_label = newid2label[u]
            fout.write(str(u) + ' ' + str(new_label) + '\n')

if __name__ == "__main__":
    if (sys.argv[1] == '-l'):
        gen_labels(sys.argv[2], sys.argv[3], 100)
    elif (sys.argv[1] == '-t'):
        trans_newid_labels(sys.argv[2], sys.argv[3], sys.argv[4])