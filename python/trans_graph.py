import sys
import time
def gen_continuous_id_graph(inputFile, outputFile, isUndirected=False):
    with open(inputFile, 'r') as fin, open(outputFile, 'w') as fout:
        cur_idx = 0
        idmap = {}
        for line in fin:
            if line.startswith('#') or line.startswith('%'):
                continue
            org_u,org_v = line.split()
            if org_u not in idmap:
                idmap[org_u] = cur_idx
                cur_idx += 1
            if org_v not in idmap:
                idmap[org_v] = cur_idx
                cur_idx += 1
            u = idmap[org_u]
            v = idmap[org_v]
            fout.write(str(u) + '\t' + str(v) + '\n')
            if isUndirected:
                fout.write(str(v) + '\t' + str(u) + '\n')
        print 'cur_idx=', cur_idx

def gen_orgorder_graph(inputFile, outputFile, isUndirected=False):
    edges = []
    org_id_map = {}
    org_id_list = []
    with open(inputFile, 'r') as fin:
        for line in fin:
            if line.startswith('#') or line.startswith('%'):
                continue
            org_u,org_v = line.split()
            u = int(org_u)
            v = int(org_v)
            if (u not in org_id_map):
                org_id_map[u] = u
                org_id_list.append(u)
            if (v not in org_id_map):
                org_id_map[v] = v
                org_id_list.append(v)
            edges.append((u, v))
    org_id_list.sort()
    for i in xrange(len(org_id_list)):
        org_id_map[org_id_list[i]] = i
    for i in xrange(len(edges)):
        u,v = org_id_map[edges[i][0]],org_id_map[edges[i][1]]
        edges[i] = (u,v)
    # edges.sort()
    with open(outputFile, 'w') as fout:
        for u,v in edges:
            fout.write(str(u) + ' ' + str(v) + '\n')
            if isUndirected:
                fout.write(str(v) + ' ' + str(u) + '\n')            

               

if __name__ == "__main__":
    # get_types(sys.argv[1], sys.argv[2])
    isUndirected = False
    if (len(sys.argv) > 3 and sys.argv[3] == '-u'):
        isUndirected = True

    time_start = time.time()
    gen_continuous_id_graph(sys.argv[1], sys.argv[2], isUndirected)
    # gen_orgorder_graph(sys.argv[1], sys.argv[2], isUndirected)

    time_end = time.time()
    time_cost = (time_end - time_start) * 100.0
    print 'time_cost = %.3fms' % time_cost 