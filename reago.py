import os
import sys
import time
import networkx as nx
import matplotlib.pyplot as plt
import pygraphviz as pgv



def get_fa(fn):
    r, r_pos, cm_pos = {}, {}, {}
    f = open(fn)
    for line in f:
        if line[0] == ">":
            read_id, m_st, m_ed, s_st, s_ed = line[1:].split()
            r[read_id]      = ""
            r_pos[read_id]  = [int(s_st), int(s_ed)]
            cm_pos[read_id] = [int(m_st), int(m_ed)]
        else:
            r[read_id]     += line[:-1]
            
    f.close()
    return r, r_pos, cm_pos




def write_fa(seq_container, fn, width):
    try:
        f = open(fn, "w")
        for elem in seq_container:
            
            if type(seq_container) == list:
                seq_id, seq = elem
            elif type(seq_container) == dict:
                seq_id, seq = elem, seq_container[elem]
            

            if seq_id[0] != ">":
                seq_id = ">" + seq_id
            
            f.write(seq_id + "\n")
            
            if width == 0:
                f.write(seq + "\n")
            else: 
                while seq:
                    f.write(seq[:width] + "\n")
                    seq = seq[width:]
        f.close()
        return True
    except:
        return False



def ts():
    return time.asctime() 



def n_read_in_node(node):
    read_list = node.split("|")
    return len(read_list)




def init_read_pos(read_dict):
    pos = {}
    for read_id in read_dict:
        pos[read_id] = 0
    return pos 




def combine_dup(read_list):
    rev = {}
    for seq_id, seq in read_list.items():
        if seq not in rev:
            rev[seq] = []
        rev[seq].append(seq_id)
    
    dup_removed = {}
    for seq in rev:
        new_id = "|".join(rev[seq]) 
        dup_removed[new_id] = seq 
            
    return dup_removed




def create_graph_using_rj(read_dict, graph_fn):
    G = nx.DiGraph()

    non_dup_fn = fn_base + "/rj/" + graph_fn + ".fasta"
    write_fa(read_dict, non_dup_fn, 0)
    os.system("gt readjoiner prefilter -q -des -readset " + rj_dir + graph_fn+ ".set -db " + rj_dir + graph_fn + ".fasta")
    os.system("gt readjoiner overlap -memlimit 100MB -q -l " + str(MIN_OVERLAP) + " -readset " + rj_dir + graph_fn + ".set")
    os.system("gt readjoiner spmtest -readset " + rj_dir + graph_fn + ".set.0 -test showlist > " + rj_dir + graph_fn + ".edge.list")

    read_map = {}
    cnt = 0
    f = open(fn_base + "/rj/" + graph_fn + ".set.des")
    for line in f:
        read_map[str(cnt)] = line[:-1]
        cnt += 1
    f.close()

    f = open(fn_base + "/rj/" + graph_fn + ".edge.list")
    for line in f:
        if "-" in line:
            continue

        read_1, read_2, over = line.split(" + ")
        read_id_1, read_id_2 = read_map[read_1], read_map[read_2]
        
        G.add_edge(read_id_1, read_id_2, overlap = int(over))
        
    f.close()
    return G




# do error correction, G is the original, uncombined graph
def do_ec(G, ratio):
    if len(G.nodes()) <= 1:
        return
        
    align   = {} 
    st_node = []
    
    visited = set([])
    # get starting nodes
    for node_str in G.nodes():
        if len(G.predecessors(node_str)) == 0:
            st_node.append(node_str) # every node represents a read
    
    for st in st_node:
        align[st] = {}
        align[st][st], max_st_pos = 0, 0
        
        que = [st] # BFS
        while que != []:
            cur = que.pop(0)
            if cur in visited:
                continue
            else:
                visited.add(cur) # ownership of cur

            succ = G.successors(cur)
            que += succ 
            
            for s in succ:
                overlap = G[cur][s]['overlap']
                align[st][s] = align[st][cur] + READ_LEN - 1 - overlap
                max_st_pos = max(max_st_pos, align[st][s])

        align_disp = []
        for a in align[st]:
            align_disp.append([" " * align[st][a] + read_dict[a], a])


        # correcting...
        for i in range(max_st_pos + READ_LEN):
            composition = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0}
            involved_read = []
            for aligned_read, read_id in align_disp: # ____ACGATGC..ACGATC 23431.CAJ.1
                if i < len(aligned_read) and aligned_read[i] != ' ':
                    composition[aligned_read[i]] += len(read_id.split("|")) + 1
                    involved_read.append(read_id)
            
            ttl_cnt = sum(composition[k] for k in composition) 

            dominant = 'X'
            for base in composition:
                base_cnt = composition[base]
                if float(base_cnt) / ((ttl_cnt - base_cnt) + 1) > ratio:
                    dominant = base
                    break
            
            if dominant == 'X': # when no dominant base
                continue
            
            for read_id in involved_read:
                orig_seq = list(read_dict[read_id])
                cur_base = orig_seq[i - align[st][read_id]]
                if float(composition[cur_base]) / ttl_cnt < 0.05:
                    orig_seq[i - align[st][read_id]] = dominant
                    read_dict[read_id] = "".join(orig_seq)
       



# do error correction, G is the original, uncombined graph
def do_ec_rev(G, ratio):
    align   = {} 
    st_node = []
    
    visited = set([])
    # get starting nodes
    for node_str in G.nodes():
        if len(G.successors(node_str)) == 0:
            st_node.append(node_str) # every node represents a read
    
    for st in st_node:
        align[st] = {}
        align[st][st], min_st_pos = 0, 0
        
        que = [st] # BFS
        while que != []:
            cur = que.pop(0)
            if cur in visited:
                continue
            else:
                visited.add(cur) # ownership of cur

            pred = G.predecessors(cur)
            que += pred
            
            for p in pred:
                overlap = G[p][cur]['overlap']
                #align[st][p] = align[st][cur] - (100 - overlap)
                align[st][p] = align[st][cur] - (READ_LEN - 1 - overlap)
                min_st_pos = min(min_st_pos, align[st][p])


        for r in align[st]:
            align[st][r] -= min_st_pos

        align_disp_orig = []
        for a in align[st]:
            align_disp_orig.append([" " * align[st][a] + "+" + read_dict[a][::-1], a])
       
        
        align_disp = []
        for seq, r_id in align_disp_orig:
            ws, seq_rev = seq.split("+")
            align_disp.append([ws + seq_rev[::-1], r_id])
                
        # correcting...
        for i in range(-min_st_pos + READ_LEN):
            composition = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            involved_read = []
            for aligned_read, read_id in align_disp: # ____ACGATGC..ACGATC 23431.CAJ.1
                if i < len(aligned_read) and aligned_read[i] != ' ':
                    composition[aligned_read[i]] += len(read_id.split("|")) + 1
                    involved_read.append(read_id)
            
            ttl_cnt = sum(composition[k] for k in composition) 

            dominant = 'X'
            for base in composition:
                base_cnt = composition[base]
                if float(base_cnt) / ((ttl_cnt - base_cnt) + 1) > ratio:
                    dominant = base
                    break
            
            if dominant == 'X': # when no dominant base
                continue
            
            for read_id in involved_read:
                orig_seq = list(read_dict[read_id])
                orig_seq[i - align[st][read_id]] = dominant
                read_dict[read_id] = "".join(orig_seq)




def collapse_graph(G, candidate):
    while True:
        node_to_combine = []
        if not candidate:
            all_node = G.nodes()
        else:
            all_node = candidate
        
        for node in all_node:
            if G.in_degree(node) == 1 and \
               G.out_degree(G.predecessors(node)[0]) == 1:
                 
                node_to_combine.append(node) 
                if candidate:
                    candidate.remove(node)

        if not node_to_combine:
            break
        
        for node in node_to_combine:
            pred        = G.predecessors(node)[0] # node is single-in-single-out
            pred_pred   = G.predecessors(pred)
            succ        = G.successors(node)
            
            # update graph
            new_node = pred + "|" + node
            new_over = G[pred][node]["overlap"]

            G.add_node(new_node)
            for p in pred_pred:
                o = G[p][pred]["overlap"]
                G.add_edge(p, new_node, overlap = o)
            
            for s in succ:
                o = G[node][s]["overlap"]
                G.add_edge(new_node, s, overlap = o)
            
            # update sequences
            offset = len(read_dict[pred]) - new_over
            for read_id in node.split("|"):
                read_pos_dict[read_id] += offset

            pred_seq = read_dict[pred]
            node_seq = read_dict[node]
            new_seq  = pred_seq + node_seq[new_over:]
            
            read_dict[new_node] = new_seq

            # clean up
            G.remove_node(node)
            G.remove_node(pred) 
            
            del read_dict[node]
            del read_dict[pred]

            if node in node_to_combine:
                node_to_combine.remove(node)
            if pred in node_to_combine:
                node_to_combine.remove(pred)
    return G




# merge a list of node to dst, direction: 0 (fork out), 1 (fork in)
# return merged node if merged, ow return None
def merge_node(src_list, dst, shared, G, direction):
    
    N_MIS = 3
    
    dst_seq  = read_dict[dst]
    dst_over = G[shared][dst]["overlap"]        if direction == 1 else G[dst][shared]["overlap"]
    dst_rem  = dst_seq[dst_over: ]              if direction == 1 else dst_seq[ :-dst_over][::-1]
    
    to_remove = []
    to_merge  = [] 
    for src in src_list:
        src_seq  = read_dict[src]
        src_over = G[shared][src]["overlap"]    if direction == 1 else G[src][shared]["overlap"]
        src_rem  = src_seq[src_over: ]          if direction == 1 else src_seq[ :-src_over][::-1]
        
        if n_read_in_node(src) >= 1.2 * n_read_in_node(dst):
            continue

        mis = 0
        for i in range(min(len(src_rem), len(dst_rem))):
            if src_rem[i] != dst_rem[i]:
                mis += 1
                if mis > N_MIS:
                    break
            
        if mis > N_MIS:
            if n_read_in_node(src) < TIP_SIZE:
                to_remove.append(src)
            continue
        
        #offset = dst_over - src_over if direction == 1 else (len(dst_seq) - len(src_seq) + dst_over - src_over)
        offset = dst_over - src_over if direction == 1 else ( (len(dst_seq) - dst_over) - (len(src_seq) - src_over) ) 

        for read_id in src.split("|"):
            read_pos_dict[read_id] += offset
        
        to_merge.append(src)
        

    if not to_remove + to_merge:
        return None
    
    for n in to_remove:
        G.remove_node(n)

    if to_merge:
        new_node            = dst + "|" + "|".join(to_merge)
        G                   = nx.relabel_nodes(G, {dst: new_node}, copy = False)
        read_dict[new_node] = read_dict.pop(dst)

        for n in to_merge:
            G.remove_node(n)
        
        return new_node

    else:
        return dst




def merge_bif(G):
    while True:
        merged = False
        # fork out
        collapse_candidate = set([])
        for node in G.nodes():
            
            if node not in G.nodes():
                continue

            succ = set(G.successors(node)) 
            if len(succ) < 2:
                continue

            tip_candidate = set([s for s in succ if G.out_degree(s) == 0])
            if len(tip_candidate) == 0:
                continue

            dst_candidate = succ - tip_candidate  
            if len(dst_candidate) == 0:
                dst_node = max([[n_read_in_node(t), t] for t in tip_candidate])[1]
                tip_candidate.remove(dst_node)
            else:
                dst_node = max([[n_read_in_node(d), d] for d in dst_candidate])[1]  # only one dst node
                dst_candidate.remove(dst_node)                                      # remove dst
                dst_candidate = [d for d in dst_candidate if G.out_degree(d) == 0]
                tip_candidate = tip_candidate.union(dst_candidate)
           
            new_node = merge_node(tip_candidate, dst_node, node, G, 1)

            if new_node:
                merged = True
                collapse_candidate.add(node)
 
        G = collapse_graph(G, list(collapse_candidate))

        # fork in
        collapse_candidate = set([])
        for node in G.nodes():
            
            if node not in G.nodes():
                continue
            
            pred = set(G.predecessors(node))
            if len(pred) < 2:
                continue

            tip_candidate = set([p for p in pred if G.in_degree(p) == 0])# and G.out_degree(p) == 1])
            if len(tip_candidate) == 0:
                continue

            dst_candidate = pred - tip_candidate  
            if len(dst_candidate) == 0:
                dst_node = max([[n_read_in_node(t), t] for t in tip_candidate])[1]
                tip_candidate.remove(dst_node)
            else:
                dst_node = max([[n_read_in_node(d), d] for d in dst_candidate])[1]  # only one dst node
                dst_candidate.remove(dst_node)                                      # remove dst
                dst_candidate = [d for d in dst_candidate if G.in_degree(d) == 0]   # and G.out_degree(d) == 1]  # only if its out-deg is 0, a node will be considered tip
                tip_candidate = tip_candidate.union(dst_candidate)

            new_node = merge_node(tip_candidate, dst_node, node, G, -1)
            if new_node:
                merged = True
                collapse_candidate.add(node)

        G = collapse_graph(G, list(collapse_candidate))

        if merged == False:
            break
        
    G = collapse_graph(G, [])
    return G




def remove_bubble(G):
    while True: 
        bubble_removed      = False
        all_node            = G.nodes()
        collapse_candidate  = set([])
        
        for node in all_node:

            if node not in G.nodes():
                continue

            succ = [s for s in G.successors(node) if G.in_degree(s) == 1 and G.out_degree(s) == 1]
            if len(succ) < 2:
                continue

            d = {}
            for s in succ:
                to_node = G.successors(s)[0] # s has only one successor
                if to_node not in d:
                    d[to_node] = []
                d[to_node].append(s)
            
            for to_node in [n for n in d if len(d[n]) > 1]:
                new_node = merge_node(d[to_node][1:], d[to_node][0], node, G, 1)
                if new_node:
                    bubble_removed = True
                    collapse_candidate.add(new_node)
        
        G = collapse_graph(G, list(collapse_candidate))
        if not bubble_removed:
            break

    return G




examined_edge = set([])
def remove_bad_link(G):
    bad_edge_removed = False
    c = cls.classifier()

    for sid_1, sid_2 in G.edges(): 
        seq_1, seq_2    = read_dict[sid_1], read_dict[sid_2]
        edge            = sid_1 + "+" + sid_2 
        if  edge not in examined_edge and len(seq_1) > 120 and len(seq_2) > 120 and len(sid_1.split("|")) > 15 and len(sid_2.split("|")) > 15 and\
            c.is_bad_edge(seq_1, seq_2):
            G.remove_edge(sid_1, sid_2)
            bad_edge_removed = True

    return G, bad_edge_removed




def remove_junk_node(G):
    for node in G.nodes():
        if  not G.in_degree(node) and not G.out_degree(node) and \
            (n_read_in_node(node) < 5 or read_dict[node] < READ_LEN * 1.05):
            G.remove_node(node)

    return G
        
def get_branching_aid(G_orig):
    G       = G_orig.reverse(copy = True)
    d       = {}
    st_node = []

    for node_str in G.nodes():
        d[node_str] = set([node_str])
        if G.in_degree(node_str) == 0:
            st_node.append(node_str)
    
    # BFS
    for st in st_node:
        que = [st]
        while que != []:
            cur     = que.pop(0) 
            succ    = G.successors(cur)
            for s in succ:
                d[s] = d[s].union(d[cur])
                if s not in que:
                    que.append(s)
    return d




# core function, DFS recursive
all_path = []
def get_all_path(G, cur_path):
    last_node   = cur_path[-1] 
    succs       = G.successors(last_node)
    
    # ending node, stop recursion. 
    if succs == []:
        all_path.append(cur_path)
        return
    else:
        
        if len(succs) > 1: 
            candidate = sorted([ [conf_inc(cur_path, s), s] for s in succs ])
            next_node = candidate[-1][1]
        else:
            next_node = succs[0]
         
        get_all_path(G, cur_path + [next_node])




def conf_inc(cur_path, next_node):
    d, n_pe         = {}, 0
    subsequent_node = branching_aid[next_node]

    for idx, node in enumerate(cur_path):
        for read_id in node.split("|"):
            try:
                base, sp, end   = read_id.split(".")
            except:
                base, end   = read_id.split(".")
            d[base]         = len(cur_path) - idx - 1
    

    for node in subsequent_node:
        for read_id in node.split("|"):
            try:
                base, sp, end = read_id.split(".")
            except:
                base, end = read_id.split(".")
                

            if base in d:
                n_pe += 1 * (CONFIDENCE_BASE  ** d[base])
                #n_pe += 1 * (2 ** d[base])

    return n_pe
                



def get_contig(path, G):
    
    contig = read_dict[path[0]]
    for idx in range(1, len(path)):
        prev, cur   = path[idx-1], path[idx]
        seq         = read_dict[cur]
        overlap     = G[prev][cur]["overlap"]
        contig     += seq[overlap:]

    return contig




def get_cm_pos(path, contig):
    min_cm_st = 9999
    max_cm_ed = 0
    
    for read_id in [r for r in"|".join(path).split("|") if read_pos_dict[r] >= 0 and (read_pos_dict[r] + READ_LEN <= len(contig))]:
        min_cm_st = min(min_cm_st, cm_pos[read_id][0])
        max_cm_ed = max(max_cm_ed, cm_pos[read_id][1])
        
    return min_cm_st, max_cm_ed  




def get_assemblie(G):
    global all_path
    global branching_aid
    
    branching_aid = get_branching_aid(G)
    
    full_gene           = []
    scaffold_candidate  = []
    for node in [n for n in G.nodes() if G.in_degree(n) == 0]:
        all_path = []
        get_all_path(G, [node])
        
        for path in all_path:
            contig = get_contig(path, G)
            if len(contig) >= FULL_LENGTH:
                
                st_pos = max([ r_pos[r][0] - read_pos_dict[r]           for r in path[0].split("|")     ])
                ed_pos = max([ len(read_dict_orig[r]) - r_pos[r][1]     for r in path[-1].split("|")    ])
                
                deflanked_contig = contig[st_pos : len(contig) - ed_pos]
                
                deflanked_contig = contig
                full_gene.append([path, deflanked_contig])

            else:

                m_st, m_ed = get_cm_pos(path, contig)
                if len(contig) > 120:
                    scaffold_candidate.append([path, m_st, m_ed, contig])

    return full_gene, scaffold_candidate



def conf_connect(path_1, path_2):
    d, n_pe = {}, 0
    for idx, node in enumerate(path_1):
        for read_id in node.split("|"):
            base, sp, end = read_id.split(".")
            d[base] = len(path_1) - idx - 1
    
    for node in path_2:
        for read_id in node.split("|"):
            base, sp, end = read_id.split(".")
            if base in d:
                n_pe += 1 * (CONFIDENCE_BASE  ** d[base])
                #n_pe += 1 * (2 ** d[base])
    return n_pe




def scaffold(scaffold_candidate):
    cont                = True
    full_gene           = []

    while cont:
        cont            = False
        candidate_next  = [] 
        n_candidate     = len(scaffold_candidate)
        conf            = [[0 for i in range(n_candidate)] for j in range(n_candidate)]
        
        for i, [path_1, m_st_1, m_ed_1, contig_1] in enumerate(scaffold_candidate):
            for j, [path_2, m_st_2, m_ed_2, contig_2] in enumerate(scaffold_candidate):
                
                if i == j or min(m_ed_1, m_ed_2) - max(m_st_1, m_st_2) < 10:# or (m_st_1 >= m_st_2 and m_ed_1 <= m_ed_2):
                    conf[i][j] = 0
                else:
                    conf[i][j] = max(conf_connect(path_1, path_2), conf_connect(path_2, path_1))
        

        used_candidate = [] 
        for i, row in enumerate(conf):
            if i in used_candidate:
                continue

            max_conf_val = max(row)
            max_conf_idx = [m for m, n in enumerate(row) if n == max_conf_val][0]

            max_conf_val     = max(conf[max_conf_idx])
            max_conf_idx_rev = [m for m, n  in enumerate(conf[max_conf_idx]) if n == max_conf_val][0]

            if i == max_conf_idx_rev and max_conf_val != 0:
                used_candidate += [i, max_conf_idx]
                
                candidate_1 = scaffold_candidate[i]
                candidate_2 = scaffold_candidate[max_conf_idx]
                path_1, m_st_1, m_ed_1, contig_1 = candidate_1
                path_2, m_st_2, m_ed_2, contig_2 = candidate_2
                
                contig_new          = connect_contig(contig_1, m_st_1, m_ed_1, contig_2, m_st_2, m_ed_2)    if m_st_1 < m_st_2 else connect_contig(contig_1, m_st_1, m_ed_1, contig_2, m_st_2, m_ed_2)
                path_new            = path_1 + path_2                                                       if m_st_1 < m_st_2 else path_2 + path_1
                m_st_new, m_ed_new  = min(m_st_1, m_st_2), max(m_ed_1, m_ed_2)

                if len(contig_new) < FULL_LENGTH:
                    candidate_next.append([path_new, m_st_new, m_ed_new, contig_new])
                else:
                    full_gene.append([path_new, contig_new])

                cont = True
            else:
                candidate_next.append(scaffold_candidate[i])
        
        if cont:
            scaffold_candidate = candidate_next
        else:
            break

    return full_gene, [[path, contig] for path, du, du, contig in scaffold_candidate]




def connect_contig(seg_1, m_st_1, m_ed_1, seg_2, m_st_2, m_ed_2):
    if m_st_1 >= m_st_2 and m_ed_1 <= m_ed_2 or m_st_2 >= m_st_1 and m_ed_2 <= m_ed_1:
        return seg_1 if len(seg_1) >= len(seg_2) else seg_2
    
    if m_st_1 > m_st_2:
        seg_1, seg_2 = seg_2, seg_1

    N_MIS = int(min(len(seg_1), len(seg_2)) * 0.08)
    overlap = 0
    for i in range(min(len(seg_1), len(seg_2)), 10, -1):
        suffix = seg_1[-i:]
        prefix = seg_2[:i]
        

        n_mis = 0
        for j in range(i):
            if suffix[j] != prefix[j]:
                n_mis += 1
            if n_mis > N_MIS:
                break

        if n_mis <= N_MIS:
            overlap = i
            break

    if overlap > 0:
        return seg_1 + seg_2[overlap:]
    else:
        return seg_1 + "....." + seg_2
                



gene_cnt = 1
def write_gene(data):
    global gene_cnt

    f = open(fn_base + "/full_gene.fasta", "w")
    for path, gene in data:
        f.write(">gene_" + str(gene_cnt) + "_len=" + str(len(gene)) + "\n")
        f.write(gene + "\n")
        gene_cnt += 1
    f.close()




frag_cnt = 1
def write_frag(data):
    global frag_cnt

    f = open(fn_base + "/fragment.fasta" , "w")
    for path, gene in data:
        f.write(">frag_" + str(frag_cnt) + "_len=" + str(len(gene)) + "\n")
        f.write(gene + "\n")
        frag_cnt += 1

    f.close()



def print_help_info():
    print "-----------------------------------------------------"
    print "Usage: python reago.py filename.fasta -l READ_LENGTH"
    print "Optional parameters:"
    print "-o OVERLAP, default 0.8"
    print "-e ERROR_CORRECTION_THRESHOLD, default 0.05"
    print "-t TIP_SIZE, default 30"
    print "-b PATH_FINDING_PARAMETER, default 10"
    print "Dependencies:"
    print "1. Readjoiner 1.2"
    print "2. Infernal 1.1.1"
    print "-----------------------------------------------------"
    sys.exit(1)




##### main procedure
arg_range = {"-l" : [1, float("inf"), None], \
                "-o" : [0.5, 1, 0.8], \
                "-e" : [1, float("inf"), 0.05],\
                "-t" : [1, float("inf"), 30],\
                "-b" : [2, 11, 10], \
                "-f" : [1, float("inf"), 1350]}

args = sys.argv
if len(args) < 2:
    print_help_info()

fn = args[1]
fn_base = args[2]

if os.path.exists(fn) == False:
    print "Error: File", "\'" + fn + "\'", "doesn't exist."
    print_help_info()

for idx in range(3, len(args) - 1, 2):
    arg, val = args[idx], float(args[idx+1])
    if arg not in arg_range:
        print "Error - Invalid arg:", arg
        print_help_info()

    min_val, max_val = arg_range[arg][:2]
    if val < min_val or val >= max_val:
        print "Error: Invalid value for", arg
        print "valid range:", "[" + str(min_val) + ", " + str(max_val) + ")"
        print_help_info()
    else:
        arg_range[arg][2] = val


if arg_range["-l"][2] == None:
    print "Error: read length is required"
    print_help_info()


MIN_OVERLAP = int(arg_range["-l"][2] * arg_range["-o"][2])
READ_LEN = int(arg_range["-l"][2])
TIP_SIZE = int(arg_range["-t"][2])
CONFIDENCE_BASE = int(arg_range["-b"][2])
ERROR_CORRECTION_THRESHOLD = int(arg_range["-e"][2])
FULL_LENGTH = int(arg_range["-f"][2])


#fn_base     = fn.split("/")[-1].split(".")[0]
graph_path  = fn_base + "/" + "graph.data"
plot_dir    = fn_base + "/plot/"
rj_dir      = fn_base + "/rj/"

if os.path.exists(fn) == False:
    print "Error: File", "\'" + fn + "\'", "doesn't exist."
    sys.exit(1)

if os.path.exists(fn_base) == False:
    os.mkdir(fn_base)

if os.path.exists(rj_dir) == False:
    os.mkdir(rj_dir)



print ts(), "REAGO started..."
print "Input file:", fn
print "Parameters:"
for arg in arg_range:
    print arg, arg_range[arg][2]

print ts(), "Reading input file..."
read_dict, r_pos, cm_pos = get_fa(fn)
read_dict_orig = dict(read_dict)
read_pos_dict = init_read_pos(read_dict)
read_dict = combine_dup(read_dict)


print ts(), "Initializing overlap graph..."
G = create_graph_using_rj(read_dict, "graph") 
subgraph = nx.weakly_connected_component_subgraphs(G)


print ts(), "Recovering 16S rRNAs..."
full_gene           = []
scaffold_candidate  = []
for subg in subgraph:
    do_ec(subg, 5)
    do_ec_rev(subg, 5)

    if len(subg.nodes()) == 0:
        continue

    read_dict_temp = {}
    for n in subg.nodes():
        read_dict_temp[n] = read_dict[n]
    
    if len(subg.nodes()) == 0:
        continue
    subg = create_graph_using_rj(read_dict_temp, "subg_temp")
    subg = collapse_graph(subg, [])
    
    while True:
        removed = False
        subg = merge_bif(subg)
        subg = remove_bubble(subg)
        subg = remove_junk_node(subg)

        #subg, removed = remove_bad_link(subg)
        subg = collapse_graph(subg, [])

        if not removed:
            break

    full, scaf = get_assemblie(subg)
    full_gene           += full
    scaffold_candidate  += scaf


print ts(), "Scaffolding contigs..."
scaf, remaining  = scaffold(scaffold_candidate)
full_gene       += scaf


write_gene(full_gene)
write_frag(remaining)




