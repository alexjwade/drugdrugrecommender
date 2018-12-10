### DRUG-DRUG SIMILARITY ###
### We use the methods from Park et al. 2016 to calculate drug-drug similarity
### via random walk with restart over protein-protein interaction network
### followed by comparison of drug target scores from the walk

"""PSEUDOCODE"""
#convert the SNAP graph of the PPI into an adjacency matrix
#   initialize to 0, iterate through edges and add to matrix

#get the two drugs we want to compare
#find their targets in the protein-protein interaction network using the
#drug-target network (limit to the small one at first)
#then perform the random walk with restart for each one until they converge
#set up ROC

#test algorithm on small debug set first before moving to large network
#test known pairs of drugs using http://biosoft.kaist.ac.kr/targetrw/ as gold standard

#play with the restart

#once we're SURE that it's perfect...run it on all the drugs and write to a pickle file

"""IMPORTS"""

#import snap #not needed after pickling
#import numpy as np
#import matplotlib.pyplot as plt
#from scipy import spatial
import collections
import pickle
import copy
import math
import time
import operator
import timeit

"""CONSTANTS"""

r = 0.7
epsilon = 1.0e-5

def dd():
    return collections.defaultdict(float)

def dotProduct(d1, d2):
    if len(d1) < len(d2):
        return dotProduct(d2, d1)
    else:
        return sum(d1.get(f, 0) * v for f, v in d2.items())

def scale(d1, n): #1d
    for f, v in d1.items():
        d1[f] = d1.get(f, 0) * n

def add(d1, d2):
    for f, v in d2.items():
        d1[f] = d1.get(f, 0) + v

def normalize1d(d1):
    dict_sum = len(d1)
    if dict_sum >= epsilon:
        scale(d1, 1.0 / dict_sum)

def normalize(d1): #2d
    for key, value in d1.items():
        #get sum
        dict_sum = len(value)
        #print dict_sum
        #print value
        if dict_sum >= epsilon:
            scale(value, 1.0 / dict_sum)
        #print value
        #divide by sum

def find_difference(d1, d2):
    sum = 0.0
    for f, v in d2.items():
        sum += abs(d1.get(f, 0) - v)
    return sum

"""LOAD NETWORKS"""


"""drug to target"""



dtt_pickle = 'drug_to_target.pkl'

"""
dtt_output = open(dtt_pickle, 'wb')

drug_to_target_filepath = 'bio-decagon-targets-all.csv'
#drug_to_target_filepath = 'test_targets.csv'
drug_to_target_filepath_2 = 'processed_' + drug_to_target_filepath

drug_ids = []
target_ids = []
drug_to_target = collections.defaultdict(dd)

#process the CID names into ints
with open(drug_to_target_filepath, 'r') as f:
    i = 0
    for line in f:
        if i == 0:
            i = 1
        else:
            line_bits = line.split(',')

            #line_bits[0] is drug and line_bits[1] is target
            drug_id = line_bits[0]
            if drug_id not in drug_ids:
                drug_ids.append(drug_id)

            target_id = int(line_bits[1])
            if target_id not in target_ids:
                target_ids.append(target_id)

            line_bits[0] = line_bits[0].replace('CID', '99')
            drug_to_target[line_bits[0]][line_bits[1].rstrip()] = 1.0
#print drug_to_target

pickle.dump(drug_to_target, dtt_output)
#ppi_output.close()
dtt_output.close()


"""

dtt_read = open(dtt_pickle, 'rb')
drug_to_target = pickle.load(dtt_read)

"""double check"""
#print len(drug_to_target.keys())
#sum = 0
#for key, value in drug_to_target.items():
#    sum += len(value)
#print sum

"""protein-protein network"""
#make adjacency matrix


ppi_filepath = 'bio-decagon-ppi.csv'
#ppi_filepath = 'test-ppi.csv'

ppi_pickle = 'ppi.pkl'
proteinlist_pickle = 'proteinlist.pkl'

"""

ppi_output = open(ppi_pickle, 'wb')
proteinlist_output = open(proteinlist_pickle, 'wb')

#protein_position = []

context = snap.TTableContext()
schema = snap.Schema()
schema.Add(snap.TStrTAttrPr("Gene 1", snap.atInt))
schema.Add(snap.TStrTAttrPr("Gene 2", snap.atInt))
table = snap.TTable.LoadSS(schema, ppi_filepath, context, ",", snap.TBool(True))
proteins1 = snap.TIntV()
proteins2 = snap.TIntV()
table.ReadIntCol("Gene 1", proteins1)
table.ReadIntCol("Gene 2", proteins2)
proteins1.Union(proteins2)
proteinlist = list(set(list(proteins1)))
proteinlist.sort()
#print len(proteinlist)





ppi = collections.defaultdict(dd)

with open(ppi_filepath, 'r') as f:
    i = 0
    for line in f:
        if i == 0:
            i = 1
        else:
            i += 1
            if i % 10000 == 0:
                pass
                #print i
            line_bits = line.split(',')
            protein1 = line_bits[0]
            protein2 = line_bits[1].rstrip()
            ppi[protein1][protein2] = 1.0
            ppi[protein2][protein1] = 1.0



pickle.dump(ppi, ppi_output)
pickle.dump(proteinlist, proteinlist_output)
ppi_output.close()
proteinlist_output.close()

"""



ppi_read = open(ppi_pickle, 'rb')
proteinlist_read = open(proteinlist_pickle, 'rb')
ppi = pickle.load(ppi_read)
proteinlist = pickle.load(proteinlist_read)

#print ppi
#print len(proteinlist)



normalize(ppi)

def set_p0(drug):
    p0 = drug_to_target[drug]
    normalize1d(p0)
    return p0

def rwr(drug):
    p0 = set_p0(drug)
    pt = copy.copy(p0)
    pt1 = copy.copy(p0)
    teleport_vector = copy.copy(p0)
    scale(teleport_vector, r)
    difference = 1.0 #sentinel
    while(difference > epsilon):
        pt = copy.copy(pt1)
        for k, v in ppi.items():
            pt1[k] = dotProduct(ppi[k], pt)
        scale(pt1, 1 - r)
        add(pt1, teleport_vector)
        difference = find_difference(pt, pt1)
    #print pt1
    return pt1

def ddiscore(d1, d2):
    if len(d1) < len(d2):
        return ddiscore (d2, d1)
    else:
        return sum(math.sqrt(d1.get(f, 0) * v) for f, v in d2.items())

#sim_a = rwr('99000001983')
#sim_gly = rwr('99000003488')
#sim_er = rwr('99000012560')
#sim_ibu =        rwr('99000003672')
#sim_ketoprofen = rwr('99000003825')
#sim_codeine =    rwr('99005284371')
#sim_naproxen =   rwr('99000156391')
#sim_dextro =     rwr('99000010100')
#sim_oxy =        rwr('99005284603')

drwr_pickle = 'drug_to_rwr.pkl'
# drwr_output = open(drwr_pickle, 'wb')
# drug_to_rwr = collections.defaultdict(dd) #data struct

# start_time = time.time()
# i = 0
# for drug in drug_to_target.keys():
#     drug_to_rwr[drug] = rwr(drug)
#     i += 1
#     if i % 10 == 0:
#         curr_time = time.time()
#     	print "processed %d of %d drugs in %f seconds" % (i, len(drug_to_target.keys()),curr_time - start_time)

# pickle.dump(drug_to_rwr, drwr_output)
# drwr_output.close()

drwr_read = open(drwr_pickle, 'rb')
drug_to_rwr = pickle.load(drwr_read)
# print drug_to_rwr


# pairs to try:
# drug_pair_dict = [('99000148192', '99000054454'), ('99000002585', '99002724385'), ('99000002585', '99000439260'), 
#     ('99002724385','99000054454'), ('99002724385','99000005578'), ('99000000681', '99000439260'), ('99000000681', '99000004782'),
#     ('99000000838', '99000004171'), ('99000000838', '99000004782'), ('99000003345', '99000004192'), ('99000003386', '99000003559'),
#     ('99000084029', '99000008223'), ('99000054454' ,'99000084029'), ('99000006167', '99000002520'), ('99000005487', '99000002764')]

# drug_pair_dict = [
# ('99000004058','99000005538'), 
# ('99000002764','99000005538'),
# ('99000004819','99000005073'),
# ('99000003878','99000004819'),
# ('99000002764','99000004058'),
# ('99000004585','99000004819'),
# ('99000004927','99000005538'),
# ('99000002771','99000004601'),
# ('99000002802','99000004819'),
# ('99000000861','99000004819')]

drug_pair_dict = [('99000004819', '99000005073'),  
('99000003878', '99000004819'), 
('99000000861', '99000004819'), 
('99000004585', '99000004819'), 
('99000000861', '99000005073'), 
('99000004158', '99000004819'), 
('99000004058', '99000005538'),
('99000004440', '99000005523'),
('99000000861', '99000004158'),
('99000004543', '99000004601')
]

drug_to_top5 = {}

def getTop5(drugA, drug_to_top5, drug_to_rwr):
    if drugA in drug_to_top5:
        return drug_to_top5[drugA]
    else:
        simA = drug_to_rwr[drugA]
        drugA_sim_scores = {}
        for drug, sim in drug_to_rwr.items():
            if drug != drugA:
                drugA_sim_scores[drug] = ddiscore(simA, sim)
        drugA_sim_scores = sorted(drugA_sim_scores.items(), key=operator.itemgetter(1), reverse = True)
        i = 0
        drugA_top5 = []
        for drug, score in drugA_sim_scores:
            if i < 5:
                drugA_top5.append((drug, score))
            i += 1

        print("drugA: {}, top5: {}".format(drugA, drugA_top5))
        drug_to_top5[drugA] = drugA_top5
        return drugA_top5

for pair in drug_pair_dict:
    drug1_to_top5 = getTop5(pair[0], drug_to_top5, drug_to_rwr)
    drug2_to_top5 = getTop5(pair[1], drug_to_top5, drug_to_rwr)

    filepath = pair[0] + "_" + pair[1] + ".txt"
    with open(filepath, 'w') as f:
        output = "%s:%s:%s:%s" % (pair[0], drug1_to_top5, pair[1], drug2_to_top5)
        f.write(output)

# # start_time = time.time()
# drugA = '99100003009'
# # print(simA)
# drugB = '99100003008'
# simB = drug_to_rwr[drugB]
# # print(simB)
# # print('time elapsed: {}'.format(time.time() - start_time))

    # drugA_sim_scores = {}
    # drugB_sim_scores = {}
    # # print('calculating sim scores')
    # # start_time = time.time()
    # for drug, sim in drug_to_rwr.items():
    #     if drug != drugA:
    #         # # time for 1 ddi score:
    #         # ddi_start = time.time()
    #         drugA_sim_scores[drug] = ddiscore(simA, sim)
    #         # print('ddi time: {}'.format(time.time() - ddi_start))
    

    #     if drug != drugB:
    #         drugB_sim_scores[drug] = ddiscore(simB, sim)
    # # print(drugA_sim_scores)
    # # print(drugB_sim_scores)

    # # print('time elapsed: {}'.format(time.time() - start_time)) 

    # # sort dicts by value (descending order)
    # drugA_sim_scores = sorted(drugA_sim_scores.items(), key=operator.itemgetter(1), reverse = True)
    # drugB_sim_scores = sorted(drugB_sim_scores.items(), key=operator.itemgetter(1), reverse = True)

    # # print('time elapsed: {}'.format(time.time() - start_time))

    # i = 0
    # drugA_top5 = []
    # for drug, score in drugA_sim_scores.items():
    #     if i < 5:
    #         drugA_top5.append((drug, score))
    #     i += 1

    # i = 0
    # drugB_top5 = []
    # for drug, score in drugB_sim_scores.items():
    #     if i < 5:
    #         drugB_top5.append((drug, score))
    #     i += 1

    # # print('time elapsed: {}'.format(time.time() - start_time))
    # print("drugA: {}, top5: {}".format(pair[0], drugA_top5))
    # print("drugB: {}, top5: {}".format(pair[1], drugB_top5))



#similarity_ibu = ddiscore(sim_a, sim_ibu)
"""similarity_nap = ddiscore(sim_a, sim_naproxen)
similarity_cod = ddiscore(sim_a, sim_codeine)
similarity_keto = ddiscore(sim_a, sim_ketoprofen)
similarity_dex = ddiscore(sim_a, sim_dextro)
similarity_oxy = ddiscore(sim_a, sim_oxy)"""
#similarity_er = ddiscore(sim_a, sim_er)
#print type(drug_to_target)
#print similarity_ibu
"""print similarity_nap
print similarity_cod
print similarity_keto
print similarity_dex
print similarity_oxy"""
#print similarity_er
#print ppi['1']
#print test_vector
#print dotProduct(test_vector, ppi['1'])

#
# row_sums = a.sum(axis=1)
# new_matrix = a / row_sums[:, numpy.newaxis]
