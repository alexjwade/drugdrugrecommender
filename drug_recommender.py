'''File: drug_recommender.py
Authors: Lea Jabbour, Marcos Torres, and Alex Wade 
---------------------------
Accepts a query of two drug CIDs and searches the CHCHSE-decagon
network for adverse interactions between them. In the event of an interaction,
it searches for similar drugs that can be substituted to minimize interactions.
'''

import snap
import argparse
import collections
import os
import ast
import operator

########################## old functions ##################################
'''
reformats a given CID string as an nID
'''
def format_cid(CID_str):
    return int("99" + CID_str[3:])

'''
loads ChChSe csv as TTable
'''
def load_chchse(path):

    #load table
    context = snap.TTableContext()
    schema = snap.Schema()
    schema.Add(snap.TStrTAttrPr("STITCH 1", snap.atStr))
    schema.Add(snap.TStrTAttrPr("STITCH 2", snap.atStr))
    schema.Add(snap.TStrTAttrPr("Polypharmacy Side Effect", snap.atStr))
    schema.Add(snap.TStrTAttrPr("Side Effect Name", snap.atStr))
    table = snap.TTable.LoadSS(schema, path, context, ",", snap.TBool(True))
 
    #reformat CIDs and seIDs as strings
    raw_cid1s = snap.TStrV()
    cid1s = snap.TIntV()
    table.ReadStrCol("STITCH 1", raw_cid1s)
    for raw_cid in raw_cid1s:
        cid = format_cid(raw_cid)
        cid1s.Add(cid)
    table.StoreIntCol("cid1", cid1s)

    raw_cid2s = snap.TStrV()
    cid2s = snap.TIntV()
    table.ReadStrCol("STITCH 2", raw_cid2s)
    for raw_cid in raw_cid2s:
        cid = format_cid(raw_cid)
        cid2s.Add(cid)
    table.StoreIntCol("cid2", cid2s)

    #save table as binary
    #cache_path = "../cache/ChChSe-decagon_table.tsv"
    #table.Save(snap.TFOut(cache_path))

    #TEST: checks the number of side effect types
    seVec = snap.TStrV()
    table.ReadStrCol("Side Effect Name", seVec)
    print len(set(list(seVec)))

    return table
############################################################################

'''
parses input
'''
# def arg_parsing():
#     #define input
#     parser = argparse.ArgumentParser(description='Return similar drug pairs that miniize interactions.')
#     parser.add_argument('CID1', type=str)
#     parser.add_argument('CID2', type=str)
#     args = parser.parse_args()
#     #return formatted CIDs
#     cid1 = args.CID1
#     cid2 = args.CID2
#     return cid1, cid2

def parse_file(file):
    file = "./output2_top_twosides_weighted/similarity_scores2/" + file
    f = open(file, "r")
    line = f.read().split(':')
    cid1 = line[0]
    cid1_top5 = ast.literal_eval(line[1]) # converts string of a list into a list
    cid2 = line[2]
    cid2_top5 = ast.literal_eval(line[3])
    f.close()
    return cid1, cid2, cid1_top5, cid2_top5

def count_interactions_twosides(path):
    counts = collections.defaultdict(int)
    with open(path, 'r') as f:
        for line in f:
            tokenized = line.split(',')
            tokenized[0] = tokenized[0].replace('CID', '99') # added this for consistency
            tokenized[1] = tokenized[1].replace('CID', '99')
            counts[(tokenized[0], tokenized[1])] = int(tokenized[2])

    return counts 


def count_interactions_decagon(path):
    counts = collections.defaultdict(int)
    with open(path, 'r') as f:
        i = 0
        for line in f:
            if i == 0:
                i = 1
            else:
                tokenized = line.split(',')
                tokenized[0] = tokenized[0].replace('CID', '99') # added this for consistency
                tokenized[1] = tokenized[1].replace('CID', '99')
                counts[(tokenized[0], tokenized[1])] += 1

    return counts 


# TODO: FIX
def get_recommendation(cid1, cid2, cid1_top5, cid2_top5, SE_counts):
    # #get current interaction score
    # curr_pair = (cid1[0], cid2[0])
    # curr_score = SE_counts[curr_pair]
    # print curr_pair, curr_score
    
    # #get best similar interaction score
    # cid1_top5 = [cid1] + cid1_top5
    # cid2_top5 = [cid2] + cid2_top5
    # for new_cid1 in cid1_top5:
    #     for new_cid2 in cid2_top5:
    #         new_pair = (new_cid1[0], new_cid2[0])
    #         new_score = SE_counts[new_pair]
    #         if new_score < curr_score:
    #             print new_score
    #             print curr_score
    #             print new_pair
    #             print curr_pair
    #             curr_pair = new_pair
    #             curr_score = new_score
    # print curr_pair, curr_score
    # return curr_pair, curr_score

    filepath = "twosides_recommendation_" + str(cid1[0]) + "_" + str(cid2[0]) + ".txt"
    # filepath = "decagon_recommendation_" + str(cid1[0]) + "_" + str(cid2[0]) + ".txt"
    with open(filepath, 'w') as f:
        curr_pair = (cid1[0], cid2[0])
        curr_score = SE_counts[curr_pair]
        
        #get scores for all combos of the similar drugs 
        cid1_top5 = [cid1] + cid1_top5
        cid2_top5 = [cid2] + cid2_top5
        for new_cid1 in cid1_top5:
            for new_cid2 in cid2_top5:
                new_pair = (new_cid1[0], new_cid2[0])
                new_score = SE_counts[new_pair]
                curr_pair = new_pair
                curr_score = new_score
                f.write(str(curr_pair) + " " + str(curr_score) + '\n')
                print curr_pair, curr_score


    return curr_pair, curr_score

def main():
    # path = "bio-decagon-combo.csv"
    # SE_counts = count_interactions_decagon(path)
    # path = "twosides_processed_weighted.csv"
    # SE_counts = count_interactions_twosides(path)

    path = "twosides_processed.csv"
    SE_counts = count_interactions_twosides(path)

    # sorted_counts = sorted(SE_counts.items(), key=operator.itemgetter(1), reverse = True)
    # print(sorted_counts[0:20])
    # i = 0
    # top10_counts = []
    # for pair, count in sorted_counts:
    #     if i < 17:
    #         top10_counts.append(pair)
    #     i += 1

    # print top10_counts
    

    # parse files to get drug ids and top5 similar list
    # files = [filename for filename in os.listdir('./old_output') if filename.startswith("99")]
    files = [filename for filename in os.listdir('./output2_top_twosides_weighted/similarity_scores2') if filename.startswith("99")]
    for file in files:
        d1, d2, d1_top5, d2_top5 = parse_file(file)
        d1 = (d1, 0)
        d2 = (d2, 0)
        recommendation = get_recommendation(d1, d2, d1_top5, d2_top5, SE_counts)
        print("drug1: {}, drug2: {}, recommendation: {}".format(d1, d2, recommendation))
        


if __name__ == "__main__":
    main()
