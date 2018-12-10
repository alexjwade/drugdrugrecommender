#side_effects = set()
import collections
import matplotlib.pyplot as plt

counts = collections.defaultdict(int)

# with open("bio-decagon-combo.csv", 'r') as f:
#     i = 0
#     for line in f:
#         if i == 0:
#             i = 1
#         else:
#             tokenized = line.split(',')
#             counts[(tokenized[0], tokenized[1])] += 1



with open("twosides_processed.csv", 'r') as f:
    for line in f:
        tokenized = line.split(',')
        counts[(tokenized[0], tokenized[1])] = int(tokenized[2])


freq = collections.defaultdict(int)
for tup, count in counts.items():
    freq[count] += 1

print freq
freq = sorted(freq.items()) # make sure it's sorted by keys (count)
x, y = zip(*freq)
plt.plot(x, y)
plt.xlabel('Number of edges between two drug nodes (drug-drug score),  confidence = 5')
plt.ylabel('Frequency')
plt.title('Distribution of number of edges between drug nodes in TWOSIDES')
plt.savefig('twosides_count_dist.png')
plt.show()

"""with open("side-effects-processed.csv", 'w') as g:
    for side_tuple in side_effects:
        g.write(','.join(side_tuple))"""


#print(side_effects)
