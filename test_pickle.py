import pickle

test_pickle = 'test.pkl'
test_output = open(test_pickle, 'wb')

test = "Hello World"

pickle.dump(test, test_output)
test_output.close()
