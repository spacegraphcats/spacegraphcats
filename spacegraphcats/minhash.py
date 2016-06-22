import string
import random
from khmer import MinHash as kMinHash

from .graph_parser import parse_minhash

class MinHash:
    """
        Class of MinHashes
    """
    def __init__(self, hash_size):
        self.mh = kMinHash(hash_size, 31)

    @staticmethod
    def from_list(hashes):
        res = MinHash(len(hashes))
        for h in hashes:
            res.mh.add_hash(h)
        return res

    def add(self, val):
        self.mh.add_hash(val)

    def merge(self, h, size):
        self.mh.merge(h.mh)

    def intersect(self, other):
        return self.mh.count_common(other.mh)

    def __contains__(self, val):
        return val in self.mh.get_mins()

    def __iter__(self):
        return iter(self.mh.get_mins())

    def __len__(self):
        return len(self.mh)

    def __str__(self):
        return ','.join(str(i) for i in self.get_mins())

    @staticmethod
    def string_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

if __name__ == "__main__":
    A = MinHash.from_list([1,2,3,5])
    B = MinHash.from_list([2,3,4,5])
    print(A.intersect(B))
