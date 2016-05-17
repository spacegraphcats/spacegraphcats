import string
import random
from parser import parse_minhash

class MinHash:
    """
        Class of MinHashes
    """
    def __init__(self, hash_size):
        self.values = set()
        self.size = hash_size

    @staticmethod
    def from_list(hashes):
        res = MinHash(len(hashes))
        res.values = set(hashes)
        return res

    def normalize(self):
        self.values = set(sorted(self.values)[:self.size])
        assert len(self.values) <= self.size

    def add(self, val):
        self.values.add(val)
        self.normalize()

    def merge(self, h, size):
        res = MinHash(size)
        self.values.update(h.values)
        res.normalize()
        assert len(res.values) <= size
        return res

    def intersect(self, other):
        return [x for x in other if x in self]

    def __contains__(self, val):
        return val in self.values

    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.values)

    def __str__(self):
        return ','.join(str(i) for i in self.values)

    @staticmethod
    def string_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

if __name__ == "__main__":
    A = MinHash.from_list([1,2,3,5])
    B = MinHash.from_list([2,3,4,5])
    print(A.intersect(B))
