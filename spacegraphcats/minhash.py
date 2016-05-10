import string
import random
from parser import parse_minhash

class MinHash:
    """
        Class of MinHashes
    """
    def __init__(self, hash_size):
        self.list = []
        self.size = hash_size

    @staticmethod
    def from_list(hashes):
        res = MinHash(len(hashes))
        res.list = list(hashes)
        return res

    def normalize(self):
        self.list = sorted(self.list)[:self.size]
        assert len(self.list) <= self.size

    def add(self, val):
        self.list.append(val)
        self.normalize()

    def merge(self, h, size):
        res = MinHash(size)
        res.list = self.list + h.list
        res.normalize()
        assert len(res.list) <= size
        return res

    def intersect(self, other):
        return [x for x in other if x in self]

    def __contains__(self, val):
        return val in self.list

    def __iter__(self):
        return iter(self.list)

    def __len__(self):
        return len(self.list)

    def __str__(self):
        return ','.join(str(i) for i in self.list)

    @staticmethod
    def string_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

if __name__ == "__main__":
    A = MinHash.from_list([1,2,3,5])
    B = MinHash.from_list([2,3,4,5])
    print A.intersect(B)
