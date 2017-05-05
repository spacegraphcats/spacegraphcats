import functools

def memoize(func):
    cache = func.cache = {}
    @functools.wraps(func)
    def memoized_func(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            val = func(*args, **kwargs)
            cache[key] = val
            return val
        return cache[key]
    return memoized_func
