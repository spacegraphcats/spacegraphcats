import pstats

if __name__ == "__main__":
    p = pstats.Stats("search_stats")
    p.strip_dirs()
    p.sort_stats("time")  # also try: time, cumulative
    p.print_stats("search")  # pass in string to limit. e.g. 'search'
