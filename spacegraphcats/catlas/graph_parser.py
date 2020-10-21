"""Parser and writer for simple graph format."""

from typing import Callable, List, Tuple, Iterable


def _parse_line(line: str, split_on=" ") -> List[int]:
    return list(map(int, map(str.strip, line.split(split_on))))


def parse(
    graph_file,
    create_vertices: Callable[[int], None] = None,
    add_edge: Callable[[int, int], None] = None,
):
    """Parse a graph and call provided methods with vertices and edges."""
    # read vertices
    num_vertices = int(graph_file.readline().strip())

    if create_vertices is not None:
        create_vertices(num_vertices)

    # read edges
    next_line = graph_file.readline()
    while len(next_line) > 1:
        parsed = _parse_line(next_line)
        add_edge(parsed[0], parsed[1])
        next_line = graph_file.readline()


def parse_minhash(minhash_file, add_minhash: Callable[[int, List[int]], None]):
    """Parse minhash (.mxt) file."""
    for line in minhash_file:
        if len(line) < 2:
            continue
        parsed = _parse_line(line)
        add_minhash(parsed[0], parsed[1:])


def write(graph_file, num_vertices: int, edges: Iterable[Tuple[int, int]]):
    """Write an edgelist into an .ext file."""
    graph_file.write("{}\n".format(num_vertices))
    for u, v in edges:
        graph_file.write("{} {}\n".format(u, v))
