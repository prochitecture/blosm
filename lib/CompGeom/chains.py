
from collections import defaultdict

# find_paths(edges) identifies valid paths in a directed acyclic graph (DAG) based on specific constraints:
#
# - Each edge is represented as a tuple (start_node, end_node, edge_id).
# - A valid path consists of edges where each intermediate node is visited exactly twice (once as an end and once as a start).
# - Paths stop at branching nodes (nodes with more than one incoming or outgoing edge) or at endpoints.
#
# The function:
# 1. Constructs an adjacency list representation of the graph.
# 2. Identifies start nodes (nodes with in-degree or out-degree ≠ 1).
# 3. Iterates through the graph to find all valid paths.
# 4. Returns a list of these paths.
#
# Example input:
#     edges = [('a','b','A'), ('b','c','B'), ('c','d','C'), ('d','e','D'),
#              ('e','f','E'), ('e','g','F'), ('h','i','G'), ('i','a','H'), ('k','a','J')]
#
# Example output:
#     [
#         [('a','b','A'), ('b','c','B'), ('c','d','C'), ('d','e','D')],
#         [('e','f','E')],
#         [('e','g','F')],
#         [('h','i','G'), ('i','a','H')],
#         [('k','a','J')]
#     ]
#
def find_paths(edges):
    graph = defaultdict(list)
    indegree = defaultdict(int)
    outdegree = defaultdict(int)
    
    # Step 1: Constructs an adjacency list representation of the graph.
    for u, v, eid in edges:
        graph[u].append((v, eid))
        indegree[v] += 1
        outdegree[u] += 1
    
    # Step 2: Identifies start nodes (nodes with in-degree or out-degree ≠ 1).
    start_nodes = set()
    for node in set(indegree) | set(outdegree):
        if indegree[node] != 1 or outdegree[node] != 1:
            start_nodes.add(node)
    
    paths = []
    visited = set()
    
    def extend_path(node):
        path = []
        while node in graph and len(graph[node]) == 1 and indegree[node] == 1:
            next_node, edge_id = graph[node][0]
            path.append((node, next_node, edge_id))
            visited.add(node)
            node = next_node
            if node in start_nodes:
                break
        return path
    
    # Step 3: Iterates through the graph to find all valid paths.
    for node in start_nodes:
        for neighbor, edge_id in graph[node]:
            if node not in visited:
                path = [(node, neighbor, edge_id)] + extend_path(neighbor)
                paths.append(path)
    
    return paths


def build_adjacency_list(segments):
    # Creates an adjacency list from the given segments.
    # Each node stores which other nodes it is connected to via a segment.
    # Additionally, a mapping of segment IDs to their endpoints is stored.

    # Args:
    #     segments (list of tuples): Each tuple consists of two endpoints
    #                                (tuples of coordinates) and a segment ID (string).

    # Returns:
    #     tuple: A dictionary representing the adjacency list and a dictionary mapping
    #            segment IDs to their endpoints.

    adjacency_list = defaultdict(list)
    segment_map = {}
    
    for e0, e1, ID in segments:
        adjacency_list[e0].append((e1, ID))
        adjacency_list[e1].append((e0, ID))
        segment_map[ID] = (e0, e1)
    
    return adjacency_list, segment_map

def find_longest_path(start, adjacency_list, used_segments):
    # Creates a maximally long path starting from 'start'.
    # The next unused segment is always chosen until no more segments are available.

    # Args:
    #     start (tuple): The starting node (coordinates as a tuple).
    #     adjacency_list (dict): The adjacency list representing the graph.
    #     used_segments (set): A set of segment IDs that have already been used.

    # Returns:
    #     list: A list of tuples representing the path, where each tuple contains
    #           (start_node, end_node, segment_id).

    path = []
    current = start
    
    while True:
        next_step = None
        for neighbor, segment_id in adjacency_list[current]:
            if segment_id not in used_segments:
                next_step = (neighbor, segment_id)
                break
        
        if not next_step:
            break  # Kein weiteres Segment verfügbar
        
        neighbor, segment_id = next_step
        path.append((current, neighbor, segment_id))
        used_segments.add(segment_id)
        current = neighbor
    
    return path

def find_all_lines(segments):
    # Decomposes the segments into maximally long, connected paths.
    # Each segment may only appear in a single path.
    # Preference is given to starting at nodes with only one connection.

    # Args:
    #     segments (list of tuples): Each tuple consists of two endpoints
    #                                (tuples of coordinates) and a segment ID (string).

    # Returns:
    #     list: A list of paths, where each path is a list of tuples
    #           (start_node, end_node, segment_id).

    adjacency_list, segment_map = build_adjacency_list(segments)
    used_segments = set()
    results = []
    
    # Startpunkte mit nur einer Verbindung bevorzugen
    endpoints = [node for node, neighbors in adjacency_list.items() if len(neighbors) == 1]
    
    # Falls keine Endpunkte existieren (Zyklen), einfach irgendwo starten
    nodes_to_explore = endpoints if endpoints else list(adjacency_list.keys())
    
    for node in nodes_to_explore:
        for neighbor, segment_id in adjacency_list[node]:
            if segment_id not in used_segments:
                path = find_longest_path(node, adjacency_list, used_segments)
                if path:
                    results.append(path)
    
    return results
