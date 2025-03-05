
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

