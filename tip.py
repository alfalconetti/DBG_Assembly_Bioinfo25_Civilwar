from statistics import median
from contig import Contig

class Tip:
    """
        This class represents a "Tip" in the De Bruijn Graph.

        A Tip is a short branch that disconnects from the main graph structure.
        Biologically, tips are usually caused by sequencing errors occurring at the
        very ends of the reads. Since the error k-mers don't match anything else,
        they form a short path that leads nowhere (or starts from nowhere).

        We distinguish two types of tips:
            - Type I (Outgoing):** Starts at a graph hub and ends in a dead end (no outgoing edges).
            - Type II (Incoming):** Starts at a dead end (no incoming edges) and merges into a hub.

        :ivar initial_node: The starting node of the tip.
        :type initial_node: Node
        :ivar final_node: The ending node of the tip.
        :type final_node: Node
        :ivar path: The list of edges forming the tip sequence.
        :type path: list['Edge']
        :ivar tip_type: The classification of the tip (1 for Outgoing, 2 for Incoming).
        :type tip_type: int
        :ivar visited: Flag to track if we have processed this bubble. Defaults to False.
        :type visited: bool
    """
    def __init__(self, initial_node: "Node", final_node: "Node", path: list["Edge"], tip_type: int):
        """
            Creates a new Tip instance and validates its topology.

            :param initial_node: The start of the branch.
            :type initial_node: Node
            :param final_node: The end of the branch.
            :type final_node: Node
            :param path: The sequence of edges.
            :type path: list[Edge]
            :param tip_type: 1 for Outgoing Tip, 2 for Incoming Tip.
            :type tip_type: int
            :raises ValueError: If the nodes do not match the expected topology for the declared tip type.
        """
        self.initial_node = initial_node
        self.final_node = final_node
        self.path = path
        self.tip_type = tip_type
        self.visited = False
        if tip_type == 1:
            if not self.initial_node.is_hub():
                raise ValueError("Type I tip must start from hub")
            if len(final_node.outgoing) != 0:
                raise ValueError("Type I tip must end in dead end")
        elif tip_type == 2: 
            if len(initial_node.incoming) != 0:
                raise ValueError("Type II tip must start from node with no incoming")
            if not final_node.is_receiving_hub():
                raise ValueError("Type II tip must end in receiving hub")
        else:
            raise ValueError(f"Invalid tip_type: {tip_type}. Must be 1 or 2")
    
    def _get_path_coverage(self, path: list["Edge"]) -> float:
        """
            Calculates the median coverage of the tip's path.

            :param path: The list of edges to analyze.
            :type path: list[Edge]
            :return: The median weight.
            :rtype: float
        """
        if not path:
            return 0
        weights = [edge.weight for edge in path]
        return median(weights)

    def _resolve(self, graph: "Dbg", cov_threshold: int) -> bool:
        """
            Removes the tip from the graph if its coverage is below the threshold.
.
            :param graph: The main graph object (needed to delete edges).
            :type graph: Dbg
            :param cov_threshold: The coverage cutoff. Tips below this are removed.
            :type cov_threshold: int
            :return: True if the tip was removed, False otherwise.
            :rtype: bool
        """
        popped = False
        cov = self._get_path_coverage(self.path)
        if cov < cov_threshold: #sotto threshold, da rimuovere
            graph._rm_path_edges(self.path) #rimozione di quel path
            popped = True
        return popped

    def _get_contig_from_path(self, k):
        """
            Reconstructs the DNA sequence for the tip (dead-end branch).

            It walks through the edges from the branching hub to the sink node,
            stitching the k-mers together (handling the k-1 overlap) and collecting
            coverage weights.

            :param k: The k-mer length used for the assembly.
            :type k: int
            :return: A fully initialized Contig instance representing the dead-end path.
            :rtype: Contig
        """
        sequence = self.path[0].seq
        path_weights = [self.path[0].weight]
        receiving_hubs_visited = set()
        if self.path[0].source.is_receiving_hub():
                receiving_hubs_visited.add(self.path[0].source.seq)
        for edge in self.path[1:]:
            sequence += edge.seq[-1]
            path_weights.append(edge.weight)
            curr = edge.source
            if curr.is_receiving_hub():
                receiving_hubs_visited.add(curr.seq)
        return Contig(seq = sequence, start_node = self.initial_node, final_node = self.final_node, walk_weights = path_weights, receiving_hubs_visited = receiving_hubs_visited, k = k, tip = True, tip_hub = self.initial_node.seq)