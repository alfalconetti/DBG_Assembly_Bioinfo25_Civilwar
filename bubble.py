from statistics import median
from contig import Contig
class Bubble:
    """
        This class represents a "Bubble" in the De Bruijn Graph.

        A bubble occurs when the graph splits into multiple paths (at a source node)
        and then rejoins later (at a destination node).
        Biologically, this usually represents a local variation in the genome,
        such as a SNP (Single Nucleotide Polymorphism), a small indel, or a sequencing error.

        :ivar source_node: The node where the divergence begins (the fork).
        :type source_node: Node
        :ivar dest: The node where the paths reconverge.
        :type dest: Node
        :ivar paths: A list containing the different routes (lists of edges) inside the bubble.
        :type paths: list[list['Edge']]
        :ivar visited: Flag to track if we have processed this bubble. Defaults to False.
        :type visited: bool
    """
    def __init__(self, source_node: "Node", dest_node: "Node", paths: list[list["Edge"]]):
        """
            Creates a new Bubble instance.

            :param source_node: The starting hub of the bubble.
            :type source_node: Node
            :param dest_node: The ending receiving hub of the bubble.
            :type dest_node: Node
            :param paths: A list of paths, where each path is a list of Edge objects.
            :type paths: list[list[Edge]]
        """
        self.source_node = source_node #nodo hub
        self.dest_node = dest_node #nodo receiving hub
        self.paths = paths #lista di liste di edges
        self.visited = False

    def get_path_coverage(self, path: list["Edge"]) -> float:
        """
            Calculates the median coverage of a specific path within the bubble.

            :param path: The path to analyze (a list of edges).
            :type path: list[Edge]
            :return: The median weight of the edges in the path.
            :rtype: float
        """
        if not path:
            return 0
        weights = [edge.weight for edge in path]
        return median(weights)

    def _resolve(self, graph: "Dbg", cov_threshold: int, ratio_threshold: float =0.1) -> bool:
        """
            Attempts to "pop" the bubble by removing low-quality paths.

            It compares the paths based on their coverage. It keeps the best path
            and removes the "losers" if they are either:
            1. Below a hard coverage threshold.
            2. Too weak compared to the best path (likely sequencing errors).

            :param graph: The main graph object (needed to physically delete edges).
            :type graph: Dbg
            :param cov_threshold: Minimum coverage required to keep a path.
            :type cov_threshold: int
            :param ratio_threshold: Minimum ratio required relative to the best path's coverage.
            :type ratio_threshold: float
            :return: True if at least one path was removed (bubble popped), False otherwise.
            :rtype: bool
        """
        if len(self.paths) < 2: #se il nodo ricevente ha meno di due path... non è una bolla
            return False

        path_stats = [] #lista di tuple, con coverage e path(lista edges)
        for p in self.paths: #itero sulla lista di liste di edges
            cov = self.get_path_coverage(p)
            path_stats.append((cov, p))

        path_stats.sort(key=lambda x: x[0], reverse=True) #ordinamento in base al cov

        best_cov, best_path = path_stats[0]
        popped = False

        for i in range(1, len(path_stats)): #dal secondo in poi, primo path sempre tenuto
            loser_cov, loser_path = path_stats[i]

            should_remove = False

            if loser_cov < cov_threshold: #sotto treshold, da rimuovere
                should_remove = True

            elif loser_cov < (best_cov * ratio_threshold): #inferiore di (10) volte rispetto alla mediana cov del path migliore
                should_remove = True

            if should_remove:
                graph._rm_path_edges(loser_path) #rimozione di quel path
                popped = True

        return popped

    def _get_contigs_from_path(self, k):
        """
            A generator that reconstructs the DNA sequences for all paths in the bubble.

            It walks through the edges of each path, stitching the k-mers together
            (handling the k-1 overlap) and collecting the coverage weights.

            :yield: A tuple containing (assembled_sequence, list_of_weights).
            :rtype: Iterator[tuple[str, list[int]]]
        """
        for path in self.paths: #generatore di sequenza e lista pesi per ciascun path(lista edge)
            sequence = path[0].seq
            path_weights = [path[0].weight]
            receiving_hubs_visited = set()
            if path[0].source.is_receiving_hub():
                receiving_hubs_visited.add(path[0].source.seq)
            for edge in path[1:]:
                sequence += edge.seq[-1]
                path_weights.append(edge.weight)
                curr = edge.source
                if curr.is_receiving_hub():
                    receiving_hubs_visited.add(curr.seq)
            yield Contig(seq = sequence, start_node = self.source_node, final_node = self.dest_node, walk_weights = path_weights, receiving_hubs_visited = receiving_hubs_visited, k = k, bubble = True)