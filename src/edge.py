class Edge:
    """
        This class represents a connection (edge) in the De Bruijn Graph.

        Basically, an edge is a k-mer found in your reads. It acts as a bridge
        connecting a prefix (source node) to a suffix (destination node).
        We use this object to keep track of the sequence, how many times it appears
        (weight), and whether we have already used it during the assembly process.

        :ivar seq: The nucleic acid sequence of the k-mer.
        :type seq: str
        :ivar source: The source node (prefix).
        :type source: Node
        :ivar destination: The destination node (suffix).
        :type destination: Node
        :ivar weight: The multiplicity of this k-mer in the dataset
        :type weight: int
        :ivar visited: Flag indicating if the edge has been traversed during assembly. Defaults to False.
        :type visited: bool
    """

    def __init__(self, seq_k: str, source: "Node", destination: "Node"):
        """
            Creates a new Edge.

            :param seq_k: The k-mer string.
            :type seq_k: str
            :param source: The node object where this edge starts.
            :type source: Node
            :param destination: The node object where this edge ends.
            :type destination: Node
        """
        self.seq = seq_k
        self.source = source
        self.destination = destination
        self.weight = 1
        self.visited = False

    def __repr__(self):
        return self.seq
