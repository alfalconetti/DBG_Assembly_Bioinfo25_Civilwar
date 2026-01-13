class Node:
    """
        This class represents a node in the De Bruijn Graph, corresponding to a (k-1)-mer.

        :ivar seq: The DNA sequence of this (k-1)-mer.
        :type seq: str
        :ivar incoming: A list of edges pointing TO this node (incoming connections).
        :type incoming: list['Edge']
        :ivar outgoing: A list of edges starting FROM this node (outgoing connections).
        :type outgoing: list['Edge']
    """

    def __init__(self, seq: str):
        """
            Creates a new Node instance.

            :param seq: The (k-1)-mer sequence string.
            :type seq: str
        """
        self.seq = seq
        self.incoming = []
        self.outgoing = []

    def is_hub(self):
        """
            Checks if this node is a branching point (bifurcation).

            Biologically, this often indicates a repeat or a variant where the
            sequence splits into multiple possible paths.

            :return: True if there is more than one outgoing edge, False otherwise.
            :rtype: bool
        """
        return len(self.outgoing) > 1

    def is_receiving_hub(self):
        """
            Checks if this node is a merge point (convergence).

            This indicates a point where multiple paths rejoin, for example
            at the end of a bubble (variant) or a repeat.

            :return: True if there is more than one incoming edge, False otherwise.
            :rtype: bool
        """
        return len(self.incoming) > 1

    def __eq__(self, other):
        return self.seq == other.seq
