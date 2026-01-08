class Contig:
    """
    """
    def __init__(self, seq, start_node, final_node, walk_weights, receiving_hubs_visited, k, bubble = False, tip = False):
        self.seq = seq
        self.start_node = start_node
        self.final_node = final_node
        self.walk_weights = walk_weights
        self.receiving_hubs_visited = receiving_hubs_visited
        self.k = k
        self.bubble = bubble
        self.tip = tip
        self.tag = ""
        self.branch_kmer = None

    def __iadd__(self, other):
        if not isinstance(other, Contig):
            raise TypeError(f"You can't sum {type(self)} object with {type(other)}")
        self.seq += other.seq[self.k-1:]
        self.start_node = self.start_node
        self.final_node = other.final_node
        self.walk_weights += other.walk_weights[self.k-1:]
        self.receiving_hubs_visited = self.receiving_hubs_visited.union(other.receiving_hubs_visited)
        return self

    def __add__(self, other):
        if not isinstance(other, Contig):
            raise TypeError(f"You can't sum {type(self)} object with {type(other)}")
        new = Contig(
            seq=self.seq + other.seq[self.k-1:],
            start_node=self.start_node,
            final_node=other.final_node,
            walk_weights=self.walk_weights + other.walk_weights[self.k-1:],
            receiving_hubs_visited=self.receiving_hubs_visited.union(other.receiving_hubs_visited),
            k=self.k
        )
        new.bubble = self.bubble or other.bubble
        new.tip = self.tip or other.tip
        return new

    def copy(self):
        new = Contig(
            seq=self.seq, 
            start_node=self.start_node, 
            final_node=self.final_node, 
            walk_weights=list(self.walk_weights), 
            receiving_hubs_visited=set(self.receiving_hubs_visited), 
            k=self.k, 
            bubble=self.bubble, 
            tip=self.tip
        )
        new.tag = self.tag  
        new.branch_kmer = self.branch_kmer  
        return new

    @classmethod
    def from_walk_linear(cls, walk_result, k, initial_node):
        seq, final_node, weights, hubs = walk_result
        return cls(seq, initial_node, final_node, weights, hubs, k)

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return f"Contig(len={len(self)}, tag={self.tag or 'untagged'})"
    
    def is_circular(self):
        return self.start_node == self.final_node

    def is_lariat(self):
        return self.final_node.seq in self.receiving_hubs_visited and not self.is_circular()