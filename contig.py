class Contig:
    """
    """
    def __init__(self, seq, start_node, final_node, walk_weights, receiving_hubs_visited, k, hub = False, bubble = False, tip = False, bubble_ends = (None, None), tip_hub = None, hub_seq = None, remaining = False, circular = False):
        self.seq = seq
        self.start_node = start_node
        self.final_node = final_node
        self.walk_weights = walk_weights
        self.receiving_hubs_visited = receiving_hubs_visited
        self.k = k
        self.bubble = bubble
        self.tip = tip
        self.hub = hub
        self.remaining = remaining
        self.circular = circular
        self.tag = ""
        self.branch_kmer = None
        self.bubble_ends = bubble_ends
        self.tip_hub = tip_hub
        self.hub_seq = hub_seq


    def __iadd__(self, other):
        if not isinstance(other, Contig):
            raise TypeError(f"You can't sum {type(self)} object with {type(other)}")
        self.seq += other.seq[self.k-1:]
        self.start_node = self.start_node
        self.final_node = other.final_node
        self.walk_weights += other.walk_weights
        self.receiving_hubs_visited = self.receiving_hubs_visited.union(other.receiving_hubs_visited)
        return self

    def __add__(self, other):
        if not isinstance(other, Contig):
            raise TypeError(f"You can't sum {type(self)} object with {type(other)}")
        new = Contig(
            seq=self.seq + other.seq[self.k-1:],
            start_node=self.start_node,
            final_node=other.final_node,
            walk_weights=self.walk_weights + other.walk_weights,
            receiving_hubs_visited=self.receiving_hubs_visited.union(other.receiving_hubs_visited),
            k=self.k,
            bubble=self.bubble or other.bubble,
            tip = self.tip or other.tip,
            hub = self.hub or other.hub,
            bubble_ends = other.bubble_ends if other.bubble_ends != (None, None) else self.bubble_ends,
            tip_hub = other.tip_hub if other.tip_hub is not None else self.tip_hub,
            hub_seq = other.hub_seq if other.hub_seq is not None else self.hub_seq
        )
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
            tip=self.tip,
            hub=self.hub,
            bubble_ends=self.bubble_ends,
            tip_hub=self.tip_hub,
            hub_seq=self.hub_seq
        )
        new.tag = self.tag  
        new.branch_kmer = self.branch_kmer  
        return new

    @classmethod
    def from_walk_linear(cls, walk_result, k, initial_node, hub = False, hub_seq = None, remaining = False, circular = False):
        seq, final_node, weights, hubs = walk_result
        return cls(seq, initial_node, final_node, weights, hubs, k, hub = hub, hub_seq = hub_seq, remaining = remaining, circular = circular)

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return f"Contig(len={len(self)}, tag={self.tag or 'untagged'})"
    
    def is_circular(self):
        return self.start_node == self.final_node

    def is_lariat(self):
        return self.final_node.seq in self.receiving_hubs_visited and not self.is_circular()

    def is_bubble_ending(self):
        return self.bubble and self.bubble_ends[1] == self.final_node.seq
    

    def get_tag(self):
        if self.tip:
            self.tag = "linear_with_tip"
        if self.is_circular():
            if self.remaining:
                self.tag = "circular_fork_hub_residual"
            elif self.circular:
                self.tag = "circular_isolated"
            elif self.is_bubble_ending():
                self.tag = "circular_bubble_resolved_ending" 
            elif self.bubble and not self.is_bubble_ending():
                self.tag = "circular_with_bubble"
            elif self.hub:
                self.tag = "circular_through_hub"
            elif not self.hub and self.final_node.is_hub():
                self.tag = "circular_ending_hub"
            elif not self.hub:
                self.tag = "circular_simple"
        elif self.is_lariat():
            if self.remaining:
                self.tag = "lariat_fork_hub_residual"
            elif self.circular:
                self.tag = "unknown_contig_type_lariat_like"
            elif self.is_bubble_ending():
                self.tag = "lariat_bubble_resolved_ending"
            elif self.bubble and not self.is_bubble_ending():
                self.tag = "lariat_with_bubble"
            elif self.hub:
                self.tag = "lariat_through_hub"
            elif not self.hub and self.final_node.is_hub():
                self.tag = "lariat_hub_branching"
            elif not self.hub:
                self.tag = "lariat_simple"
        else:
            if self.remaining:
                self.tag = "linear_fork_hub_residual"
            elif self.circular:
                self.tag = "unknown_contig_type"
            elif self.is_bubble_ending():
                self.tag = "linear_bubble_resolved_ending"
            elif self.bubble and not self.is_bubble_ending():
                self.tag = "linear_bubble_resolved_extended"
            elif not self.bubble and not self.tip and self.hub:
                self.tag = "linear_through_hub"
            elif not self.hub and self.final_node.is_hub():
                self.tag = "linear_ending_hub"
            elif not self.hub:
                self.tag = "linear_simple"
