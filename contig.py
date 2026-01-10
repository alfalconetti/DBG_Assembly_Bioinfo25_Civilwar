class Contig:
    """
    Represents an assembled contig with all the necessary attributes
    for statistics and for topological classification.
    
    A Contig is the result of walking through the De Bruijn Graph,
    potentially traversing bubbles, tips, and hubs.
    
    :ivar seq: Assembled DNA sequence
    :type seq: str
    :ivar start_node: Node where the contig begins
    :type start_node: Node
    :ivar final_node: Node where the contig ends
    :type final_node: Node
    :ivar walk_weights: Coverage weights for each k-mer in the path
    :type walk_weights: list[int]
    :ivar receiving_hubs_visited: Set of receiving hub k-mers encountered
    :type receiving_hubs_visited: set[str]
    :ivar k: K-mer size
    :type k: int
    :ivar hub: Whether the contig traversed through a hub
    :type hub: bool
    :ivar bubble: Whether the contig includes a resolved bubble
    :type bubble: bool
    :ivar tip: Whether the contig includes an assembled tip
    :type tip: bool
    :ivar remaining: Whether from fork_in_hub residual processing
    :type remaining: bool
    :ivar from_catch_all: Whether from final circular/isolated search
    :type from_catch_all: bool
    :ivar tag: Topological classification tag
    :type tag: str
    :ivar branch_kmer: K-mer where lariat branches (for lariat structures)
    :type branch_kmer: str or None
    :ivar bubble_ends: (source_seq, dest_seq) of bubble if present
    :type bubble_ends: tuple[str, str]
    :ivar tip_hub: Hub k-mer where tip originated
    :type tip_hub: str or None
    :ivar hub_seq: K-mer of hub traversed
    :type hub_seq: str or None
    """
    def __init__(
        self,
        seq,
        start_node,
        final_node,
        walk_weights,
        receiving_hubs_visited,
        k,
        hub=False,
        bubble=False,
        tip=False,
        bubble_ends=(None, None),
        tip_hub=None,
        hub_seq=None,
        remaining=False,
        from_catch_all=False
    ):
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
        self.from_catch_all = from_catch_all
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
    def from_walk_linear(cls, walk_result, k, initial_node, hub = False, hub_seq = None, remaining = False, from_catch_all = False):
        seq, final_node, weights, hubs = walk_result
        return cls(seq, initial_node, final_node, weights, hubs, k, hub = hub, hub_seq = hub_seq, remaining = remaining, from_catch_all = from_catch_all)

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return f"Contig(len={len(self)}, tag={self.tag or 'untagged'})"
    
    def is_circular(self):
        return self.start_node == self.final_node

    def is_lariat(self):
        if self.is_circular():
            return False
        if self.final_node.seq in self.receiving_hubs_visited and self.final_node.seq in self.seq[:-1]:
            self.branch_kmer = self.final_node.seq
            return True
        else:
            return False

    def is_bubble_ending(self):
        return self.bubble and self.bubble_ends[1] == self.final_node.seq
    
    def get_fasta_header(self, contig_id):
        from statistics import median
        # All k-mer sequences are truncated to 10 bp
        median_cov = median(self.walk_weights) if self.walk_weights else 0
        header_parts = [
            f"contig_{contig_id}",
            f"len={len(self)}",
            f"tag={self.tag}",
            f"cov={median_cov:.1f}"
        ]

        # 2. Metadati topologici condizionali (solo se attinenti)
        if self.branch_kmer:
            header_parts.append(f"branch={self.branch_kmer}")
        
        if self.bubble_ends != (None, None):
            header_parts.append(f"bubble_ends={self.bubble_ends[0][:10]}..-{self.bubble_ends[1][:10]}..")
        
        if self.tip_hub:
            header_parts.append(f"tip_from={self.tip_hub[:10]}..")
            
        if self.hub_seq:
            header_parts.append(f"hub_origin={self.hub_seq[:10]}..")

        # Uniamo tutto con il pipe
        return ">" + "|".join(header_parts)

    def get_tag(self):
        # ALL CIRCULAR NOT ISOLATED START FROM A HUB AND END IN HUB. RESIDUAL TUTTI START IN HUB
        if self.start_node.is_hub():
            hub_starting = "_hub_starting"
        else:
            hub_starting = ""
        if self.final_node.is_hub():
            hub_ending = "_hub_ending"
        else:
            hub_ending = ""
        if self.tip:
            self.tag = "linear_with_tip"
            return # Se tip non servono altri controlli, solo lineare
        if self.is_circular(): # ALL CIRCULAR NOT ISOLATED START FROM A HUB
            if self.remaining:
                self.tag = "circular_fork_hub_residual"
            elif self.from_catch_all:
                self.tag = "circular_isolated"
            elif self.is_bubble_ending():
                self.tag = "circular_bubble_resolved_ending"
            elif self.bubble and not self.is_bubble_ending():
                if self.hub:
                    self.tag = "circular_with_bubble_through_hub"
                else:
                    self.tag = "circular_with_bubble"
            elif self.hub:
                self.tag = "circular_through_hub" 
            elif not self.hub and self.final_node.is_hub():
                self.tag = "circular_self_loop_hub" # parte e ritorna da hub
            else:
                self.tag = "unknown_circular_like"
        elif self.is_lariat():
            if self.remaining:
                self.tag = "lariat_fork_hub_residual" + hub_ending
            elif self.from_catch_all:
                self.tag = "unknown_lariat_like_from_catch_all" 
            elif self.is_bubble_ending():
                self.tag = "lariat_bubble_resolved_ending" + hub_starting
            elif self.bubble and not self.is_bubble_ending():
                if self.hub:
                    self.tag = "lariat_with_bubble_through_hub" + hub_starting + hub_ending
                else:
                    self.tag = "lariat_with_bubble" + hub_starting + hub_ending
            elif self.hub:
                self.tag = "lariat_through_hub" + hub_starting
            elif not self.hub and self.final_node.is_hub():
                self.tag = "lariat_branch_at_hub" + hub_starting
            elif not self.hub:
                self.tag = "lariat_simple" + hub_starting
            else:
                self.tag = "unknown_lariat_like" + hub_starting + hub_ending
        else:
            if self.remaining:
                self.tag = "linear_fork_hub_residual" + hub_ending
            elif self.from_catch_all:
                self.tag = "unknown_linear_like_from_catch_all"
            elif self.is_bubble_ending():
                self.tag = "linear_bubble_resolved_ending" + hub_starting + hub_ending
            elif self.bubble and not self.is_bubble_ending():
                if self.hub:
                    self.tag = "linear_bubble_resolved_extended_through_hub" + hub_starting + hub_ending
                else:
                    self.tag = "linear_bubble_resolved_extended" + hub_starting + hub_ending
            elif not self.bubble and not self.tip and self.hub:
                self.tag = "linear_through_hub" + hub_starting + hub_ending
            elif not self.hub and self.final_node.is_hub():
                self.tag = "linear_hub_ending" + hub_starting
            elif not self.hub:
                self.tag = "linear_simple" + hub_starting
            else:
                self.tag = "unknown_linear_like" + hub_starting + hub_ending