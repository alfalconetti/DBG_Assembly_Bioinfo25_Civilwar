from statistics import median

class Tip:
    def __init__(self, initial_node, final_node, path, tip_type):
        self.initial_node = initial_node
        self.final_node = final_node
        self.path = path
        self.tip_type = tip_type
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
    
    def _get_path_coverage(self, path):
        if not path:
            return 0
        weights = [edge.weight for edge in path]
        return median(weights)    

    def _resolve(self, graph, cov_threshold):
        popped = False
        cov = self._get_path_coverage(self.path)
        if cov < cov_threshold: #sotto treshold, da rimuovere
            graph._rm_path_edges(self.path) #rimozione di quel path
            popped = True
        return popped

    def _get_seq_from_path(self):
        sequence = path[0].seq
        for edge in path[1:]:
            sequence += edge.seq[-1]
        return sequence