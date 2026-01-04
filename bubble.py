class Bubble:
    def __init__(self, source_node, dest_node, paths):
        self.source = source_node #nodo hub
        self.dest = dest_node #nodo receiving hub
        self.paths = paths #lista di liste di edges