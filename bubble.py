from statistics import median

class Bubble:
    def __init__(self, source_node, dest_node, paths):
        self.source = source_node #nodo hub
        self.dest = dest_node #nodo receiving hub
        self.paths = paths #lista di liste di edges

    def get_path_coverage(self, path):
        if not path:
            return 0
        weights = [edge.weight for edge in path]
        return median(weights)

    def _resolve(self, graph, cov_threshold, ratio_threshold=0.1):
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

    def _get_seq_from_path(self):
        for path in self.paths:
            sequence = path[0].seq
            path_weights = [path[0].weight]
            for edge in path[1:]:
                sequence += edge.seq[-1]
                path_weights.append(edge.weight)
        yield sequence, path_weights