class Contig:
	"""
	"""
	def __init__(seq, final_node, walk_weights, receiving_hubs_visited):
		self.seq = seq
		self.final_node = final_node
		self.walk_weights = walk_weights
		self.receiving_hubs_visited

	def __iadd__(self, other):
		
