class Node:
    def __init__(self, seq):
        self.seq = seq
        self.incoming = []
        self.outgoing = []

    def is_hub(self):
        return len(self.outgoing) > 1
