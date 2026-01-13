Modules documentation
=====================

High-Level Assembly Pipeline
----------------------------
This is the main entry point for the assembler. Use this function for
automated workflows.

.. automodule:: src.pipeline_function
   :members:
   :undoc-members:
   :show-inheritance:

Graph Engine & Topology
-----------------------
These modules handle the core De Bruijn Graph logic and topological structures.

.. automodule:: src.dbg
   :members:
   :show-inheritance:

.. automodule:: src.bubble
   :members:
   :show-inheritance:

.. automodule:: src.tip
   :members:
   :show-inheritance:


.. automodule:: src.contig
   :members:
   :show-inheritance:

Low-Level Components
--------------------
Internal graph elements (Node, Edge)

.. automodule:: src.node
   :members:

.. automodule:: src.edge
   :members:



Data Parsing & Utilities
------------------------
Low-level components for file handling and DNA sequence processing.

.. automodule:: src.parser
   :members:
   :show-inheritance:

.. automodule:: src.kmer
   :members:

.. automodule:: src.utils
   :members:


