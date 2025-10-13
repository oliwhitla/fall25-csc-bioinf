AI Model: ChatGPT 5

I followed the assignment’s suggested prompts (e.g., “Port UPGMA to Codon with homogeneous lists,” “Use Codon’s numpy import style”) as the starting point for my implementation. I applied them to my files (upgma_codon.py, tree_codon.py) and adjusted type annotations and array handling to satisfy Codon’s ahead-of-time type checking.

1) I ported the original code to plain Python (CPython)

You are an expert in Codon / typed-Python AOT compilation. Port the following CPython module to typed Python (TPython/Codon). Keep behavior identical but enforce Codon’s rules:

All collections must be homogeneous; annotate as needed (e.g., List[int], Dict[str, float]).

Avoid Any; use concrete types. Mark nullable fields as Optional[T].

Classes must declare fields explicitly (dataclass-style).

Ensure float/int operations are type-consistent (e.g., use 2.0 for float division).
