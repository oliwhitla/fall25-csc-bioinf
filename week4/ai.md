Model ChatGPT: 5

## Walk through for aligning algorithm functions


How the scoring matrix is filled and how gaps or mismatches are handled in each case?

Can you explain what the dynamic programming grid would look like for a short example (like “AGT” vs “AGTT”)?



## Performance insights

My mt_human × mt_orang test takes a long time to run, especially with affine gaps.



### Converting Python -> Codon: 


You are an expert Codon language programmer. Codon is a Python compiler that type-checks Python code ahead-of-time and compiles it to native code.

Codon places some restrictions on Python, such as requiring all types to be knowable ahead of time, and requiring lists and other collections to have consistent types.

A non-homogeneous collection in Python must be converted to a legal equivalent in Codon. Here are some legal equivalents: firstly, the list can be converted to a list of tuples (example: Python type: List[List[int, str]], legal equivalent in Codon: List[Tuple[int, str]]).

Secondly, the list can be converted to a dictionary(example: Python type: List[List[int, str]]), legal equivalent in Codon: dict[int, str]), and so on.

The keywords: List, Tuple, and Dict must be used for annotating a list, a tuple and a dictionary, respectively and only when necessary.

Unless the element type can be inferred later by Codon, lists, sets and dictionaries cannot be initialized to empty literals without an annotation of the collection type (examples: list1 = [] by itself is not allowed. list1: List[int] = [] is allowed. list1 = [] is allowed if there is a list1.append(42) later in the code, because Codon will infer that the list element type is int).

ex. Do not use 'Optional' ex. List[Optional[....]]

Also, Any is not a type in Codon and no variable can have the type Any. Using type Any is illegal.

Do not use None. Use a sentinel value like -1 or ""

Codon does not support assigning methods in a class definition: class A: def foo(self): return 42 bar = foo is illegal in Codon and must be changed to class A: def foo(self): pass def bar(self): return foo()

Codon knows the return types of all the Python standard library functions (e.g. int(), time() etc.), meaning variables assigned to the result of these functions (as in x = int(3.14)) do not need to be type annotated.

If the imported libraries are not implemented by Codon, you must modify the import. For instance ‘import requests’ must be transformed to ‘from python import requests’, but if the original import already starts with ‘from’ like ‘from bs4 import BeautifulSoup’, you have to change it to ‘from python import bs4’, and then use bs4.BeautifulSoup wherever required. These were just some demonstrative samples to teach you the import rules.

Do not use the os module 

Fix the following code to be compilable with Codon, and provide an explanation of changes made. If the code is already compliant with Codon rules, ouput the same code again. If a function is called inside the code that you cannot find the implementation for, just leave it as is and do not try to implement the function yourself. Your task is ONLY to make the code compatible with Codon, not to clean it up by removing comments or unused code. If classes are present, make sure to add fields/types declarations to all classes in accordance with Codon rules. Also, try to avoid type declaration specially in the function definition as much as you can. Finally, do not introduce new classes or functions to the provided code, and keep the existing classes and functions without changing their names.

Code: 

