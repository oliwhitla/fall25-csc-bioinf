AI Model: Chat GPT 5

Prompts: 

Context: Report.md was in madarin 

Prompt: Translate the attached report from mandarin to english. Keep technical term accurate. Don't summarize.


Context: Task was to re-implement a Python module so it would run under codon.

Prompt: Convert the following Python module to a Codon-compatible version. Avoid libraries Codon may not support. No None in fields or returns. 

Notes: 
- codon rejects Optional int|None, str|None fields/returns
- chat added extra type hints which was sometimes overkill
- replace None with sentinel values ex. -1
- Debugging involved resolving errors like:
'>contig_%d\n' % … not compatible
Node[NoneType,…] does not match expected type …


Context: Creating evaluate.sh

Prompt: write a helper function to format seconds - 0:00:00

Prompt: I need to loop datasets 1, 2, 3, 4. Get the runtime for each dataset and format it using my helper function. Data 4 needs ulimit -s 524288 to allow for more stack space. Ignore writing a function for N50 for now. Output formation should look like the following: 

Dataset	Language 	Runtime 	N50
-------------------------------------------------------------------------------------------------------
data1	python		0:20:00		9118
data1	codon		0:10:00		9118
...

Context: Didn't know what N50 was 

Prompt: How to calculate N50? 

Prompt: Convert N50 function from python to bash 

Context: could set stack size at 8192000 since os limit are smaller (~64 MB), also asked a couple follow up questions after the original prompt because I was confused

Prompt: what does this mean: (base) oliviawhitelaw@Olivias-Air genome-assembly % ulimit -s 8192000
ulimit: value exceeds hard limit



