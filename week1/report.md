
## Steps Taken
### 1. Cloning the Repository

Getting the repo set up was straightforward. The only hiccup at the start was figuring out which values I was supposed to reproduce from the original implementation

### 2. Converting to Codon

The only tricky part was just tracking down what was actually causing compilation errors. 
ex. non comptible libraries, NoneType issues, etc 
After I isolated those issues, the changes themselves were small and Codon ran smoothly.

### 3. Writing the evaluate.sh Script

Probably the most confusing par since I’d never written a Bash script before. 
On top of that, I didn’t know how to calculate N50 in Bash at first, so I had to spend extra time figuring that part out.


## Gotchas: 

macOS Stack Limit

My biggest surprise was that macOS hard-caps the stack size at 64 MB. Every time I tried to run the recommended command:

ulimit -s 8192000


I just got back:

ulimit: value exceeds hard limit


At first this completely confused me — I thought I was doing something wrong. Eventually I realized it wasn’t a mistake on my end, it’s just a restriction of macOS. On Linux or GitHub CI runners this command works fine, but locally on macOS you’ll always hit that error.

Slightly Different Results for Data 1

Another weird thing I noticed: the N50 result for data1 didn’t match exactly between Python and Codon. The difference was small, and the other datasets matched perfectly. My guess is that it’s either randomness in the assembler or small floating-point differences between Python and Codon. Either way, it didn’t seem like a major problem, but it was worth noting.


## Bonus BLAST 

Used: https://blast.ncbi.nlm.nih.gov/Blast.cgi


