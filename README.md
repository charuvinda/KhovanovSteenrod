# KhovanovSteenrod
This is a tool to compute the function St_1, St_2, St_3, and so on, on the families of even and odd Khovanov stable homotopy types X_0(L), X_1(L), X_2(L), and so on (St_i gives information about the spectra X_i). Really, the function St_i only depends on i mod 4, since the first and second Steenrod squares on X_i onldy depend on i mod 4. If you only want the second Steenrod square on the spectra X_i(L), i = 0 mod 4, this program by Robert Lipshitz and Sucharit Sarkar also does the job: https://github.com/sucharit/KhovanovSteenrod. There are 3 files of interest:
F2Algebra.py
EvenKhovanov.py
Even Kh loader.py
GenericKhovanov.py
St loader.py

F2Algebra.py simply contains classes and functions needed to do linear algebra in the later python files, such as vector spaces, sparse matrices, reduced row echelon form, kernel, image, rank and nullity.

EvenKhovanov.py computes the formula for St_0 for a MorseLink presentation of a link diagram $L$. This contains the information on Sq^1 and Sq^2 on the even Khovanov spectrum X_0(L). This program is simpler than GenericKhovanov.py, but it cannot compute St_i for other i.

GenericKhovanov computes St_i for a MorseLink presentation of a link diagram $L$, for any i. (Again, it only depends on i mod 4).

Even Kh loader.py computes St_0 for a list of MorseLink presentations. This is useful if you have multiple links you want to evaluate St on.

St loader.py computes St_i, for any fixed i, for a list of MorseLink presentations.


INSTRUCTIONS FOR USE:
CASE 1a. You would like to compute St_0 for a MorseLink presentation.
Step 1. Open EvenKhovanov.py
Step 2. 

CASE 1b. You would like to compute St_i for a MorseLink presentation.

CASE 2a. You would like to compute St_0 for a list of MorseLink presentations.
CASE 2b. You would like to compute St_i for a list of MorseLink presentations.

Remember that if 
1. r_1 is the rank of the map Sq^2 on gradings (i,j),
2. r_2 is the rank of Sq^2 restricted to the kernel of Sq^1 on gradings (i,j),
3. r_3 is the dimension of the intersection of the image of Sq^2 and image of Sq^1 in gradings (i+2,j),
4. r_4 is the dimension of the intersection of the image of Sq^1 and the image of (Sq^2 restricted to the kernel of Sq^1) in grading (i+2,j),
then St(i,j) = (r_2-r_4,r_1-r_2-r_3+r_4,r_4,r_3-r_4).
