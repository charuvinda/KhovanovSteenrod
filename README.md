# KhovanovSteenrod
This is a tool to compute the function St_1, St_2, St_3, and so on, on the families of even and odd Khovanov stable homotopy types X_0(L), X_1(L), X_2(L), and so on (St_i gives information about the spectra X_i). Really, the function St_i only depends on i mod 4, since the first and second Steenrod squares on X_i onldy depend on i mod 4. If you only want the second Steenrod square on the spectra X_i(L), i = 0 mod 4, this program by Robert Lipshitz and Sucharit Sarkar also does the job: https://github.com/sucharit/KhovanovSteenrod. There are 3 files of interest:
F2Algebra.py
EvenKhovanov.py
Even Kh loader.py
GenericKhovanov.py
St loader.py

F2Algebra.py simply contains classes and functions needed to do linear algebra in the later python files, such as vector spaces, sparse matrices, reduced row echelon form, kernel, image, rank and nullity.

EvenKhovanov.py contains the formula for St_0, which contains the information on Sq^1 and Sq^2 on the even Khovanov spectrum X_0(L). 
<img src="https://latex.codecogs.com/png.latex?O_t=\text { Onset event at time bin } t " /> 
