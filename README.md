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
Step 2. Enter your MorseLink presentation into the prompt "linkDict = ..." towards the bottom of the code (instructions below on MorseLink presentations)
Step 3. Run the code. Your output should be a set of tuples of the form ((i,j),(a,b,c,d)). Each tuple like this means that St(L)(i,j) = (a,b,c,d)

CASE 1b. You would like to compute St_i for a MorseLink presentation.
Step 1. Open EvenKhovanov.py. 
Step 2. Replace the text "khovanovParameter = 0" with "khovanovParameter = i."
Step 3. Enter your MorseLink presentation into the prompt "linkDict = ..." towards the bottom of the code (instructions below on MorseLink presentations)
Step 4. Run the code. Your output should be a set of tuples of the form ((i,j),(a,b,c,d)). Each tuple like this means that St(L)(i,j) = (a,b,c,d)


CASE 2a. You would like to compute St_0 for a list of MorseLink presentations.
Step 1. Create a .txt file with the list of links you would like to compute. (There are some examples in the folder MorseLink tables for lists of n-crossing knots and links, and their mirrors.)
Step 2. (Optional, if you want faster runtime) If you know the (i,j) gradings you want to compute for each link, write them down in order, matching the order of the MorseLinks. This should be written as a Python set variable (For examplle, write {(0, 1), (-1, -1), (-2, -3)} in the nth row if you want to compute St(L) for these gradings (i,j) or set() if you want to compute St(L) for none of the gradings).
Step 3. Open Even Kh loader.py
Step 4. If you have not done Step 2,

CASE 2b. You would like to compute St_i for a list of MorseLink presentations.

MorseLink Presentations:
In this section, we give a recipe for how to make your own Morselink presentation you can input
We encode an oriented link diagram as a MorseLink, which, by our definition, is a sequence of Morse moves and crossings moving left to right. It is written in the form (Move_1,Move_2,Move_3,...,Move_n), where
1. Move_i is of the form ((a,b),X),
2. X is either 0, 1, '+', or '-'
3. a,b are distinct integers, indicating the two rows that involve this move. In the case X is not '0', a must be less than b
4. In the case that X is '0', then a may be greater than b. However, the orientation of the component containing this critical point is determined by the direction traveling from row a to row b.

X=0 indicates an index 0 critical point, or a birth. X=1 indicates an index 1 critical point, or a death. For example, the move ((2,3),0) would indicate the birth of two new strands, moving from row 2 to row 3.
X = '+' indicates that the strand starting at higher index will move over the strand starting at lower index during the crossing.
Likewise, X = '-' indicates that the strand starting at higher index will move over the strand starting at lower index during the crossing.
Here's an example of a MorseLink diagram of an oriented negative trefoil: (((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
<img width="1072" height="411" alt="image" src="https://github.com/user-attachments/assets/affa27a0-bcac-4f97-8b71-f373026d76c8" />

Here's an example of another MorseLink diagram for an oriented L6n1: (((1,4),0), ((6,5), 0), ((3,2), 0), ((4,5),'-'), ((1,2),'+'), ((3,4),'+'), ((2,3),'+'), ((1,2),'-'), ((3,4),'-'), ((2,3),1), ((4,5),1), ((1,6),1))
<img width="1245" height="455" alt="image" src="https://github.com/user-attachments/assets/a6172a96-443d-424d-b00c-fbb4ef968442" />
Note how the first three moves are index 0 critical points, which orient the three link components. We first orient going from row 1 to row 4, row 6 to row 5, and row 3 to row 2.  

Remember that if 
1. r_1 is the rank of the map Sq^2 on gradings (i,j),
2. r_2 is the rank of Sq^2 restricted to the kernel of Sq^1 on gradings (i,j),
3. r_3 is the dimension of the intersection of the image of Sq^2 and image of Sq^1 in gradings (i+2,j),
4. r_4 is the dimension of the intersection of the image of Sq^1 and the image of (Sq^2 restricted to the kernel of Sq^1) in grading (i+2,j),
then St(i,j) = (r_2-r_4,r_1-r_2-r_3+r_4,r_4,r_3-r_4).
