from EvenKhovanov import *
import ast

def KhDataPrint(MorseLinkFilePath):

  try:
    with open(MorseLinkFilePath, 'r') as file:
        MorseLinkContent = file.read()
        MorseLinkList = MorseLinkContent.splitlines()
    newContent = ""
    with open("data.txt", 'w') as file:
      pass
    for i in range(len(MorseLinkList)):
      linkDict = ast.literal_eval(MorseLinkList[i])
      linkMorseLink = morseLink(linkDict)
      link = khovanovBasis(linkMorseLink)
      hqToHomologyLoader(link,'all')
      validGradings = set()
      for grading in link.hqToHomology.keys():
        if link.hqToHomology[grading].sparseEchelon != dict():
          print(link.hqToHomology[grading] )
          validGradings.add(grading)
      newLine = str(tuple(validGradings))
      newContent += newLine + "\n"
      with open("data.txt", 'r+') as file:
        file.seek(0, 2)
        file.write(newLine + "\n")
  except FileNotFoundError:
    print(f"Error: File '{file_path}' not found.")

KhDataPrint('MorseLink L10 Table.txt')

"""
linkDict = (((1, 6), 0), ((7, 10), 0), ((6, 7), "+"), ((2, 3), 0), ((1, 2), "-"), ((1, 2), "-"), ((9, 8), 0), ((7, 8), "-"), ((8, 9), "+"), ((4, 5), 0), ((3, 4), "+"), ((3, 4), "+"), ((5, 6), "+"), ((5, 6), "+"), ((6, 7), "-"), ((5, 6), "+"), ((5, 4), 1), ((6, 7), 1), ((3, 8), 1), ((9, 2), 1), ((10, 1), 1))
linkMorseLink = morseLink(linkDict)
gradingSet = {(-2, -5), (-7, -15), (-4, -9), (-6, -13), (-3, -7), (-5, -11)}

x = ast.literal_eval('(((1, 6), 0), ((7, 10), 0), ((6, 7), "+"), ((2, 3), 0), ((1, 2), "-"), ((1, 2), "-"), ((9, 8), 0), ((7, 8), "-"), ((8, 9), "+"), ((4, 5), 0), ((3, 4), "+"), ((3, 4), "+"), ((5, 6), "+"), ((5, 6), "+"), ((6, 7), "-"), ((5, 6), "+"), ((5, 4), 1), ((6, 7), 1), ((3, 8), 1), ((9, 2), 1), ((10, 1), 1))')
print(x)
print(StList(linkMorseLink, gradingSet))
"""

