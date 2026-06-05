from GenericKhovanov import *
#from EvenKhovanov import *
#UNCOMMENT THE SECOND LINE AND COMMENT THE FIRST LINE IF YOU ONLY WANT TO COMPUTE THE ST FNCTION FOR THE EVEN HOMOTOPY TYPE X_0
import ast

def StDataPrint(MorseLinkFilePath, gradingsFilePath):
  try:
    with open('MorseLink tables/' + MorseLinkFilePath, 'r') as file:
        MorseLinkContent = file.read()
        MorseLinkList = MorseLinkContent.splitlines()
    with open('Grading list/' + gradingsFilePath, 'r') as file:
        gradingsContent = file.read()
        gradingsList = gradingsContent.splitlines()
    newContent = ""
    with open("data.txt", 'w') as file:
      pass
    for i in range(len(MorseLinkList)):
      linkDict = ast.literal_eval(MorseLinkList[i])
      linkMorseLink = morseLink(linkDict)
      print(ast.literal_eval(gradingsList[i]))
      newLine = str(StList(linkMorseLink, set(ast.literal_eval(gradingsList[i]))))
      newContent += newLine + "\n"
      with open("data.txt", 'r+') as file:
        file.seek(0, 2)
        file.write(newLine + "\n")
  except FileNotFoundError:
    print(f"Error: File '{file_path}' not found.")

def StDataNoGradingsPrint(MorseLinkFilePath):
  try:
    with open('MorseLink tables/' + MorseLinkFilePath, 'r') as file:
        MorseLinkContent = file.read()
        MorseLinkList = MorseLinkContent.splitlines()
    newContent = ""
    with open("data.txt", 'w') as file:
      pass
    for i in range(len(MorseLinkList)):
      linkDict = ast.literal_eval(MorseLinkList[i])
      linkMorseLink = morseLink(linkDict)
      newLine = str(StList(linkMorseLink, 'all'))
      newContent += newLine + "\n"
      with open("data.txt", 'r+') as file:
        file.seek(0, 2)
        file.write(newLine + "\n")
  except FileNotFoundError:
    print(f"Error: File '{file_path}' not found.")


StDataNoGradingsPrint('MorseLink 8 Table.txt')
#StDataPrint('MorseLink 10 Table.txt', 'Grading List 10.txt')

