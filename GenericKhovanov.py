from itertools import chain, combinations
from F2Algebra import *
khovanovParameter = 2
def powerSet(iterable):
    "powerset([1,2,3]) → () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = iterable
    #s = list(iterable)
    S = chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    return {frozenset(x) for x in S}

class morseLink:
    def __init__(self, moveTuple):#input (positiveCrossingNo, negativeCrossingNo)
        moveNo = len(moveTuple)
        self.moveTuple = moveTuple
        crossingTuple = ()
        crossingNos = set()
        for n in range(len(moveTuple)):
            if moveTuple[n][1] == '+' or moveTuple[n][1] == '-':
                crossingNos.add(n)
                crossingTuple = crossingTuple + (moveTuple[n][1],)
        self.crossingNos = crossingNos
        self.crossingTuple = crossingTuple
        indexToCrossingNo = list(crossingNos)
        indexToCrossingNo.sort()
        self.indexToCrossingNo = indexToCrossingNo
        noCrossingNos = len(indexToCrossingNo)
        self.indexToCrossingNo = indexToCrossingNo
        self.noCrossingNos = noCrossingNos
        #Now we find the number of components in the link
        cupGroup = set()
        runningDict = dict({})
        for n in range(moveNo):
            activeHeight = n
            activeTuple = self.moveTuple[activeHeight][0]
            activeMove = self.moveTuple[activeHeight][1]
            if activeMove == 0:
                cupGroup.add(frozenset({activeHeight}))
                for i in activeTuple:
                    runningDict[i] = activeHeight
            if activeMove == 1:
                setToBeAdded = frozenset()
                setsToBeRemoved = frozenset()
                for i in activeTuple:
                    curveComponentNo = runningDict[i]
                    for x in cupGroup:
                        if curveComponentNo in x:
                            setToBeAdded = setToBeAdded | x
                            setsToBeRemoved = setsToBeRemoved | frozenset({x})
                cupGroup = cupGroup - setsToBeRemoved | frozenset({setToBeAdded})
            if activeMove == '+' or activeMove == '-':
                previousRunningDictEntries = (runningDict[activeTuple[0]], runningDict[activeTuple[1]])
                for i in activeTuple:
                    runningDict[activeTuple[0]] = previousRunningDictEntries[1]
                    runningDict[activeTuple[1]] = previousRunningDictEntries[0]
        componentMins = {min(x) for x in cupGroup}
        self.componentNo = len(componentMins)
        #print('components:', self.componentNo)
        #Now we figure out how many crossings are positive and negative
        crossingDirectionsDict = dict([[x,[0,0]] for x in crossingNos])
        widths = set()
        for n in range(moveNo):
            for i in self.moveTuple[n][0]:
                widths.add(i)
        movesWidthI = dict()
        for i in widths:
            movesWidthI[i] = set()
        for n in range(moveNo):
            for i in moveTuple[n][0]:
                movesWidthI[i].add(n)
        #print('widths:', movesWidthI)
        for k in componentMins:
            #print('k:', k)
            currentComponent = k
            #print('startingComponent:', currentComponent)
            currentStrandWidth = moveTuple[k][0][0]
            nextStrandWidth = moveTuple[k][0][1]
            nextComponent = min(movesWidthI[nextStrandWidth] & set(range(k+1, moveNo + 1)))
            #print('nextComponent:', nextComponent)
            direction = 'up'
            counter = 0
            while nextComponent != k:
                currentComponent = nextComponent
                currentStrandWidth = nextStrandWidth
                activeTuple = moveTuple[currentComponent][0]
                #print(currentComponent, currentStrandWidth, activeTuple)
                if currentStrandWidth == activeTuple[0]:
                    nextStrandWidth = activeTuple[1]
                elif currentStrandWidth == activeTuple[1]:
                    nextStrandWidth = activeTuple[0]
                else:
                    print('help!')
                if moveTuple[currentComponent][1] == 0:
                    nextComponent = min(movesWidthI[nextStrandWidth] & set(range(currentComponent + 1, moveNo + 1)))
                    direction = 'up'
                elif moveTuple[currentComponent][1] == 1:
                    nextComponent = max(movesWidthI[nextStrandWidth] & set(range(currentComponent)))
                    direction = 'down'
                elif moveTuple[currentComponent][1] == '+' or moveTuple[currentComponent][1] == '-':
                    if direction  == 'up':
                        if currentStrandWidth == moveTuple[currentComponent][0][0]:
                            crossingDirectionsDict[currentComponent][0] = 'up'
                        elif currentStrandWidth == moveTuple[currentComponent][0][1]:
                            crossingDirectionsDict[currentComponent][1] = 'up'
                    elif direction == 'down':
                        if currentStrandWidth == moveTuple[currentComponent][0][0]:
                            crossingDirectionsDict[currentComponent][1] = 'down'
                        elif currentStrandWidth == moveTuple[currentComponent][0][1]:
                            crossingDirectionsDict[currentComponent][0] = 'down'
                    if direction == 'up':
                        nextComponent = min(movesWidthI[nextStrandWidth] & set(range(currentComponent + 1, moveNo + 1)))
                    elif direction == 'down':
                        nextComponent = max(movesWidthI[nextStrandWidth] & set(range(currentComponent)))
        positiveCrossingNo = 0
        negativeCrossingNo = 0
        for x in crossingNos:
            if moveTuple[x][1] == '+':
                if (crossingDirectionsDict[x][0] == 'up' and crossingDirectionsDict[x][1] == 'up') or (crossingDirectionsDict[x][0] == 'down' and crossingDirectionsDict[x][1] == 'down'):
                    negativeCrossingNo += 1
                else:
                    positiveCrossingNo += 1
            elif moveTuple[x][1] == '-':
                if (crossingDirectionsDict[x][0] == 'up' and crossingDirectionsDict[x][1] == 'up') or (crossingDirectionsDict[x][0] == 'down' and crossingDirectionsDict[x][1] == 'down'):
                    positiveCrossingNo += 1
                else:
                    negativeCrossingNo += 1
            else:
                print('help')
        self.positiveCrossingNo = positiveCrossingNo
        self.negativeCrossingNo = negativeCrossingNo
            
class resolvedMorseLink: #data is a dictionary where the entries are of the form (n,((i,j), move))
    def __init__(self, resolutionList, diagram):
        crossingTuple = diagram.crossingTuple
        noCrossingNos = diagram.noCrossingNos
        moveTuple = diagram.moveTuple
        self.moveTuple = moveTuple
        crossingNos = diagram.crossingNos
        moveNo = len(moveTuple)
        cupGroup = set()
        #cupGroup is the set of cups in a given component
        runningDict = dict({})
        crossingToComponentsDict = dict({})
        crossingToComponentsList = []
        crossingsPassed = 0
        for n in range(moveNo):
            activeHeight = n
            activeTuple = self.moveTuple[activeHeight][0]
            activeMove = self.moveTuple[activeHeight][1]
            if n in crossingNos:
                if resolutionList[crossingsPassed] == 'stay':
                    crossingToComponentsList.append((runningDict[activeTuple[0]], runningDict[activeTuple[1]]))
                if resolutionList[crossingsPassed] == '10':
                    setToBeAdded = frozenset()
                    setsToBeRemoved = frozenset()
                    for i in activeTuple:
                        curveComponentNo = runningDict[i]
                        for x in cupGroup:
                            if curveComponentNo in x:
                                setToBeAdded = setToBeAdded | x
                                setsToBeRemoved = setsToBeRemoved | frozenset({x})
                                cupGroup = cupGroup - setsToBeRemoved | frozenset({setToBeAdded})
                                belowHeight = runningDict[activeTuple[0]]
                                #ends the capping off part, starts the cup part
                    cupGroup.add(frozenset({activeHeight}))
                    for i in activeTuple:
                        runningDict[i] = activeHeight
                    crossingToComponentsList.append((belowHeight,activeHeight))
                crossingsPassed += 1
            elif activeMove == 0:
                cupGroup.add(frozenset({activeHeight}))
                for i in activeTuple:
                    runningDict[i] = activeHeight
            elif activeMove == 1:
                setToBeAdded = frozenset()
                setsToBeRemoved = frozenset()
                for i in activeTuple:
                    curveComponentNo = runningDict[i]
                    for x in cupGroup:
                        if curveComponentNo in x:
                            setToBeAdded = setToBeAdded | x
                            setsToBeRemoved = setsToBeRemoved | frozenset({x})
                cupGroup = cupGroup - setsToBeRemoved | frozenset({setToBeAdded})
        self.components = sorted((min(x) for x in cupGroup))
        #components are indexed by their minimum point (which must be index 0)
        self.componentNo = len(cupGroup)
            #This tuple takes the crossing index and outputs the lowest cup
        for n in range(noCrossingNos):
            componentList = []
            for x in cupGroup:
                if crossingToComponentsList[n][0] in x and crossingToComponentsList[n][1] in x:
                    componentList = [min(x)]
                    break
                elif crossingToComponentsList[n][0] in x:
                    componentList.insert(0, min(x))
                elif crossingToComponentsList[n][1] in x:
                    componentList.insert(1, min(x))
            crossingToComponentsList[n] = tuple(componentList)
        self.crossingToComponentsTuple = tuple(crossingToComponentsList)
        #print(moveTuple)
        #print(crossingToComponentsDict)
        #print(self.crossingToComponentsTuple)
        #print(crossingToComponentsDict)
        
def resolution(diagram, u):
    #u has to be in the power set of {0,1,...,n-1}, where n is the number of crossings
    crossingTuple = diagram.crossingTuple
    moveList = list(diagram.moveTuple)
    noCrossingNos = diagram.noCrossingNos
    indexToCrossingNo = diagram.indexToCrossingNo
    resolutionList = []
    for n in range(noCrossingNos):
        activeTuple = moveList[indexToCrossingNo[n]][0]
        if crossingTuple[n] == '-':
            if n in u:
                resolutionList.append('10')
            else:
                resolutionList.append('stay')
        elif crossingTuple[n] == '+':
            if n in u:
                resolutionList.append('stay')
            else:
                resolutionList.append('10')
    return resolvedMorseLink(resolutionList, diagram)

#crossings are always oriented so that the arrow either points to the right (sequentially) or down (in ascending order of crossings)
class xFaceAssignment: #returns the set of faces which have sign 1 (not 0)
    def __init__(self, diagram):
        noCrossingNos = diagram.noCrossingNos
        setOfUValues = powerSet(list(n for n in range(noCrossingNos)))
        crossingSet = {j for j in range(noCrossingNos)}
        faceSet = set()
        for u in setOfUValues:
            uResolution = resolution(diagram,u)
            for i in crossingSet - u:
                for j in crossingSet - u & set(range(i+1,noCrossingNos+1)): #looking at pairs i,j where j>i
                    w = u | {i,j}
                    uResolution = resolution(diagram,u)
                    wResolution = resolution(diagram,w)
                    if wResolution.componentNo == uResolution.componentNo + 2:
                        faceSet.add((u,i,j))
                    elif wResolution.componentNo == uResolution.componentNo:
                        if uResolution.crossingToComponentsTuple[i] == uResolution.crossingToComponentsTuple[j] and len(uResolution.crossingToComponentsTuple[i]) == 2:
                            faceSet.add((u,i,j)) #detects the face of type A that is two disjoint circles with the arrows pointing from the left circle to the right circle
                        elif uResolution.crossingToComponentsTuple[i] == uResolution.crossingToComponentsTuple[j] and len(uResolution.crossingToComponentsTuple[i]) == 1:
                            v = u | {i}
                            iCrossing = diagram.crossingTuple[i]
                            vResolution = resolution(diagram,v)
                            vResITuple = vResolution.crossingToComponentsTuple[i]
                            vResJTuple = vResolution.crossingToComponentsTuple[j]
                            if (iCrossing == '+' and vResITuple != vResJTuple) or (iCrossing == '-' and vResITuple == vResJTuple):
                                faceSet.add((u,i,j))
        self.faceSet = faceSet
        #print((frozenset({1,3}),2,5) in faceSet)
        for u in setOfUValues:
            for i in crossingSet - u:
                for j in crossingSet - u & set(range(i+1,noCrossingNos+1)): #looking at pairs i,j where j>i
                    for k in crossingSet - u & set(range(j+1,noCrossingNos+1)):
                        if (((u,i,j) in faceSet) + ((u,i,k) in faceSet) + ((u,j,k) in faceSet) + ((u | {k},i,j) in faceSet) + ((u | {j} ,i,k) in faceSet) + ((u | {i},j,k) in faceSet)) % 2 != 0:
                            print('faceError:', (u,i,j,k))
                            print('faceserrors:', ((u,i,j) in faceSet) , ((u,i,k) in faceSet) , ((u,j,k) in faceSet) , ((u | {k},i,j) in faceSet) , ((u | {j} ,i,k) in faceSet) , ((u | {i},j,k) in faceSet) )
        
    def __repr__(self):
        return "xFaceAssignment(%s)" %self.faceSet


class xEdgeAssignment:
    def __init__(self, diagram):
        faceAssignment = xFaceAssignment(diagram)
        noCrossingNos = diagram.noCrossingNos
        setOfUValues = powerSet(list(n for n in range(noCrossingNos)))
        faceSet = faceAssignment.faceSet
        edgeSet = set()
        crossingSet = {x for x in range(noCrossingNos)}
        for u in setOfUValues:
            for i in crossingSet - u: #u and u|{i} are the vertices of the edge
                squareCount = 0
                for k in u & set(range(i)): #for k<i in u we are going to remove {0,...,k} from u and u|{i}
                    u1 = u - set(range(k+1))
                    if (u1,k,i) in faceSet:
                        squareCount += 1
                if squareCount % 2 == 1:
                    edgeSet.add((u,i))
        self.edgeSet = edgeSet
    def __repr__(self):
        return "xEdgeAssignment(%s)" %self.edgeSet
     
class khovanovBasis:
    def __init__(self, diagram):
        noCrossingNos = diagram.noCrossingNos
        negativeCrossingNo = diagram.negativeCrossingNo
        positiveCrossingNo = diagram.positiveCrossingNo
        minHomGrading = -1 - negativeCrossingNo
        maxHomGrading = noCrossingNos - negativeCrossingNo
        minQuantumGrading = positiveCrossingNo - 2*negativeCrossingNo - 2*noCrossingNos
        #print('(minHomGrading,maxHomGrading):', minHomGrading,maxHomGrading)
        maxQuantumGrading = positiveCrossingNo - 2*negativeCrossingNo + noCrossingNos + 2*noCrossingNos
        edgeAssignment = xEdgeAssignment(diagram)
        self.edgeSet = edgeAssignment.edgeSet
        self.minHomGrading = minHomGrading
        self.maxHomGrading = maxHomGrading
        self.minQuantumGrading = minQuantumGrading
        self.maxQuantumGrading = maxQuantumGrading
        #for the signed part of Odd Khovanov homology. Some generators are going to be mapped to another generator via a minus sign
        differentialNegativeSignsDict = dict({})
        powerSetDict = tuple(powerSet(list((i for i in range(m)))) for m in range(2*noCrossingNos + 1))
        setOfUValues = powerSet(list(n for n in range(noCrossingNos))) #locations on the Kauffman cube
        indexToComponentDict = dict({}) #going to be a dictionary that takes a cube position u and outputs a sequencing of components
        componentToIndexDict = dict({})
        gradingsToGeneratorsDict = dict({}) #This dict will take in a quantum grading and homological grading and output the generators
        componentNoDict = dict({})
        hqToSparseMatrixSetup = dict({})
        crossingToComponentsTupleDict = dict({})
        
        qSet = set()
        if diagram.componentNo%2 == 1:
            for q in range(minQuantumGrading, maxQuantumGrading + 1):
                if q%2 == 1:
                    qSet.add(q)
        if diagram.componentNo%2 == 0:
            for q in range(minQuantumGrading, maxQuantumGrading + 1):
                if q%2 == 0:
                    qSet.add(q)
        self.qSet = qSet
        for q in qSet:
            for n in range(minHomGrading, maxHomGrading + 1):
                gradingsToGeneratorsDict[(n,q)] = set()
                hqToSparseMatrixSetup[(n,q)] = set()
        qGradingDict = dict({})
        resolutionDict = dict({})
        for u in setOfUValues:
            uResolution = resolution(diagram,u)
            resolutionDict[u] = uResolution
            crossingToComponentsTupleDict[u] = uResolution.crossingToComponentsTuple
            componentNoDict[u] = uResolution.componentNo
            componentTuple = uResolution.components
            indexToComponentDict[u] = dict({(k, componentTuple[k]) for k in range(uResolution.componentNo)})
            componentToIndexDict[u] = dict({(componentTuple[k], k) for k in range(uResolution.componentNo)})
            homGradingU = len(u) - negativeCrossingNo
            for S in powerSetDict[componentNoDict[u]]: #S represents the components that are labeled with x
                qGrading = positiveCrossingNo - 2 * negativeCrossingNo + len(u) + (uResolution.componentNo - len(S)) - (len(S))
                gradingsToGeneratorsDict[(homGradingU,qGrading)].add((u,S))
                qGradingDict[(u,S)] = qGrading
        self.crossingToComponentsTupleDict = crossingToComponentsTupleDict
        hqToIndexToGenerator = dict({})
        hqToGeneratorToIndex = dict({})
        hqToGeneratorNo = dict({})
        for q in qSet:
            for n in range(minHomGrading, maxHomGrading + 1):
                generatorTuple = tuple(gradingsToGeneratorsDict[(n, q)])
                hqToGeneratorNo[(n,q)] = len(generatorTuple)
                indexToGenerator = dict({(i, generatorTuple[i]) for i in range(hqToGeneratorNo[(n,q)])})
                generatorToIndex = dict({(generatorTuple[i], i) for i in range(hqToGeneratorNo[(n,q)])})
                hqToIndexToGenerator[(n,q)] = indexToGenerator
                hqToGeneratorToIndex[(n,q)] = generatorToIndex
                differentialNegativeSignsDict[(n,q)] = set()
        self.hqToGeneratorNo = hqToGeneratorNo
        self.hqToIndexToGenerator = hqToIndexToGenerator
        self.hqToGeneratorToIndex = hqToGeneratorToIndex
        self.componentToIndexDict = componentToIndexDict
        self.indexToComponentDict = indexToComponentDict
        crossingSet = {j for j in range(noCrossingNos)}
        #Now is the Khovanov differential part
        for u in setOfUValues:
            homGradingU = len(u) - negativeCrossingNo
            for i in crossingSet - u:
                v = u | {i}
                homGradingV = homGradingU + 1
                uSpecialComponents = crossingToComponentsTupleDict[u][i]
                vSpecialComponents = crossingToComponentsTupleDict[v][i]
                if componentNoDict[u] < componentNoDict[v]: #if the differential is along a split
                    for uS in powerSetDict[componentNoDict[u]]:
                        uLabeling = {indexToComponentDict[u][s] for s in uS}
                        vLabelingOfOtherComponents = uLabeling - set(uSpecialComponents)
                        qGradingU = positiveCrossingNo - 2 * negativeCrossingNo + len(u) + (componentNoDict[u] - len(uS)) - (len(uS))
                        qGradingV = qGradingU
                        jCoordinate = hqToGeneratorToIndex[(homGradingU,qGradingU)][(u,uS)]
                        if diagram.crossingTuple[i] == '+':
                            a1 = vSpecialComponents[0]
                            a2 = vSpecialComponents[1]
                        elif diagram.crossingTuple[i] == '-':
                            a1 = vSpecialComponents[1]
                            a2 = vSpecialComponents[0]
                        if min(uSpecialComponents) in uLabeling: #if the circle that is being split is labeled x
                            generatorComponentForm = vLabelingOfOtherComponents | set(vSpecialComponents)
                            vS = frozenset(componentToIndexDict[v][x] for x in generatorComponentForm)
                            iCoordinate = hqToGeneratorToIndex[(homGradingV,qGradingV)][(v, vS)]
                            hqToSparseMatrixSetup[(homGradingU, qGradingU)].add((iCoordinate, jCoordinate))
                            if a2<a1:
                                if len(uLabeling & set(range(a1))) % 2 == 1:
                                    differentialNegativeSignsDict[(homGradingU,qGradingU)].add((jCoordinate,iCoordinate))
                            elif a1<a2:
                                if len(uLabeling & set(range(a2))) % 2 == 0:
                                    differentialNegativeSignsDict[(homGradingU,qGradingU)].add((jCoordinate,iCoordinate))
                        else: #if the circle that is being split is labled 1
                            generatorsComponentForm = (vLabelingOfOtherComponents | {a1}, vLabelingOfOtherComponents | {a2})
                            vS1 = frozenset(componentToIndexDict[v][x] for x in generatorsComponentForm[0])
                            vS2 = frozenset(componentToIndexDict[v][x] for x in generatorsComponentForm[1])
                            iCoordinate1 = hqToGeneratorToIndex[(homGradingV, qGradingV)][(v, vS1)]
                            iCoordinate2 = hqToGeneratorToIndex[(homGradingV, qGradingV)][(v, vS2)]
                            hqToSparseMatrixSetup[(homGradingU, qGradingU)].update(((iCoordinate1, jCoordinate,), (iCoordinate2, jCoordinate)))
                            #we move the a1 circle up all the way to its position, and then we move the -a2 circle up all the way to its position
                            if len(uLabeling & set(range(a1))) % 2 == 1:
                                differentialNegativeSignsDict[(homGradingU,qGradingU)].add((jCoordinate,iCoordinate1))
                            if len(uLabeling & set(range(a2))) % 2 == 0:
                                differentialNegativeSignsDict[(homGradingU,qGradingU)].add((jCoordinate,iCoordinate2))                      
                if componentNoDict[u] > componentNoDict[v]: #if the differential is along a merge
                    for uS in powerSetDict[componentNoDict[u]]:
                        uLabeling = {indexToComponentDict[u][s] for s in uS}
                        vLabelingOfOtherComponents = uLabeling - set(uSpecialComponents)
                        qGradingU = positiveCrossingNo - 2 * negativeCrossingNo + len(u) + (componentNoDict[u] - len(uS)) - (len(uS))
                        qGradingV = qGradingU
                        jCoordinate = hqToGeneratorToIndex[(homGradingU, qGradingU)][(u,uS)]
                        a1 = uSpecialComponents[0]
                        a2 = uSpecialComponents[1]
                        if (a1 in uLabeling) ^ (a2 in uLabeling):
                            generatorComponentForm = vLabelingOfOtherComponents | set(vSpecialComponents)
                            vS = frozenset(componentToIndexDict[v][x] for x in generatorComponentForm)
                            iCoordinate = hqToGeneratorToIndex[(homGradingV, qGradingV)][(v, vS)]
                            hqToSparseMatrixSetup[(homGradingU, qGradingU)].add((iCoordinate,jCoordinate))
                            if a1 in uLabeling and a2 < a1: #if the merge decreases the lowest point of the circle, the label a1 should "move down," at the expense switching signs every time you move past a circle
                                if len(set(range(a2 + 1, a1)) & uLabeling) % 2 == 1:
                                    differentialNegativeSignsDict[(homGradingU,qGradingU)].add((jCoordinate,iCoordinate))
                            if a2 in uLabeling and a1 < a2: #similar story here
                                if len(set(range(a1 + 1, a2)) & uLabeling) % 2 == 1:
                                    differentialNegativeSignsDict[(homGradingU,qGradingU)].add((jCoordinate,iCoordinate))
                        elif (a1 not in uLabeling) and (a2 not in uLabeling):
                            generatorComponentForm = vLabelingOfOtherComponents
                            vS = frozenset(componentToIndexDict[v][x] for x in generatorComponentForm)
                            iCoordinate = hqToGeneratorToIndex[(homGradingV, qGradingV)][(v, vS)]
                            hqToSparseMatrixSetup[(homGradingU, qGradingU)].add((iCoordinate,jCoordinate))
        self.differentialNegativeSignsDict = differentialNegativeSignsDict
        #print(differentialNegativeSignsDict)
        hqToSparseMatrix = dict()
        for q in qSet:
            for n in range(minHomGrading, maxHomGrading + 1):
                sparseMatrix = F2sparseMatrix(hqToSparseMatrixSetup[(n,q)])
                sparseMatrix.officialColumnRange = (0, hqToGeneratorNo[(n,q)] - 1)
                hqToSparseMatrix[(n,q)] = sparseMatrix
        self.hqToSparseMatrix = hqToSparseMatrix
        self.hqToHomology = dict()
def hqToHomologyLoader(khovanovGenerators, gradings):
    qSet = khovanovGenerators.qSet
    for n,q in gradings:
        if (n,q) not in set(khovanovGenerators.hqToHomology.keys()):
            khovanovGenerators.hqToHomology[(n,q)] = homology(khovanovGenerators.hqToSparseMatrix[(n,q)], khovanovGenerators.hqToSparseMatrix[(n-1,q)])
                
class flowCategory:
    def __init__(self, morseLink, khovanovGenerators, gradings):
        qGrading = gradings[1]
        bottomHGrading = gradings[0]
        componentToIndexDict = khovanovGenerators.componentToIndexDict
        crossingToComponentsTupleDict = khovanovGenerators.crossingToComponentsTupleDict
        crossingTuple = morseLink.crossingTuple
        moveTuple = morseLink.moveTuple
        indexToComponentDict = khovanovGenerators.indexToComponentDict
        componentToIndexDict = khovanovGenerators.componentToIndexDict
        indexToGenerator = khovanovGenerators.hqToIndexToGenerator[(bottomHGrading + 1, qGrading)]
        generatorToIndex = khovanovGenerators.hqToGeneratorToIndex[(bottomHGrading + 1, qGrading)]
        nGeneratorDict = khovanovGenerators.hqToIndexToGenerator[(bottomHGrading, qGrading)]
        nPlusOneGeneratorDict = khovanovGenerators.hqToIndexToGenerator[(bottomHGrading + 1, qGrading)]
        nPlusTwoGeneratorDict = khovanovGenerators.hqToIndexToGenerator[(bottomHGrading + 2, qGrading)]
        self.nGeneratorDict = nGeneratorDict 
        self.nPlusOneGeneratorDict = nPlusOneGeneratorDict
        self.nPlusTwoGeneratorDict = nPlusTwoGeneratorDict
        nMinusOneMatrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading - 1, qGrading)]
        nMatrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading, qGrading)]
        nPlusOneMatrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading + 1, qGrading)]
        nPlusTwoMatrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading + 2, qGrading)]
        pontrjaginThomBoundaries = dict({})
        pontrjaginThoms = dict({})
        boundariesY = dict({})
        pontrjaginThomAdjacency = dict({})
        zToXyPairs = dict({})
        for x in nMatrix.columns:
            for y in nMatrix.entriesInColumn[nMatrix.columnDict[x]]:
                if y in nPlusOneMatrix.columns:
                    for z in nPlusOneMatrix.entriesInColumn[nPlusOneMatrix.columnDict[y]]:
                        yzDifference = min(nPlusTwoGeneratorDict[z][0]-nPlusOneGeneratorDict[y][0])
                        pontrjaginThomBoundaries[(x,z)] = set()
                        boundariesY[y] = set()
                        zToXyPairs[z] = set()
                        pontrjaginThomAdjacency[z] = dict({})
                        
        for x in nMatrix.columns:
            for y in nMatrix.entriesInColumn[nMatrix.columnDict[x]]:
                if y in nPlusOneMatrix.columns:
                    for z in nPlusOneMatrix.entriesInColumn[nPlusOneMatrix.columnDict[y]]:
                        pontrjaginThomBoundaries[(x,z)].add(y)
                        boundariesY[y].add(x)
                        zToXyPairs[z].add((x,y))
        self.pontrjaginThomBoundaries = pontrjaginThomBoundaries
        for xz in pontrjaginThomBoundaries.keys():
            x = xz[0]
            if len(pontrjaginThomBoundaries[xz]) == 2:
                pontrjaginThoms[xz] = (tuple(pontrjaginThomBoundaries[xz]),)
            elif len(pontrjaginThomBoundaries[xz]) == 4:
                
                fourTuple = tuple(pontrjaginThomBoundaries[xz])
                fourTupleGenerators = tuple(nPlusOneGeneratorDict[fourTuple[i]] for i in range(4))
                ABCrossings = tuple({min(fourTupleGenerators[i][0] - nGeneratorDict[x][0]) for i in range(4)}) #is a tuple that has the crossing 2 indices that get from where x is on the cube to where z is
                #print(ABCrossings)
                crossingAChange = ABCrossings[0]
                crossingBChange = ABCrossings[1]
                zeroAtVAAndVB = nGeneratorDict[x][0]
                oneAtVA = nGeneratorDict[x][0] | {crossingAChange} #The "u" of the (u,uS) paring that is a generator
                oneAtVB = nGeneratorDict[x][0] | {crossingBChange}
                oneAtVAComponents = crossingToComponentsTupleDict[oneAtVA][crossingAChange] #if +,left, then right at v1
                oneAtVBComponents = crossingToComponentsTupleDict[oneAtVB][crossingAChange] #if +, bottom, then top at v1
                oneAtVAComponentIndices = (componentToIndexDict[oneAtVA][oneAtVAComponents[0]], componentToIndexDict[oneAtVA][oneAtVAComponents[1]])
                oneAtVBComponentIndices = (componentToIndexDict[oneAtVB][oneAtVBComponents[0]], componentToIndexDict[oneAtVB][oneAtVBComponents[1]]) 
                u = nGeneratorDict[x][0]
                uS = nGeneratorDict[x][1]
                uLabeling = {indexToComponentDict[u][s] for s in uS} #set of positions of components labeled 'x'
                labelingWithoutLadybug = nGeneratorDict[x][1]
                #Because in a ladybug matching, the ladybug component at the 00 resoultion is always labeled with a 1 
                labelingWithoutLadybugComponentForm = frozenset({indexToComponentDict[u][i] for i in labelingWithoutLadybug})
                #positions in diagram of components of diagram labeled 'x'. Note that the ladybug isn't labeled 'x' in the first place or else it wouldn't be a ladybug matching. So no need to exclude the ladybug
                labelingWithoutLadybugIndexForm1A0B = frozenset({componentToIndexDict[oneAtVA][s] for s in labelingWithoutLadybugComponentForm})
                #indices of components labeled 'x' in the diagram where the 'A' vertex is resolved, with the ladybug components not included
                labelingWithoutLadybugIndexForm0A1B = frozenset({componentToIndexDict[oneAtVB][s] for s in labelingWithoutLadybugComponentForm})
                #indices of components labeled 'x' in the diagram where the 'B' vertex is resolved, with the ladybug components not included
                labeling1A0BIndexForm = (labelingWithoutLadybugIndexForm1A0B | frozenset({oneAtVAComponentIndices[0]}), labelingWithoutLadybugIndexForm1A0B | frozenset({oneAtVAComponentIndices[1]})) #labeling for when V1 is 0 and V2 is 1
                #A pair where each is an intermediate diagram. In the first, the bottom (or left) circle to the 'A' crossing is labeled 'x' and in the second, the top (or right) circle to the 'A' crossing is lableled 'x'
                labeling0A1BIndexForm = (labelingWithoutLadybugIndexForm0A1B | frozenset({oneAtVBComponentIndices[0]}), labelingWithoutLadybugIndexForm0A1B | frozenset({oneAtVBComponentIndices[1]}))
                #Similar to the above
                if crossingTuple[crossingAChange] == '+':
                    #0-resolution at vertex1 looks like cup and cap
                    pair1 = (generatorToIndex[(oneAtVA,labeling1A0BIndexForm[0])],generatorToIndex[(oneAtVB,labeling0A1BIndexForm[0])])
                    pair2 = (generatorToIndex[(oneAtVA,labeling1A0BIndexForm[1])],generatorToIndex[(oneAtVB,labeling0A1BIndexForm[1])])
                    
                elif crossingTuple[crossingAChange] == '-':
                    #0-resoultion at vertex1 looks like a continuation
                    pair1 = (generatorToIndex[(oneAtVA,labeling1A0BIndexForm[0])],generatorToIndex[(oneAtVB,labeling0A1BIndexForm[1])])
                    pair2 = (generatorToIndex[(oneAtVA,labeling1A0BIndexForm[1])],generatorToIndex[(oneAtVB,labeling0A1BIndexForm[0])])
                else:
                    raise ValueError('crossings!')
                pontrjaginThoms[xz] = (pair1, pair2)
                #print(pontrjaginThoms[xz])
                
            else: print('Error!')
        self.pontrjaginThoms = pontrjaginThoms
        self.boundariesY = boundariesY
        self.zToXyPairs = zToXyPairs
        #print('zToXyPairs:', zToXyPairs)
        for xz in pontrjaginThoms.keys():
            for yPair in pontrjaginThoms[xz]:
                pontrjaginThomAdjacency[xz[1]][(xz[0], yPair[0])] = (xz[0], yPair[1])
                pontrjaginThomAdjacency[xz[1]][(xz[0], yPair[1])] = (xz[0], yPair[0])
        self.pontrjaginThomAdjacency = pontrjaginThomAdjacency
        #print(pontrjaginThomAdjacency[1][(3,4)])
        #print('ladybugNumbers:', ladybug)
        
def Sq1(khovanovGenerators, gradings, kernel):
    qGrading = gradings[1]
    bottomHGrading = gradings[0]
    nMatrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading, qGrading)]
    nIndexToGenerator = khovanovGenerators.hqToIndexToGenerator[(bottomHGrading, qGrading)]
    nPlus1IndexToGenerator = khovanovGenerators.hqToIndexToGenerator[(bottomHGrading + 1, qGrading)]
    edgeAssignmentSet = khovanovGenerators.edgeSet
    oddArrowXYSet = khovanovGenerators.differentialNegativeSignsDict[gradings]
    minusSignMatrixSet = set()
    plusSignMatrixSet = set()
    for ij in nMatrix.sparseMatrix:
        i = ij[0]
        j = ij[1]
        u = nIndexToGenerator[j][0]
        v = nPlus1IndexToGenerator[i][0]
        differencePoint = min(v-u)
        #The first thing being added is the usual Kh sign infusement that makes a cube that commutes on 2D faces now anticommute
        #The second thing being added is the edge assignment that switches the signs of the A and X faces
        #The third thing being added is not an edge assignment, but the coefficients of the chosen Odd Khovanov generator i in the image of the Odd Khovanov generator j
        
        if khovanovParameter % 2 == 0:
            totalDifferentialSign = (len(set(u) & set(range(differencePoint))))%2
        if khovanovParameter % 2 == 1:
            totalDifferentialSign = (len(set(u) & set(range(differencePoint))) + ((u,differencePoint) in edgeAssignmentSet) + ((j,i) in oddArrowXYSet))%2
        if totalDifferentialSign == 1:
            minusSignMatrixSet.add(ij)
        elif totalDifferentialSign == 0:
            plusSignMatrixSet.add(ij)
    minusSignMatrix = F2sparseMatrix(minusSignMatrixSet)
    plusSignMatrix = F2sparseMatrix(plusSignMatrixSet)
    #print(minusSignMatrix+plusSignMatrix+nMatrix)
    output = set()
    for i in nMatrix.rows:
        if i in minusSignMatrix.rows:
            #Takes the dot product of the minus sign matrix with kernel
            minusMatrixDotProduct = len(minusSignMatrix.entriesInRow[minusSignMatrix.rowDict[i]] & kernel) % 4
        else:
            minusMatrixDotProduct = 0
        if i in plusSignMatrix.rows:
            #Takes the dot product of the plus sign matrix with kernel
            plusMatrixDotProduct = len(plusSignMatrix.entriesInRow[plusSignMatrix.rowDict[i]] & kernel) % 4
        else:
            plusMatrixDotProduct = 0
        #print(minusMatrixDotProduct + plusMatrixDotProduct)
        #print(plusMatrixDotProduct - minusMatrixDotProduct)
        if (plusMatrixDotProduct - minusMatrixDotProduct) % 4 == 2:
            output.add(i)
        elif (plusMatrixDotProduct - minusMatrixDotProduct) % 2 == 1:
            raise ValueError('Bockstein not working!')
    #print('Sq1:', gradings, kernel, khovanovGenerators.hqToSparseMatrix[(gradings[0]+1,gradings[1])]*output)
    return output
    
class matchings:
    def __init__(self, flowCategory, kernel):
        pontrjaginThomBoundaries = flowCategory.pontrjaginThomBoundaries
        pontrjaginThoms = flowCategory.pontrjaginThoms
        boundariesY = flowCategory.boundariesY
        activeBoundariesY = dict({})
        boundaryMatchingAdjacency = dict({})
        zToActiveXyPairs = dict({})
        for z in flowCategory.zToXyPairs.keys():
            zToActiveXyPairs[z] = set()
            for xy in flowCategory.zToXyPairs[z]:
                if xy[0] in kernel:
                    zToActiveXyPairs[z].add(xy)
        self.zToActiveXyPairs = zToActiveXyPairs
        #print('zToActiveXyPairs:', zToActiveXyPairs)
        for y in boundariesY.keys():
            if (boundariesY[y] & kernel != set()):
                activeBoundariesY[y] = boundariesY[y] & kernel
                #print(activeBoundariesY[y])
        #print(activeBoundaries)
        boundaryMatchings = dict({})
        for y in activeBoundariesY.keys():
            yBoundaryTuple = tuple(activeBoundariesY[y])
            #print(yBoundaryTuple)
            yBoundaryNo = len(yBoundaryTuple)
            halfNo = yBoundaryNo // 2
            subTuple1 = tuple(yBoundaryTuple[n] for n in range(halfNo))
            subTuple2 = tuple(yBoundaryTuple[n] for n in range(halfNo, yBoundaryNo))
            boundaryMatchings[y] = {(subTuple1[i], subTuple2[i]) for i in range(halfNo)}
        self.boundaryMatchings = boundaryMatchings
        for y in boundaryMatchings.keys():
            for xPair in boundaryMatchings[y]:
                boundaryMatchingAdjacency[(xPair[0],y)] = (xPair[1],y)
                boundaryMatchingAdjacency[(xPair[1],y)] = (xPair[0],y)
        self.boundaryMatchingAdjacency = boundaryMatchingAdjacency
        #print('boundarymatchings:', boundaryMatchings)
        #print(self.boundaryMatchings)

def zCoefficient(khovanovGenerators, gradings, signedFlowCategory, matchings, z):
    hGrading = gradings[0]
    qGrading = gradings[1]
    #takes the position on the cube of the generator z
    zPosition = khovanovGenerators.hqToIndexToGenerator[(hGrading + 2, qGrading)][z][0]
    pontrjaginThomAdjacency = signedFlowCategory.pontrjaginThomAdjacency
    boundaryAdjacencies = matchings.boundaryMatchingAdjacency
    edgeAssignmentSet = khovanovGenerators.edgeSet
    oddArrowYZSet = khovanovGenerators.differentialNegativeSignsDict[(hGrading+1,qGrading)]
    oddArrowXYSet = khovanovGenerators.differentialNegativeSignsDict[(hGrading,qGrading)]
    #print('activeXyPairs:', matchings.zToActiveXyPairs)
    xyPairs = matchings.zToActiveXyPairs[z]
    leftOverXyPairList = list(xyPairs)
    
    cycleSet = set()
    while leftOverXyPairList != []:
        boundaryAdjacency1 = leftOverXyPairList[0]
        boundaryAdjacency2 = boundaryAdjacencies[boundaryAdjacency1]
        leftOverXyPairList.remove(boundaryAdjacency1)
        leftOverXyPairList.remove(boundaryAdjacency2)
        cycleToAdd = [(boundaryAdjacency1,boundaryAdjacency2)]
        x = pontrjaginThomAdjacency[z][boundaryAdjacency2]
        while pontrjaginThomAdjacency[z][boundaryAdjacency2] != (cycleToAdd[0][0]):
            boundaryAdjacency1 = pontrjaginThomAdjacency[z][boundaryAdjacency2]
            boundaryAdjacency2 = boundaryAdjacencies[boundaryAdjacency1]
            cycleToAdd.append((boundaryAdjacency1,boundaryAdjacency2))
            leftOverXyPairList.remove(boundaryAdjacency1)
            leftOverXyPairList.remove(boundaryAdjacency2)
        cycleSet.add(tuple(cycleToAdd))

    Z2Element = 0
    for K in cycleSet:
        #print('K:', K)
        yCycle = tuple(K[i][0][1] for i in range(len(K)))
        yPositions = tuple(signedFlowCategory.nPlusOneGeneratorDict[yCycle[i]][0] for i in range(len(K))) #actual coordinates of the y's on the cube
        #print(yPositions)
        c = tuple(len(yPositions[i] & set(range(min(zPosition-yPositions[i])+1))) for i in range(len(K)))
        #c = tuple(yzPositionToYNoDict[(K[i][0][1], z)] for i in range(len(K)))
        #print('c:', c)
        cLength = len(c)
        KElement = 0
        switchBackNo = 0
        doubleSignedTwists = 0
        for i in range(cLength): #Adds the products of adjacencies
            KElement += (c[i%cLength] * c[(i+1)%cLength])
            KElement += max(c[i%cLength], c[(i+1)%cLength])
            KElement += c[i]#Adds each vertex
            triple = (c[i], c[(i+1)%cLength], c[(i+2)%cLength])
            #Odd Kh part: accounts for flipping
            y = K[(i+1)%cLength][0][1]
            selectedY = yPositions[(i+1)%cLength]
            yP1 = K[(i+2)%cLength][0][1]
            selectedYP1 = yPositions[(i+2)%cLength]
            selectedXPair = (K[(i+1)%cLength][0][0],K[(i+1)%cLength][1][0])
            selectedX = K[(i+1)%cLength][1][0]
            selectedXP1Pair = (K[(i+2)%cLength][0][0],K[(i+2)%cLength][1][0])
            x0Position = signedFlowCategory.nGeneratorDict[selectedXPair[0]][0]
            x1Position = signedFlowCategory.nGeneratorDict[selectedXPair[1]][0]
            selectedXPosition = signedFlowCategory.nGeneratorDict[selectedX][0]
            #print(selectedXPairPositions)
            bottomTentEdges = ((x0Position ,min(selectedY - x0Position)),(x1Position,min(selectedY - x1Position)))
            #print(selectedY, selectedXPairPositions)
            selectedEdge = (selectedY, min(zPosition - selectedY))
            selectedEdgeP1 = (selectedYP1, min(zPosition - selectedYP1))
            #cubeEdgeSum = (bottomDiamondCubeEdges[0] in edgeAssignmentSet) + (bottomDiamondCubeEdges[1] in edgeAssignmentSet) + (selectedEdge in edgeAssignmentSet) + (selectedEdgeP1 in edgeAssignmentSet)
            #alternatingAlgebraSum = ((selectedX,y) in oddArrowXYSet) + ((selectedX,yP1) in oddArrowXYSet) + ((y,z) in oddArrowYZSet) + ((yP1,z) in oddArrowYZSet)
            #diamondSign = cubeEdgeSum + alternatingAlgebraSum
            #print('diamondSign:', diamondSign, selectedXPosition, zPosition, cubeEdgeSum)
            #print(selectedEdge)
            #print(K[(i+1)%cLength])
            xPair = (K[(i+1)%cLength][0][0],K[(i+1)%cLength][1][0])
            selectedXYSignDiff = ((bottomTentEdges[0] in edgeAssignmentSet) + (bottomTentEdges[1] in edgeAssignmentSet) + ((xPair[0],y) in oddArrowXYSet) + ((xPair[1],y) in oddArrowXYSet)) % 2
            #print((selectedXYEdges[0] in edgeAssignmentSet) , (selectedXYEdges[1] in edgeAssignmentSet) , ((xPair[0],y) in oddArrowXYSet) , ((xPair[1],y) in oddArrowXYSet))
            #print(selectedXYEdges[0])
            #print(selectedEdge)
            #print(selectedXYSignDiff)
            if selectedXYSignDiff == 1:
                if triple[0] < triple[2]:
                    doubleSignedTwists += 1
                elif triple[0] > triple[2]:
                    doubleSignedTwists -= 1
                else:
                    #print('what')
                    middleBoundaryArc = K[(i+1)%cLength]
                    if (middleBoundaryArc[0][0], middleBoundaryArc[1][0]) not in matchings.boundaryMatchings[middleBoundaryArc[0][1]]:
                    #print(i, 'reversed')
                        doubleSignedTwists -= 1
                    elif (middleBoundaryArc[0][0], middleBoundaryArc[1][0]) in matchings.boundaryMatchings[middleBoundaryArc[0][1]]:
                        doubleSignedTwists += 1
            if sorted(triple)[1] == triple[1]:
                if (triple[2]-triple[0]) % 2 == 1 and (((selectedEdge in edgeAssignmentSet) + ((y,z) in oddArrowYZSet)) % 2 == 1):
                    if khovanovParameter % 2 == 1:
                        KElement += 1 #adds 1 if we are in the odd or second odd realization
            else:
                if (triple[2]-triple[0]) % 2 == 0 and (((selectedEdge in edgeAssignmentSet) + ((y,z) in oddArrowYZSet)) % 2 == 1):
                    if khovanovParameter % 2 == 1:
                        KElement += 1 #adds 1 if we are in the odd or second odd realization
            KElement += sorted(triple)[1]#Adds the middle of the value set of  each triple
            if sorted(triple)[1] != triple[1]:#Counts the number of times the path turns around (should be even)
                switchBackNo +=  1
            if triple[0] == triple[2]:
                if triple[1] < triple[0]:
                    KElement += 1
                elif triple[1] > triple[2]:
                    KElement += 0
            elif sorted(triple)[1] == triple[0]:
                KElement += 1
            if triple[0] == triple[2]:
                #print('hello!')
                middleBoundaryArc = K[(i+1)%cLength]
                if (middleBoundaryArc[0][0], middleBoundaryArc[1][0]) not in matchings.boundaryMatchings[middleBoundaryArc[0][1]]:
                    #print(i, 'reversed')
                    KElement += 1
                elif (middleBoundaryArc[0][0], middleBoundaryArc[1][0]) in matchings.boundaryMatchings[middleBoundaryArc[0][1]]:
                    KElement += 0
                else:
                    print('help!')
        #print('doubleSignedTwists:', doubleSignedTwists)
        if  (khovanovParameter % 4) in {2,3}:
            KElement += (doubleSignedTwists // 2)
        #print(switchBackNo)
        KElement += (switchBackNo // 2) #Adds half the number of times the path turns around
        KElement += 1 #adds 1 at the end
        Z2Element += KElement
    return Z2Element%2

def Sq2(khovanovGenerators, gradings, signedFlowCategory, kernel):
    boundaryMatchings = matchings(signedFlowCategory, kernel)
    Sq2Set = set({})
    for z in boundaryMatchings.zToActiveXyPairs.keys():
        #print('zValue:', z, '.')
        #print(matchings.zToActiveXyPairs[z])
        if zCoefficient(khovanovGenerators, gradings, signedFlowCategory, boundaryMatchings, z) == 1:
            Sq2Set.add(z)
    #print('Sq2:', gradings, kernel, khovanovGenerators.hqToSparseMatrix[(gradings[0]+2,gradings[1])]*Sq2Set)
    #print('check:', gradings, khovanovGenerators.hqToSparseMatrix[(gradings[0]+2,gradings[1])]*Sq2Set)
    if khovanovGenerators.hqToSparseMatrix[(gradings[0]+2,gradings[1])]*Sq2Set != set():
        print('help!')
    return Sq2Set

def Sq1Matrix(khovanovGenerators, gradings):
    bottomHGrading = gradings[0]
    qGrading = gradings[1]
    nMinus1Matrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading - 1, qGrading)]
    nMatrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading, qGrading)]
    nPlus1Matrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading + 1, qGrading)]
    hqToHomologyLoader(khovanovGenerators, {(bottomHGrading, qGrading), (bottomHGrading + 1, qGrading)})

    nHomology = khovanovGenerators.hqToHomology[(bottomHGrading, qGrading)]
    nPlus1Homology = khovanovGenerators.hqToHomology[(bottomHGrading + 1, qGrading)]
    B_N = image(nMatrix)
    nHomologyDim = len(nHomology.sortedPivots)
    Sq1SparseMatrix = set()
    for j in range(nHomologyDim):
        inputElement = nHomology.orderedSparseEchelonSet[j]
        outputElement = Sq1(khovanovGenerators, gradings, inputElement)
        column = homologyProjection(nPlus1Homology, B_N, outputElement)
        for i in column:
            Sq1SparseMatrix.add((i,j))
    Sq1SparseMatrix = F2sparseMatrix(Sq1SparseMatrix)
    Sq1SparseMatrix.officialColumnRange = (0,nHomologyDim - 1)
    return Sq1SparseMatrix

def Sq2Matrix(khovanovGenerators, flowCategory, gradings):
    bottomHGrading = gradings[0]
    qGrading = gradings[1]
    nMinus1Matrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading - 1, qGrading)]
    nMatrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading, qGrading)]
    nPlus1Matrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading + 1, qGrading)]
    nPlus2Matrix = khovanovGenerators.hqToSparseMatrix[(bottomHGrading + 2, qGrading)]
    hqToHomologyLoader(khovanovGenerators, {(bottomHGrading, qGrading), (bottomHGrading + 2, qGrading)})
    nHomology = khovanovGenerators.hqToHomology[(bottomHGrading, qGrading)]
    nPlus2Homology = khovanovGenerators.hqToHomology[(bottomHGrading + 2, qGrading)]
    B_NPlus2 = image(nPlus1Matrix)
    nHomologyDim = len(nHomology.sortedPivots)
    Sq2SparseMatrix = set()
    for j in range(nHomologyDim):
        inputElement = nHomology.orderedSparseEchelonSet[j]
        outputElement = Sq2(khovanovGenerators, gradings, flowCategory, inputElement)
        column = homologyProjection(nPlus2Homology, B_NPlus2, outputElement)
        for i in column:
            Sq2SparseMatrix.add((i,j))
    Sq2SparseMatrix = F2sparseMatrix(Sq2SparseMatrix)
    Sq2SparseMatrix.officialColumnRange = (0, nHomologyDim - 1)
    return Sq2SparseMatrix

def St(flowCategory, khovanovGenerators, gradings):
    bottomHGrading = gradings[0]
    qGrading = gradings[1]
    if khovanovGenerators.hqToSparseMatrix[gradings].sparseMatrix == set():
        return (0,0,0,0)
    #print(1)
    matrixSq2 = Sq2Matrix(khovanovGenerators, flowCategory, gradings)
    #print(2)
    matrixSq1iGrading = Sq1Matrix(khovanovGenerators, gradings)
    #print(3)
    matrixSq1iPlus1Grading = Sq1Matrix(khovanovGenerators,(bottomHGrading + 1, qGrading))
    #print(4)
    kernelMatrixSq1iGrading = kernel(matrixSq1iGrading)
    #print(5)    
    Sq2RestrictedToKernelColumns = {frozenset(matrixSq2 * v) for v in kernelMatrixSq1iGrading.sparseEchelonSet}
    r1 = rank(matrixSq2)
    Sq2RestrictedToKernelImage = reduceRowSparse(Sq2RestrictedToKernelColumns)
    Sq2RestrictedToKernelRank = len(Sq2RestrictedToKernelImage.pivots)
    r2 = Sq2RestrictedToKernelRank
    r3 = intersectionDim(image(matrixSq2),image(matrixSq1iPlus1Grading))
    r4 = intersectionDim(image(matrixSq1iPlus1Grading),Sq2RestrictedToKernelImage)
    return (r2-r4,r1-r2-r3+r4,r4,r3-r4)

def StList(diagram, homologies):
    if homologies == set():
        return tuple()
    else:
        khovanovGenerators = khovanovBasis(diagram)
        hValues = range(khovanovGenerators.minHomGrading + 1, khovanovGenerators.maxHomGrading - 2 + 1) #the minimum hom grading is actually one less than the lowest nontrivial chain complex grading. This is for the purposes of taking homology: we have to reach 1 lower to take the nth homology
        qSet = khovanovGenerators.qSet
        #print('done with Khovanov Generators')
        if isinstance(homologies, set):
            l = []
            homologiesPlusMirrorHomologies = set()
            for n,q in homologies:
                if n in hValues and q in qSet:
                    homologiesPlusMirrorHomologies.add((n,q))
                if -n-2 in hValues and -q in qSet:
                    homologiesPlusMirrorHomologies.add((-n-2,-q))
            for n,q in homologiesPlusMirrorHomologies:
                hqFlowCategory = flowCategory(diagram, khovanovGenerators, (n,q))
                StOutput = St(hqFlowCategory, khovanovGenerators, (n,q))
                if StOutput != (0,0,0,0):
                    l.append(((n,q), StOutput))
        elif homologies == 'all':
            l = []
            for n in hValues:
                for q in qSet:
                    #print(n,q)
                    hqToHomologyLoader(khovanovGenerators, {(n,q)})
            print('done loading homologies')
            for n in hValues:
                for q in qSet:
                    hqFlowCategory = flowCategory(diagram, khovanovGenerators, (n,q))
                    StOutput = St(hqFlowCategory, khovanovGenerators, (n,q))
                    if StOutput != (0,0,0,0):
                        l.append(((n,q), StOutput))
    return tuple(l)

knot4_1Dict = (((1,2),0), ((3,4),0), ((2,3),'-'), ((2,3),'-'), ((3,4),'+'), ((3,4),'+'), ((2,3),1),((1,4),1))
connectSumKnot4_1Dict = (((1,2),0), ((3,4),0), ((2,3),'-'), ((2,3),'-'), ((3,4),'+'), ((3,4),'+'), ((2,3),1), ((5,6),0), ((4,5),'-'), ((4,5),'-'), ((5,6),'+'), ((5,6),'+'), ((4,5),1),((1,6),1))
doubleKnot4_1Dict = (((1,2),0), ((3,4),0), ((2,3),'-'), ((2,3),'-'), ((3,4),'+'), ((3,4),'+'), ((2,3),1),((1,4),1), ((1,2),0), ((3,4),0), ((2,3),'-'), ((2,3),'-'), ((3,4),'+'), ((3,4),'+'), ((2,3),1),((1,4),1))
linkL6n1Dict = (((1,4),0), ((6,5), 0), ((2,3), 0), ((4,5),'-'), ((1,2),'+'), ((3,4),'+'), ((2,3),'+'), ((1,2),'-'), ((3,4),'-'), ((2,3),1), ((4,5),1), ((1,6),1))
linkL6n1MirrorDict = (((1,4),0), ((6,5), 0), ((2,3), 0), ((4,5),'+'), ((1,2),'-'), ((3,4),'-'), ((2,3),'-'), ((1,2),'+'), ((3,4),'+'), ((2,3),1), ((4,5),1), ((1,6),1))
knot8_19Dict = (((1, 6), 0), ((7, 8), 0), ((6, 7), "+"), ((2, 3), 0), ((1, 2), "-"), ((1, 2), "-"), ((4, 5), 0), ((3, 4), "+"), ((2, 3), "-"), ((2, 3), "-"), ((4, 5), "-"), ((5, 6), "+"), ((4, 5), 1), ((6, 3), 1), ((7, 2), 1), ((8, 1), 1))
knot8_19MirrorDict = (((1, 6), 0), ((7, 8), 0), ((6, 7), "-"), ((2, 3), 0), ((1, 2), "+"), ((1, 2), "+"), ((4, 5), 0), ((3, 4), "-"), ((2, 3), "+"), ((2, 3), "+"), ((4, 5), "+"), ((5, 6), "-"), ((4, 5), 1), ((6, 3), 1), ((7, 2), 1), ((8, 1), 1))
knot10_124Dict = (((1,2), 0), ((3,8),0), ((2,3), '-'), ((4,5), 0), ((6,7), 0), ((5,6), '-'), ((7,8), '+'), ((4,5), '+'), ((7,8), '+'), ((4,5), '+'), ((6,7), '+'), ((4,5), '+'), ((6,7), '+'), ((3,4), '-'), ((4,5), 1), ((3,6), 1), ((2,7), 1), ((1,8), 1))
otherKnot10_124Dict = (((1,9),0),((2,8),0),((1,2),'-'),((3,4),0),((5,7),0),((4,5),'-'),((7,8),'+'), ((7,8),'+'),((3,4),'+'),((8,9),1),((3,4),'+'),((5,7),'+'),((5,7),'+'),((3,4),'+'),((2,3),'-'),((3,4),1),((2,5),1),((1,7),1))
knot10_124mDict = (((1,2), 0), ((3,8),0), ((2,3), '+'), ((4,5), 0), ((6,7), 0), ((5,6), '+'), ((7,8), '-'), ((4,5), '-'), ((7,8), '-'), ((4,5), '-'), ((6,7), '-'), ((4,5), '-'), ((6,7), '-'), ((3,4), '+'), ((4,5), 1), ((3,6), 1), ((2,7), 1), ((1,8), 1))
knot10_132Dict = (((1,4), 0), ((5,8), 0), ((2,3), 0), ((4,5), '-'), ((6,7), 0), ((3,4), '-'), ((7,8), '+'), ((4,5), '+'), ((7,8), '+'), ((4,5), '+'), ((4,5), '+'), ((5,6), '-'), ((5,6), '-'), ((4,5), 1), ((3,6), '-'), ((2,3), 1), ((6,7), 1), ((1,8), 1))
knot10_132mDict = (((1,4), 0), ((5,8), 0), ((2,3), 0), ((4,5), '+'), ((6,7), 0), ((3,4), '+'), ((7,8), '-'), ((4,5), '-'), ((7,8), '-'), ((4,5), '-'), ((4,5), '-'), ((5,6), '+'), ((5,6), '+'), ((4,5), 1), ((3,6), '+'), ((2,3), 1), ((6,7), 1), ((1,8), 1))
knot10_152Dict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((5, 6), 0), ((4, 5), "-"), ((3, 4), "+"), ((4, 5), 1), ((7, 8), 0), ((6, 7), "+"), ((3, 6), "-"), ((2, 3), "+"), ((3, 6), 1), ((1, 2), "-"), ((7, 8), "-"), ((7, 2), 1), ((8, 1), 1))
knot10_152mDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((5, 6), 0), ((4, 5), "+"), ((3, 4), "-"), ((4, 5), 1), ((7, 8), 0), ((6, 7), "-"), ((3, 6), "+"), ((2, 3), "-"), ((3, 6), 1), ((1, 2), "+"), ((7, 8), "+"), ((7, 2), 1), ((8, 1), 1))
knot10_154Dict = (((1, 6), 0), ((7, 10), 0), ((6, 7), "-"), ((2, 3), 0), ((1, 2), "+"), ((1, 2), "+"), ((9, 8), 0), ((7, 8), "+"), ((8, 9), "-"), ((4, 5), 0), ((3, 4), "-"), ((5, 6), "-"), ((6, 7), "+"), ((5, 6), "-"), ((7, 6), 1), ((5, 8), "-"), ((4, 5), 1), ((8, 3), 1), ((9, 2), 1), ((10, 1), 1))
knot10_154mDict = (((1, 6), 0), ((7, 10), 0), ((6, 7), "+"), ((2, 3), 0), ((1, 2), "-"), ((1, 2), "-"), ((9, 8), 0), ((7, 8), "-"), ((8, 9), "+"), ((4, 5), 0), ((3, 4), "+"), ((5, 6), "+"), ((6, 7), "-"), ((5, 6), "+"), ((7, 6), 1), ((5, 8), "+"), ((4, 5), 1), ((8, 3), 1), ((9, 2), 1), ((10, 1), 1))
knot10_161Dict = (((1, 2), 0), ((6, 5), 0), ((2, 5), "-"), ((7, 8), 0), ((6, 7), "-"), ((5, 6), "+"), ((6, 7), 1), ((3, 4), 0), ((2, 3), "+"), ((4, 5), "+"), ((4, 5), "+"), ((3, 4), "-"), ((5, 8), "+"), ((4, 5), "-"), ((5, 8), 1), ((3, 4), "-"), ((3, 2), 1), ((4, 1), 1))
knot10_161mDict = (((1, 2), 0), ((6, 5), 0), ((2, 5), "+"), ((7, 8), 0), ((6, 7), "+"), ((5, 6), "-"), ((6, 7), 1), ((3, 4), 0), ((2, 3), "-"), ((4, 5), "-"), ((4, 5), "-"), ((3, 4), "+"), ((5, 8), "-"), ((4, 5), "+"), ((5, 8), 1), ((3, 4), "+"), ((3, 2), 1), ((4, 1), 1))
knot10_162Dict = (((1,2),0),((5,6),0),((7,8),0),((2,5),'+'),((6,7),'+'),((3,4),0),((2,3),'+'),((4,5),'+'),((2,3),'+'),((4,5),'+'),((1,2),'-'),((4,5),'+'),((3,4),'+'),((5,6),1),((3,4),'+'),((2,3),1),((4,7),1),((1,8),1))
knot11_n6Dict = (((1, 12), 0), ((13, 14), 0), ((12, 13), "+"), ((7, 2), 0), ((1, 2), "-"), ((1, 2), "-"), ((3, 6), 0), ((2, 3), "+"), ((8, 9), 0), ((7, 8), "-"), ((6, 7), "+"), ((7, 8), 1), ((5, 4), 0), ((3, 4), "+"), ((4, 5), "-"), ((5, 6), 1), ((11, 10), 0), ((9, 10), "+"), ((11, 12), "+"), ((9, 10), "+"), ((9, 4), 1), ((10, 11), 1), ((3, 12), 1), ((13, 2), 1), ((14, 1), 1))
conwayKnotDict = (((5,10), 0), ((1,4), 0), ((2,3), 0), ((4,5), '-'), ((6,7), 0), ((8,9), 0), ((3,4), '+'), ((5,6), '-'), ((7,8), '-'), ((9,10), '+'), ((2,3), '-'), ((4,5), '+'), ((6,7), '+'), ((9,10), '+'), ((3,4), '-'), ((5,6), 1), ((3,4), '-'), ((4,7), 1), ((3,8), 1), ((2,9), 1), ((1,10), 1))
modifiedDisjointTrefoilDict = (((1,6), 0), ((2, 3), 0), ((5, 4), 0), ((3, 4), "+"), ((3, 4), "+"), ((3, 4), "+"), ((3, 2), 1), ((4, 5), 1), ((5, 4), 0), ((1, 4), "+"), ((1, 4), "+"), ((1, 4), "+"), ((4, 5), 1), ((1,6), 1))
disjointTrefoilDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
connectSumTrefoilDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((1, 2), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
disjointTrefoilMirrorDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((2, 1), 1), ((3, 4), 1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((2, 1), 1), ((3, 4), 1))
trefoilPlus8_1Dict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((1, 2), "-"), ((1, 2), "-"), ((1, 2), "-"), ((1, 2), "-"), ((1, 2), "-"), ((2, 3), "+"), ((2, 3), "+"), ((1, 2), 1), ((3, 4), 1))
trefoilPlus8_1MirrorDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((2, 1), 1), ((3, 4), 1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((1, 2), "+"), ((1, 2), "+"), ((1, 2), "+"), ((1, 2), "+"), ((1, 2), "+"), ((2, 3), "-"), ((2, 3), "-"), ((1, 2), 1), ((3, 4), 1))
trefoilDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
figure8TrefoilDict = (((1,2),0), ((3,4),0), ((2,3),'-'), ((2,3),'-'), ((3,4),'+'), ((3,4),'+'), ((2,3),1),((1,4),1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
figure8TrefoilMirrorDict = (((1,2),0), ((3,4),0), ((2,3),'-'), ((2,3),'-'), ((3,4),'+'), ((3,4),'+'), ((2,3),1),((1,4),1),    ((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((2, 1), 1), ((3, 4), 1))
figure8MirrorTrefoilDict = (((1,2),0), ((3,4),0), ((2,3),'+'), ((2,3),'+'), ((3,4),'-'), ((3,4),'-'), ((2,3),1),((1,4),1),    ((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
randomDisjointUnionDict = (((1,4),0), ((6,5), 0), ((2,3), 0), ((4,5),'-'), ((1,2),'+'), ((3,4),'+'), ((2,3),'+'), ((1,2),'-'), ((3,4),'-'), ((2,3),1), ((4,5),1), ((1,6),1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
randomDisjointUnionMirrorDict = (((1,4),0), ((6,5), 0), ((2,3), 0), ((4,5),'+'), ((1,2),'-'), ((3,4),'-'), ((2,3),'-'), ((1,2),'+'), ((3,4),'+'), ((2,3),1), ((4,5),1), ((1,6),1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((2, 1), 1), ((3, 4), 1))
trefoil3ConnectSumDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((1, 2), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((1, 2), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
trefoil3xDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "+"), ((2, 3), "+"), ((2, 3), "+"), ((2, 1), 1), ((3, 4), 1))
trefoilMirror3xDict = (((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((2, 1), 1), ((3, 4), 1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((2, 1), 1), ((3, 4), 1), ((1, 2), 0), ((4, 3), 0), ((2, 3), "-"), ((2, 3), "-"), ((2, 3), "-"), ((2, 1), 1), ((3, 4), 1))
linkDict = trefoil3ConnectSumDict
linkMorseLink = morseLink(linkDict)
print(linkMorseLink.componentNo)
print('componentNo:', linkMorseLink.componentNo)
#link = khovanovBasis(linkMorseLink)
Q = -17
#mh = link.hqToSparseMatrix[(-8,Q)]
#mg = link.hqToSparseMatrix[(-7,Q)]
#mf = link.hqToSparseMatrix[(-6,Q)]
#me = link.hqToSparseMatrix[(-5,Q)]
#md = link.hqToSparseMatrix[(-4,Q)]
#mc = link.hqToSparseMatrix[(-3,Q)]
#mb = link.hqToSparseMatrix[(-2,Q)]
#ma = link.hqToSparseMatrix[(-1, Q)]
#a = link.hqToSparseMatrix[(0, Q)]
#b = link.hqToSparseMatrix[(1, Q)]
#c = link.hqToSparseMatrix[(2, Q)]
#d = link.hqToSparseMatrix[(3, Q)]
#e = link.hqToSparseMatrix[(4, Q)]
#f = link.hqToSparseMatrix[(5, Q)]
#g = link.hqToSparseMatrix[(6, Q)]
#h = link.hqToSparseMatrix[(7, Q)]
#i = link.hqToSparseMatrix[(8, Q)]
#j = link.hqToSparseMatrix[(9, Q)]

#print(a)
#print(b)
#print(kernel(b))
#print(image(a))
#print(b*a)
#print(ma)
#print(mb)
#print(ma*mb)
#print(c)
#print(d)
#print(e)
#print(b*a)
#print(c*b)
#print(d*c)
#print(e*d)

"""
print('-7Homology:', homology(mg,mh))
print('-6Homology:', homology(mf,mg))
print('-5Homology:', homology(me,mf))
print('-4Homology:', homology(md,me))
print('-3Homology:', homology(mc,md))
print('-2Homology:', homology(mb,mc))
print('-1Homology:', homology(ma,mb))
print('0Homology:', homology(a,ma))
print('1Homology:', homology(b,a))
print('2Homology:', homology(c,b))
print('3Homology:', homology(d,c))
print('4Homology:', homology(e,d))
print('5Homology:', homology(f,e))
print('6Homology:', homology(g,f))
print('7Homology:', homology(h,g))
print('8Homology:', homology(i,h))
print('9Homology:', homology(j,i))
"""

#fC = flowCategory(linkMorseLink, link, (5, 19))
#kernel for (-4. -13) on 10_124
#kernelElement = {768, 770, 772, 773, 774, 775, 776, 777, 778, 786, 788, 795, 797, 799, 803, 804, 815, 562, 819, 565, 823, 569, 575, 832, 584, 596, 599, 601, 604, 610, 613, 627, 628, 632, 635, 636, 638, 639, 640, 642, 645, 649, 650, 654, 660, 662, 666, 670, 677, 684, 690, 693, 695, 697, 700, 701, 445, 708, 709, 710, 715, 719, 738, 740, 742, 756, 759, 762, 764}
#kernelElement = {0, 1, 4, 6, 7, 11, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 27, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 48, 49, 50, 52, 53, 56, 57, 58, 59, 60, 61, 62, 65, 66, 67, 69, 70, 71, 72, 74, 75, 77, 78, 79, 80, 81, 85, 86, 87, 89, 90, 92, 93, 94, 95, 97, 98, 101, 103, 104, 105, 109, 111, 112, 113, 114, 115, 116, 117, 119, 120, 121, 122, 125, 127, 128, 130, 132, 134, 135, 136, 138, 139, 140, 141, 142, 144, 145, 147, 148, 149, 150, 151, 152, 153}
#kernelElement = {2}
#bM = matchings(fC, kernelElement)
#square = Sq2(link, (5, 19), fC, kernelElement)
#print('square:', square)
#print('Sq2 of element times n+2 matrix:', h * square)
#print(xFaceAssignment(linkMorseLink))
#print(xEdgeAssignment(linkMorseLink))

#print(St(fC, link, (3, 11)))
print(StList(linkMorseLink, 'all'))
