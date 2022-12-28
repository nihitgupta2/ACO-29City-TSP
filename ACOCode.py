# Nihit Gupta
import math
import random
import matplotlib.pyplot as plt

def euclidean_distance_matrix(nodes):
    #creating a matrix containing all the distance from one node to another node
    totalNodes = len(nodes)
    rows, cols = (totalNodes, totalNodes)
    nodeDistances = [[0 for c in range(cols)] for r in range(rows)]
    for i in range(totalNodes-1):
        for j in range(i+1,totalNodes):
            distance = 0
            distance = math.sqrt(pow(nodes[j][0]-nodes[i][0], 2)+pow(nodes[j][1]-nodes[i][1], 2))
            nodeDistances[i][j]=distance
            nodeDistances[j][i]=distance
    return nodeDistances


def ACO_Algorithm(nodes,iterations,m_ants,distanceMatrix,q):
    # ACO algorithm where all ants start from the same starting point 
    globalBestDistance = 1000000
    globalBestRoute = None
    listBestPerIteration = []
    listAveragePerIteration = []
    totalNodes = len(nodes)
    rows, cols = (totalNodes, totalNodes)
    pheromoneTable = [[0.95 for c in range(cols)] for r in range(rows)] # global pheromone table 
    nodeList = nodes[:]
    decay_p = 0.965 # p E (0,1] quatity left after decay
    alpha = 1
    beta = 6 # 6 was good
    
    iter=0
    while(iter<iterations):
        antcounter = 0
        recordSolutionList = []
        pheromoneTableDeposit = [[0 for c in range(cols)] for r in range(rows)]
        while(antcounter<m_ants):
            currentNode = 0   #Node number 0 of our list is the starting node number
            nodesVisited = []
            nodeNumberVisted = []
            i = 0
            while(len(nodesVisited)<totalNodes):
                probabilityTracker = [-1 for items in range(totalNodes)]
                nodesVisited.append(nodeList[currentNode])
                nodeNumberVisted.append(currentNode)
                if len(nodeNumberVisted)==totalNodes:
                    break

                for z in range(totalNodes):
                    if nodeList[z] not in nodesVisited:
                        pk = (pow(pheromoneTable[currentNode][z],alpha))/(pow(distanceMatrix[currentNode][z],beta)) #probability calculation
                        probabilityTracker[z] = pk

                sum = 0
                for element in probabilityTracker:
                    if element>0:
                        sum += element 
                for elementCount in range(len(probabilityTracker)):
                    if probabilityTracker[elementCount] > 0:
                        probabilityTracker[elementCount]=probabilityTracker[elementCount]/sum
                listNo_NodeProbability = []
                listVal_NodeProbability = []
                for element in range(len(probabilityTracker)):
                    if probabilityTracker[element]>0:
                        listNo_NodeProbability.append(element)
                        listVal_NodeProbability.append(probabilityTracker[element])
                selectedNode = random.choices(listNo_NodeProbability,weights=listVal_NodeProbability,k=1) # selecting the next node based probability
                currentNode = selectedNode[0]           #making next node current node

            #nodesVisited example [1,5,3,8,7,9,2,4,6]
            totalCurrentDistance = 0
            for eachNode in range(len(nodeNumberVisted)-1):
                temp1 = nodeNumberVisted[eachNode]
                temp2 = nodeNumberVisted[eachNode+1]
                totalCurrentDistance=totalCurrentDistance+distanceMatrix[temp1][temp2]
            
            for eachNode in range(len(nodeNumberVisted)-1):
                temp1 = nodeNumberVisted[eachNode]
                temp2 = nodeNumberVisted[eachNode+1]
                pheromoneTableDeposit[temp1][temp2] = pheromoneTableDeposit[temp1][temp2] + (q/totalCurrentDistance)
                pheromoneTableDeposit[temp2][temp1] = pheromoneTableDeposit[temp2][temp1] + (q/totalCurrentDistance)

            if totalCurrentDistance<globalBestDistance:
                globalBestDistance=totalCurrentDistance
                globalBestRoute = nodesVisited
            recordSolutionList.append(totalCurrentDistance)
            antcounter+=1
        #artificial pheromone evaporation
        pheromoneTable = [[decay_p*pheromoneTable[r][c] for c in range(cols)] for r in range(rows)]
        #depositing the pheromons according to the explored paths of the ants
        for icounter in range(rows):
            for jcounter in range(cols):
                pheromoneTable[icounter][jcounter] = pheromoneTable[icounter][jcounter] + pheromoneTableDeposit[icounter][jcounter]
        recordSolutionListSorted = sorted(recordSolutionList)
        listBestPerIteration.append(recordSolutionListSorted[0])
        sumAverage = 0
        sumAverage = math.fsum(recordSolutionList)/m_ants 
        listAveragePerIteration.append(sumAverage)
        iter+=1
    return globalBestRoute,globalBestDistance,listBestPerIteration,listAveragePerIteration


def main():
    pointNodes = [[1150, 1760], 
                  [630, 1660],
                  [40, 2090],
                  [750, 1100],
                  [750, 2030],
                  [1030, 2070],
                  [1650, 650],
                  [1490, 1630],
                  [790, 2260],
                  [710, 1310],
                  [840, 550],
                  [1170, 2300],
                  [970, 1340],
                  [510, 700],
                  [750, 900],
                  [1280, 1200],
                  [230, 590],
                  [460, 860],
                  [1040, 950],
                  [590, 1390],
                  [830, 1770],
                  [490, 500],
                  [1840, 1240],
                  [1260, 1500],
                  [1280, 790],
                  [490, 2130],
                  [1460, 1420],
                  [1260, 1910],
                  [360, 1980]]

    n2nDistances = euclidean_distance_matrix(pointNodes) #node to node distance matrix
    iterations = 75
    m_ants = 150
    q = 10
    route,distances,bestPerIter,avgPerIter = ACO_Algorithm(pointNodes,iterations,m_ants,n2nDistances,q)
    print (route)
    print(distances)
    #For generating plots
    iterXComponenet = [i+1 for i in range(iterations)]
    fig,(ax1,ax2) = plt.subplots(1,2)
    fig.suptitle("Graphs")
    ax1.plot(iterXComponenet,bestPerIter)
    ax1.set_title("Best Solution Vs Iteration")
    ax2.plot(iterXComponenet,avgPerIter)
    ax2.set_title("Average value of solution Vs Iteration")
    plt.show()
    return True


if __name__ == '__main__':
    main()
