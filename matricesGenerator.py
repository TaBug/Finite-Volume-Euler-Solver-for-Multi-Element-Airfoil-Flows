import numpy as np


def readgri(fname):
    f = open(fname, 'r')
    Nn, Ne, dim = [int(s) for s in f.readline().split()]
    # read vertices
    V = np.array([[float(s) for s in f.readline().split()] for n in range(Nn)])
    # read boundaries
    NB = int(f.readline())
    B = np.array([[]])
    Bname = np.array([])
    for i in range(NB):
        s = f.readline().split()
        Nb = int(s[0])
        Bname = np.append(Bname, [i for s in range(Nb)])
        Bi = np.array([[int(s) for s in f.readline().split()] for n in range(Nb)])
        B = Bi if (B.size == 0) else np.append(B, Bi, axis=0)

    # read elements
    Ne0 = 0
    E = np.array([[]])
    while Ne0 < Ne:
        s = f.readline().split()
        ne = int(s[0])
        Ei = np.array([[int(s) for s in f.readline().split()] for n in range(ne)])
        E = Ei if (Ne0 == 0) else np.concatenate((E, Ei), axis=0)
        Ne0 += ne
    f.close()
    Mesh = {'V': V, 'E': E, 'B': B, 'Bname': Bname}
    return V, E, B, Bname


def getI2E(fnameInput, fnameOutput):
    nodes, NE, NB, NBName = readgri(fnameInput)
    print(NB, NB.shape)
    output = np.array([[]])
    faces = np.array([[]])
    with open(fnameOutput, 'w') as f:
        for i, elem in enumerate(NE):
            for e in range(3):
                if e == 0:
                    node1 = 1
                    node2 = 2
                elif e == 1:
                    node1 = 2
                    node2 = 0
                else:
                    node1 = 1
                    node2 = 2

                newFace = np.array([[elem[node1], elem[node2]]])
                if np.isin(NB, newFace).all(axis=1).any():
                    print(np.isin(NB, newFace).all(axis=1))
                    continue

                if faces.size == 0:
                    faces = newFace
                    output = np.array([[i + 1, e + 1, 0, 0]])
                else:
                    isin = np.isin(faces, newFace).all(axis=1)
                    if isin.any():
                        iFace = np.where(isin)[0][0]
                        output[iFace][2] = i + 1
                        output[iFace][3] = e + 1
                    else:
                        faces = np.append(faces, newFace, axis=0)
                        face = np.array([[i + 1, e + 1, 0, 0]])
                        output = np.append(output, face, axis=0)

        for i in range(len(output)):
            f.write(f'{int(output[i][0])} {int(output[i][1])} {int(output[i][2])} {int(output[i][3])}\n')
        f.close()
    return output


def getB2E(fnameInput):
    nodes, NE, NB, NBName = readgri(fnameInput)
    elems = []
    faces = []
    nBGroups = []

    for ib, nb in enumerate(NB):
        iElem = np.where(np.isin(NE, nb).all(axis=1))
        node1 = nb[0]
        node2 = nb[1]
        elem = NE[iElem]
        inode1 = np.where(elem, node1)
        inode2 = np.where(elem, node2)
        iface = 3 - inode1 - inode2

        elems.append(elem)
        faces.append(iface)
        nBGroups.append(NBName[ib])

    output = np.zeros([len(elems), 3])
    with open('B2E.txt', 'w') as f:
        for i in range(len(elems)):
            f.write(f'{int(elems[i])} {int(faces[i])} {int(nBGroups[i])}\n')
            output[i][0] = elems[i]
            output[i][1] = faces[i]
            output[i][2] = nBGroups[i]
        f.close()
    return output


def edgehash(fnameInput, fnameOutput):
    nodes, NE, NB, NBName = readgri(fnameInput)
    I2E = getI2E(fnameInput, fnameOutput)
    B2E = getB2E(fnameInput)
    In = []
    Bn = []
    faces = []
    lIn = []
    lBn = []
    with open('In.txt', 'w') as f:
        for iface, face in enumerate(I2E):
            elemL = face[0]
            faceL = face[1]
            elemR = face[2]
            faceR = face[3]
            if faceL == 0:
                node1 = 1
                node2 = 2
            elif faceL == 1:
                node1 = 2
                node2 = 0
            else:
                node1 = 0
                node2 = 1
            node1Global = NE[elemL][node1]
            node2Global = NE[elemL][node2]
            node1c = nodes[node1Global]
            node2c = nodes[node2Global]
            l = np.sqrt((node2c[0] - node1c[0]) ** 2 + (node2c[1] - node2c[0]) ** 2)
            n = np.array([(node2c[1] - node1c[1]) / l, -(node2c[0] - node1c[0]) / l])
            In.append(n)
            lIn.append(l)
            f.write(f'{n[0]} {n[1]}\n')

    with open(f'Bn.txt', 'w') as f:
        for iface, face in enumerate(B2E):
            elem = face[0]
            face = face[1]
            if face == 0:
                node1 = 1
                node2 = 2
            elif face == 1:
                node1 = 2
                node2 = 0
            else:
                node1 = 0
                node2 = 1
            node1Global = NE[elem][node1]
            node2Global = NE[elem][node2]
            node1c = nodes[node1Global]
            node2c = nodes[node2Global]
            l = np.sqrt((node2c[0] - node1c[0]) ** 2 + (node2c[1] - node2c[0]) ** 2)
            n = np.array([(node2c[1] - node1c[1]) / l, -(node2c[0] - node1c[0]) / l])
            Bn.append(n)
            lBn.append(l)
            f.write(f'{n[0]} {n[1]}\n')

    return In, Bn, lIn, lBn


# input: element matrix, node coordinate matrix
# output: element area matrix (index = element index)
def area(fnameInput, fnameOutput):
    nodes, NE, NB, NBName = readgri(fnameInput)
    with open(fnameOutput, 'w') as f:
        for i, ne in enumerate(NE):
            coor0 = nodes[int(ne[0]) - 1]
            coor1 = nodes[int(ne[1]) - 1]
            coor2 = nodes[int(ne[2]) - 1]
            area = 1 / 2 * (coor0[0] * (coor1[1] - coor2[1]) + coor1[0] * (coor2[1] - coor0[1]) + coor2[0] * (
                    coor0[1] - coor1[1]))
            f.write(f'area\n')


def main():
    getI2E('all.gri', 'I2E.txt')
    # B2E('all.gri', 'B2E.txt')


if __name__ == "__main__":
    main()
