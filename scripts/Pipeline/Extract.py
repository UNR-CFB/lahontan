#!/usr/bin/python3

import glob
import os
import shutil
from pprint import pprint
import subprocess

def getDir(read):
    GroupNum = os.path.basename(read).split('_')[0][4]
    Sample = os.path.basename(read).split('_')[0][-1]
    Lane = os.path.basename(read).split('_')[2][-1]
    if GroupNum == '2':
        GroupNum = '24'
    Dir = '/data/ceph/alberton/HExp' + '/Group{}Stretch{}Lane{}'.format(GroupNum,Sample,Lane)
    return Dir

def getDir2(read,AoB):
    GroupNum = os.path.basename(read).split('_')[0][4]
    Sample = os.path.basename(read).split('_')[0][-1]
    Lane = os.path.basename(read).split('_')[2][-1]
    if GroupNum == '2':
        GroupNum = '24'
    Dir1 = '/data/ceph/alberton/HExp' + '/Group{1}Stretch{2}-{0}'.format(AoB,GroupNum,Sample)
    Dir2 = Dir1 + '/Lane{0}'.format(Lane)
    return Dir1,Dir2

def getDir3(read,AoB):
    GroupNum = os.path.basename(read).split('_')[0][4]
    Sample = os.path.basename(read).split('_')[0][-1]
    Lane = os.path.basename(read).split('_')[2][-1]
    if GroupNum == '2':
        GroupNum = '24'
    Dir = '/data/ceph/alberton/HExp' + '/Group{1}Stretch{2}Lane{3}-{0}'.format(AoB,GroupNum,Sample,Lane)
    return Dir

def getName(read,AoB):
    a = os.path.basename(read).split('_')
    a.insert(4,str(AoB))
    if AoB == 'B' and a[1][-1] == '7':
        a[1] = 'S8'
    newname = '_'.join(a)
    return newname

def makeLink(read,Dir,name):
    try:
        os.mkdir(Dir)
    except OSError:
        pass
    try:
        os.symlink(read, os.path.join(Dir,name))
    except FileExistsError:
        pass

def makeLink2(read,Dir1,Dir2,name):
    try:
        os.mkdir(Dir1)
    except OSError:
        pass
    try:
        os.mkdir(Dir2)
    except OSError:
        pass
    try:
        os.symlink(read, os.path.join(Dir2,name))
    except FileExistsError:
        pass

def gatherStuff1(data,AoB):
    for group in data:
        Groups = sorted(glob.glob(group + '/*'))
        for read in Groups:
            Dir = getDir(read)
            newname = getName(read,AoB)
            makeLink(read,Dir,newname)

def gatherStuff2(data,AoB):
    for group in data:
        Groups = sorted(glob.glob(group + '/*'))
        for read in Groups:
            Dir1,Dir2 = getDir2(read,AoB)
            newname = getName(read,AoB)
            makeLink2(read,Dir1,Dir2,newname)

def gatherStuff3(data,AoB):
    for group in data:
        Groups = sorted(glob.glob(group + '/*'))
        for read in Groups:
            Dir = getDir3(read,AoB)
            newname = getName(read,AoB)
            makeLink(read,Dir,newname)

A = sorted(glob.glob('/data/ceph/alberton/HBurkin-A/*'))
B = sorted(glob.glob('/data/ceph/alberton/HBurkin-B/*'))

try:
    os.mkdir('/data/ceph/alberton/HExp')
except OSError:
    pass

#gatherStuff1(A,'A')
#gatherStuff1(B,'B')
#gatherStuff2(A,'A')
#gatherStuff2(B,'B')
#gatherStuff3(A,'A')
#gatherStuff3(B,'B')

def extract(data,path):
    try:
        os.mkdir(path)
    except OSError:
        pass
    Experiment = glob.glob(data + '/*')
    for sample in Experiment:
        reads = glob.glob(sample + '/*')
        for read in reads:
            os.symlink(read, os.path.join(path,os.path.basename(read)))

#extract('/data/ceph/alberton/HExpA','/data/ceph/alberton/Data_HB_A')
#extract('/data/ceph/alberton/HExpB','/data/ceph/alberton/Data_HB_B')

def extract2(data,path):
    try:
        os.mkdir(path)
    except OSError:
        pass
    Experiment = glob.glob(data + '/*')
    for sample in Experiment:
        newDir = os.path.join(path,os.path.basename(sample))
        try:
            os.mkdir(newDir)
        except OSError:
            pass
        R1 = glob.glob(os.path.join(sample,'*R1*'))
        R2 = glob.glob(os.path.join(sample,'*R2*'))
        R1 = sorted(R1,key=lambda x: x.split('_')[-2])
        R2 = sorted(R2,key=lambda x: x.split('_')[-2])
        N1 = os.path.basename(R1[0]).split('_')
        N1[4] = 'AB'
        del N1[1]
        name1 = os.path.join(newDir,'_'.join(N1))
        N2 = os.path.basename(R2[0]).split('_')
        N2[4] = 'AB'
        del N2[1]
        name2 = os.path.join(newDir,'_'.join(N2))
        command1 = r'cat {} {} > {}'.format(os.path.realpath(R1[0]),os.path.realpath(R1[1]),name1)
        command2 = r'cat {} {} > {}'.format(os.path.realpath(R2[0]),os.path.realpath(R2[1]),name2)
        subprocess.run(command1,shell=True,check=True)
        subprocess.run(command2,shell=True,check=True)

def extract3(data,path):
    try:
        os.mkdir(path)
    except OSError:
        pass
    letter = 'B'
    Experiment = glob.glob(data + '/*-{}'.format(letter))
    for sample in Experiment:
        newDir = os.path.join(path,os.path.basename(sample))
        try:
            os.mkdir(newDir)
        except OSError:
            pass
        R1,R2 = [],[]
        for lane in glob.glob(sample + '/*'):
            R1.append(glob.glob(os.path.join(lane,'*R1*'))[0])
            R2.append(glob.glob(os.path.join(lane,'*R2*'))[0])
        R1 = sorted(R1)
        R2 = sorted(R2)
        N1 = os.path.basename(R1[0]).split('_')
        N1[4] = letter
        N1[2] = 'alllanes'
        del N1[1]
        name1 = os.path.join(newDir,'_'.join(N1))
        N2 = os.path.basename(R2[0]).split('_')
        N2[4] = letter
        N2[2] = 'alllanes'
        del N2[1]
        name2 = os.path.join(newDir,'_'.join(N2))
        command1 = r'cat {} {} {} {} > {}'.format(os.path.realpath(R1[0]),os.path.realpath(R1[1]),os.path.realpath(R1[2]),os.path.realpath(R1[3]),name1)
        command2 = r'cat {} {} {} {} > {}'.format(os.path.realpath(R2[0]),os.path.realpath(R2[1]),os.path.realpath(R2[2]),os.path.realpath(R2[3]),name2)
        subprocess.run(command1,shell=True,check=True)
        subprocess.run(command2,shell=True,check=True)

def combine(a,b,path):
    try:
        os.mkdir(path)
    except OSError:
        pass
    A = sorted(glob.glob(os.path.join(a,'*')))
    B = sorted(glob.glob(os.path.join(b,'*')))
    for read in A:
        N1 = os.path.basename(read).split('_')
        N1[3] = 'AB'
        name1 = os.path.join(path,'_'.join(N1))
        command1 = r'cat {} {} > {}'.format(os.path.realpath(A[A.index(read)]),os.path.realpath(B[A.index(read)]),name1)
        subprocess.run(command1,shell=True,check=True)




#extract2('/data/ceph/alberton/HExpAB','/data/ceph/alberton/HBurkin-AB')
#extract('/data/ceph/alberton/HBurkin-AB','/data/ceph/alberton/Data_HB_AB')
#extract3('/data/ceph/alberton/HExp2','/data/ceph/alberton/Data_HB_A-lanes')
#extract3('/data/ceph/alberton/HExp2','/data/ceph/alberton/Data_HB_B-lanes')
#extract('/data/ceph/alberton/HExpB-lanes','/data/ceph/alberton/Data_HB_B-lanes')
combine('/data/ceph/alberton/Data_HB_A-lanes','/data/ceph/alberton/Data_HB_B-lanes','/data/ceph/alberton/Data_HB_AB-lanes')






















