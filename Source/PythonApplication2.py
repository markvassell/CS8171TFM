seq=open("../Data/T0880/T0880.fasta","r")
pos=open("../Data/T0880/t000_.200.9mers","r")
seq.readline()
se=seq.readline()
length=len(se)-1
L1=list(se);
xunhuan=length//9+1

L2=[]
L3=[]

for i in range(0, xunhuan):
    #print(i)
    tp=pos.readline()
    tp2=tp.split(" ");
    while tp2[0]!= 'position:':
        tp=pos.readline()
        tp2=tp.split(" ");
    pos.readline()
    for j in range(0,9):
        tp=pos.readline()
        tp2=tp.split(" ");
        while '' in tp2:
            tp2.remove('')
        #print((tp2))
        #break
        L2.append(tp2[5])
        L3.append(tp2[6])
L2=L2[0:length]
L3=L3[0:length]
L=[L1,L2,L3]

L=[list(i) for i in zip(*L)]
print(L)

seq.close()
pos.close()