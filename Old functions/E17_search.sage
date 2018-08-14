bound1,bound2,bound3,bound4= 0,3,-2,3

input1 = [[a,b,c,d,e,f] for a in srange(bound1,bound2,1) for b in srange(bound1,bound2,1) for c in srange(bound1,bound2,1) for d in srange(bound1,bound2,1) for e in srange(bound1,bound2,1) for f in srange(bound1,bound2,1)]
input2 = [[a,b,c,d,e,f] for a in srange(bound3,bound4,1) for b in srange(bound3,bound4,1) for c in srange(bound3,bound4,1) for d in srange(bound3,bound4,1) for e in srange(bound3,bound4,1) for f in srange(bound3,bound4,1)]

for aa,bb,cc,dd,ee,ff in input1:
    F = [aa,bb,cc,dd,ee,ff]
    if F == [0,0,0,0,0,0]:
        continue
    for gg,hh,ii,jj,kk,ll in input2:
        G=[gg,hh,ii,jj,kk,ll]
        if F==G or G == [0,0,0,0,0,0] or proj_eq(vector(F),vector(G)):
            continue
        for xx in [0,1]:
            for yy in [0,1]:
                for xxx in [-1,0,1]:
                    for yyy in [-1,0,1]:
                        print F,G,xx,yy,  xxx, yyy
                        flag = 0
                        for kA in range(6):
                            if tame_symb(F,G, kA, xx, yy, xxx,yyy) != 1 and tame_symb(F,G, kA, xx, yy,xxx,yyy) != -1:
                                    flag = 1
                                    break
                        if flag == 0:
                            if not check_Steinberg(F,G,xx,yy,xxx,yyy):
                                print'Found one!'
                                fi.write(str([F,G,xx,yy,xxx,yyy])+'\n')
