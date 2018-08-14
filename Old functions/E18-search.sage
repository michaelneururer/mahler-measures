bound1,bound2,bound3,bound4= 0,2,-1,2

input1 = [[a,b,c,d,e,f,g] for a in srange(bound1,bound2,1) for b in srange(bound1,bound2,1) for c in srange(bound1,bound2,1) for d in srange(bound1,bound2,1) for e in srange(bound1,bound2,1) for f in srange(bound1,bound2,1) for g in srange(bound1,bound2,1)]
input2 = [[a,b,c,d,e,f,g] for a in srange(bound3,bound4,1) for b in srange(bound3,bound4,1) for c in srange(bound3,bound4,1) for d in srange(bound3,bound4,1) for e in srange(bound3,bound4,1) for f in srange(bound3,bound4,1) for g in srange(bound3,bound4,1)]

for aa,bb,cc,dd,ee,ff,gg in input1:
    F = [aa,bb,cc,dd,ee,ff,gg]
    if F == [0,0,0,0,0,0,0]:
        continue
    for gg,hh,ii,jj,kk,ll,mm in input2:
        G=[gg,hh,ii,jj,kk,ll,mm]
        if F==G or G == [0,0,0,0,0,0,0] or proj_eq(vector(F),vector(G)):
            continue
        for xx in [0,1]:
            for yy in [0,1]:
                for zz in [0,1]:
                    for xxx in [-1,0,1]:
                        for yyy in [-1,0,1]:
                            for zzz in [-1,0,1]:
                                print F,G,xx,yy,zz,  xxx, yyy,zzz
                                flag = 0
                                for kA in range(7):
                                    if tame_symb(F,G, kA, xx, yy, zz, xxx,yyy,zzz) != 1 and tame_symb(F,G, kA, xx, yy,zz, xxx,yyy,zzz) != -1:
                                        flag = 1
                                        break
                                if flag == 0:
                                    if not check_Steinberg(F,G,xx,yy,zz,xxx,yyy,zzz):
                                        print'Found one!'
                                        fi.write(str([F,G,xx,yy,zz,xxx,yyy,zzz])+'\n')
