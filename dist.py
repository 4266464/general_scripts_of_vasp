#!/home/phys/qif/anaconda3/bin/python
#!/share/home/ckduan/anaconda3/bin/python

import scc_lib
import poscar
import numpy as np
import argparse

# better distort
# Cartesian


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("poscar_name",default="CONTCAR",help="the name of POSCAR type file, eg POSCAR CONTCAR")
    parser.add_argument("-a","--atom_number",default=0,type=int,help="the number specify which atom is manipulated, eg 0 1 2, default 0")
    parser.add_argument("-m","--manipulate",nargs="+",default=["show"],help="detailed manipulation to atoms, eg random 6*0.1, default show")
    args = parser.parse_args()
    print(args.manipulate)

    p=poscar.Poscar()
    p.read(args.poscar_name)

    def distance_string(aposcar):
        string = ""
        for i in ["atom", "distance", "x", "y", "z"]:
            string += f"{i:10}"
        string += "\n"
        for i in aposcar.distance_from_point_to_line(aposcar.label[args.atom_number], )[::-1]:
            string += f"{i[0][0]}#{i[0][1]:3}"
            string += f"{aposcar.label.index(i[0]):4}"
            string += f"{np.linalg.norm(i[1]):10.5f}"
            for j in i[1]:
                if abs(j) < 1e-5:
                    string += f"{0:10}"
                else:
                    string += f"{j:10.5f}"
            string += "\n"
        return string
    
    if args.manipulate==["show"]:
        # for distance display

        print(distance_string(p),end="")
        
    elif args.manipulate==["random"]:
        scc_lib.linux_command("mv "+args.poscar_name+" "+args.poscar_name+".old")
        p.distort_randomly().write(args.poscar_name)

    elif args.manipulate==["inter"]:
        #scc_lib.linux_command("mv "+args.poscar_name+" "+args.poscar_name+".old")
        poscar.local_inter(args.poscar_name)
        #p.make_inter_file(args.poscar_name)

    elif args.manipulate==["env"]:
        string=""
        for i in p.env(p.label[args.atom_number], 0.1):
            string += f"{i[0]} "
            for j in i[1]:
                string+=f"{j} "
            string+="\n"
        print(string,end="")

    elif "move" in args.manipulate:
        p.move([float(i) for i in args.manipulate[1:] if i])
        scc_lib.linux_command("mv "+args.poscar_name+" "+args.poscar_name+".old")
        p.write(args.poscar_name)

    else:
        distort_list=[]
        for i in args.manipulate:
            if "*" in i:
                t,d=i.split("*")
                t=int(t)
                d=float(d)
                for j in range(t):
                    distort_list.append(d)
            else:
                distort_list.append(float(i))
        scc_lib.linux_command("mv "+args.poscar_name+" "+args.poscar_name+".old")
        q= p.distort_sphere(distort_list, args.atom_number)
        q.write(args.poscar_name)
        print(distance_string(q),end="")
