import matplotlib.pyplot as plt

def absorb(filename):
    temp_list=[]
    with open(filename,"r") as f:
        for i in f.readlines()[1:]:
            temp_list.append([float(j) for j in i.split(" ") if j][:2])
    temp_list=list(zip(*temp_list))

    plt.plot(temp_list[0],temp_list[1])
    maxy=[]
    for i in range(len(temp_list[1]))[1:-1]:
        if temp_list[1][i-1]<temp_list[1][i]>temp_list[1][i+1]:
            maxy.append(i)
    for i in maxy:
        if temp_list[0][i]<7:
            plt.scatter(temp_list[0][i],temp_list[1][i],color="red")
            print(temp_list[0][i],temp_list[1][i])
            plt.annotate(f"{temp_list[0][i]:.2f}",xy=(temp_list[0][i],temp_list[1][i]),xytext=(-10, 10), textcoords='offset points')
    plt.xlim(0, 10)
    #plt.ylim(0, 4000)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Absorption (cm$^{-1}$)")
    #plt.title(f"peaking at {temp_list[0][(temp_list[1].index(max(temp_list[1])))]:.2f} eV")
    plt.show()
    #plt.savefig("absorb_nb.png", dpi=300)
    #plt.close('all')

absorb("ABSORB_nb.dat")


