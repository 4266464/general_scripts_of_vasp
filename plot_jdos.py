import matplotlib.pyplot as plt

def jdos(filename):
    temp_list=[]
    with open(filename,"r") as f:
        for i in f.readlines()[1:]:
            temp_list.append([float(j) for j in i.split(" ") if j])
    temp_list=list(zip(*temp_list))
    maxy = []
    for i in range(len(temp_list[1]))[1:-1]:
        if temp_list[1][i - 1] < temp_list[1][i] > temp_list[1][i + 1]:
            maxy.append(i)
    for i in maxy:
        if temp_list[0][i] < 7:
            plt.scatter(temp_list[0][i], temp_list[1][i], color="red")
            print(temp_list[0][i], temp_list[1][i])
            plt.annotate(f"{temp_list[0][i]:.2f}", xy=(temp_list[0][i], temp_list[1][i]), xytext=(-10, 10),
                         textcoords='offset points')

    plt.plot(temp_list[0],temp_list[1])
    #plt.xlim(0, 7)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Joint DOS")
    #plt.title(f"peaking at {temp_list[0][(temp_list[1].index(max(temp_list[1])))]:.2f} eV")
    plt.savefig("jdos_nb.png", dpi=300)
    plt.close('all')

jdos("JDOS_nb.dat")


