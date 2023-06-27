import os
import numpy as np
import random



class poscar:
    def __init__(self):
        pass

    def read(self, poscar_path):
        lines = []
        with open(poscar_path, "r") as poscar_file:
            for i in poscar_file.readlines():
                lines.append([tofloat(j) for j in i.replace("\n", "").split(" ") if j])

        self.comment = " ".join(lines[0])
        self.scale = lines[1][0]
        self.axis = np.array([i for i in lines[2:5]])
        self.element = lines[5]
        self.number = [int(i) for i in lines[6]]
        self.label = []
        # in case 相同元素未在一起
        for i, j in zip(self.element, self.number):
            for k in range(int(j) + 1)[1:]:
                while (i, k) in self.label:
                    k += 1
                self.label.append((i, k))

        self.position = []
        self.additional = []
        if lines[7][0][0] in ["S", "s"]:
            self.selective = True
            self.mode = lines[8][0]
            next_line = 9
        else:
            self.selective = False
            self.mode = lines[7][0]
            next_line = 8

        for i in lines[next_line:]:
            if i:
                self.position.append(i[:3])
                self.additional.append(i[3:])
            else:
                break

        next_line += len(self.position)
        self.position = np.array(self.position)
        # read CHGCAR
        start_read_chg = False
        for i in lines[next_line:]:
            if not start_read_chg:
                # 存在 CHGCAR
                if len(i) == 3 and all(i):
                    start_read_chg = True
                    self.chg = []
                    self.chg_diff = i
                    count = 0
            else:
                if i:
                    for j in i:
                        x = count % self.chg_diff[0]
                        y = count // self.chg_diff[0]
                        z = y // self.chg_diff[1]
                        y = y % self.chg_diff[1]
                        self.chg.append([x, y, z, j])
                        count += 1
        if not start_read_chg:
            self.chg = []
            self.chg_diff = []

        return None

    def d2c(self):
        aposcar = self.copy()
        if aposcar.mode[0] in ["D", "d"]:
            for i in range(len(aposcar.position)):
                aposcar.position[i] = aposcar.position[i].dot(aposcar.axis)
            aposcar.mode = "C"

        return aposcar

    def c2d(self):
        aposcar = self.copy()
        if aposcar.mode[0] in ["C", "c"]:
            taxis = np.linalg.inv(aposcar.axis)
            for i in range(len(aposcar.position)):
                aposcar.position[i] = aposcar.position[i].dot(taxis)
            aposcar.mode = "D"
        return aposcar

    def write(self, path):
        print(f"#atom {len(self.position)}")
        string = f"{self.comment}\n{self.scale}\n"
        for i in self.axis:
            for j in i:
                string += f"{j} "
            string += "\n"
        for i in self.element:
            string += f"{i} "
        string += "\n"
        for i in self.number:
            string += f"{i} "
        string += "\n"
        if self.selective:
            string += "Selective\n"
        string += f"{self.mode}\n"

        for i in range(len(self.position)):
            for k in self.position[i]:
                string += f"{k} "
            if i < len(self.additional):
                for k in self.additional[i]:
                    string += f"{k} "
            string += "\n"
        string += "\n"

        # for i, j in zip(self.position, self.tag):
        #     for k in i:
        #         string += f"{k} "
        #     for k in j:
        #         string += f"{k} "
        #     string += "\n"
        # string += "\n"

        if self.chg:
            for i in self.chg_diff:
                string += f"{int(i)} "
            for i in range(len(self.chg)):
                if i % 5 == 0:
                    string += "\n"
                string += f"{self.chg[i][-1]} "

        with open(path, "w") as f:
            f.write(string)
        return None

    def copy(self):
        # p = poscar()
        # p.comment = self.comment
        # p.scale = self.scale
        # p.lattice = self.lattice.copy()
        # p.element = self.element[:]
        # p.number = self.number[:]
        # p.label = self.label[:]
        # p.type = self.type
        # p.selective = self.selective
        # p.tag = self.tag[:]
        # p.chg_diff = self.chg_diff[:]
        # p.position = self.position.copy()
        # p.chg = self.chg[:]
        aposcar=poscar()
        aposcar.__dict__=self.__dict__.copy()
        return aposcar

    def __add__(self, other_poscar):
        aposcar = self.copy()
        aposcar.position = [i + j for i, j in zip(self.position, other_poscar.position)]
        return aposcar

    def add_chg(self, other_poscar):
        aposcar = self.copy()
        for i in range(len(aposcar.chg)):
            aposcar.chg[i][-1] += other_poscar.chg[i][-1]
        return aposcar

    def __sub__(self, other_poscar):
        aposcar = self.copy()
        aposcar.position = []
        for i in range(len(self.position)):
            diff = self.position[i] - other_poscar.position[i]
            for j in range(len(diff)):
                if diff[j] > 0.5:
                    diff[j] -= 1
                elif diff[j] < -0.5:
                    diff[j] += 1
            aposcar.position.append(diff)
        aposcar.position = np.array(aposcar.position)
        return aposcar

    def __mul__(self, multiplier):
        aposcar = self.copy()
        aposcar.position = multiplier * self.position
        return aposcar

    def add_dummy(self, a, b, c):
        app = [[i / a, j / b, k / c] for i in range(a) for j in range(b) for k in range(c)]
        self.element.append("H")
        self.number.append(len(app))
        self.position = np.append(self.position, np.array(app), 0)
        return None

    def move(self, vector):
        for i in range(len(self.position)):
            self.position[i] += vector
            for j in range(3):
                if self.position[i][j] >= 1:
                    self.position[i][j] -= 1
                elif self.position[i][j] < 0:
                    self.position[i][j] += 1

        if self.chg:
            vector = [round(i * j) for i, j in zip(vector, self.chg_diff)]
            for i in range(len(self.chg)):
                for j in range(3):
                    self.chg[i][j] += vector[j]
                    if self.chg[i][j] >= self.chg_diff[j]:
                        self.chg[i][j] -= self.chg_diff[j]
                    elif self.chg[i][j] < 0:
                        self.chg[i][j] += self.chg_diff[j]
        self.chg.sort(key=lambda x: x[0] + (x[1] + x[2] * self.chg_diff[1]) * self.chg_diff[0])

        return None

    def is_direct(self):
        if self.mode[0] in ["D", "d"]:
            return True
        else:
            return False

    def distance(self, site, full=False):
        if len(site) == 3:
            # 输入坐标, 相对!!
            site = np.array(site)
        else:
            site = self.position[self.label.index(site)]

        distance_list = []
        for i, j in zip(self.label, self.position):
            # 距离
            if full:
                distance_list.append([i, np.array(
                    [k if abs(k) < 0.5 else (k + 1 if k < 0 else k - 1) for k in (j - site)]).dot(self.axis)])
            else:
                distance_list.append([i, np.linalg.norm(
                    np.array([k if abs(k) < 0.5 else (k + 1 if k < 0 else k - 1) for k in (j - site)]).dot(self.axis))])

        distance_list.sort(key=lambda x: x[0])
        distance_list.sort(key=lambda x: np.linalg.norm(x[1]))
        return distance_list

    def distort_sphere(self, distort_list, atom_number=0):
        p = self.copy()
        for i, j in zip(self.distance(self.label[atom_number], False)[1:], distort_list):
            temp_position = []
            for k in self.position[self.label.index(i[0])] - self.position[atom_number]:
                if k > 0.5:
                    k -= 1
                elif k < -0.5:
                    k += 1
                k *= 1 + j / i[1]
                if k > 0.5:
                    k -= 1
                elif k < -0.5:
                    k += 1
                temp_position.append(k)
            p.position[self.label.index(i[0])] = np.array(temp_position) + self.position[atom_number]
        return p

    # def distort_randomly(self):
    #     p = poscar()
    #     p.position = [[np.array((random.random() - .5) * 5e-3) for j in i] for i in self.position]
    #     return self + p

    def random_distort(self, displace=0.1):
        p = poscar()
        # displace /= np.sqrt(3)
        p.position = [np.array([(2 * random.random() - 1) * displace for j in i]) for i in self.position]
        if self.mode[0] in ["D", "d"]:
            return self.d2c() + p
        else:
            return self + p

    def atom(self, element):
        aposcar = self.copy()
        aposcar.element = [element]
        aposcar.number = [1]
        aposcar.position = np.array([[0.5, 0.5, 0.5]])
        return aposcar

    # configuration coordinate curves
    def cc(self, vib, ini, fin, step):
        poscar_list = []
        for i in range(step):
            factor = round(ini + i * (fin - ini) / (step - 1), 2)
            poscar_list.append((factor, self + vib * factor))
        return poscar_list

    # def cc2(self):

    #     z2x=self.copy()
    #     z2y=self.copy()
    #     for i in self.position:
    #         z2x.position.append(np.array([i[2],i[0],i[1]]))
    #         z2y.position.append(np.array([i[1],i[2],i[0]]))
    #     temp_pos=[]
    #     for i,j in zip(z2x.position,z2x.label):
    #         temp_pos.append(z2y.position[z2y.label.index(z2y.distance(i,False)[0][0])])
    #     z2y.position=temp_pos
    #     return z2x-z2y

    def avg(self):
        self_avg = {}
        for i, j in zip(self.label, self.position):
            temp_position = []
            for k in range(3):
                if j[k] > 0.9:
                    # temp_position.append(j[k]-0.5)
                    temp_position.append(j[k] - 1)
                # elif j[k]<0.1:
                #     temp_position.append(j[k]+0.5)
                else:
                    temp_position.append(j[k])

            temp_position = np.array(temp_position)
            print(temp_position)
            if i[0] in self_avg:
                self_avg[i[0]] = [self_avg[i[0]][0] + temp_position, self_avg[i[0]][1] + 1]
            else:
                self_avg[i[0]] = [temp_position, 1]

        return self_avg

    # def mimic(self,ref_poscar):

    def cc_xtick(self, mode, num):
        distance_list = []
        for i in self.distance(self.label[0], False)[1:]:
            if len(distance_list) < num:
                distance_list.append(i[1])

        avg = sum(distance_list) / num
        if mode == "avg":
            return avg
        else:
            if distance_list[0] - avg > -1e-3:
                return 0
            elif (distance_list[-1] - avg) / (distance_list[0] - avg) < -1:
                return distance_list[-1] - avg
            else:
                return distance_list[0] - avg

    def env(self, site, prec, number=True):
        distance_list = self.distance(site, False)
        env_dict = {}
        if number:
            for i in distance_list:
                i[1] = round(round(i[1] / prec) * prec, 3)
                if i[1] in env_dict:
                    if i[0][0] in env_dict[i[1]]:
                        env_dict[i[1]][i[0][0]] += 1
                    else:
                        env_dict[i[1]][i[0][0]] = 1
                else:
                    env_dict[i[1]] = {i[0][0]: 1}
        else:
            for i in distance_list:
                i[1] = round(round(i[1] / prec) * prec, 3)
                if i[1] in env_dict:
                    env_dict[i[1]].append(i[0])
                else:
                    env_dict[i[1]] = [i[0]]

        # prev = False
        # united = {}
        # for i in env_dict:
        #     if prev and (i - prev) < prec:
        #         if prev in united:
        #             united[i] = united[prev]
        #             united[i].append(i)
        #             del united[prev]
        #         else:
        #             united[i] = [prev, i]
        #     prev = i
        # # print(united)
        # for i in united:
        #     temp_dict = {}
        #     for j in united[i]:
        #         for k in env_dict[j]:
        #             if k in temp_dict:
        #                 temp_dict[k] += env_dict[j][k]
        #             else:
        #                 temp_dict[k] = env_dict[j][k]
        #         del env_dict[j]
        #     env_dict[round(round(sum(united[i]) / len(united[i]) / prec) * prec, 5)] = temp_dict
        # # return env_dict
        return [[i, set(env_dict[i].items())] for i in sorted(env_dict.keys())] if number else env_dict

    def sites(self, site=False, nei=3, prec=0.01):
        site_list = []
        for i in self.label:
            if site and not site == i[0]:
                continue
            site_env = self.env(i, prec)[:nei]
            if site_list:
                site_saved = False
                for j in site_list:
                    if j[0][0][0] == i[0] and self.same_env(j[1], site_env, prec):
                        j[0].append(i)
                        site_saved = True
                        break
                if not site_saved:
                    site_list.append([[i], site_env])
            else:
                site_list.append([[i], site_env])
        for i in range(len(site_list)):
            site_list[i] = site_list[i][0]
        return site_list

    def same_env(self, a, b, prec):
        issame = True
        for i, j in list(zip(a, b))[:3]:
            if i[1] == j[1]:
                if abs(i[0] - j[0]) < 1.1 * prec:
                    pass
                else:
                    issame = False
                    break
            else:
                issame = False
                break
        return issame

    # def inter(self, new_site, prec=0.1, loop=1000):
    #     inter_list = []

    #     for i in range(loop):
    #         print(i)
    #         # 生成出生点
    #         current_site = np.array([random.random() for _ in range(3)])
    #         count = 0
    #         while True:
    #             neighbours = self.distance(current_site)
    #             nearest_distance=neighbours[0][1]
    #             for neighbour in neighbours[1:]:
    #                 if neighbour[1] - nearest_distance > prec:





    #             # 移动远离最近邻
    #             while True:
    #                 bpos = current_site + 0.02 * random.random() * np.array([random.random() - 0.5 for _ in range(3)])
    #                 if all([1 > xyz > 0 for xyz in bpos]):
    #                     break
    #             b_nearest = self.distance(bpos)[0][1]
    #             if b_nearest > nearest_distance:
    #                 current_site = bpos
    #                 count = 0
    #             else:
    #                 count += 1
    #             if count > 100:
    #                 break

    #         site_env = self.env(current_site, prec)[:1]
    #         # 是否是重复的格位
    #         if inter_list:
    #             site_saved = False
    #             for j in inter_list:
    #                 if self.same_env(j[0], site_env, prec):
    #                     site_saved = True
    #                     j[1] += 1
    #                     break
    #             if not site_saved:
    #                 aposcar = self.copy()
    #                 aposcar.element.append(new_site)
    #                 aposcar.number.append(1)
    #                 aposcar.position = np.insert(aposcar.position, len(aposcar.position), current_site, axis=0)
    #                 inter_list.append([site_env, 1, aposcar])
    #         else:
    #             aposcar = self.copy()
    #             aposcar.element.append(new_site)
    #             aposcar.number.append(1)
    #             aposcar.position = np.insert(aposcar.position, len(aposcar.position), current_site, axis=0)
    #             inter_list.append([site_env, 1, aposcar])

    #     inter_list = [i for i in inter_list if i[1] > 10]
    #     inter_list.sort(key=lambda x: x[0][0][0], reverse=True)

    #     for i in inter_list:
    #         print(i[0], i[1])
    #         i[0] = "".join(filter(str.isalnum, f"{i[0]}{i[1]}"))
    #         del i[1]
    #     return inter_list

    def sub(self, origin_site, new_site):
        poscar_list = []
        for i in self.sites(site=origin_site, nei=3):
            if i[0][0] == origin_site:
                p = self.copy()
                origin_index = p.element.index(origin_site)
                if p.number[origin_index] == 1:
                    del p.element[origin_index]
                    del p.number[origin_index]
                else:
                    p.number[origin_index] -= 1
                origin_position = p.position[p.label.index(i[0])]
                p.position = np.delete(p.position, p.label.index(i[0]), 0)
                # del p.position[p.label.index(i[0])]
                if not new_site == "vac":
                    if not new_site in p.element:
                        p.element.insert(0, new_site)
                        p.number.insert(0, 1)
                        p.position = np.insert(p.position, 0, origin_position, axis=0)
                        # p.position.insert(0, origin_position)
                        p.move(0.5 - origin_position)
                    else:
                        new_index = p.element.index(new_site)
                        p.number[new_index] += 1

                        p.position = np.insert(p.position, sum(p.number[:new_index]), origin_position, axis=0)
                        # p.position.insert(sum(p.number[:new_index]), origin_position)
                p.label = [(i, k + 1) for i, j in zip(p.element, p.number) for k in range(j)]
                poscar_list.append(["".join(filter(str.isalnum, f"{i[0]}")), p])
        return poscar_list

    def slab(self, vac_length, multi=5):

        # z方向扩展2倍
        self.number = [multi * i for i in self.number]
        c_length = np.linalg.norm(self.axis[-1])
        ratio = (c_length * multi + vac_length) / c_length
        for i in range(len(self.axis[-1])):
            self.axis[-1][i] *= ratio
        new_pos = []
        for i in range(len(self.position)):
            for k in range(multi):
                new_pos.append(
                    np.array([(self.position[i][j] + k) / ratio if j == 2 else self.position[i][j] for j in range(3)]))
        self.position = new_pos
        return None

    def slab2(self, vac_length):

        self.number = [2 * i for i in self.number]
        c_length = np.linalg.norm(self.axis[-1])
        ratio = (c_length * 2 + vac_length) / c_length
        for i in range(len(self.axis[-1])):
            self.axis[-1][i] *= ratio
        move_list = []
        for i in self.layer[-1]:
            if abs(i[2][0] - i[2][1]) < 0.01:
                move_list.append((i[0], i[1]))

        new_pos = []
        for i, j in zip(self.position, self.label):
            new_pos.append(np.array([i[j] / ratio if j == 2 else i[j] for j in range(3)]))
            if j in move_list:
                new_pos.append(np.array([1 + (i[j] - 1) / ratio if j == 2 else i[j] for j in range(3)]))
            else:
                new_pos.append(np.array([(i[j] + 1) / ratio if j == 2 else i[j] for j in range(3)]))
        self.position = new_pos
        return None

    def expand(self, tm):
        # transform_matrix

        self.axis = np.dot(tm, self.axis)
        # self.label=[(i,k+1) for i,j in zip(self.element,self.number) for k in range(j)]
        for i in range(len(self.number)):
            self.number[i] = 0
        new_positions = []
        for i, j in zip(self.label, self.position):
            for x in range(-3, 4):
                for y in range(-3, 4):
                    for z in range(-3, 4):
                        a_new_position = (j + np.array([x, y, z])).dot(np.linalg.inv(tm))
                        if max(a_new_position) < 1 and min(a_new_position) >= 0:
                            new_positions.append(a_new_position)
                            self.number[self.element.index(i[0])] += 1
        self.position = new_positions
        return None

    def ecp(self, center_atom, atom_number):
        # for crystalprep
        # read xyz file

        # # 以xyz坐标范围的平均值为中心
        # resize_position=list(zip(*self.position))
        #
        # # 找到符合条件的中心
        # for i in self.distance(np.array([(max(i)+min(i))*.5 for i in resize_position])):
        #     if i[0][0]==center_atom:
        #         center_atom=i[0]
        #         break
        #
        # # 初始化各个分层的原子列表
        # atom_list=[[] for i in range(len(atom_number))]
        #
        # centered_env=self.env(center_atom,prec=0.1,number=False)
        # for distance in centered_env:
        #     for j in range(len(atom_list)):
        #         if len(atom_list[j])>=atom_number[j]:
        #             continue
        #         else:
        #             for label in centered_env[distance]:
        #                 atom_list[j].append(self.label.index(label))
        #             break
        #
        # # 排序并输出
        # for i in range(len(atom_list)):
        #     print(len(atom_list[i]))
        #     atom_list[i].sort()
        #     print(", ".join([str(j) for j in atom_list[i]]))
        # return atom_list

        # 直接输出到inp

        aposcar=self.c2d()
        aposcar.move(np.array([0.5, 0.5, 0.5]) - aposcar.position[aposcar.label.index((center_atom, 1))])
        aposcar=aposcar.d2c()
        centered_env = aposcar.env((center_atom, 1), prec=0.1, number=False)
        actual_number = False
        sphere_position = []
        sphere_element = []
        sphere_number = []

        for radius in centered_env:
            for atom_label in centered_env[radius]:
                sphere_position.append(aposcar.position[aposcar.label.index(atom_label)])
                if not sphere_element or not atom_label[0] == sphere_element[-1]:
                    sphere_element.append(atom_label[0])
                    sphere_number.append(1)
                else:
                    sphere_number[-1] += 1

            if not actual_number and len(sphere_position) > atom_number:
                actual_number = len(sphere_position)
                actual_radius = radius
                break
        aposcar.position = np.array(sphere_position)
        aposcar.element = sphere_element
        aposcar.number = sphere_number

        return aposcar, actual_number, actual_radius


def local_inter(poscar_path):
    if os.path.isdir(f"inter1{poscar_path}"):
        pass
    else:
        p = poscar()
        p.read(poscar_path)
        os.mkdir(f"inter1{poscar_path}")
        for i in p.inter("inter_atom"):
            i[1].write(f"inter1{poscar_path}/{i[0]}")
    return None


def make_inter_file(poscar_path, inter_atom_list):
    p = poscar()
    if not os.path.isdir(f"inter2{poscar_path}"):
        os.mkdir(f"inter2{poscar_path}")
    for i in inter_atom_list:
        if os.path.isdir(f"inter2{poscar_path}/i{i}"):
            break
        else:
            os.mkdir(f"inter2{poscar_path}/i{i}")
            for j in os.listdir(f"inter1{poscar_path}"):
                p.read(f"inter1{poscar_path}/{j}")
                p.element[p.element.index("inter_atom")] = i
                os.mkdir(f"inter2{poscar_path}/i{i}/{j}")
                p.write(f"inter2{poscar_path}/i{i}/{j}/CONTCAR")

    return None


if __name__ == '__main__':
    # p = poscar()
    # p.read("CONTCAR_Hf")
    # print(p.sites(site="Cl", nei=2))
    # print(p.sub("Cl","vac"))
    # print(p.inter("Cl",prec=0.2))
    # local_inter("CONTCAR_cs3bicl6")
    # make_inter_file("CONTCAR_cs3bicl6", ["Cs","Bi","Cl"])
    # local_inter("CONTCAR_cs3bi2cl9")
    # make_inter_file("CONTCAR_cs3bi2cl9", ["Cs","Bi","Cl"])
    # local_inter("CONTCAR_cs2nabicl6")
    # make_inter_file("CONTCAR_cs2nabicl6", ["Cs","Na","Bi","Cl"])
    # local_inter("CONTCAR_small_cs2nabicl6")
    # make_inter_file("CONTCAR_small_cs2nabicl6", ["Cs","Na","Bi","Cl"])
    # local_inter("POSCAR678")

    a = poscar()
    path = r"e:\Desktop\orca\Cs2NaYCl6\222\Cs2NaYCl6_666.vasp"
    a.read(path)
    # print(a.chg)
    sphere_poscar, actual_number, actual_radius = a.ecp("Y", 7)

    sphere_poscar.write(r"e:\Desktop\orca\Cs2NaYCl6\222\Cs2NaYCl6_666_sphere.vasp")
    print(actual_number)
    print(actual_radius)

    # __sub__ test
    # a1 = poscar()
    # a2 = poscar()
    # a1.read("CONTCAR1")
    # a2.read("CONTCAR2")
    # print(a1.position[0])
    # print(a2.position[0])
    # print((a1 - a2).position[0])
