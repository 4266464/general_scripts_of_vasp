import copy
import os

import json_tricks
import numpy as np
import random
import json

# set numpy print options
np.set_printoptions(threshold=np.inf)


# remove \n, split with " ", and filter empty string
def split_and_filter(string):
    return list(filter(lambda x: x != "", string.strip().split(" ")))


class Poscar:
    def __init__(self, poscar_path=None, json_path=None):
        if poscar_path:
            self.read(poscar_path)
        elif json_path:
            self.load(json_path)

    def read(self, poscar_path):
        assert type(poscar_path) == str and os.path.isfile(poscar_path), \
            "The parameter should be path of POSCAR-class file."

        '''
        Ca4 O4                   # comment
        1.0                      # scaling
        9.678532 0.0 0.0         # lattice vector 1
        0.0 9.678532 0.0         # lattice vector 2
        0.0 0.0 9.678532         # lattice vector 3
        Ca                       # element list
        1                        # element number list
        Direct# or Cartesian     # coordinate type
        0.5 0.5 0.5              # coordinate list
        ...
        '''

        with open(poscar_path, "r") as poscar:
            self.comment = poscar.readline().strip()

            scaling = float(poscar.readline().strip())
            lattice = []
            for i in range(3):
                lattice.append(split_and_filter(poscar.readline()))
            lattice = np.array(lattice)
            self.lattice = lattice.reshape(3, 3).astype(float) * scaling
            self.reciprocal = np.linalg.inv(self.lattice)

            self.element = np.array(split_and_filter(poscar.readline()))
            self.number = np.array(split_and_filter(poscar.readline())).astype(int)
            # label for every atom (element of the atom, index of the atom in the element list)
            self.label = np.array(list(zip(np.repeat(self.element, self.number),
                                           sum([list(range(_number)) for _number in self.number], []))))
            self.name = "".join([_element + ("" if _number == 1 else str(_number))
                                 for _element, _number in zip(self.element, self.number)])

            if self.name not in self.comment:
                self.comment += " " + self.name

            # read selective dynamics and type = Direct or Cartesian
            selective = poscar.readline().strip()
            if selective[0] in ["S", "s"]:
                self.selective = True
                self.coordinate = poscar.readline().strip()
            else:
                self.selective = False
                self.coordinate = selective
            self.coordinate = self.coordinate[0].upper()
            assert self.coordinate in ["D", "C"]

            # read atomic position
            position = []
            addition = []

            while True:
                line = poscar.readline().strip()
                if line:
                    line = list(filter(None, line.split(" ")))
                    position.append(line[:3])
                    addition.append(" ".join(line[3:]))
                else:
                    break

            self.position = np.array(position).astype(float)
            self.addition = np.array(addition) if addition[0] else None

    def dump(self, json_path):
        with open(json_path, "w") as json_file:
            json_tricks.dump(self.__dict__, json_file)

    def load(self, json_path):
        with open(json_path, "r") as json_file:
            self.__dict__ = json_tricks.load(json_file)

    def write(self, poscar_path):
        poscar_string = self.comment + "\n1.\n"
        poscar_string += self.lattice.__str__().replace("[", "").replace("]", "") + "\n"
        poscar_string += " ".join(self.element) + "\n"
        poscar_string += " ".join(np.array(self.number).astype("str")) + "\n"

        if self.selective:
            poscar_string += "S\n"
        poscar_string += f"{self.coordinate}\n"

        if self.addition is not None:
            # the shape of self.position is (natom, 3), and self.addition is (natom,)
            for atom_str in np.concatenate((self.position, self.addition[:, np.newaxis]), axis=1):
                poscar_string += " ".join(atom_str) + "\n"
            poscar_string += "\n"
        else:
            for atom_str in self.position.astype("str"):
                poscar_string += " ".join(atom_str) + "\n"
            poscar_string += "\n"

        with open(poscar_path, "w") as f:
            f.write(poscar_string)

    # make a copy from self
    def copy(self):
        return copy.deepcopy(self)

    # covert Direct coordinate to Cartesian coordinate
    def d2c(self):
        if self.coordinate == "D":
            self.position = self.position.dot(self.lattice)
            self.coordinate = "C"

    # covert Cartesian coordinate to Direct coordinate
    def c2d(self):
        if self.coordinate == "C":
            self.position = self.position.dot(self.reciprocal)
            self.coordinate = "D"

    # move the atoms and make their positions between 0 and 1 or -0.5 and 0.5
    def std(self, diff=False):
        # self.coordinate = "D" required
        final_d2c = self.coordinate == "C"
        self.c2d()
        self.position = np.where(self.position >= (0.5 if diff else 1), self.position - 1, self.position)
        self.position = np.where(self.position < (-0.5 if diff else 0), self.position + 1, self.position)

        if final_d2c:
            self.d2c()

    def move(self, vector, diff=False):
        assert isinstance(vector, (np.ndarray, list, tuple)), \
            f"the type of vector should be np.ndarray, list or tuple, but {type(vector)} is given."
        vector = np.array(vector)
        assert vector.shape == (3,), \
            f"the shape of the array should be (3,), but {vector.shape} is given."

        poscar_copy = self.copy()
        poscar_copy.position += vector
        poscar_copy.std(diff)

        return poscar_copy

    def __add__(self, poscar):
        assert isinstance(poscar, Poscar), "the parameter should be a Poscar-class."
        poscar_copy = self.copy()
        assert np.allclose(poscar_copy.lattice, poscar.lattice), "the lattice of poscars are not match."
        poscar_copy.position += poscar.position
        poscar_copy.std()

        return poscar_copy

    def __sub__(self, poscar):
        assert isinstance(poscar, Poscar), "the parameter should be a Poscar-class object."
        poscar_copy = self.copy()
        assert np.allclose(poscar_copy.lattice, poscar.lattice), "the lattice of poscars are not match."
        poscar_copy.position -= poscar.position
        poscar_copy.std(diff=True)

        return poscar_copy

    # only for diff
    def __mul__(self, scale):
        assert isinstance(scale, (int, float)), "the parameter should be a real number"
        poscar_copy = self.copy()
        poscar_copy.position *= scale

        return poscar_copy

    def __rmul__(self, scale):
        return self * scale

    def distance(self, position):
        assert isinstance(position, (list, tuple, np.ndarray)), \
            f"parameter should be a coordinate, but {type(position)} is given."
        position = np.array(position)
        assert position.shape == (3,), \
            f"the shape of the array should be (3,), but {position.shape} is given."

        # set position as (0, 0, 0) and calculate the distance from the position
        distance = np.linalg.norm(self.move(-position, True).position.dot(self.lattice), axis=1)
        # sort by distance
        argsort = distance.argsort()
        distance = distance[argsort]
        label = self.label[argsort]
        return distance, label

    # TODO: need tags
    def environment(self, site, prec=0.1):
        env_tuple = namedtuple(
            "env_tuple", ["distance", "counter", "specific_site"])
        env_list = []
        for atom in self.distance(site):
            distance = round(round(atom.distance_from_point_to_line / prec) * prec, 3)
            for env in env_list:
                if distance == env.distance:
                    if atom.label[0] in env.counter:
                        env.counter[atom.label[0]] += 1
                    else:
                        env.counter[atom.label[0]] = 1
                    env.specific_site.append(atom.label)
                    break
            else:
                env_list.append(env_tuple(distance=distance, counter={
                    atom.label[0]: 1}, specific_site=[atom.label]))

        return env_list

    # TODO: need cluster
    def same_site(self, env1, env2):
        for aenv1, aenv2 in list(zip(env1, env2)):
            if aenv1.distance_from_point_to_line == aenv2.distance_from_point_to_line and aenv1.counter == aenv2.counter:
                continue
            else:
                return False
        return True

    def sites(self, element, nei=3, prec=0.1):
        site_tuple = namedtuple("site_tuple", ["env", "specific_site"])
        site_list = []
        for label in self.label:
            if element == label[0]:
                env = self.environment(label, prec)[:nei]
                for site in site_list:
                    if self.same_site(env, site.environment):
                        site.specific_site.append(label)
                        break
                else:
                    site_list.append(site_tuple(env=env, specific_site=[label]))

        return site_list

    def distort_sphere(self, distort_list=[0.1], site=0):
        aposcar = Poscar(self)
        site = self.every2pos(site)
        for atom, distort in zip(self.distance(site)[1:], distort_list):
            aposcar.position[self.label.index(atom.label)] = atom.position.dot(self.reciprocal) * (
                    1 + distort / atom.distance_from_point_to_line) + site
        return aposcar

    def random_distort(self, displace=0.1, site=0, nei=3):
        aposcar = Poscar(self)
        aposcar.position = np.zeros((len(aposcar.position), 3))
        aposcar.coordinate = "C"
        env = self.environment(site, 0.1)
        for neighbours in env[:nei]:
            for label in neighbours.specific_site:
                vector = np.array([random.random() for _ in range(3)]) - .5
                vector /= np.linalg.norm(vector)
                aposcar.position[self.label.index(label)] = displace * random.random() * vector

        return self + aposcar.c2d()

    def atom(self, element):
        aposcar = Poscar(self)
        aposcar.element = [element]
        aposcar.number = [1]
        aposcar.position = np.array([[0.5, 0.5, 0.5]])
        return aposcar

    # configuration coordinate curves
    def cc(self, other, ini, fin, step):
        poscar_list = []
        for i in range(step):
            factor = round(ini + i * (fin - ini) / (step - 1), 2)
            poscar_list.append([factor, self + other * factor])
        return poscar_list

    # under construct
    def inter(self, prec=0.1, loop=1000):
        inter_tuple = namedtuple("inter_tuple", ["env", "counter", "poscar"])
        inter_list = []

        for i in range(loop):
            # 生成出生点
            current_site = np.array([random.random() for _ in range(3)])
            count = 1
            neighbours = self.distance(current_site)

            while True:
                # 移动远离最近邻
                new_site = current_site + 0.01 / (count / 100 if count > 100 else 1) * (
                        np.array([random.random() for _ in range(3)]) - .5)
                for i in range(3):
                    if new_site[i] < 0:
                        new_site[i] += 1
                    elif new_site[i] >= 1:
                        new_site[i] -= 1

                new_neighbours = self.distance(new_site)

                if new_neighbours[0].distance_from_point_to_line > neighbours[0].distance_from_point_to_line:
                    current_site = new_site
                    count = 1
                    neighbours = new_neighbours
                else:
                    count += 1

                if count > 500:
                    break

            current_env = self.environment(current_site, prec)[:2]
            # 是否是重复的格位

            for inter_site in inter_list:
                if self.same_site(inter_site.environment, current_env):
                    inter_site.counter[0] += 1
                    break
            else:
                aposcar = Poscar(self)
                aposcar.element.insert(0, "i")
                aposcar.number.insert(0, 1)
                aposcar.position = np.insert(aposcar.position, 0, current_site, axis=0)
                aposcar.tag.insert(0, [])
                print(current_site, current_env[0])
                inter_list.append(inter_tuple(env=current_env, counter=[1], poscar=aposcar))
        return inter_list

    def sub(self, origin_site, new_site):
        sub_tuple = namedtuple("sub_tuple", ["site", "poscar"])
        sub_list = []
        for site in self.sites(element=origin_site, nei=3):
            aposcar = Poscar(self)
            origin_element_index = aposcar.element.index(origin_site)
            if aposcar.number[origin_element_index] == 1:
                del aposcar.element[origin_element_index]
                del aposcar.number[origin_element_index]
            else:
                aposcar.number[origin_element_index] -= 1
            index = aposcar.label.index(site.specific_site[0])
            origin_position = aposcar.position[index]
            aposcar.position = np.delete(aposcar.position, index, 0)
            del aposcar.tag[index]
            if not new_site == "vac":
                aposcar.element.insert(0, new_site)
                aposcar.number.insert(0, 1)
                aposcar.position = np.insert(aposcar.position, 0, origin_position, axis=0)
                aposcar.tag.insert(0, [])
                aposcar = aposcar + (0.5 - origin_position)

            aposcar.label = []
            for element, element_number in zip(aposcar.element, aposcar.number):
                for number in range(element_number):
                    # in case element is separated
                    while (element, number) in aposcar.label:
                        number += 1
                    aposcar.label.append((element, number))
            sub_list.append(sub_tuple(site=site, poscar=aposcar))
        return sub_list

    # under constuct
    def slab(self, vac_length, multi=5):
        # z方向扩展2倍
        self.number = [multi * i for i in self.number]
        c_length = np.linalg.norm(self.lattice[-1])
        ratio = (c_length * multi + vac_length) / c_length
        for i in range(len(self.lattice[-1])):
            self.lattice[-1][i] *= ratio
        new_pos = []
        for i in range(len(self.position)):
            for k in range(multi):
                new_pos.append(np.array([(self.position[i][j] + k) / ratio if j == 2 else self.position[i][j] for j in
                                         range(3)]))
        self.position = new_pos
        return None

    def slab2(self, vac_length):
        self.number = [2 * i for i in self.number]
        c_length = np.linalg.norm(self.lattice[-1])
        ratio = (c_length * 2 + vac_length) / c_length
        for i in range(len(self.lattice[-1])):
            self.lattice[-1][i] *= ratio
        move_list = []
        for i in self.layer[-1]:
            if abs(i[2][0] - i[2][1]) < 0.01:
                move_list.append((i[0], i[1]))

        new_pos = []
        for i, j in zip(self.position, self.label):
            new_pos.append(
                np.array([i[j] / ratio if j == 2 else i[j] for j in range(3)]))
            if j in move_list:
                new_pos.append(np.array([1 + (i[j] - 1) / ratio if j == 2 else i[j] for j in range(3)]))
            else:
                new_pos.append(np.array([(i[j] + 1) / ratio if j == 2 else i[j] for j in range(3)]))
        self.position = new_pos
        return None

    def expand(self, tm):
        # transform_matrix

        self.lattice = np.dot(tm, self.lattice)
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
        # 直接输出到inp
        aposcar = self.c2d()
        aposcar = aposcar + \
                  (0.5 - aposcar.position[aposcar.label.index((center_atom, 1))])
        aposcar = aposcar.d2c()
        centered_env = aposcar.environment((center_atom, 1), prec=0.1)
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


if __name__ == '__main__':
    p = Poscar(poscar_path="poscar_test.vasp")
    # test read and write, dump and load
    # p.write("poscar_write_test")
    # p_write = Poscar(poscar_path="poscar_write_test")
    # for i in p.__dict__:
    #     if isinstance(p.__dict__[i], np.ndarray):
    #         assert np.all(p.__dict__[i] == p_write.__dict__[i]), f"{i} not match."
    #     elif not p.__dict__[i] == p_write.__dict__[i]:
    #         f"{i} not match, which are {p.__dict__[i]} and {p_write.__dict__[i]}."
    #
    # p.dump("poscar_dump_test")
    # p_load = Poscar(json_path="poscar_dump_test")
    # for i in p.__dict__:
    #     if isinstance(p.__dict__[i], np.ndarray):
    #         assert np.all(p.__dict__[i] == p_load.__dict__[i]), f"{i} not match."
    #     elif not p.__dict__[i] == p_load.__dict__[i]:
    #         f"{i} not match, which are {p.__dict__[i]} and {p_load.__dict__[i]}."

    # test move
    # p_move = p.move([0.1, 0.2, 0.3])
    # p_move.write("poscar_move_test.vasp")

    # test +, - and *
    # shape = p.position.shape
    # random_matrix = np.random.rand(*shape)
    # p_random = p.copy()
    # p_random.position = random_matrix
    # p_add = p + p_random
    # print(p_add.position - (p.position + p_random.position))
    #
    # p_sub = p_add - p
    # print(p_sub.position - (p.position - p_random.position))
    #
    # p_mul = p_sub * 1.3
    # print(p_mul.position - (random_matrix * 1.3))

    # test distance
    distance, label = p.distance(p.position[0])
    print(distance)
    print(label)

    # print((p + [0.4, 0, 0]).__dict__)

    # print(p.position[:5])
    # print(p.c2d().position[:5])
    # print(p.d2c().position[:5])
    # print(p.d2c().c2d().position[:5])
    # print(p.distance(p.label[0]))
    # print(p.distance([0, 0, 0]))
    # print(p.sites(site="Cl", nei=2))
    # print(p.sub("Cl","vac"))
    # print(p.inter("Cl",prec=0.2))

    # __sub__ test
    # a1 = Poscar("CONTCAR1")
    # a2 = Poscar("CONTCAR2")
    # print(a1.position[0])
    # print(a2.position[0])
    # print((a1 - a2).position[0])
