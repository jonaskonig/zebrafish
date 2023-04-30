import datetime
import os
import random
import copy
import numpy as np
from PIL import Image, ImageOps


class Chromatophores:

    def __init__(self, chromatophores_type, position):
        self.possible = ["M", "VL", "VR", "HT", "HB"]
        if position not in self.possible:
            raise NotImplementedError
        # xant 0 Men 1
        self.chromatophores_type = chromatophores_type
        self.age = 0
        self.__position = position


    def age_cell(self):
        self.age += 1

    def set_position(self, position):
        if position not in self.possible:
            raise NotImplementedError
        self.__position = position

    def get_position(self):
        return self.__position


class Lattice:

    def __init__(self, size: list, stem_cell_density, elastic_energy_weight: float, adhesion_coe_matrix: list,
                 mean_cell_life: int, lattice_groth_interval: int, runtime: int, y_growth_factor:int|bool =False,chromatophore_birth_rate:int = 3, paper_move_method:bool = False):
        self.size = size
        # self.lattice = np.full((size[0], size[1]), None)
        self.lattice = [[None] * size[1]] * size[0]
        self.stem_cell_density = stem_cell_density
        self.elastic_energy_weight = elastic_energy_weight
        self.adhesion_coe_matrix = adhesion_coe_matrix
        self.mean_cell_life = mean_cell_life
        self.chromatophore_birth_rate = chromatophore_birth_rate
        self.lattice_groth_interval = lattice_groth_interval
        self.runtime = runtime
        self.create_start_stripes()
        self.y_growth_factor =y_growth_factor
        self.grow_down = True
        self.grow_right = True
        self.paper_move_method = paper_move_method


    def create_start_stripes(self):
        stripe = []
        for x in range(len(self.lattice[0])):
            stripe.append([Chromatophores(1, "M")])
        # stripe = np.asarray(stripe)
        s1 = int(0.15 * len(self.lattice))
        s2 = int(0.4 * len(self.lattice))
        s3 = int(0.6 * len(self.lattice))
        s4 = int(0.85 * len(self.lattice))
        self.lattice = np.asarray(self.lattice)
        self.lattice[s1 - 1] = copy.deepcopy(stripe)
        self.lattice[s1] = copy.deepcopy(stripe)
       # self.lattice[s1 + 1] = copy.deepcopy(stripe)
        self.lattice[s2 - 1] = copy.deepcopy(stripe)
        self.lattice[s2] = copy.deepcopy(stripe)
        #self.lattice[s2 + 1] = copy.deepcopy(stripe)
        #self.lattice[s3 - 1] = copy.deepcopy(stripe)
        self.lattice[s3] = copy.deepcopy(stripe)
        self.lattice[s3 + 1] = copy.deepcopy(stripe)
        #self.lattice[s4 - 1] = copy.deepcopy(stripe)
        self.lattice[s4] = copy.deepcopy(stripe)
        self.lattice[s4 + 1] = copy.deepcopy(stripe)

    def get_rand_position(self):
        possible = {0: "VL", 1: "VR", 2: "HT", 3: "HB"}
        randnum = random.randint(0,3)
        if randnum == 0:
            return possible[0], possible[1]
        if randnum == 1:
            return possible[1], possible[0]
        if randnum == 2:
            return possible[2], possible[3]
        if randnum == 3:
            return possible[3], possible[2]

    def distributeStemCells(self):
        dense = int(round(self.size[0] * self.size[1] * self.stem_cell_density))#*self.chromatophore_birth_rate
        for k in range(dense):
            posx = random.randint(0,len(self.lattice)-1)
            posy = random.randint(0,len(self.lattice[0])-1)
            if self.lattice[posx][posy] is None:
                if self.get_neighbor_count(posx, posy) > 2:
                    self.lattice[posx][posy] = [Chromatophores(0, "M")]
                else:
                    self.lattice[posx][posy] = [Chromatophores(1, "M")]
            elif len(self.lattice[posx][posy]) < 2:
                pos1, pos2 = self.get_rand_position()
                if self.get_neighbor_count(posx, posy, pos2) + 1 > 2:
                    if self.lattice[posx][posy][0].chromatophores_type == 0:
                        self.lattice[posx][posy][0].set_position(pos1)
                        self.lattice[posx][posy].append(Chromatophores(0, pos2))
                        self.lattice[posx][posy] = self.order_cells(self.lattice[posx][posy])
                else:
                    if self.lattice[posx][posy][0].chromatophores_type == 0:
                        self.lattice[posx][posy][0].set_position(pos1)
                        self.lattice[posx][posy].append(Chromatophores(1, pos2))
                        self.lattice[posx][posy] = self.order_cells(self.lattice[posx][posy])
                    else:
                        self.lattice[posx][posy][0].set_position(pos1)
                        self.lattice[posx][posy].append(Chromatophores(0, pos2))
                        self.lattice[posx][posy] = self.order_cells(self.lattice[posx][posy])

    def order_cells(self, cells: list) -> list:
        if len(cells) == 1:
            cells[0].set_position("M")
            return [cells[0]]
        if cells[0].get_position == "VL" or cells[0].get_position == "HT":
            return cells
        else:
            return [cells[1], cells[0]]

    def age_cells(self):
        for x in self.lattice:
            for y in x:
                if y:
                    y[0].age_cell()
                    if len(y) > 1:
                        y[1].age_cell()

    def cell_rearangement(self):
        for x in range(len(self.lattice)):
            for y in range(len(self.lattice[0])):
                if self.lattice[x][y] is None:
                    continue
                if len(self.lattice[x][y]) > 1:
                    number = random.randint(0, 1)
                    self.cell_deth_and_move(x, y, number)
                else:
                    self.cell_deth_and_move(x, y, 0)

    def grow_lattice(self, round):
        ap = np.full(len(self.lattice[0]), None)
        if round < 0.45:
            insertposition = random.randint(int((len(self.lattice)-1)*0.3),int((len(self.lattice)-1)*0.7))
        elif random.randint(0,1) == 0:
            insertposition = random.randint(0, int((len(self.lattice)-1)*0.3))
        else:
            insertposition = random.randint(int((len(self.lattice) - 1) * 0.7),len(self.lattice) )
        #if self.insertposition:
        #    self.lattice = np.append(self.lattice, [ap], axis=0)
        self.lattice = np.insert(self.lattice,insertposition,[ap],axis=0)
        #else:
        #    self.lattice = np.concatenate(([ap],self.lattice), axis=0)
        #self.grow_down = not self.grow_down
        self.size[0] += 1

    def grow_lattice_x(self,round):
        #ap = np.full((len(self.lattice),1), None)
        if round < 0.45:
            insertposition = random.randint(int((len(self.lattice[0]) - 1) * 0.3), int((len(self.lattice[0]) - 1) * 0.7))
        elif random.randint(0, 1) == 0:
            insertposition = random.randint(0, int((len(self.lattice[0]) - 1) * 0.3))
        else:
            insertposition = random.randint(int((len(self.lattice) - 1) * 0.7), len(self.lattice[0]))

        self.lattice = np.insert(self.lattice, insertposition, None, axis=1)
        self.size[1] += 1

    def cell_deth_and_move(self, x, y, num):
        vs, vd = self.get_cellneigbors_vd_vs(x, y, num, self.lattice)
        if self.lattice[x][y][num].age > self.mean_cell_life and vs > vd:
            #print(f"Cell age is: {self.lattice[x][y][num].age}")
            if len(self.lattice[x][y]) > 1:
                self.lattice[x][y].pop(num)
                self.lattice[x][y][0].set_position("M")
            else:
                self.lattice[x][y] = None
        else:
            dir = random.randint(0,3)
            eb = 0
            for n in range(vs):
                eb += self.adhesion_coe_matrix[self.lattice[x][y][num].chromatophores_type] + len(
                    self.lattice[x][y]) * self.elastic_energy_weight
            for n in range(vd):
                eb += self.adhesion_coe_matrix[2] + len(self.lattice[x][y]) * self.elastic_energy_weight
            newlatice, newnum, newx, newy = self.move_cell(self.lattice, x, y, num, dir, orgi=self.paper_move_method)
            if len(newlatice) > 0:
                ea = 0
                vs, vd = self.get_cellneigbors_vd_vs(newx, newy, newnum, newlatice)
                for n in range(vs):
                    ea += self.adhesion_coe_matrix[newlatice[newx][newy][newnum].chromatophores_type] + len(
                        newlatice[newx][newy]) * self.elastic_energy_weight
                for n in range(vd):
                    ea += self.adhesion_coe_matrix[2] + len(newlatice[newx][newy]) * self.elastic_energy_weight
                if ea > eb:
                    self.lattice = newlatice


    def get_new_xy(self, x, y, dir):
        if dir == 0:
            return x + 1, y
        if dir == 1:
            return x, y + 1
        if dir == 2:
            return x - 1, y
        if dir == 3:
            return x, y - 1

    def move_cell(self, lattice, x, y, num, dir, orgi=True):
        newx = 0
        newy = 0
        num2 = 0 if num == 1 else 1
        if len(lattice[x][y]) == 1 and lattice[x][y][0].get_position() != "M":
            lattice[x][y][0].set_position("M")
        if lattice[x][y][num].get_position() == "VR" and dir == 3:
            lattice[x][y] = [lattice[x][y][num], lattice[x][y][num2]]
            lattice[x][y][0].set_position("VL")
            lattice[x][y][1].set_position("VR")
            lattice[x][y] = self.order_cells(lattice[x][y])
            return lattice, 0, x, y
        if lattice[x][y][num].get_position() == "VL" and dir == 1:
            lattice[x][y] = [lattice[x][y][num2], lattice[x][y][num]]
            lattice[x][y][0].set_position("VL")
            lattice[x][y][1].set_position("VR")
            lattice[x][y] = self.order_cells(lattice[x][y])
            return lattice, 1, x, y
        if lattice[x][y][num].get_position() == "HT" and dir == 2:
            lattice[x][y] = [lattice[x][y][num2], lattice[x][y][num]]
            lattice[x][y][1].set_position("HB")
            lattice[x][y][0].set_position("HT")
            lattice[x][y] = self.order_cells(lattice[x][y])
            return lattice, 1, x, y
        if lattice[x][y][num].get_position() == "HB" and dir == 1:
            lattice[x][y] = [lattice[x][y][num], lattice[x][y][num2]]
            lattice[x][y][0].set_position("HT")
            lattice[x][y][1].set_position("HB")
            lattice[x][y] = self.order_cells(lattice[x][y])
            return lattice, 0, x, y
        if dir == 0:
            if x + 1 < len(lattice):
                newx = x + 1
            else:
                return [], 0, 0, 0
        elif dir == 2:
            if not (x - 1 < 0):
                newx = x - 1
            else:
                return [], 0, 0, 0
        elif dir == 1:
            if y + 1 < len(lattice[0]):
                newy = y + 1
            else:
                return [], 0, 0, 0
        elif dir == 3:
            if not (y - 1 < 0):
                newy = y - 1
            else:
                return [], 0, 0, 0
        else:
            raise Exception(f"Direction out of bound for dir = {dir}")
        if lattice[newx][newy] is None:
            lattice[newx][newy] = [lattice[x][y][num]]
            lattice[newx][newy][0].set_position("M")
            if len(lattice[x][y]) > 1:
                lattice[x][y].pop(num)
                lattice[x][y][0].set_position("M")
            else:
                lattice[x][y] = None
            return lattice, 0, newx, newy
        if len(lattice[newx][newy]) == 1:
            pos1, pos2 = self.get_rand_position()
            lattice[newx][newy][0].set_position(pos1)
            lattice[newx][newy].append(lattice[x][y][num])
            lattice[newx][newy][1].set_position(pos2)
            lattice[newx][newy] = self.order_cells(lattice[newx][newy])
            if len(lattice[x][y]) > 1:
                lattice[x][y].pop(num)
                lattice[x][y][0].set_position("M")
            else:
                lattice[x][y] = None
            return lattice, 0, newx, newy
        else:
            if orgi:
                return [], 0, 0,0
            if (lattice[x][y][num].get_position() == "VL" or lattice[x][y][num].get_position() == "VR") and \
                    lattice[newx][newy][0] == "VL":
                # should not happen with pos = VR
                if dir == 3:
                    mov = lattice[x][y][num]
                    lattice[x][y][num] = lattice[newx][newy][1]
                    lattice[x][y][num].set_position("VL")
                    lattice[newx][newy][1] = mov
                    lattice[newx][newy][1].set_position("VR")
                # should not happen with pos = VL
                if dir == 1:
                    mov = lattice[x][y][num]
                    lattice[x][y][num] = lattice[newx][newy][0]
                    lattice[x][y][num].set_position("VR")
                    lattice[newx][newy][0] = mov
                    lattice[newx][newy][1].set_position("VL")
                if dir == 0 or dir == 2:
                    temp = lattice[x][y][num]
                    lattice[x][y][num] = lattice[newx][newy][num]
                    lattice[newx][newy][num] = temp
                return lattice, num, newx, newy
            if (lattice[x][y][num].get_position() == "HT" or lattice[x][y][num].get_position() == "HB") and \
                    lattice[newx][newy][0] == "HT":
                if dir == 1 or dir == 3:
                    temp = lattice[x][y][num]
                    lattice[x][y][num] = lattice[newx][newy][num]
                    lattice[newx][newy][num] = temp
                if dir == 0:
                    mov = lattice[x][y][num]
                    lattice[x][y][num] = lattice[newx][newy][1]
                    lattice[x][y][num].set_position("HT")
                    lattice[newx][newy][1] = mov
                    lattice[newx][newy][1].set_position("HB")
                if dir == 2:
                    mov = lattice[x][y][num]
                    lattice[x][y][num] = lattice[newx][newy][0]
                    lattice[x][y][num].set_position("HB")
                    lattice[newx][newy][0] = mov
                    lattice[newx][newy][0].set_position("HT")
                return lattice, num, newx, newy
            else:
                which = random.randint(0,1)
                temp = lattice[x][y][num]
                lattice[x][y][num] = lattice[newx][newy][which]
                lattice[x][y][num].set_position(temp.get_position())
                newp = lattice[newx][newy][which].get_position()
                lattice[newx][newy][which] = temp
                lattice[newx][newy][which].set_position(newp)
                return lattice, which, newx, newy

    def get_cellneigbors_vd_vs(self, x, y, num, lattice):
        vd = 0
        vs = 0
        num2 = 0 if num == 1 else 1
        if len(lattice[x][y]) > 1:
            if lattice[x][y][num2].chromatophores_type == lattice[x][y][num].chromatophores_type:
                vs += 1
            else:
                vd += 1
        for pos in range(4):
            if (pos == 0 and lattice[x][y][num].get_position() == "HB") or (
                    pos == 1 and lattice[x][y][num].get_position() == "VL") or (
                    pos == 2 and lattice[x][y][num].get_position() == "HT") or (
                    pos == 3 and lattice[x][y][num].get_position() == "VR"):
                continue
            newx = x
            newy = y
            if dir == 0:
                if x + 1 < len(lattice):
                    newx = x + 1
                else:
                    continue
            if dir == 2:
                if not (x - 1 < 0):
                    newx = x - 1
                else:
                    continue
            if dir == 1:
                if y + 1 < len(lattice[0]):
                    newy = y + 1
                else:
                    continue
            if dir == 3:
                if not (y - 1 < 0):
                    newy = y - 1
                else:
                    continue
            if lattice[newx][newy] is None:
                continue
            if len(lattice[newx][newy]) == 1:
                if pos == 0 and not lattice[x][y][num].get_position() == "HB":
                    if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                        vs += 1
                    else:
                        vd += 1

                if pos == 1 and not lattice[x][y][num].get_position() == "VR":
                    if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                        vs += 1
                    else:
                        vd += 1

                if pos == 2 and not lattice[x][y][num].get_position() == "HT":
                    if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                        vs += 1
                    else:
                        vd += 1

                if pos == 3 and not lattice[x][y][num].get_position() == "VL":
                    if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                        vs += 1
                    else:
                        vd += 1
                continue

            if lattice[x][y][num].get_position() == "VL" or lattice[x][y][num].get_position() == "VR":
                if lattice[newx][newy][0].get_position() == "HT" and pos == 2:
                    if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                        vs += 1
                    else:
                        vd += 1
                if lattice[newx][newy][1].get_position() == "HB" and pos == 0:
                    if lattice[newx][newy][1].chromatophores_type == lattice[x][y][num].chromatophores_type:
                        vs += 1
                    else:
                        vd += 1
                if lattice[x][y][num].get_position() == "VL":
                    if lattice[newx][newy][0].get_position() == "VL" and pos != 3:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                    if lattice[newx][newy][1].get_position() == "VR" and pos == 3:
                        if lattice[newx][newy][1].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                    if lattice[newx][newy][1].get_position() == "HB" and pos == 3:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                        if lattice[newx][newy][1].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                if lattice[x][y][num].get_position() == "VR":
                    if lattice[newx][newy][1].get_position() == "VR" and pos != 1:
                        if lattice[newx][newy][1].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                    if lattice[newx][newy][0].get_position() == "VL" and pos == 1:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                    if lattice[newx][newy][1].get_position() == "HB" and pos == 1:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                        if lattice[newx][newy][1].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
            if lattice[x][y][num].get_position() == "HT" or lattice[x][y][num].get_position() == "HB":
                if lattice[newx][newy][0].get_position() == "VL" and pos == 1:
                    if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                        vs += 1
                    else:
                        vd += 1
                if lattice[newx][newy][1].get_position() == "VR" and pos == 3:
                    if lattice[newx][newy][1].chromatophores_type == lattice[x][y][num].chromatophores_type:
                        vs += 1
                    else:
                        vd += 1
                if lattice[x][y][num].get_position() == "HT":
                    if lattice[newx][newy][1].get_position() == "HB" and pos == 0:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                    if lattice[newx][newy][0].get_position() == "VL" and pos == 0:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                        if lattice[newx][newy][1].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                if lattice[x][y][num].get_position() == "HB":
                    if lattice[newx][newy][0].get_position() == "HT" and pos == 3:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                    if lattice[newx][newy][1].get_position() == "HB" and pos != 3:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                    if lattice[newx][newy][0].get_position() == "VL" and pos == 3:
                        if lattice[newx][newy][0].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
                        if lattice[newx][newy][1].chromatophores_type == lattice[x][y][num].chromatophores_type:
                            vs += 1
                        else:
                            vd += 1
        return vs, vd

    def get_neighbor_count(self, x, y, position="M"):
        count = 0
        if position == "M":
            for pos in range(4):
                if pos == 0 and x + 1 < len(self.lattice):
                    count = count if self.lattice[x + 1][y] is None else count + 1
                    continue
                if pos == 1 and y + 1 < len(self.lattice[0]):
                    count = count if self.lattice[x][y + 1] is None else count + 1
                    continue
                if pos == 2 and x - 1 >= 0:
                    count = count if self.lattice[x - 1][y] is None else count + 1
                    continue
                if pos == 3 and y - 1 >= 0:
                    count = count if self.lattice[x][y - 1] is None else count + 1
                    continue
        else:
            for pos in range(4):
                if pos == 0 and x + 1 < len(self.lattice):
                    if position == "HB":
                        count += 1
                        continue
                    if self.lattice[x + 1][y] is None:
                        continue
                    if ((position == "VL" or position == "VR") and self.lattice[x + 1][y][
                        0].get_position() == "HT") or (
                            position == "HT" and self.lattice[x + 1][y][0].get_position() == "VL"):
                        count += 2
                    else:
                        count += 1
                    continue
                if pos == 1 and y + 1 < len(self.lattice[0]):
                    if position == "VL":
                        count += 1
                        continue
                    if self.lattice[x][y + 1] is None:
                        continue
                    if ((position == "HT" or position == "HB") and self.lattice[x][y + 1][
                        0].get_position() == "VL") or (
                            position == "VR" and self.lattice[x][y + 1][0].get_position() == "HT"):
                        count += 2
                    else:
                        count += 1
                    continue
                if pos == 2 and x - 1 >= 0:
                    if position == "HT":
                        count += 1
                        continue
                    if self.lattice[x - 1][y] is None:
                        continue
                    if ((position == "VL" or position == "VR") and self.lattice[x - 1][y][
                        0].get_position() == "HT") or (
                            position == "HB" and self.lattice[x - 1][y][0].get_position() == "VL"):
                        count += 2
                    else:
                        count += 1
                    continue
                if pos == 3 and y - 1 >= 0:
                    if position == "VR":
                        count += 1
                        continue
                    if self.lattice[x][y - 1] is None:
                        continue
                    if ((position == "HT" or position == "HB") and self.lattice[x][y - 1][
                        0].get_position() == "VL") or (
                            position == "VL" and self.lattice[x][y - 1][0].get_position() == "HT"):
                        count += 2
                    else:
                        count += 1
                    continue
        return count

    def run_sim(self):
        lattices = []
        # lc = copy.deepcopy(self.lattice)
        lattices.append(copy.deepcopy(self.lattice))
        for d in range(self.runtime):
            if d % self.lattice_groth_interval == 0:
                self.grow_lattice(d/self.runtime)
            if self.y_growth_factor and d % self.y_growth_factor == 0:
                self.grow_lattice_x(d/self.runtime)
            if d % self.chromatophore_birth_rate == 0:
                self.distributeStemCells()
                #print(f"D is {d}")
            self.age_cells()
            self.cell_rearangement()
            lattices.append(copy.deepcopy(self.lattice))
        return lattices

    def to_picture(self, lattice):
        #print(len(lattice))
        size = 4
        space = 4
        sizex = (size + space) * len(lattice) + space
        sizey = (size + space) * len(lattice[0]) + space
        pic = np.full((sizex, sizey, 3), [122, 192, 234])
        color = [[240, 234, 122], [64, 64, 59]]
        picx = space
        picy = space
        for x in range(len(lattice)):
            picy = space
            for y in range(len(lattice[0])):
                if lattice[x][y] is None:
                    picy += size + space
                    continue
                if len(lattice[x][y]) == 1:
                    for s1 in range(size):
                        for s2 in range(size):
                            pic[picx + s1][picy + s2] = color[lattice[x][y][0].chromatophores_type]
                else:
                    if lattice[x][y][0].get_position == "VL":
                        for s1 in range(size):
                            for s2 in range(size):
                                if s1 % 2 == 0:
                                    pic[picx + s1][picy + s2] = color[lattice[x][y][0].chromatophores_type]
                                else:
                                    pic[picx + s1][picy + s2] = color[lattice[x][y][1].chromatophores_type]
                    else:
                        for s1 in range(size):
                            for s2 in range(size):
                                if s2 % 2 == 0:
                                    pic[picx + s1][picy + s2] = color[lattice[x][y][0].chromatophores_type]
                                else:
                                    pic[picx + s1][picy + s2] = color[lattice[x][y][1].chromatophores_type]
                picy += (size + space)
            picx += (size + space)
        return Image.fromarray(pic.astype('uint8'), "RGB")

    def create_animated_gif(self, lattices, name= str(datetime.datetime.now())):
        frames = [self.to_picture(lattice) for lattice in lattices]
        #print(len(frames))
        if not os.path.exists(name):
            os.makedirs(name)
        maxwidth , maxheight = frames[-1].size

        for f in range(len(frames)):
            frames[f].save(f"{name}/{f}.png")
            width, height = frames[f].size
            borderh = int((maxheight - height) / 2)
            borderw = int((maxwidth-width) / 2)
            frames[f] = ImageOps.expand(frames[f], border=(borderw, borderh, borderw, borderh), fill=(122, 192, 234))
        frames[0].save(f'{name}.gif',
                       save_all=True, append_images=frames[1:],
                       optimize=False, duration=10)
