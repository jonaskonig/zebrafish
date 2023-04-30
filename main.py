from Fish import Lattice


if __name__ == '__main__':

    #Matrix respresanting the adhesion coefficients Jxx, Jmm,Jxm
    matrix = [4,1,-1]
    #stammzellen dichtte
    dense = 0.6
    # elastic energy weight
    elastic_energy_weight = 7
    # mean cell lofe
    mean_cell_life = 30000
    # simmulation duratuin
    duration = 42
    #lattice size
    lattice_size = [20,100]
    # growth interval y
    growth_interval_y = 1
    # groth intervall x
    growth_interval_x = 1
    #special cell move
    special_cell_move = False
    # Cell birth rate
    cell_birth_rate =3
    fish = Lattice(lattice_size,dense,elastic_energy_weight,matrix,mean_cell_life,growth_interval_y,duration,growth_interval_x,cell_birth_rate,special_cell_move)
    fish.create_animated_gif(fish.run_sim())
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
