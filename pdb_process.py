# pdb_process.py
# process .pdb file for use with spectralab
#
# suggested method is to inspect structure in pymol and apply anyd
#   desired symmetry operations. any file saved from pymol changes the
#   order of atomic coordinates, however, this will remedied below.


def pdb_process():
    """modify protein data base file for use with spectralab"""
    

    from tempfile import mkstemp
    from shutil import move
    from os import remove, close

    # for each chlorin molecule, the resName and (atom) names
    # chl(orophyll) a, chl b, pheophytin a,
    # bacteriochlorophyll a, bacteriopheophytin a

    atom_order = {'CLA':['G  ', 'CHA', 'CHB', 'CHC', 'CHD', 'NA ', 'C1A',
                         'C2A', 'C3A', 'C4A', 'CMA', 'CAA', 'CBA', 'CGA',
                         'O1A', 'O2A', 'NB ', 'C1B', 'C2B', 'C3B', 'C4B',
                         'CMB', 'CAB', 'CBB', 'NC ', 'C1C', 'C2C', 'C3C',
                         'C4C', 'CMC', 'CAC', 'CBC', 'ND ', 'C1D', 'C2D',
                         'C3D', 'C4D', 'CMD', 'CAD', 'OBD', 'CBD', 'CGD',
                         'O1D', 'O2D', 'CED', 'C1 ', 'C2 ', 'C3 ', 'C4 ',
                         'C5 ', 'C6 ', 'C7 ', 'C8 ', 'C9 ', 'C10', 'C11',
                         'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18',
                         'C19', 'C20'],
                  'CHL':['G  ', 'CHA', 'CHB', 'CHC', 'CHD', 'NA ', 'C1A',
                         'C2A', 'C3A', 'C4A', 'CMA', 'CAA', 'CBA', 'CGA',
                         'O1A', 'O2A', 'NB ', 'C1B', 'C2B', 'C3B', 'C4B',
                         'CMB', 'CAB', 'CBB', 'NC ', 'C1C', 'C2C', 'C3C',
                         'C4C', 'CMC', 'OMC', 'CAC', 'CBC', 'ND ', 'C1D',
                         'C2D', 'C3D', 'C4D', 'CMD', 'CAD', 'OBD', 'CBD',
                         'CGD', 'O1D', 'O2D', 'CED', 'C1 ', 'C2 ', 'C3 ',
                         'C4 ', 'C5 ', 'C6 ', 'C7 ', 'C8 ', 'C9 ', 'C10',
                         'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17',
                         'C18', 'C19', 'C20'],
                  'PHO':['CHA', 'CHB', 'CHC', 'CHD', 'NA ', 'C1A', 'C2A',
                         'C3A', 'C4A', 'CMA', 'CAA', 'CBA', 'CGA', 'O1A',
                         'O2A', 'NB ', 'C1B', 'C2B', 'C3B', 'C4B', 'CMB',
                         'CAB', 'CBB', 'NC ', 'C1C', 'C2C', 'C3C', 'C4C',
                         'CMC', 'CAC', 'CBC', 'ND ', 'C1D', 'C2D', 'C3D',
                         'C4D', 'CMD', 'CAD', 'OBD', 'CBD', 'CGD', 'O1D',
                         'O2D', 'CED', 'C1 ', 'C2 ', 'C3 ', 'C4 ', 'C5 ',
                         'C6 ', 'C7 ', 'C8 ', 'C9 ', 'C10', 'C11', 'C12',
                         'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19',
                         'C20'],
                  'BCL':['G  ', 'CHA', 'CHB', 'CHC', 'CHD', 'NA ', 'C1A',
                         'C2A', 'C3A', 'C4A', 'CMA', 'CAA', 'CBA', 'CGA',
                         'O1A', 'O2A', 'NB ', 'C1B', 'C2B', 'C3B', 'C4B',
                         'CMB', 'CAB', 'OBB', 'CBB', 'NC ', 'C1C', 'C2C',
                         'C3C', 'C4C', 'CMC', 'CAC', 'CBC', 'ND ', 'C1D',
                         'C2D', 'C3D', 'C4D', 'CMD', 'CAD', 'OBD', 'CBD',
                         'CGD', 'O1D', 'O2D', 'CED', 'C1 ', 'C2 ', 'C3 ',
                         'C4 ', 'C5 ', 'C6 ', 'C7 ', 'C8 ', 'C9 ', 'C10',
                         'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17',
                         'C18', 'C19', 'C20'],
                  'BPH':['CHA', 'CHB', 'CHC', 'CHD', 'NA ', 'C1A', 'C2A',
                         'C3A', 'C4A', 'CMA', 'CAA', 'CBA', 'CGA', 'O1A',
                         'O2A', 'NB ', 'C1B', 'C2B', 'C3B', 'C4B', 'CMB',
                         'CAB', 'OBB', 'CBB', 'NC ', 'C1C', 'C2C', 'C3C',
                         'C4C', 'CMC', 'CAC', 'CBC', 'ND ', 'C1D', 'C2D',
                         'C3D', 'C4D', 'CMD', 'CAD', 'OBD', 'CBD', 'CGD',
                         'O1D', 'O2D', 'CED', 'C1 ', 'C2 ', 'C3 ', 'C4 ',
                         'C5 ', 'C6 ', 'C7 ', 'C8 ', 'C9 ', 'C10', 'C11',
                         'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18',
                         'C19', 'C20']}

    file_name = input('enter filename: ')
    print('\nextracting data...')

    # pull only relevant data from file
    # store resName in pigments and resSeq in names (keys are resSeq)
    with open(file_name, 'r') as old_file:
        coords = [line for line in old_file if line[:6] == 'HETATM'
                  and line[17:20] in ('CLA', 'CHL', 'PHO', 'BCL', 'BPH')]
    names = {}
    pigments = {}
    for line in coords:
        if line[22:26] not in names:
            names[line[22:26]] = line[22:26]
            pigments[line[22:26]] = line[17:20]
    print(f'pigment names: {list(names.values())}\n')
    
    # resSeq numbers are not useful for identification
    # change names as desired
    cmnd = input('rename a pigment? <y> or <n> ')
    
    if cmnd == 'y':
        while True:
            old = input('replace: (<enter> to stop) ')
            # check for stop
            if old == '':
                break
            # add whitespace if needed
            length = len(old)
            if length < 4:
                old = old.rjust(4)
            # check if available to change
            if old in names.values():
                while True:
                    new = input('with: ')
                    length = len(new)
                    # check length, if too long try again
                    if length > 4:
                        print('keep names to four or less characters')
                        continue
                    elif length <= 4:
                        new = new.rjust(4)
                        key = list(names.keys())[list(names.values())
                                                 .index(old)]
                        names[key] = new
                        print(f'pigment names: {list(names.values())}\n')
                        break
            # if not available name, try again until stop
            else:
                print('name not found')
    # reorder pigments if desired
    cmnd = input('\nreorder pigments? <y> or <n> ')
    if cmnd == 'y':
        while True:
            new_order = input('new order: (<enter> for example) ')
            # explain input style
            if new_order == '':
                print(('\ngiven names are ["Chlb", "Chla", "Chlc"]'
                       '\nto switch order input [2, 1, 3]'
                       '\nand result will be ["Chla", "Chlb", "Chlc"]\n'))
                continue
            # convert input string to list
            new_order = eval(new_order)
            # check 
            if len(new_order) != len(names):
                print(f'there are {len(names)} pigments')
                continue
            # list of keys reordered
            new_order = [list(names.keys())[ind]
                         for ind in [int(num) - 1 for num in new_order]]
            # names dictionary reordered
            names = {key:names[key] for key in new_order}
            break
    print('\nplease wait...')
    # create temp file 
    fh, abs_path = mkstemp()
    for i, line in enumerate(coords):
        # apply renaming
        if names[line[22:26]] != line[22:26]:
            coords[i] = line[:22] + names[line[22:26]] + line[26:]
    # open temp file for writing
    with open(abs_path, 'a') as new_file: 
        for key in names:
            for atom in atom_order[pigments[key]]:
                for line in coords:
                    if line[13:16] == atom and line[22:26] == names[key]:
                        new_file.write(line)
    close(fh)
    # remove original file
    remove(file_name)
    # move new file
    move(abs_path, file_name)
    print('done')