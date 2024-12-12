from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
from Matrices import blosum62, pam250


def escoger_loc_glob():
    metodo = input('Local (L) or Global (G):')
    if metodo == 'L':
        return False
    return True
def escoger_mat():
    metodo = input('Match/mismatch (MM) or Matrix (MT):')
    if metodo == 'MM':
        return False
    return True

def need_wunch_MM(seq_one, seq_two, match_score, mismatch_score, gap_score):

    #Creamos listas de letras sequencias
    lis_seq_one = [' ']
    lis_seq_two = [' ']
    for aa in seq_one:
        lis_seq_one.append(aa)
    for aa in seq_two:
        lis_seq_two.append(aa)
    #matriz de 0s del tamaño necesario (vamos a hacer que la seq1 este en vertical y la seq2 en hor
    matriz = np.matrix(np.zeros((len(lis_seq_one), len(lis_seq_two))))
    #seteamos primera fila y columna a los gaps
    set_gap_vertical = []
    set_gap_hor = []
    for i in range(len(lis_seq_one)):
        set_gap_vertical.append(i * gap_score)
    for i in range(len(lis_seq_two)):
        set_gap_hor.append(i * gap_score)

    set_gap_verticall = np.array(set_gap_vertical)
    set_gap_verticalll = set_gap_verticall.reshape(-1, 1)
    matriz[0, :] = set_gap_hor
    matriz[: , 0] = set_gap_verticalll
    dict_de_donde_viene = {}
    #iteramos para score matrix
    for i in range(1, matriz.shape[0]):
        for j in range(1, matriz.shape[1]):
            lista_temp_opciones = []
            #diagonal
            if lis_seq_one[i] == lis_seq_two[j]:
                lista_temp_opciones.append(match_score + matriz[i - 1, j - 1])
            elif lis_seq_one[i] != lis_seq_two[j]:
                lista_temp_opciones.append(mismatch_score + matriz[i - 1, j - 1])
            #vertical-horizontal
            lista_temp_opciones.append(matriz[i - 1, j] + gap_score)
            lista_temp_opciones.append(matriz[i, j - 1] + gap_score)
            matriz[i, j] = max(lista_temp_opciones)
            #ahora tenemos ya la matrix_score, pero falta saber de donde viene cada score
            #en este diccionario se ve de donde viene cada posicion
            dict_de_donde_viene[(i, j)] = lista_temp_opciones.index(max(lista_temp_opciones))
    #aliniacion
    seq_alin_one = []
    seq_alin_two = []
    #empezamos por abajo derecha
    a = len(lis_seq_one) -1
    b = len(lis_seq_two) -1
    while True:
        #diagonal
        if dict_de_donde_viene[a, b] == 0:
            seq_alin_one.insert(0, lis_seq_one[-1])
            lis_seq_one.pop()
            seq_alin_two.insert(0, lis_seq_two[-1])
            lis_seq_two.pop()
            a = a - 1
            b = b - 1
        #vertical
        elif dict_de_donde_viene[a, b] == 1:
            seq_alin_one.insert(0, lis_seq_one[-1])
            lis_seq_one.pop()
            seq_alin_two.insert(0, '-')
            a = a -1
        #horizontal
        elif dict_de_donde_viene[a, b] == 2:
            seq_alin_one.insert(0, '-')
            seq_alin_two.insert(0, lis_seq_two[-1])
            lis_seq_two.pop()
            b = b - 1

        #verificamos el bucle
        if a == 0 or b == 0 and (a, b) not in dict_de_donde_viene:
            break
    score = matriz[len(seq_one), len(seq_two)]
    print(seq_alin_one)
    print(seq_alin_two)
    print('Aligment score:', score)
    print(' ')
    print(' ')
    print(matriz)
    with open('outputalin.txt', 'w') as out:
        out.write(''.join(seq_alin_one))
        out.write('\n')
        out.write(''.join(seq_alin_two))
        out.write('\n')
        out.write('Aligment score: ')
        out.write(str(score))
    return seq_alin_one, seq_alin_two, score
def need_wunch_MT(seq_one, seq_two, matrixx, gap_score):

    #Creamos listas de letras sequencias
    lis_seq_one = [' ']
    lis_seq_two = [' ']
    for aa in seq_one:
        lis_seq_one.append(aa)
    for aa in seq_two:
        lis_seq_two.append(aa)
    #matriz de 0s del tamaño necesario (vamos a hacer que la seq1 este en vertical y la seq2 en hor
    matriz = np.matrix(np.zeros((len(lis_seq_one), len(lis_seq_two))))
    #seteamos primera fila y columna a los gaps
    set_gap_vertical = []
    set_gap_hor = []
    for i in range(len(lis_seq_one)):
        set_gap_vertical.append(i * gap_score)
    for i in range(len(lis_seq_two)):
        set_gap_hor.append(i * gap_score)

    set_gap_verticall = np.array(set_gap_vertical)
    set_gap_verticalll = set_gap_verticall.reshape(-1, 1)
    matriz[0, :] = set_gap_hor
    matriz[: , 0] = set_gap_verticalll
    dict_de_donde_viene = {}
    #iteramos para score matrix
    for i in range(1, matriz.shape[0]):
        for j in range(1, matriz.shape[1]):
            lista_temp_opciones = []
            #diagonal
            lista_temp_opciones.append(matrixx[lis_seq_one[i]][lis_seq_two[j]] + matriz[i - 1, j - 1])
            #vertical-horizontal
            lista_temp_opciones.append(matriz[i - 1, j] + gap_score)
            lista_temp_opciones.append(matriz[i, j - 1] + gap_score)
            matriz[i, j] = max(lista_temp_opciones)
            #ahora tenemos ya la matrix_score, pero falta saber de donde viene cada score
            #en este diccionario se ve de donde viene cada posicion
            dict_de_donde_viene[(i, j)] = lista_temp_opciones.index(max(lista_temp_opciones))
    #aliniacion
    seq_alin_one = []
    seq_alin_two = []
    #empezamos por abajo derecha
    a = len(lis_seq_one) -1
    b = len(lis_seq_two) -1
    while True:
        #diagonal
        if dict_de_donde_viene[a, b] == 0:
            seq_alin_one.insert(0, lis_seq_one[-1])
            lis_seq_one.pop()
            seq_alin_two.insert(0, lis_seq_two[-1])
            lis_seq_two.pop()
            a = a - 1
            b = b - 1
        #vertical
        elif dict_de_donde_viene[a, b] == 1:
            seq_alin_one.insert(0, lis_seq_one[-1])
            lis_seq_one.pop()
            seq_alin_two.insert(0, '-')
            a = a -1
        #horizontal
        elif dict_de_donde_viene[a, b] == 2:
            seq_alin_one.insert(0, '-')
            seq_alin_two.insert(0, lis_seq_two[-1])
            lis_seq_two.pop()
            b = b - 1

        #verificamos el bucle
        if a == 0 and b == 0:
            break

    score = matriz[len(seq_one), len(seq_two)]
    print(seq_alin_one)
    print(seq_alin_two)
    print('Aligment score:', score)
    print(' ')
    print(' ')
    print(matriz)
    with open('outputalin.txt', 'w') as out:
        out.write(''.join(seq_alin_one))
        out.write('\n')
        out.write(''.join(seq_alin_two))
        out.write('\n')
        out.write('Aligment score: ')
        out.write(str(score))
    return seq_alin_one, seq_alin_two, score



def smith_waterman_MM(seq_one, seq_two, match_score, mismatch_score, gap_score, num_out):
    # Creamos listas de letras sequencias
    lis_seq_one = [' ']
    lis_seq_two = [' ']
    for aa in seq_one:
        lis_seq_one.append(aa)
    for aa in seq_two:
        lis_seq_two.append(aa)
    # matriz de 0s del tamaño necesario (vamos a hacer que la seq1 este en vertical y la seq2 en hor
    matriz = np.matrix(np.zeros((len(lis_seq_one), len(lis_seq_two))))
    # iteramos para score matrix
    dict_de_donde_viene = {}
    for i in range(1, matriz.shape[0]):
        for j in range(1, matriz.shape[1]):
            lista_temp_opciones = []
            # diagonal
            if lis_seq_one[i] == lis_seq_two[j]:
                lista_temp_opciones.append(match_score + matriz[i - 1, j - 1])
            elif lis_seq_one[i] != lis_seq_two[j]:
                lista_temp_opciones.append(mismatch_score + matriz[i - 1, j - 1])
            # vertical-horizontal
            lista_temp_opciones.append(matriz[i - 1, j] + gap_score)
            lista_temp_opciones.append(matriz[i, j - 1] + gap_score)
            if max(lista_temp_opciones) >= 0:
                matriz[i, j] = max(lista_temp_opciones)
            else:
                matriz[i, j] = 0
            # ahora tenemos ya la matrix_score, pero falta saber de donde viene cada score
            # en este diccionario se ve de donde viene cada posicion
            dict_de_donde_viene[(i, j)] = lista_temp_opciones.index(max(lista_temp_opciones))
    #print(matriz)
    #print(dict_de_donde_viene)
    #voy a ver en la matriz, que numeros son mayor a 0
    lista_pos_mas_0 = []
    for i in range(1, matriz.shape[0]):
        for j in range(1, matriz.shape[1]):
            if matriz[i, j] > 0:
                lista_pos_mas_0.append([i, j, matriz[i, j]])
    #oredenamos de mayor a menor de score, la lista tiene los indices y como tercer elemento el score
    lista_pos_mas_0.sort(reverse= True, key = lambda x: x[-1])
    #vamos a ver cada aliniacion
    par_de_cadena = []
    for i, j, sco in lista_pos_mas_0:
        cad1 = [i]
        cad2 = [j]

        while True:
            if dict_de_donde_viene[i, j] == 0:
                if matriz[i - 1, j - 1] == 0:
                    cad1.insert(0, lis_seq_one[i])
                    cad1.insert(0, i)
                    cad2.insert(0, lis_seq_two[j])
                    cad2.insert(0, j)
                    break
                else:
                    cad1.insert(0, lis_seq_one[i])
                    cad2.insert(0, lis_seq_two[j])
                    i = i -1
                    j = j -1
            elif dict_de_donde_viene[i, j] == 1:
                if matriz[i - 1, j] == 0:
                    cad1.insert(0, lis_seq_one[i])
                    cad2.insert(0, '_')
                    cad1.insert(0, i)
                    cad2.insert(0, j)
                    break
                else:
                    cad1.insert(0, lis_seq_one[i])
                    cad2.insert(0, '_')
                    i = i -1
            elif dict_de_donde_viene[i, j] == 2:
                if matriz[i, j - 1] == 0:
                    cad1.insert(0, '_')
                    cad2.insert(0, lis_seq_two[j])
                    cad1.insert(0, i)
                    cad2.insert(0, j)
                    break
                else:
                    cad1.insert(0, '_')
                    cad2.insert(0, lis_seq_two[j])
                    j = j - 1
#solo insertamos
        # si la cadena tiene algo
        letras = 'ACDEFGHIKLMNPQRSTVWY'
        if any(letra in cad1 for letra in letras) and any(letra in cad2 for letra in letras):
            par_de_cadena.insert(0, cad2)
            par_de_cadena.insert(0, cad1)
    if num_out > len(par_de_cadena):
        a = []
        b = []
        align_score = []
        for i in range(0, len(par_de_cadena), 2):
            a.append(''.join(map(str, par_de_cadena[i])))
            print(''.join(map(str, par_de_cadena[i])))
            b.append(''.join(map(str, par_de_cadena[i + 1])))
            print(''.join(map(str, par_de_cadena[i + 1])))
            align_score.append(int(lista_pos_mas_0[int(i/2)][2]))
            print('Aligment Score:', int(lista_pos_mas_0[int(i/2)][2]))
            print()
        return a, b, align_score
    else:
        c = []
        d = []
        aligmentscore = []
        for i in range(0, num_out * 2, 2):
            c.append(''.join(map(str, par_de_cadena[i])))
            print(''.join(map(str, par_de_cadena[i])))
            d.append(''.join(map(str, par_de_cadena[i + 1])))
            print(''.join(map(str, par_de_cadena[i + 1])))
            aligmentscore.append(int(lista_pos_mas_0[int(i / 2)][2]))
            print('Aligment Score:', int(lista_pos_mas_0[int(i / 2)][2]))
            print()
        return c, d, aligmentscore




def smith_waterman_MT(seq_one, seq_two, matrixx, gap_score, num_out):
    # Creamos listas de letras sequencias
    lis_seq_one = [' ']
    lis_seq_two = [' ']
    for aa in seq_one:
        lis_seq_one.append(aa)
    for aa in seq_two:
        lis_seq_two.append(aa)
    # matriz de 0s del tamaño necesario (vamos a hacer que la seq1 este en vertical y la seq2 en hor
    matriz = np.matrix(np.zeros((len(lis_seq_one), len(lis_seq_two))))
    # iteramos para score matrix
    dict_de_donde_viene = {}
    for i in range(1, matriz.shape[0]):
        for j in range(1, matriz.shape[1]):
            lista_temp_opciones = []
            # diagonal
            lista_temp_opciones.append(matrixx[lis_seq_one[i]][lis_seq_two[j]] + matriz[i - 1, j - 1])
            # vertical-horizontal
            lista_temp_opciones.append(matriz[i - 1, j] + gap_score)
            lista_temp_opciones.append(matriz[i, j - 1] + gap_score)
            if max(lista_temp_opciones) >= 0:
                matriz[i, j] = max(lista_temp_opciones)
            else:
                matriz[i, j] = 0
            # ahora tenemos ya la matrix_score, pero falta saber de donde viene cada score
            # en este diccionario se ve de donde viene cada posicion
            dict_de_donde_viene[(i, j)] = lista_temp_opciones.index(max(lista_temp_opciones))
    print(matriz)
    print(dict_de_donde_viene)
    #voy a ver en la matriz, que numeros son mayor a 0
    lista_pos_mas_0 = []
    for i in range(1, matriz.shape[0]):
        for j in range(1, matriz.shape[1]):
            if matriz[i, j] > 0:
                lista_pos_mas_0.append([i, j, matriz[i, j]])
    #oredenamos de mayor a menor de score, la lista tiene los indices y como tercer elemento el score
    lista_pos_mas_0.sort(reverse= True, key = lambda x: x[-1])
    #vamos a ver cada aliniacion
    par_de_cadena = []
    for i, j, sco in lista_pos_mas_0:
        cad1 = [i]
        cad2 = [j]

        while True:
            if dict_de_donde_viene[i, j] == 0:
                if matriz[i - 1, j - 1] == 0:
                    cad1.insert(0, lis_seq_one[i])
                    cad1.insert(0, i)
                    cad2.insert(0, lis_seq_two[j])
                    cad2.insert(0, j)
                    break
                else:
                    cad1.insert(0, lis_seq_one[i])
                    cad2.insert(0, lis_seq_two[j])
                    i = i -1
                    j = j -1
            elif dict_de_donde_viene[i, j] == 1:
                if matriz[i - 1, j] == 0:
                    cad1.insert(0, lis_seq_one[i])
                    cad2.insert(0, '_')
                    cad1.insert(0, i)
                    cad2.insert(0, j)
                    break
                else:
                    cad1.insert(0, lis_seq_one[i])
                    cad2.insert(0, '_')
                    i = i -1
            elif dict_de_donde_viene[i, j] == 2:
                if matriz[i, j - 1] == 0:
                    cad1.insert(0, '_')
                    cad2.insert(0, lis_seq_two[j])
                    cad1.insert(0, i)
                    cad2.insert(0, j)
                    break
                else:
                    cad1.insert(0, '_')
                    cad2.insert(0, lis_seq_two[j])
                    j = j - 1
#solo insertamos
        # si la cadena tiene algo
        letras = 'ACDEFGHIKLMNPQRSTVWY'
        if any(letra in cad1 for letra in letras) and any(letra in cad2 for letra in letras):
            par_de_cadena.insert(0, cad2)
            par_de_cadena.insert(0, cad1)
    if num_out > len(par_de_cadena):
        a = []
        b = []
        alig_score = []
        for i in range(0, len(par_de_cadena), 2):
            a.append(''.join(map(str, par_de_cadena[i])))
            print(''.join(map(str, par_de_cadena[i])))
            b.append(''.join(map(str, par_de_cadena[i + 1])))
            print(''.join(map(str, par_de_cadena[i + 1])))
            alig_score.append(int(lista_pos_mas_0[int(i/2)][2]))
            print('Aligment Score:', int(lista_pos_mas_0[int(i/2)][2]))
            print()
        return a, b, alig_score
    else:
        c = []
        d = []
        aligg_score = []
        for i in range(0, num_out * 2, 2):
            c.append(''.join(map(str, par_de_cadena[i])))
            print(''.join(map(str, par_de_cadena[i])))
            d.append(''.join(map(str, par_de_cadena[i + 1])))
            print(''.join(map(str, par_de_cadena[i + 1])))
            aligg_score.append(int(lista_pos_mas_0[int(i / 2)][2]))
            print('Aligment Score:', int(lista_pos_mas_0[int(i / 2)][2]))
            print()
        return c, d, aligg_score

if __name__ == "__main__":
    ruta = 'Sequen_alin.fasta'

    seq = []
    id = []

    for sec in SeqIO.parse(ruta, "fasta"):
        seq.append(sec.seq)
        id.append(sec.id)
    seq_one = seq[0]
    seq_two = seq[1]

#True es global
    if escoger_loc_glob() == True:
        #Falso es MM
        if escoger_mat() == False:
            match_score = int(input('Match Score:'))
            mismatch_score = int(input('Mismatch score:'))
            gap_score = int(input('Gap score:'))
            need_wunch_MM(seq_one, seq_two, match_score, mismatch_score, gap_score)
        else:
            matr_q = input('BL or PAM:')
            if matr_q == 'BL':
                gap_score = int(input('Gap score:'))
                need_wunch_MT(seq_one, seq_two, blosum62, gap_score)
            elif matr_q == 'PAM':
                gap_score = int(input('Gap score:'))
                need_wunch_MT(seq_one, seq_two, pam250, gap_score)
    else:
        if escoger_mat() == False:
            match_score = int(input('Match Score:'))
            mismatch_score = int(input('Mismatch score:'))
            gap_score = int(input('Gap score:'))
            seq_out = int(input('Cuantas seq quiere como resultado?:'))
            smith_waterman_MM(seq_one, seq_two, match_score, mismatch_score, gap_score, seq_out)
        else:
            matr_que = input('BL or PAM:')
            if matr_que == 'BL':
                gap_score = int(input('Gap score:'))
                seq_out = int(input('Cuantas seq quiere como resultado?:'))
                smith_waterman_MT(seq_one, seq_two, blosum62, gap_score, seq_out)
            elif matr_que == 'PAM':
                gap_score = int(input('Gap score:'))
                seq_out = int(input('Cuantas seq quiere como resultado?:'))
                smith_waterman_MT(seq_one, seq_two, pam250, gap_score, seq_out)







