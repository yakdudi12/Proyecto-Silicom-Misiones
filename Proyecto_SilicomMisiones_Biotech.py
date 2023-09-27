import re
import math

#Diccionario species chromosomes

species_chromosomes =   {"Especie 1": 46 , "Especie 2" : 48, "Especie 3": 42, "Especie 4": 40}
print("El número de cromosomas de las especies son:")
print(species_chromosomes)
print()
#Abrir y almacenar las enzimas de corte en un diccionario clave-valor
enzimas_corte = r"C:\Users\Usuario\Desktop\Programacion Biologica\ProyectoSilicom\enzimas_corte.txt"
dic_enzimas = {}
with open(enzimas_corte, "r") as enzimas:
    for a in enzimas:
        clave, valor = a.strip().split(" ")
        dic_enzimas[clave] = valor

#Convertí las ubicaciones en int, porque me saltaba error más adelante en el código
dic_enzimas_int = {k: int(v) for k, v in dic_enzimas.items()}
print("Las enzimas cortan en:")
print(dic_enzimas_int)

puntosde_corte = list(dic_enzimas_int.values())
#print(puntosde_corte)
#Aca hice una pequeña modificación al .txt de las secuencias de dna, para facilitar almacenar esa información
print()
secuencias_adn = r"C:\Users\Usuario\Desktop\Programacion Biologica\ProyectoSilicom\secuencias_adn.txt"
dic_secuencias_dna = {}
with open(secuencias_adn, "r") as dnasec:
    for b in dnasec:
        clave2, valor2 = b.strip().split("=")
        dic_secuencias_dna[clave2] = valor2

#Me di cuenta que no todos los valores del .txt estan en mayúsculas
dna_seq = {k.upper(): v.upper() for k, v in dic_secuencias_dna.items()}
print("Las cadenas de ADN son:")
print(dna_seq)

#Supongo que la resolución que deberia hacer es esta
#Post-revision me di cuenta que podria resumir esto en menos lineas de código, y que lo repita automaticamente hasta que se quede sin valores. Para hacerlo más universal.
print("Las cadenas cortadas con las enzimas resultan en:")
def cortar_secuencias():
    for i in dna_seq.keys():
        corte1 = dna_seq[i][0:puntosde_corte[0]]
        corte2 = dna_seq[i][puntosde_corte[0]:puntosde_corte[1]]
        corte3 = dna_seq[i][puntosde_corte[1]:puntosde_corte[2]]
        corte4 = dna_seq[i][puntosde_corte[2]:puntosde_corte[3]]
        cortes = (corte1, corte2, corte3, corte4)
        print(("Secuencia:" + i))
        print(cortes)
        print()

cortar_secuencias()

#Calcular porcentaje de aminoácidos
proteinas = r"C:\Users\Usuario\Desktop\Programacion Biologica\ProyectoSilicom\secuencias_proteinas.txt"
prot_list = []
print("Las secuencias de las proteínas son:")
#Aca busqué en internet otra manera de eliminar el \n, porque luego tenia problemas con usarlo como lista.
with open(proteinas, "r") as archivo:
    for a in archivo:
        lista = a.strip() #la solución al final fue unicamente poner .strip , no eran necesario "\n"
        prot_list.append(lista)

print(prot_list)

aa_code = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","Q","S","T","W","Y","V"]
#Prácticamente lo mismo que ya estaba en el colab, solo que modificado para que hagalo que necesitaba
def calcular_porcentage_aminoacido(prot_list):
    for proteina in prot_list:
        longitud = len(proteina)
        aa_count = {}
        pp_ = {}
        for aa_search in aa_code:
            conteo = proteina.count(aa_search)
            if conteo > 0:
                aa_count[aa_search] = conteo
            percentage = (conteo / longitud)*100
            if percentage > 0:
                pp_[aa_search] = percentage

        print(f"Proteína: {proteina}")
        print(f"Longitud: {longitud}")
        print(f"La cantidad de aminoácidos es: {aa_count}")
        print(f"Y los porcentajes son: {pp_}")
        print()


calcular_porcentage_aminoacido(prot_list)

#Ahora para el final lo más sencillo
list_dna = list(dna_seq.values())
#Bases_complementarias= []
def generar_complementaria(bases):
    traductor = str.maketrans("ATCG", "TAGC")
    complementaridad = bases.translate(traductor)
    return complementaridad
bases = str(list_dna)
traduccion_adn = generar_complementaria(bases)
print(f"Las cadenas complementarias son: {traduccion_adn}")

#Primero debo traducir adn to arn, realmento no es necesario si lo hacemos del punto de vista de computadoras, pero lo hice
#Lo hice para seguir el orden biológico
print()
list_adn = list(dna_seq.values())
arn_list = []
arn_list_app = []
def traducir_a_arn(list_adn):
    global arn_list
    for secuencia in list_adn:
        transcripcion = secuencia.replace("T","U")
        arn_list.append(transcripcion)
        arn_list_app = transcripcion
        print(f"La traducción a ARN de las secuencia resultó en:")
        print(arn_list_app)
        print()
    return arn_list


traducir_a_arn(list_adn)

def dividir_en_grupos_de_tres(string):
    return [string[i:i+3] for i in range(0, len(string), 3)]
resultados = [] #Muy importante guardarlo en lista, porque más adelante me saltó error en la traducción a aminos.
for elemento in arn_list:
    grupos = dividir_en_grupos_de_tres(elemento)
    resultados.append(grupos)

for i, grupos in enumerate(resultados, 1):
    print(f"Secuencia {i}: {' '.join(grupos)}")
print()
Secuencia1 = resultados[0] #Probé si funcionaba bien.
codones = resultados
#print("Codones de la Secuencia 1:", Secuencia1)
print()
diccionario_codones = {"L": ["CUU", "CUC","CUA", "CUG","UUA","UUG"], "F": ["UUU","UUC"],"I":["AUU","AUC","AUA"],"M":"AUG","V":["GUU","GUC","GUA","GUG"],"S":["UCU","UCC","UCA","UCG","AGU","AGC"],"P":["CCU","CCC","CCA","CCG"],"T":["ACU","ACC","ACA","ACG"],"A":["GCU","GCC","GCA","GCG"],"Y":["UAU","UAC"],"H":["CAU","CAC"],"Q":["CAA","CAG"],"K":["AAA","AAG"],"E":["GAA","GAG"],"C":["UGU","UGC"],"W":"UGG","R":["CGU","CGC","CGA","CGG","AGA","AGG"],"G":["GGU","GGC","GGA","GGG"],"*":["UAA","UAG","UGA"],"N":["AAU","AAC"],"D":["GAU","GAC"]}
#El diccionario lo hice a mano, desconozco si la cantidad de "desconocidos" sea debido a algún missclick
lista_aminoacidos = []
all_resultados = str(resultados)
# Función para traducir una lista de codones a aminoácidos
def traducir_a_aminoacidos(codones, diccionario_codones):
    aminoacidos = []
    for codon in codones:
        codon_str = str(codon)  #Tuve problemas con el "codon", asi que tuve que cambiarlo a str
        for clave, valores in diccionario_codones.items():
            if codon_str in valores:
                aminoacidos.append(clave)
                break
        else:
            aminoacidos.append("-")
    return aminoacidos


lista_aminos = []
for codones in resultados:
    traducciones = traducir_a_aminoacidos(codones,diccionario_codones)
    if traducciones:
        lista_aminos.extend(traducciones)
    else:
        lista_aminos.append(codones)
#print(lista_aminos)
#Acá quedó medio feo, pero era porque ya era tarde y solo queria que quede bien en la consola
print("Las secuencias Traduciadas a Aminoacidos resultaron en:")
amino1 = traducir_a_aminoacidos(resultados[0],diccionario_codones)
amino2 =  traducir_a_aminoacidos(resultados[1],diccionario_codones)
amino3 = traducir_a_aminoacidos(resultados[2],diccionario_codones)
amino4 = traducir_a_aminoacidos(resultados[3],diccionario_codones)
secuencia1 = "".join(amino1)
secuencia2 = "".join(amino2)
secuencia3 = "".join(amino3)
secuencia4 = "".join(amino4)
print(f"Secuencia 1: {secuencia1}")
print(f"Secuencia 2: {secuencia2}")
print(f"Secuencia 3: {secuencia3}")
print(f"Secuencia 4: {secuencia4}")

