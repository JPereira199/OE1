#!/usr/bin/env python3

import csv
from collections import defaultdict
from itertools import combinations, islice

# Cambia 'blast_output.tsv' por el nombre de tu archivo de BLAST
input_file = 'blast_output.tsv'

# Diccionario para almacenar los mapeos por read
reads = defaultdict(list)

# Leer el archivo de BLAST
with open(input_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        # Extraer la información necesaria
        read_id = row[0]
        block = row[1]
        start = int(row[6])
        end = int(row[7])

        # Agregar la información al diccionario
        reads[read_id].append({'block': block, 'start': start, 'end': end})

# Diccionarios para contar las combinaciones
pair_counts = defaultdict(int)
triplet_counts = defaultdict(int)

# Procesar cada read
for read_id, mappings in reads.items():
    # Ordenar los mapeos por la posición de inicio en el read
    sorted_mappings = sorted(mappings, key=lambda x: x['start'])

    # Obtener la lista ordenada de bloques
    ordered_blocks = [m['block'] for m in sorted_mappings]

    # Generar las combinaciones de a dos (pares) consecutivas
    for i in range(len(ordered_blocks) - 1):
        pair = (ordered_blocks[i], ordered_blocks[i + 1])
        pair_counts[pair] += 1

    # Generar las combinaciones de a tres (tripletes) consecutivas
    for i in range(len(ordered_blocks) - 2):
        triplet = (ordered_blocks[i], ordered_blocks[i + 1], ordered_blocks[i + 2])
        triplet_counts[triplet] += 1

# Escribir los resultados en archivos de salida
with open('pair_counts.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Block1', 'Block2', 'Count'])
    for pair, count in pair_counts.items():
        writer.writerow([pair[0], pair[1], count])

with open('triplet_counts.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Block1', 'Block2', 'Block3', 'Count'])
    for triplet, count in triplet_counts.items():
        writer.writerow([triplet[0], triplet[1], triplet[2], count])

print("Análisis completado. Los resultados están en 'pair_counts.tsv' y 'triplet_counts.tsv'.")
