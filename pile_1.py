from typing import NamedTuple
import requests
import re
import sys


class Protein(NamedTuple):
    start: int
    end: int
    length: int
    sequence: str


CODON_TAB = {
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "AGA": "R", "AGG": "R", "AGT": "S", "AGC": "S", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L"
}

URL = "https://raw.githubusercontent.com/Anscharius-anq/PILE/refs/heads/main/tarea_1.txt"

MIN_RESIDUES = 30 # puedes ajustar este valor según lo que consideres un gen significativo


def url_content(url: str) -> str:
    """
    Obtiene el contenido de una URL.

    Args:
        url (str): La URL de la que obtener el contenido.

    Returns:
        str: El contenido de la URL.
    """
    # se hace la solicitud HTTP y se maneja cualquier error que pueda ocurrir
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.text
    except requests.exceptions.HTTPError as e:
        raise RuntimeError(
            f"Error HTTP {e.response.status_code}: {url}"
        ) from e
    except requests.exceptions.RequestException as e:
        raise RuntimeError(
            f"No se pudo conectar con el servidor: {url}"
        ) from e


def get_dna(text: str) -> str:
    '''
    Extrae y devuelve la secuencia de ADN del texto dado en el argumento `text`. 
    Si no se encuentra una secuencia de ADN válida, se lanza un error.
    '''
    dna_lines = []
    # se corta el texto en líneas a una lista y se procesa cada línea
    for line in text.splitlines():
        clean = re.sub(r"\s", "", line).upper() # se eliminan espacios y se convierten a mayúsculas
        if clean and re.fullmatch(r"[ATCG]+", clean): # se extraen secuencias ATCG dentro de la línea
            dna_lines.append(clean)

    if not dna_lines:
        raise ValueError("No se encontró una secuencia de ADN válida en el texto.")

    return "".join(dna_lines)


def dna_indices(dna: str, /, min_residues: int = MIN_RESIDUES) -> list:
    '''
    Devuelve una lista de tuplas `[(start, end), ...]` con los indices de inicio (ATG) y termino
    (TAA, TAG o TGA) de cada gen encontrado en la secuencia de ADN dada en el argumento `dna`.

    El argumento `min_residues` se utiliza para filtrar genes que tengan una longitud menor a ese número
    de residuos. 

    '''

    pattern = r"(?=(ATG(?:[ACGT]{3})*?T(?:GA|AA|AG)))"
    matches = re.finditer(pattern, dna) 
    seen_stops = set() # evita genes solapados con el mismo codón de termino
    genes = [] # lista de tuplas con indices de inicio y termino de cada gen encontrado

    for m in matches:
        start, stop = m.span(1) 
        if stop in seen_stops: 
            continue

        length_aa = ((stop - start) // 3 ) - 1
        if length_aa >= min_residues:
            seen_stops.add(stop) # marca el codon de termino ocupado
            genes.append((start, stop))

    return genes


def top_genes(genes, top: int = 3):
    length = lambda x: x[1] - x[0] # función para calcular la longitud de un gen a partir de sus indices
    largest_genes = sorted(genes, key=length, reverse=True) # ordena los genes por longitud de mayor a menor

    return largest_genes[:top]


def translation(dna: str, genes: list):
    proteins = []

    for start, end in genes:
        protein_seq = []
        for j in range(start, end, 3):
            codon = dna[j : j + 3]
            amino_acid = CODON_TAB.get(codon, "?")

            if amino_acid == "*":
                break

            protein_seq.append(amino_acid)

        proteins.append(
            Protein(
                start=start,
                end=end,
                length=len(protein_seq),
                sequence="".join(protein_seq),
            )
        )
    return proteins


def main():

    
    try:
        content = url_content(URL)
        dna = get_dna(content)
    except (RuntimeError, ValueError) as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    genes = dna_indices(dna)
    proteins = translation(dna, top_genes(genes))
    
    # se muestra un resumen de los genes encontrados y sus proteínas traducidas    
    print(
        f"Se muestran {len(proteins)} de {len(genes)} proteínas encontradas\n"
    )

    for i, p in enumerate(proteins, start=1):
        print(f"Proteína {i} [pos:{p.start}]: {p.length} aa\n"
              f"{p.sequence}\n")


if __name__ == "__main__":
    main()

