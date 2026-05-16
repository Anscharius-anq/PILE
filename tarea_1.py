"""
Tarea 1: Extracción de proteínas candidatas a partir de una secuencia de ADN en una URL

El script "Tarea_1.py" implementa un flujo secuencial para obtener proteínas candidatas 
a partir de una secuencia de ADN alojada en una URL. 

Nota: Aunque existen bibliotecas especializadas como Biopython que realizan este proceso
mucho más eficiente, este desarrollo se realiza con fines didacticos.

para un futuro, instalar biopython:
        $ python -m pip install biopython

    Documentación:
        https://biopython.org/wiki/Documentation
"""

import requests
import sys
import re



CODON_TAB = {
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "AGA": "R", "AGG": "R", "AGU": "S", "AGC": "S", "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L"
}

def url_content(url: str) -> str:
    # si requests no es suficiente usar trafilatura
    try:
        response = requests.get(url, timeout=10)  
        response.raise_for_status()               
        return response.text
    except requests.exceptions.HTTPError as e:
        raise RuntimeError(f"Error HTTP {e.response.status_code}: {url}") from e
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"No se pudo conectar con el servidor: {url}") from e


def get_dna(text: str) -> str:
    dna_lines = []
    # se corta el texto en líneas a una lista y se procesa cada línea
    for line in text.splitlines():
        # se eliminan espacios y se convierten a mayúsculas
        clean = re.sub(r"\s", "", line).upper()
        # se extraen secuencias ATCG dentro de la línea
        if clean and re.fullmatch(r"[ATCG]+", clean):
            dna_lines.append(clean)

    if not dna_lines:
        raise ValueError("No se encontró una secuencia de ADN válida en el texto.")

    return "".join(dna_lines)


def transcription(dna: str) -> str:
    return dna.replace("T", "U")


def translation(mrna: str, min_length: int = 2) -> list:
    # se buscan todos los índices de inicio de codones AUG en el mRNA
    start = [m.start() for m in re.finditer("AUG", mrna)]
    if not start:
        return []

    proteins = []
    # para cada índice de inicio, se traduce la secuencia hasta el codon de termino:
    for i in start:
        protein_seq = []
        stop_index = None  
        # Se itera cada  3 nucleotidos desde el inicio
        for j in range(i, len(mrna) - 2, 3):
            codon = mrna[j: j + 3]
            amino_acid = CODON_TAB.get(codon, "?") # se marca ? si no es un codon valido
            if amino_acid == "*": # indica codon de termino
                stop_index = j
                break
            # se agrega el aminoacido a la secuencia de la proteina
            protein_seq.append(amino_acid)
        # si no se encontro un codon de termino, se ignora esta secuencia
        if stop_index is None:  
            continue
        #
        if len(protein_seq) >= min_length:
            # se verifica que no haya una proteina contenido con el mismo codon de termino 
            # antes de agregar la proteina a la lista de proteinas
            is_nested = any(p["termino"] == stop_index for p in proteins) 
            if not is_nested:
                proteins.append({
                    "secuencia": "".join(protein_seq),
                    "inicio": i,
                    "termino": stop_index,
                    "longitud": len(protein_seq)
                })

    return proteins

# su unico uso es para formatear la salida de las proteinas
def format_proteins(proteins: list, top: int = None, sort: bool = False) -> str:
    if sort:
        proteins = sorted(proteins, key=lambda p: p["longitud"], reverse=True)
    if top:
        proteins = proteins[:top]

    result = []
    for protein in proteins:
        result.append(
            f"Inicio  : {protein['inicio']}\n"
            f"Longitud: {protein['longitud']} aa\n"
            f"Secuencia: {protein['secuencia']}\n"
            f"{'-' * 60}"
        )

    return "\n".join(result)




def get_proteins_from_url(url: str, top: int = None, sort: bool = False) -> list:
    """
    Obtiene proteínas candidatas a partir de una secuencia de ADN en una URL.

    Ejecuta el flujo completo: descarga el contenido, extrae la secuencia de ADN,
    la transcribe a mRNA y la traduce a proteínas candidatas.

    Args:
        url (str): URL del archivo de texto que contiene la secuencia de ADN.
        top (int, optional): Límite de proteínas a retornar. Si es None, retorna todas.
        sort (bool, optional): Si es True, ordena las proteínas de mayor a menor longitud.
                               Por defecto es False.

    Returns:
        list[dict]: Lista de proteínas encontradas. Cada proteína es un diccionario con:
            - secuencia (str): Cadena de aminoácidos en código de una letra.
            - inicio (int):    Índice del codón AUG en el mRNA.
            - termino (int):   Índice del codón de término en el mRNA.
            - longitud (int):  Número de aminoácidos de la proteína.

    Raises:
        RuntimeError: Si no se puede conectar al servidor o la respuesta HTTP es un error.
        ValueError:   Si no se encuentra una secuencia de ADN válida en el contenido.

    Example:
        >>> proteins = get_proteins_from_url("https://ejemplo.com/gen.txt", top=3, sort=True)
        >>> for p in proteins:
        ...     print(p["longitud"], p["secuencia"])
    """
    content = url_content(url)
    dna = get_dna(content)
    mrna = transcription(dna)
    proteins = translation(mrna)

    if sort:
        proteins = sorted(proteins, key=lambda p: p["longitud"], reverse=True)
    if top:
        proteins = proteins[:top]

    return proteins




def main():
    
    url = "https://raw.githubusercontent.com/Anscharius-anq/PILE/refs/heads/main/tarea_1.txt"

    try:
        content = url_content(url)
        dna = get_dna(content)
    except (RuntimeError, ValueError) as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    mrna = transcription(dna)
    proteins = translation(mrna) 

    print(format_proteins(proteins, 3, sort=True))


if __name__ == "__main__":
    main()