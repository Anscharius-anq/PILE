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


def url_content(url: str) -> str:
    try:
        response = requests.get(url, timeout=10)
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
    dna_lines = []
    # se corta el texto en líneas a una lista y se procesa cada línea
    for line in text.splitlines():
        # se eliminan espacios y se convierten a mayúsculas
        clean = re.sub(r"\s", "", line).upper()
        # se extraen secuencias ATCG dentro de la línea
        if clean and re.fullmatch(r"[ATCG]+", clean):
            dna_lines.append(clean)

    if not dna_lines:
        raise ValueError(
            "No se encontró una secuencia de ADN válida en el texto."
        )

    return "".join(dna_lines)


def dna_indices(dna: str) -> list:
    pattern = r"(?=(ATG(?:[ACGT]{3})*?T(?:GA|AA|AG)))"

    matches = re.finditer(pattern, dna)

    seen_stops = set()
    genes = []

    for m in matches:
        stop = m.end(1)

        if stop in seen_stops:
            continue

        seen_stops.add(stop)
        genes.append(m.span(1))

    return genes


def top_genes(genes, top: int = 3):
    length = lambda x: x[1] - x[0]
    largest_genes = sorted(genes, key=length, reverse=True)[:top]

    return largest_genes


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


def display_proteins(proteins: list, total_genes: int):
    print(
        f"Se muestran {len(proteins)} de {total_genes} proteínas encontradas\n"
    )

    for i, p in enumerate(proteins, start=1):
        print(
            f"Proteína {i}:\n"
            f"[{p.start}-{p.end - 3}] - "
            f"{p.length} aminoácidos\n"
            f"{p.sequence}\n"
        )


def get_proteins(dna, top, sort):
    pass


def main():

    URL = "https://raw.githubusercontent.com/Anscharius-anq/PILE/refs/heads/main/tarea_1.txt"

    try:
        content = url_content(URL)
        dna = get_dna(content)
    except (RuntimeError, ValueError) as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    genes = dna_indices(dna)
    proteins = translation(dna, top_genes(genes))
    display_proteins(proteins, len(genes))


if __name__ == "__main__":
    main()

