import time
import xml.etree.ElementTree as ET
from typing import List, Tuple
from collections import deque, defaultdict

# Wczytywanie danych z pliku XML: oddzielne listy dla zwykłych i walidujących chipów,
# długość oczekiwanej sekwencji oraz sekwencja początkowa
def parse_xml(file_path: str) -> Tuple[List[str], List[str], int, str]:
    tree = ET.parse(file_path)
    root = tree.getroot()
    normal_chips = []
    validation_chips = []
    start_seq = root.attrib.get("start", "")
    length = int(root.attrib.get("length", "0"))
    for probe in root.findall(".//probe"):
        pattern = probe.attrib.get('pattern', '')
        cells = [cell.text.strip() for cell in probe.findall("cell") if cell.text]
        if pattern.endswith('NN'):  # chipy walidujące mają wzorzec kończący się na NN
            validation_chips.extend(cells)
        else:  # pozostałe to chipy normalne
            normal_chips.extend(cells)
    return normal_chips, validation_chips, length, start_seq

# Sprawdza, czy dany chip pasuje do sekwencji na wskazanej pozycji
def is_match(chip: str, seq: str, pos: int, is_normal: bool) -> bool:
    step = 2 if is_normal else 1  # co drugi znak dla chipów normalnych
    for i in range(0, len(chip), step):
        if pos + i >= len(seq):
            return False
        x = chip[i]
        y = seq[pos + i]
        if x == 'X' or y == 'X':  # X oznacza nieznany znak, traktowany jako pasujący
            continue
        if x != y:
            return False
    return True

# Wprowadza chip do sekwencji, nadpisując nieznane znaki
def apply_chip(seq: str, chip: str, pos: int) -> str:
    seq = list(seq)
    for i, c in enumerate(chip):
        if c != 'X':
            seq[pos + i] = c
    return ''.join(seq)

# Tworzy indeks chipów według znanych liter z ich początku (domyślnie 3 pierwsze)
def build_index(chips: List[str], known_len: int = 3) -> dict:
    index = defaultdict(list)
    for chip in chips:
        key = ''.join([c for c in chip if c != 'X'])[:known_len]
        if key:
            index[key].append(chip)
        else:
            index[''].append(chip)
    return index

# Generuje klucz z fragmentu sekwencji — używany do indeksowania chipów
def generate_keys(seq: str, pos: int, window: int = 10, key_len: int = 3) -> List[str]:
    window_seq = seq[pos:pos + window]
    key = ''
    for i in range(0, len(window_seq), 2):
        if window_seq[i] != 'X':
            key += window_seq[i]
        if len(key) == key_len:
            break
    return [key] if key else ['']

# Główna funkcja odtwarzająca sekwencję
def reconstruct_sequence(start_seq: str, normal_chips: List[str], validation_chips: List[str],
                         length: int, max_gaps=10):

    # Uzupełniamy brakujące znaki 'X' do końcowej długości
    initial_seq = start_seq + 'X' * (length - len(start_seq))
    chip_lengths = (len(normal_chips[0]), len(validation_chips[0]))

    # Indeksujemy chipy
    normal_index = build_index(normal_chips, known_len=2)
    validation_index = build_index(validation_chips, known_len=2)

    queue = deque()
    visited = set()

    # parity = 1 -> następny chip to normalny, parity = 0 -> walidujący
    queue.append((initial_seq, 1, 1, 0))  # (sekwencja, pozycja, parzystość, liczba pominięć)

    best_seq = initial_seq
    best_score = sum(1 for c in best_seq if c != 'X')

    while queue:
        seq, pos, parity, gaps = queue.popleft()

        # Przerwij jeśli przekroczono dozwoloną liczbę luk lub koniec sekwencji
        if gaps > max_gaps or pos >= length - chip_lengths[0] + 2:
            continue

        # Aktualizacja najlepszego wyniku
        current_score = sum(1 for c in seq if c != 'X')
        if current_score > best_score:
            best_seq = seq
            best_score = current_score

        # Jeżeli sekwencja została w pełni odtworzona
        if 'X' not in seq:
            return seq[:length]

        visited_key = (pos, parity, gaps, seq[pos:pos+chip_lengths[0]])
        if visited_key in visited:
            continue
        visited.add(visited_key)

        chip_index = normal_index if parity == 1 else validation_index
        chip_len = chip_lengths[0] if parity == 1 else chip_lengths[1]

        # Generujemy klucz do szybkiego wyszukiwania potencjalnych chipów
        keys = generate_keys(seq, pos, window=4, key_len=3)
        matched = False

        for key in keys:
            for chip in chip_index.get(key, []):
                if pos + chip_len > length:
                    continue
                if is_match(chip, seq, pos, parity == 1):
                    matched = True
                    new_seq = apply_chip(seq, chip, pos)
                    # Po dopasowaniu przełączamy parzystość (czyli typ chipu)
                    if parity == 1:
                        queue.append((new_seq, pos + 1, 1 - parity, gaps))
                    else:
                        queue.append((new_seq, pos, 1 - parity, gaps))
                    break
            if matched:
                break

        # Jeżeli nie znaleziono dopasowania, dopuszczamy lukę
        if not matched and gaps < max_gaps:
            queue.append((seq, pos + 1, 1 - parity, gaps + 1))

    return best_seq[:length]

# Zlicza liczbę unikalnych chipów, które zostały wykorzystane do odtworzenia sekwencji
def count_used_chips(seq: str, chips: List[str], is_normal: bool) -> int:
    used = set()
    chip_len = len(chips[0])
    for i in range(len(seq) - chip_len + 1):
        for chip in chips:
            if chip not in used:
                if is_match(chip, seq, i, is_normal):
                    used.add(chip)
                    break
    return len(used)

# Przykładowe uruchomienie algorytmu
if __name__ == "__main__":
    file_path = "daneOLD/450A.xml"
    normal_chips, validation_chips, length, start_seq = parse_xml(file_path)
    max_G = length // 4  # maksymalna liczba dopuszczalnych luk

    start_time = time.time()
    result = reconstruct_sequence(start_seq, normal_chips, validation_chips, length, max_gaps=max_G)
    end_time = time.time()

    print("Odtworzona sekwencja:" if result != "RECONSTRUCTION_FAILED" else "Nie udało się odtworzyć sekwencji.")
    print(result)

    # Statystyki dotyczące jakości odtworzenia
    num_X = result.count('X')
    found_len = len(result) - num_X
    percentage_found = (found_len / length) * 100

    used_normal = count_used_chips(result, normal_chips, is_normal=True)
    used_validation = count_used_chips(result, validation_chips, is_normal=False)
    total_chips = len(normal_chips) + len(validation_chips)
    used_chips = used_normal + used_validation
    usage_percent = (used_chips / total_chips) * 100

    print(f"\nUżyte fragmenty: {used_chips}/{total_chips} ({usage_percent:.2f}%)")
    print(f"\nDługość sekwencji: {len(result)}")
    print(f"Liczba 'X': {num_X}")
    print(f"Odnaleziony % sekwencji: {percentage_found:.2f}%")
    print(f"\nCzas działania algorytmu: {end_time - start_time:.2f} sekund")
