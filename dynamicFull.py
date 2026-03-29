import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import List, Tuple, Dict
from functools import lru_cache


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
        if pattern.endswith('NN'):
            validation_chips.extend(cells)
        else:
            normal_chips.extend(cells)
    return normal_chips, validation_chips, length, start_seq


def is_match(chip: str, seq: List[str], pos: int, is_normal: bool) -> bool:
    step = 2 if is_normal else 1
    for i in range(0, len(chip), step):
        if pos + i >= len(seq):
            return False
        x = chip[i]
        y = seq[pos + i]
        if x == 'X' or y == 'X':
            continue
        if x != y:
            return False
    return True


def apply_chip(seq: str, chip: str, pos: int) -> str:
    seq = list(seq)
    for i, c in enumerate(chip):
        if c != 'X':
            seq[pos + i] = c
    return ''.join(seq)


def known_letters(chip: str, seq: str, pos: int, is_normal: bool) -> int:
    score = 0
    step = 2 if is_normal else 1
    for i in range(0, len(chip), step):
        if pos + i < len(seq) and seq[pos + i] != 'X' and chip[i] == seq[pos + i]:
            score += 1
    return score


def build_position_index(chips, length, is_normal):
    index = defaultdict(list)
    chip_len = len(chips[0]) if chips else 0
    for idx, chip in enumerate(chips):
        for pos in range(length - chip_len + 1):
            # jeśli chip pasuje do startowej sekwencji w tej pozycji
            # albo jeśli dopuszczamy luki, to zapisujemy tę pozycję
            # (można tu zrobić dokładne dopasowanie do 'X' na start_seq)
            index[pos].append(idx)
    return index

def reconstruct_sequence(start_seq: str, normal_chips: List[str], validation_chips: List[str],
                         length: int, max_gaps=10):

    # Uzupełnij start_seq do długości length znakami 'X'
    initial_seq = start_seq + 'X' * (length - len(start_seq))

    chip_len_normal = len(normal_chips[0]) if normal_chips else 0
    chip_len_validation = len(validation_chips[0]) if validation_chips else 0

    # Indeksujemy chipy po pozycji startu (przygotowujemy słownik: pozycja -> lista indeksów chipów pasujących na tej pozycji)
    normal_index = defaultdict(list)
    for idx, chip in enumerate(normal_chips):
        for pos in range(length - chip_len_normal + 1):
            if is_match(chip, initial_seq, pos, is_normal=True):
                normal_index[pos].append(idx)

    validation_index = defaultdict(list)
    for idx, chip in enumerate(validation_chips):
        for pos in range(length - chip_len_validation + 1):
            if is_match(chip, initial_seq, pos, is_normal=False):
                validation_index[pos].append(idx)

    best_result = {'score': 0, 'seq': initial_seq}

    @lru_cache(maxsize=None)
    def dfs(seq: str, pos: int, parity: int, gaps: int,
            used_normal_mask: int, used_validation_mask: int) -> int:
        # warunki zakończenia
        if gaps > max_gaps or pos >= length:
            known_letters = sum(1 for c in seq if c != 'X')
            if known_letters > best_result['score']:
                best_result['score'] = known_letters
                best_result['seq'] = seq
            return known_letters

        chip_list = normal_chips if parity == 1 else validation_chips
        chip_len = chip_len_normal if parity == 1 else chip_len_validation
        pos_index = normal_index if parity == 1 else validation_index

        max_score = -1
        best_seq_local = None

        candidates = pos_index.get(pos, [])

        # Heurystyka: faworyzuj chipy dopasowujące najwięcej znanych liter
        def known_letters_count(chip_idx):
            chip = chip_list[chip_idx]
            count = 0
            step = 2 if parity == 1 else 1
            for i in range(0, len(chip), step):
                p = pos + i
                if p < length and seq[p] != 'X' and chip[i] == seq[p]:
                    count += 1
            return -count  # minus dla sortowania rosnącego

        candidates = sorted(candidates, key=known_letters_count)

        for chip_idx in candidates:
            # Sprawdź, czy chip już użyty
            if parity == 1 and (used_normal_mask & (1 << chip_idx)) != 0:
                continue
            if parity == 0 and (used_validation_mask & (1 << chip_idx)) != 0:
                continue

            chip = chip_list[chip_idx]
            if pos + chip_len > length:
                continue
            if is_match(chip, seq, pos, parity == 1):
                new_seq = apply_chip(seq, chip, pos)
                new_used_normal_mask = used_normal_mask
                new_used_validation_mask = used_validation_mask
                if parity == 1:
                    new_used_normal_mask |= (1 << chip_idx)
                else:
                    new_used_validation_mask |= (1 << chip_idx)

                # Po zastosowaniu chipa przesuwamy pozycję o długość chipa
                score = dfs(new_seq, pos + chip_len, 1 - parity, gaps,
                            new_used_normal_mask, new_used_validation_mask)
                if score > max_score:
                    max_score = score
                    best_seq_local = new_seq

        # Jeśli nie znaleziono dopasowania, to dodaj lukę (przesuń o 1, zwiększ gaps)
        if max_score == -1 and gaps < max_gaps:
            score = dfs(seq, pos + 1, 1 - parity, gaps + 1, used_normal_mask, used_validation_mask)
            if score > max_score:
                max_score = score
                best_seq_local = seq

        if max_score > best_result['score']:
            best_result['score'] = max_score
            if best_seq_local is not None:
                best_result['seq'] = best_seq_local

        return max_score

    # Startujemy od pozycji 0, parity 1 (normal chips), 0 gaps, i mask 0 (nieużyte chipy)
    dfs(initial_seq, 0, 1, 0, 0, 0)

    return best_result['seq'][:length], best_result['score']


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

if __name__ == "__main__":
    file_path = "daneOLD/dane2.xml"
    normal_chips, validation_chips, length, start_seq = parse_xml(file_path)
    max_G = length // 3
    chip_len_n = len(normal_chips[0])
    chip_len_v = len(validation_chips[0])

    initial_seq = list(start_seq + 'X' * (length - len(start_seq)))

    start_time = time.time()
    memo = dict()
    result_seq,score = reconstruct_sequence(start_seq, normal_chips, validation_chips, length,max_gaps=max_G)
    end_time = time.time()

    result_str = ''.join(result_seq)
    print("Odtworzona sekwencja:" if 'X' not in result_str else "Nie udało się odtworzyć całej sekwencji.")
    print(result_str)
    print(f"Długość sekwencji: {len(result_str)}")
    print(f"Liczba 'X': {result_str.count('X')}")
    print(f"Odnaleziony % sekwencji: {(length - result_str.count('X')) / length * 100:.2f}%")
    print(f"\nCzas działania algorytmu: {end_time - start_time:.2f} sekund")

    used_normal = count_used_chips(result_str, normal_chips, is_normal=True)
    used_validation = count_used_chips(result_str, validation_chips, is_normal=False)
    total_chips = len(normal_chips) + len(validation_chips)
    used_chips = used_normal + used_validation
    usage_percent = (used_chips / total_chips) * 100

    print(f"\nUżyte fragmenty: {used_chips}/{total_chips} ({usage_percent:.2f}%)")