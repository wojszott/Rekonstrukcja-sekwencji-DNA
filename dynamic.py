import time
import xml.etree.ElementTree as ET
from typing import List, Tuple
import sys
sys.setrecursionlimit(3000)  # Zwiększenie limitu rekurencji – potrzebne przy dłuższych sekwencjach

# Funkcja wczytująca dane z pliku XML
def parse_xml(file_path: str) -> Tuple[List[str], List[str], int, str]:
    tree = ET.parse(file_path)
    root = tree.getroot()
    normal_chips = []        # Lista chipów zwykłych
    validation_chips = []    # Lista chipów walidujących
    start_seq = root.attrib.get("start", "")              # Sekwencja początkowa
    length = int(root.attrib.get("length", "0"))          # Docelowa długość sekwencji
    for probe in root.findall(".//probe"):
        pattern = probe.attrib.get('pattern', '')         # Wzorzec fragmentu (np. NXNXN lub NXNN)
        cells = [cell.text.strip() for cell in probe.findall("cell") if cell.text]
        if pattern.endswith('NN'):                        # Sprawdzenie, czy to chip walidujący
            validation_chips.extend(cells)
        else:
            normal_chips.extend(cells)
    return normal_chips, validation_chips, length, start_seq

# Funkcja sprawdzająca, czy dany chip pasuje do sekwencji na podanej pozycji
def is_match(chip: str, seq: Tuple[str], pos: int, is_normal: bool) -> bool:
    step = 2 if is_normal else 1  # Dla chipów zwykłych sprawdzamy co drugi znak (NXNX...), dla walidujących każdy (NXN → NN)
    for i in range(0, len(chip), step):
        if pos + i >= len(seq):  # Wyjście poza zakres sekwencji
            return False
        x = chip[i]
        y = seq[pos + i]
        if x == 'X' or y == 'X':  # 'X' traktujemy jako miejsce nieznane, pasuje do wszystkiego
            continue
        if x != y:
            return False
    return True

# Funkcja wstawiająca chip do sekwencji – nadpisuje tylko znaki różne od 'X'
def apply_chip(seq: Tuple[str], chip: str, pos: int) -> Tuple[str]:
    seq_list = list(seq)
    for i, c in enumerate(chip):
        if c != 'X':  # Ignorujemy 'X', nie nadpisujemy nieznanych pozycji
            seq_list[pos + i] = c
    return tuple(seq_list)

# Główna funkcja rekonstrukcji sekwencji (zachłanna z programowaniem dynamicznym)
def reconstruct_sequence_dp(start_seq: str, normal_chips: List[str], validation_chips: List[str],
                            length: int, max_gaps=10) -> Tuple[str, int]:
    # Uzupełniamy sekwencję do wymaganej długości znakami 'X'
    initial_seq = tuple(start_seq + 'X' * (length - len(start_seq)))

    chip_len_normal = len(normal_chips[0])
    chip_len_validation = len(validation_chips[0])
    memo = {}  # Słownik do zapamiętywania wyników dla tych samych stanów (memoizacja)

    # Funkcja licząca liczbę znanych (odkrytych) znaków w sekwencji
    def score(seq: Tuple[str]) -> int:
        return sum(1 for c in seq if c != 'X')

    # Rekurencyjna funkcja dynamiczna
    def dp(pos: int, parity: int, gaps: int, seq: Tuple[str]) -> Tuple[int, Tuple[str]]:
        if gaps > max_gaps:  # Zbyt dużo luk – zakończ
            return -1, seq
        if pos >= length:  # Koniec sekwencji – zwracamy wynik
            return score(seq), seq

        chip_len = chip_len_normal if parity == 1 else chip_len_validation
        key = (pos, parity, gaps, seq[pos:pos+chip_len])  # Klucz do memoizacji

        if key in memo:  # Jeśli wynik już policzony – zwracamy
            return memo[key]

        best_score = score(seq)  # Początkowy wynik (jeśli nic nie dopasujemy)
        best_seq = seq
        chips = normal_chips if parity == 1 else validation_chips
        matched_any = False  # Flaga informująca, czy coś dopasowaliśmy

        # Próbujemy dopasować każdy chip w danej grupie
        for chip in chips:
            if pos + chip_len > length:  # Chip nie zmieści się – pomijamy
                continue
            if is_match(chip, seq, pos, parity == 1):
                new_seq = apply_chip(seq, chip, pos)
                # Dla chipów zwykłych – przesuwamy się o 1
                # Dla walidujących – nie przesuwamy pozycji
                if parity == 1:
                    s, res_seq = dp(pos + 1, 1 - parity, gaps, new_seq)
                else:
                    s, res_seq = dp(pos, 1 - parity, gaps, new_seq)
                # Jeśli wynik lepszy – zapamiętujemy
                if s > best_score:
                    best_score = s
                    best_seq = res_seq
                matched_any = True

        # Jeśli nie udało się nic dopasować – próbujemy wykonać lukę
        if not matched_any and gaps < max_gaps:
            if parity == 1:
                s, res_seq = dp(pos + 1, 1 - parity, gaps + 1, seq)
            else:
                s, res_seq = dp(pos, 1 - parity, gaps + 1, seq)
            if s > best_score:
                best_score = s
                best_seq = res_seq

        # Zapamiętujemy najlepszy wynik
        memo[key] = (best_score, best_seq)
        return best_score, best_seq

    # Wywołanie funkcji dynamicznej z pozycji 0 i zaczynając od chipu zwykłego (parity = 1)
    best_score, best_seq = dp(0, 1, 0, initial_seq)
    return ''.join(best_seq), best_score

# Funkcja licząca, ile różnych chipów zostało użytych w odtworzonej sekwencji
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

# Główna część programu – odczyt danych, wywołanie algorytmu, pomiar czasu
if __name__ == "__main__":
    file_path = "daneOLD/450A.xml"  # <- Ścieżka do pliku z danymi (należy zmienić na własną)
    normal_chips, validation_chips, length, start_seq = parse_xml(file_path)
    max_G = length // 4  # Maksymalna liczba luk (heurystycznie: 1/4 długości)

    start_time = time.time()
    result, best_score = reconstruct_sequence_dp(start_seq, normal_chips, validation_chips, length, max_gaps=max_G)
    end_time = time.time()

    print("Odtworzona sekwencja:")
    print(result)

    num_X = result.count('X')
    found_len = len(result) - num_X
    percentage_found = (found_len / length) * 100

    print(f"Liczba odkrytych znaków: {best_score}")
    print(f"Długość sekwencji: {len(result)}")
    print(f"Liczba 'X': {num_X}")
    print(f"Odnaleziony % sekwencji: {percentage_found:.2f}%")
    print(f"Czas działania algorytmu: {end_time - start_time:.2f} sekund")

    # Statystyki użycia chipów
    used_normal = count_used_chips(result, normal_chips, is_normal=True)
    used_validation = count_used_chips(result, validation_chips, is_normal=False)
    total_chips = len(normal_chips) + len(validation_chips)
    used_chips = used_normal + used_validation
    usage_percent = (used_chips / total_chips) * 100

    print(f"\nUżyte fragmenty: {used_chips}/{total_chips} ({usage_percent:.2f}%)")
