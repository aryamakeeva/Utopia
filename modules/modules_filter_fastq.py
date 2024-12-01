def gc_filter(seq: str, gc_bounds: tuple):
    # Проверяем, вдруг в gc_bounds подано одно значение int или float
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
        # , то устанавливаем нижнюю границу = 0, а верхнюю - поданное значение
    gc = (seq.count("C") + seq.count("G")) / len(seq) * 100
    # Если последовательность не удовлетворяет условиям - возвращаем EXCLUDE
    if gc < gc_bounds[0] or gc > gc_bounds[1]:
        return "EXCLUDE"


def length_filter(seq: str, length_bounds: tuple):
    # Проверяем, вдруг в length_bounds подано одно значение int или float
    if isinstance(length_bounds, (int, float)):
        length_bounds = (
            0,
            length_bounds,
        )  # , то устанавливаем нижнюю границу = 0,
        # а верхнюю - поданное значение
    if len(seq) < length_bounds[0] or len(seq) > length_bounds[1]:
        return "EXCLUDE"


def quality_filter(quality_string: str, threshold: int):
    unique_values = set(quality_string)
    total_quality = 0

    # считаем, сколько раз в строке
    # используется каждый символ из unique_values,
    # умножаем это количество на quality score согласно
    # шкале phred33, прибавляем в total_quality
    for i in unique_values:
        total_quality += quality_string.count(i) * (ord(i) - 33)

    # разделив суммарное качество на длину
    # последовательности, получим среднее качество
    total_quality = total_quality / len(quality_string)
    if total_quality < threshold:
        return "EXCLUDE"
